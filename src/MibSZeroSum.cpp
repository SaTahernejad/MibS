/*===========================================================================*/
/* This file is part of a Mixed Integer Bilevel Solver                       */
/* developed using the BiCePS Linear Integer Solver (BLIS).                  */
/*                                                                           */
/* Authors: Scott DeNegre, Lehigh University                                 */
/*          Ted Ralphs, Lehigh University                                    */
/*          Sahar Tahernajad, Lehigh University                              */
/*                                                                           */
/* Copyright (C) 2007-2017 Lehigh University, Scott DeNegre, and Ted Ralphs. */
/* All Rights Reserved.                                                      */
/*                                                                           */
/* This software is licensed under the Eclipse Public License. Please see    */
/* accompanying file for terms.                                              */
/*===========================================================================*/

#include "MibSModel.hpp"

#if  COIN_HAS_MPI
#include "AlpsKnowledgeBrokerMPI.h"
#else
#include "AlpsKnowledgeBrokerSerial.h"
#endif

//#############################################################################
MibSZeroSum::MibSZeroSum()
{
    model_ = 0;
    lSolver_ = 0;
    isAlgStarted_ = false;
    isFirstPhaseFinished_ = false;
}

//#############################################################################
MibSZeroSum::~MibSZeroSum()
{
    if(lSolver_) delete lSolver_;
}

//#############################################################################
void
MibSZeroSum::solveZeroSum(MibSModel *mibs, double *sol)
{
    
    model_ = mibs;
    isAlgStarted_ = true;

    int numCols(model_->getNumOrigVars());
    int lCols(model_->getLoweDim());
    bool isUnbounded(false), foundOptimal(false);

    double *lColLb = new double[lCols];
    double *lColUb = new double[lCols];
    CoinZeroN(lColLb, lCols);
    CoinZeroN(lColUb, lCols);

    double *optimalSol = new double[numCols];
    CoinZeroN(optimalSol, numCols);

    double *LLSol = new double[lCols];
    memcpy(LLSol, sol, sizeof(double) * lCols);
    
    // first-phase of the algorithm: finding a bound on lower-level
    //variables so that the lower-level problem will be feasible for
    //all feasible first-level variables 
    doFirstPhase(LLSol, lColLb, lColUb, optimalSol, isUnbounded, foundOptimal);

    delete [] lColLb;
    delete [] lColUb;
    delete [] optimalSol;
    delete [] LLSol;

}

//#############################################################################
void
MibSZeroSum::doFirstPhase(double *LLSol, double *lColLb, double *lColUb,
			  double *optimalSol, bool &isUnbounded, bool &foundOptimal)
{
    OsiSolverInterface * oSolver = model_->getSolver();

    bool isFirstPhaseFinished(false), hasPositiveArtCol(false);
    int i;
    int colIndex(0);
    int numAddedArtCols(0), numAddedRows(0), numAddedLrows(0), numAdddedURows(0);
    double infinity(oSolver->getInfinity());
    int numCols(model_->getNumOrigVars());
    int numRows(model_->getNumOrigCons());
    int lCols(model_->getLoweDim());
    int lRows(model_->getLowerRowNum());
    int newUCols(model_->getUpperDim);
    int newLCols(lCols + lRows);
    int initURows(model_->getOrigUpperRowNum());
    int newLRows(2 * lRows);
    int newNumCols(numCols + lRows);
    int initNumRows(initURows + newLRows);
    double *rowLB(model_->getOrigRowLb());
    double *rowUB(model_->getOrigRowUb());
    const char *rowSense(model_->getOrigRowSense());
    int *lColInd(model_->getLowerColInd());
    int *lRowInd(model_->getLowerRowInd());
    CoinPackedMatrix matrix = *model_->getOrigConstCoefMatrix();
    /** matrix of upper-level vars in lower-level problem **/
    CoinPackedMatrix *matrixA2;
    
    //finding the bounds on the lower-level variables
    //for the first step
    /*double *lColLb = new double[lCols];
    double *lColUb = new double[lCols];
    CoinZeroN(lColLb, lCols);
    CoinZeroN(lColUb, lCols);*/
    findBoundsOnLLCols(lColLb, lColUb, LLSol);

    //find upper bound of the artificial variables
    //double *artColsUB = findUBOnArtCols();

    //set a reasonable value for big M
    double bigM = findBigM();

    //setting the components of first phase model which are
    //the same for all steps

    //Set up lp solver
    OsiClpSolverInterface lpSolver;
    lpSolver.getModelPtr()->setDualBound(1.0e10);
    lpSolver.messageHandler()->setLogLevel(0);

    //set col type
    char * newColType = new char[newNumCols];
    memcpy(newColType, model_->colType_, numCols);
    for(i = numCols; i < newNumCols; i++){
	newColType[i] = 'C';
    }

    //set initial matrix, row bounds, row senses and lower and upper row indices
    CoinShallowPackedVector origVec;
    CoinPackedVector modifiedVec;
    CoinPackedVector appendedVec;
    CoinPackedMatrix *initMatrix = new CoinPackedMatrix(false, 0, 0);
    double *initRowLB = new double[initNumRows];
    double *initRowUB = new double[initNumRows];
    char *initRowSense = new char[initNumRows];
    int *newLRowInd = new int[newLRows];
    int *initURowInd = new int[initUROws];
    
    matrix.reverseOrdering();
    initMatrix->setDimensions(0, newNumCols);
    numAddedArtCols = 0;
    numAddedRows = 0;
    numAddedLrows = 0;
    numAdddedURows = 0;
    for(i = 0; i < numRows; i++){
	origVec = matrix.getVector(i);
	if(binarySearch(0, lRows - 1, i, lRowInd) < 0){
	    initMatrix->appendRow(origVec);
	    initRowLB[numAddedRows] = rowLB[i];
	    initRowUB[numAddedRows] = rowUB[i];
	    initRowSense[numAddedRows] = rowSense[i];
	    initURowInd[numAdddedURows] = numAddedRows;
	    numAddedRows++;
	    numAdddedURows++;
	    origVec.clear();
	}
	else{
	    if(rowSense[i] == 'L'){
		appendedVec.insert(numCols + numAddedArtCols, -1);
	    }
	    else if(origRowSense[i] == 'G'){
		appendedVec.insert(numCols + numAddedArtCols, 1);
	    }
	    else{
		assert(0);
	    }
	    modifiedVec = origVec;
	    modifiedVec.append(appendedVec);
	    initMatrix->appendRow(modifiedVec);
	    initRowLB[numAddedRows] = rowLB[i];
	    initRowUB[numAddedRows] = rowUB[i];
	    initRowSense[numAddedRows] = rowSense[i];
	    newLRowInd[numAddedLrows] = numAddedRows;
	    numAddedRows++;
	    numAddedLrows++;
	    initMatrix->appendRow(modifiedVec);
	    if(rowSense[i] == 'L'){
		initRowLB[numAddedRows] = rowUB[i];
		initRowUB[numAddedRows] = infinity;
		initRowSense[numAddedRows] = 'G';
	    }
	    else{
		initRowUB[numAddedRows] = rowLB[i];
		initRowLB[numAddedRows] = -1 * infinity;
		initRowSense[numAddedRows] = 'L';
	    }
	    newLRowInd[numAddedLrows] = numAddedRows;
	    numAddedRows++;
	    numAddedLrows++;
	    origVec.clear();
	    appendedVec.clear();
	    modifiedVec.clear();
	    numAddedArtCols++;
	}
    }
    newMatrix->reverseOrdering();

    //set col bounds
    double *newColLB = new double[newNumCols];
    double *newColUB = new double[newNumCols];
    memcpy(newColLB, model_->getOrigColLb(), sizeof(double) * numCols);
    memcpy(newColUB, model_->getOrigColUb(), sizeof(double) * numCols);
    for(i = numOrigCols; i < numCols; i++){
	newColLB[i] = 0;
	//sahar:modify it
	newColUB[i] = infinity;
    }

    //set the upper- and lower-level objectives
    double *newUObjCoeffs = new double[newNumCols];
    double *newLObjCoeffs = new double[newLCols];
    memcpy(newObjCoeffs, oSolver->getObjCoefficients(), sizeof(double) * numCols);
    memcpy(newLObjCoeffs, oSolver->getLowerObjCoeffs(), sizeof(double) * lCols);
    CoinFill(&newUObjCoeffs[numCols], &newUObjCoeffs[newNumCols], -1 * bigM);
    CoinFill(&lObjCoeffs[lCols], &lObjCoeffs[newLCols], bigM);

    //set upper and lower column indices
    int *newUColInd = new int[newUCols];
    int *newLColInd = new int[newLCols];
    memcpy(newUColInd, model_->getUpperColInd(), sizeof(int) * newUcols);
    memcpy(newLColInd, lColInd, sizeof(int) * lCols);
    CoinIotaN(newLColInd, lRows, numCols);

    //set structure row indices
    int *initStructRowInd = new int[initNumRows];
    CoinIotaN(initStructRowInd, initnumRows, 0);

    int newURows(0), newNumRows(0), newStructRowNum(0);
    double objVal(0.0), realLObj(0.0), lObj(0.0);
    bool isFirstStep(true);
    bool hasValidBound(false);
    int argc = 1;
    char **argv = new char *[1];
    argv[0] = (char *) "mibs";
    
    while(!isFirstPhaseFinished){
	for(i = 0; i < lCols; i++){
	    colIndex = lColInd[i];
	    newColLB[colIndex] = lColLb[i];
	    newColUB[colIndex] = lColUb[i];
	}

	MibSModel *firstPhaseModel = new MibSModel();
	firstPhaseModel->setSolver(&lpSolver);
	
	if(!hasValidBound){
	    firstPhaseModel->loadAuxiliaryData(newLCols, newLRows,
					       newLColInd, newLRowInd, 1,
					       newLObjCoeffs, newUCols, initURows,
					       newUColInd, initURowInd, initRowNum,
					       initStructRowInd, 0, NULL);
	    firstPhaseModel->loadProblemData(initMatrix, newColLB, newColUB, newUObjCoeffs,
					     initRowLB, initRowUB, newColType, 1, infinity,
					     initRowSense);
	}
	else{
	    newURows = initUrows + 1;
	    newNumRows = initNumRows + 1;
	    int *newURowInd = new int[newURows];
	    memcpy(newURowInd, initURowInd, sizeof(int) * initURows);
	    newURowInd[newURows] = newNumRows;
	    int *newRowLB = new double[newURows];
	    int *newRowUB = new double[newURows];
	    char *newRowSense = new char[newURows];
	    memcpy(newURowLB, initRowLB, sizeof(double) * initNumRows);
	    memcpy(newURowUB, initRowUB, sizeof(double) * initNumRows);
	    memcpy(newURowSense, initRowSense, initNumRows);
	    newRowLB[newNumRows] = objVal;
	    newRowUB[newNumRows] = infinity;
	    newRowSense[newNumRows] = 'G';
	    CoinPackedMatrix newMatrix = *new CoinPackedMatrix(false, 0, 0);
	    initMatrix->setDimensions(0, newNumCols);
	    newMatrix = *initMatrix;
	    for(i = 0; i < newNumCols; i++){
		appendedVec.insert(i, newUObjCoeffs[i]);
	    }
	    newMatrix.appendRow(appendedVec);
	    delete appendedVec;
	    int *newStructRowInd = new int[newNumRows];
	    memcpy(newStructRowInd, initStructRowInd, sizeof(int) * initNumRows);
	    newStructRowInd[newNumRows] = newNumRows;
	    firstPhaseModel->loadAuxiliaryData(newLCols, newLRows,
					       newLColInd, newLRowInd, 1,
					       newLObjCoeffs, newUCols, newURows,
					       newUColInd, newURowInd, newRowNum,
					       newStructRowInd, 0, NULL);
	    firstPhaseModel->loadProblemData(newMatrix, newColLB, newColUB, newUObjCoeffs,
					     newRowLB, newRowUB, newColType, 1, infinity,
					     newRowSense);
	}

#ifdef  COIN_HAS_MPI
	AlpsKnowledgeBrokerMPI broker(argc, argv, *firstPhaseModel);
#else
	AlpsKnowledgeBrokerSerial broker(argc, argv, *firstPhaseModel);
#endif
	broker.search(firstPhaseModel);

	if(firstPhaseModel->getNumSolutions() > 0){
	    double *firstPhaseSol = firstPhaseModel->incumbent();
	    objVal = broker.getBestNode()->getQuality();
	}
	else{
	    assert(0);
	}

	hasPositiveArtCol = false;
	for(i = 0; i < lRows; i++){
	    if(firstPhaseSol[i + numCols] > etol){
		hasPositiveArtCol = true;
		break;
	    }
	}

	if(!lSolver_){
	    initialSetUpLowerSolver(matrixA2);
	}

	//correct the lower-level solver
	modifyLowerSolver(firstPhaseSol, matrixA2);

	//solve the lower-level problem
	solveLowerProblem();

	if(hasPositiveArtCols){
	    if(!lSolver_->isProvenOptimal()){
		isUnbounded = true;
		isFirstPhaseFinished = true;
	    }
	    else{
		findLowerSol(firstPhaseSol, LLSol);
		findBoundsOnLLCols(lColLb, lColUb, LLSol);
	    }
	}
	else{
	    isFirstPhaseFinished = true;
	    if(lSolver_->isProvenOptimal()){
		realLObj = lSolver_->getObjValue();
		lObj = model_->bS_->getLowerObj(firstPhaseSol, 1);
		if(fabs(realLObj - lObj) <= etol){
		    foundOptimal = true;
		}
		else{
		    indLowerSol(firstPhaseSol, LLSol);
		    findBoundsOnLLCols(lColLb, lColUb, LLSol);
		}
	    }
	    else{
		assert(0);
	    }
		    
	}
    }

}

//#############################################################################
/*double*
MibSZeroSum::findUBOnArtCols()
{
    std::string feasCheckSolver =
	model_->MibSPar_->entry(MibSParams::feasCheckSolver);
    int i, j;
    int index(0);
    int numRows(model_->getNumOrigCons());
    int numCols(model_->getNumOrigVars());
    int lRows(model_->getLowerRowNum());
    int numArtCols(lRows);
    int *lowerRowInd = model_->getLowerRowInd();
    //int *lowerColInd = model_->getLowerColInd();
    //double *origRowLb = model_->getOrigRowLb();
    //double *origRowUb = model_->getOrigRowUb();
    //double *origColLb =model_->getOrigColLb();
    //double *origColUb = model_->getOrigColUb();

    double *upperBound = new double[numArtCols];

    OsiSolverInterface * nSolver;


    if (feasCheckSolver == "Cbc"){
	nSolver = new OsiCbcSolverInterface();
    }else if (feasCheckSolver == "SYMPHONY"){
	#ifdef COIN_HAS_SYMPHONY
	nSolver = new OsiSymSolverInterface();
	#else
	throw CoinError("SYMPHONY chosen as solver, but it has not been enabled",
			"findUBOnArtCols", "MibSZeroSum");
	#endif
    }else if (feasCheckSolver == "CPLEX"){
	#ifdef COIN_HAS_CPLEX
	nSolver = new OsiCpxSolverInterface();
	#else
	throw CoinError("CPLEX chosen as solver, but it has not been enabled",
			"findUBOnArtCols", "MibSZeroSum");
	#endif
    }else{
	throw CoinError("Unknown solver chosen",
			"findUBOnArtCols", "MibSZeroSum");
    }

    double objSense(1);
    //set integer variables
    for(i = 0; i < numCols; i++){
	if(oSolver->isInteger(i)){
	    nSolver->isInteger(i);
	}
    }

    double *objCoeffs = new double[numCols];

    const CoinPackedMatrix *rowMatrix = moel_->origConstCoefMatrix_->reverseOrdering();
    const double * matElements = rowMatrix.getElements();
    const int * matLength = rowMatrix.getVectorLengths();
    const int * matIndices = rowMatrix.getIndices();
    const int * matStarts = rowMatrix.getVectorStarts();

    //solve the original problem while upper level objective is
    //set to one of the lower-level rows
    for(i = 0; i < lRows; i++){
	//set the objective
	CoinZeroN(objCoeffs, numCols);
	index = lowerRowInd[i];
	start = matStarts[index];
	for(j = start; j < start + matLength[index]; j++){
	    objCoeffs[matIndices[j]] = matElements[j];
	    
	    
	    
    
    

	    }*/

//#############################################################################
double
MibSZeroSum::findBigM()
{
    std::string feasCheckSolver(model_->MibSPar_->
				entry(MibSParams::feasCheckSolver));
    int maxThreadsLL(model_->MibSPar_->entry
		     (MibSParams::maxThreadsLL));
    int whichCutsLL(model_->MibSPar_->entry
		    (MibSParams::whichCutsLL));

    int i;
    int intCnt(0);
    int numCols(model_->getNumOrigVars());
    OsiSolverInterface *oSolver = model_->solver();
    
    int * integerVars = new int[numCols];
    OsiSolverInterface *nSolver;

    CoinFillN(integerVars, numCols, 0);

    /** Fill in array of integer variables **/
    for(i = 0; i < numCols; i++){
	if(oSolver->isInteger(i)){
	    integerVars[intCnt] = i;
	    intCnt ++;
	}
    }

    if (feasCheckSolver == "Cbc"){
	nSolver = new OsiCbcSolverInterface();
    }else if (feasCheckSolver == "SYMPHONY"){
	        #ifdef COIN_HAS_SYMPHONY
	nSolver = new OsiSymSolverInterface();
	        #else
	throw CoinError("SYMPHONY chosen as solver, but it has not been enabled",
			"findBigM", "MibSZeroSum");
	        #endif
    }else if (feasCheckSolver == "CPLEX"){
	        #ifdef COIN_HAS_CPLEX
	nSolver = new OsiCpxSolverInterface();
	        #else
	throw CoinError("CPLEX chosen as solver, but it has not been enabled",
			"findBigM", "MibSZeroSum");
	        #endif
    }else{
	throw CoinError("Unknown solver chosen",
			"findBigM", "MibSZeroSum");
    }

    const CoinPackedMatrix *matrix = model_->getOrigConstCoefMatrix();

    nSolver->loadProblem(*matrix, model_->getOrigColLb(), model_->getOrigColUb(),
			 model_->getObjCoefficients(), model_->getOrigRowLb(),
			 model_->getOrigRowUb());

    nSolver->setObjSense(1);
    for(i = 0; i < intCnt; i++){
	nSolver->setInteger(integerVars[i]);
    }
    
    nSolver->setHintParam(OsiDoReducePrint, true, OsiHintDo);

    if (feasCheckSolver == "Cbc"){
	dynamic_cast<OsiCbcSolverInterface *>
	    (nSolver)->getModelPtr()->messageHandler()->setLogLevel(0);
    }else if (feasCheckSolver == "SYMPHONY"){
#if COIN_HAS_SYMPHONY
	sym_environment *env = dynamic_cast<OsiSymSolverInterface *>
	    (nSolver)->getSymphonyEnvironment();
	sym_set_int_param(env, "do_primal_heuristic", FALSE);
	sym_set_int_param(env, "verbosity", -2);
	sym_set_int_param(env, "prep_level", -1);
	sym_set_int_param(env, "max_active_nodes", maxThreadsLL);
	sym_set_int_param(env, "tighten_root_bounds", FALSE);
	sym_set_int_param(env, "max_sp_size", 100);
	sym_set_int_param(env, "do_reduced_cost_fixing", FALSE);
	if (whichCutsLL == 0){
	    sym_set_int_param(env, "generate_cgl_cuts", FALSE);
	}else{
	    sym_set_int_param(env, "generate_cgl_gomory_cuts", GENERATE_DEFAULT);
	}
	if (whichCutsLL == 1){
	    sym_set_int_param(env, "generate_cgl_knapsack_cuts",
			      DO_NOT_GENERATE);
	    sym_set_int_param(env, "generate_cgl_probing_cuts",
			      DO_NOT_GENERATE);
	    sym_set_int_param(env, "generate_cgl_clique_cuts",
			      DO_NOT_GENERATE);
	    sym_set_int_param(env, "generate_cgl_twomir_cuts",
			      DO_NOT_GENERATE);
	    sym_set_int_param(env, "generate_cgl_flowcover_cuts",
			      DO_NOT_GENERATE);
	}
#endif
    }else if (feasCheckSolver == "CPLEX"){
#ifdef COIN_HAS_CPLEX
	nSolver->setHintParam(OsiDoReducePrint);
	nSolver->messageHandler()->setLogLevel(0);
	CPXENVptr cpxEnv =
	    dynamic_cast<OsiCpxSolverInterface*>(lSolver)->getEnvironmentPtr();
	assert(cpxEnv);
	CPXsetintparam(cpxEnv, CPX_PARAM_SCRIND, CPX_OFF);
	CPXsetintparam(cpxEnv, CPX_PARAM_THREADS, maxThreadsLL);
#endif
    }

    nSolver->branchAndBound();

    if(nSolver->isProvenOptimal()){
	bigM = -2 * nSolver->getObjValue();
    }
    //this should not happen because we assume that the relaxation
    //is bounded and feasible
    else{
	assert(0);
    }

    return bigM;
}

//#############################################################################   
void
MibSZeroSum::findBoundsOnLLCols(double *newLb, double *newUb, double *LLSol)
{
    
    int i;
    double value(0.0);
    int colIndex(0);
    double etol(model_->etol_);
    int lCols(model_->getLoweDim());
    const double *origColLb(model_->getOrigColLb());
    const double *origColUb(model_->getOrigColUb());
    int *lColInd(model_->getLowerColInd());
    
    for(i = 0; i < lCols; i++){
	value = LLSol[i];
	colIndex = lColInd[i];
	if(fabs(floor(value + 0.5) - value) <= etol){
	    value = floor(value + 0.5);
	    if(1 - value + origColLb[colIndex] <= 0){
		newLb[i] = value - 1;
	    }
	    else{
		newLb[i] = origColLb[colIndex];
	    }
	    if(value + 1 - origColUb[colIndex] <= 0){
		newUb[i] = value + 1;
	    }
	    else{
		newUb[i] = origColUb[colIndex];
	    }
	}
	else{
	    newLb[i] = floor(value);
	    newUb[i] = ceil(value);
	}
    }
}

//#############################################################################
void
MibSZeroSum::initialSetUpLowerSolver(CoinPackedMatrix *A2Matrix)
{
    OsiSolverInterface * oSolver = model_->solver();
    int i;
    int colIndex(0), rowIndex(0), numElements(0), pos(0), cntInt(0);
    int uCols(model_->getUpperDim());
    int lCols(model_->getLowerDim());
    int lRows(model_->getLowerRowNum());
    double *origColLb(model_->getOrigColLb());
    double *origColUb(model_->getOrigColUb());
    double *origRowLb(model_->getOrigRowLb());
    double *origRowUb(model_->getOrigRowUb());
    int *uColInd(model_->getUpperColInd());
    int *lColInd(model_->getLowerColInd());
    int *lRowInd(model_->getLowerrowInd());

    CoinShallowPackedVector origVec;
    CoinPackedVector vecG2;
    CoinPackedVector vecA2;
    
    
    CoinPackedMatrix matrix = *model_->getOrigConstCoefMatrix();
    CoinPackedMatrix *matrixG2 = new CoinPackedMatrix(false, 0, 0);

    double *colLb = new double[lCols];
    double *colUb = new double[lCols];
    CoinZeroN(colLb, lCols);
    CoinZeroN(colUb, lCols);

    double *rowLb = new double[lRows];
    double *rowUb = new double[lRows];
    CoinZeroN(rowLb, lRows);
    CoinZeroN(rowUb, lRows);

    int *integerVars = new int[lCols];
    CoinZeroN(integerVars, lCols);
    
    matrix.reverseOrdering();
    matrixG2->setDimensions(0, lCols);
    matrixA2->setDimensions(0, uCols);

    for(i = 0; i < lRows; i++){
	rowIndex = lRowInd[i];
	origVec = matrix.getVector(rowIndex);
	numElements = origVec.getNumElements();
	const int *indices = origVec.getIndices();
	const double *elements = origVec.getElements();
	for(j = 0; j < numElements; j++){
	    pos = model_->binarySearch(0, lCols -1, indices[j], lColInd);
	    if(pos >= 0){
		vecG2.insert(pos, elements[j]);
	    }
	    else{
		pos = model_->binarySearch(0, uCols -1, indices[j], uColInd);
		vecA2.insert(pos, elements[j]);
	    }
	}
	matrixG2->append(vecG2);
	matrixA2->append(vecA2);
	origVec.clear();
	vecA2.clear();
	vecG2.clear();
    }

    for(i = 0; i < lCols; i++){
	colIndex = lColInd[i];
	colLb[i] = origColLb[colIndex];
	colUb[i] = origColUb[colIndex];
    }

    for(i = 0; i < lRows; i++){
	rowIndex = lRowInd[i];
	rowLb[i] = origRowLb[rowIndex];
	rowUb[i] = origRowUb[rowIndex];
    }

    cntInt = 0;
    for(i = 0; i < lCols; i++){
	colIndex = lColInd[i];
	if(oSolver->isInteger(colIndex)){
	    integerVars[cntInt] = i;
	    cntInt++;
	}
    }

    lSolver_->loadProblem(*matrixG2, colLb, colUb,
			  model_->getLowerObjCoeffs(), rowLb, rowUb);

    lSolver_->setInteger(integerVars, cntInt);
    lSolver_->setObjSense(model_->getLowerObjSense()); //1 min; -1 max
    lSolver_->setHintParam(OsiDoReducePrint, true, OsiHintDo);

    delete matrixG2;
    delete [] colLb;
    delete [] colUb;
    delete [] rowLb;
    delete [] rowUb;
    delete [] integerVars;
}

//#############################################################################
void
MibSZeroSum::modifyLowerSolver(double *sol, CoinPackedMatrix *matrixA2)
{
    int i;
    int colIndex(0);
    int uCols(model_->getUpperDim());
    int *uColInd(model_->getUpperColInd());

    double *mult = new double[lRows];
    CoinZeroN(mult, lRows);
    
    double *uSol = new double[uCols];
    CoinZeroN(uSol, uCols);

    for(i = 0; i < uCols; i++){
	colIndex = uColInd[i];
	uSol[i] = sol[colIndex];
    }

    matrixA2->times(uSol, mult);

    for(i = 0; i < lRows; i++){
	lower = lSolver_getRowLower()[i];
	upper =lSolver_getRowUpper()[i];
	value = mult[i];
	lSolver_->setRowBounds(i, lower-value, upper-value);

    }

    delete [] mult;
    delete [] uSol;
}

//#############################################################################
void
MibSZeroSum::solveLowerProblem()
{

    if(feasCheckSolver == "Cbc"){
	            dynamic_cast<OsiCbcSolverInterface *>
			(lSolver_)->getModelPtr()->messageHandler()->setLogLevel(0);
    }else if(feasCheckSolver == "SYMPHONY"){
#if COIN_HAS_SYMPHONY
	sym_environment *env = dynamic_cast<OsiSymSolverInterface *>
	    (lSolver_)->getSymphonyEnvironment();
	sym_set_int_param(env, "do_primal_heuristic", FALSE);
	sym_set_int_param(env, "verbosity", -2);
	sym_set_int_param(env, "prep_level", -1);
        sym_set_int_param(env, "max_active_nodes", maxThreadsLL);
	sym_set_int_param(env, "tighten_root_bounds", FALSE);
	sym_set_int_param(env, "max_sp_size", 100);
	sym_set_int_param(env, "do_reduced_cost_fixing", FALSE);
	if (whichCutsLL == 0){
	    sym_set_int_param(env, "generate_cgl_cuts", FALSE);
	}else{
	    sym_set_int_param(env, "generate_cgl_gomory_cuts", GENERATE_DEFAULT);
	}
	if(whichCutsLL == 1){
	    sym_set_int_param(env, "generate_cgl_knapsack_cuts",
			      DO_NOT_GENERATE);
	    sym_set_int_param(env, "generate_cgl_probing_cuts",
			      DO_NOT_GENERATE);
	    sym_set_int_param(env, "generate_cgl_clique_cuts",
			      DO_NOT_GENERATE);
	    sym_set_int_param(env, "generate_cgl_twomir_cuts",
			      DO_NOT_GENERATE);
	    sym_set_int_param(env, "generate_cgl_flowcover_cuts",
			      DO_NOT_GENERATE);
	}
#endif
    }else if(feasCheckSolver == "CPLEX"){
#ifdef COIN_HAS_CPLEX
	lSolver_->setHintParam(OsiDoReducePrint);
	lSolver_->messageHandler()->setLogLevel(0);
	CPXENVptr cpxEnv =
	    dynamic_cast<OsiCpxSolverInterface*>(lSolver_)->getEnvironmentPtr();
	assert(cpxEnv);
	CPXsetintparam(cpxEnv, CPX_PARAM_SCRIND, CPX_OFF);
	CPXsetintparam(cpxEnv, CPX_PARAM_THREADS, maxThreadsLL);
#endif
    }

}

//#############################################################################
void
MibSZeroSum::findLowerSol(double *sol, double *LLSol)
{
    for(i = 0; i < lCols; i++){
	colIndex = lColInd[i];
	LLSol[i] = sol[colIndex];
    }
}
				    
	
	

	
		
	
	    
    
	
	
