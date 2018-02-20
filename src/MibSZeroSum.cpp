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

#include "MibSZeroSum.hpp"
#include "MibSModel.hpp"
#include "MibSConfig.hpp"

#if  COIN_HAS_MPI
#include "AlpsKnowledgeBrokerMPI.h"
#else
#include "AlpsKnowledgeBrokerSerial.h"
#endif

#include "OsiCbcSolverInterface.hpp"

#include "OsiCbcSolverInterface.hpp"

#ifdef COIN_HAS_SYMPHONY
#include "symphony.h"
#include "SymConfig.h"
#include "OsiSymSolverInterface.hpp"
#endif
#ifdef COIN_HAS_CPLEX
#include "cplex.h"
#include "OsiCpxSolverInterface.hpp"
#endif


//#############################################################################
MibSZeroSum::MibSZeroSum()
{
    model_ = 0;
    isAlgStarted_ = false;
    isUnbounded_ = false;
    foundOptimal_ = false;
    returnedUpperBound_ = false;
    optimalSol_ = NULL;
}

//#############################################################################
MibSZeroSum::~MibSZeroSum()
{
    if(optimalSol_) delete [] optimalSol_;
}

//#############################################################################
void
MibSZeroSum::solveZeroSum(MibSModel *mibs, double *sol)
{
    
    model_ = mibs;
    isAlgStarted_ = true;

    int numCols(model_->getNumOrigVars());
    int lCols(model_->getLowerDim());
  
    // first-phase of the algorithm: finding a bound on lower-level
    //variables so that the lower-level problem will be feasible for
    //all feasible first-level variables
    doFirstPhase();

}

//#############################################################################
void
    MibSZeroSum::doFirstPhase()
{
    OsiSolverInterface * oSolver = model_->getSolver();

    std::string feasCheckSolver(model_->MibSPar_->
				entry(MibSParams::feasCheckSolver));

    double timeLimit(model_->AlpsPar()->entry(AlpsParams::timeLimit));

    bool isFirstPhaseFinished(false), hasPositiveArtCol(false);
    int i;
    int tmpIndex(0), colIndex(0);
    int numAddedArtCols(0), numAddedRows(0), numAddedLrows(0), numAdddedURows(0);
    double artColLower(0.0), artColUpper(0.0);
    double remainingTime(0.0);
    double etol(model_->etol_);
    double infinity(oSolver->getInfinity());
    int numCols(model_->getNumOrigVars());
    int numRows(model_->getNumOrigCons());
    int lCols(model_->getLowerDim());
    int lRows(model_->getLowerRowNum());
    int newUCols(model_->getUpperDim());
    int newLCols(lCols + 3 * lRows);
    int newURows(model_->getOrigUpperRowNum()); 
    int newLRows(5 * lRows);
    int newNumCols(newUCols + newLCols);
    int newNumRows(newURows + newLRows);
    double *rowLB(model_->getOrigRowLb());
    double *rowUB(model_->getOrigRowUb());
    const char *rowSense(model_->getOrigRowSense());
    int *lColInd(model_->getLowerColInd());
    int *lRowInd(model_->getLowerRowInd());
    double *lObjCoeffs(model_->getLowerObjCoeffs());
    CoinPackedMatrix matrix = *model_->getOrigConstCoefMatrix();
    
    //set a reasonable value for big M objective
    double bigMObj = findBigMObj();

    //set a reasonable value for big M constraints
    double bigMConst = findBigMConst(artColLower, artColUpper);

    //Set up lp solver
    OsiCpxSolverInterface lpSolver;
    lpSolver.messageHandler()->setLogLevel(0);
    CPXENVptr cpxEnv = lpSolver.getEnvironmentPtr();
    assert(cpxEnv);
    CPXsetintparam(cpxEnv, CPX_PARAM_SCRIND, CPX_OFF);
    CPXsetintparam(cpxEnv, CPX_PARAM_THREADS, 1);

    //set col type
    char * newColType = new char[newNumCols];
    memcpy(newColType, model_->colType_, numCols);
    CoinFillN(newColType + numCols, 2 * lRows, 'C');
    CoinFillN(newColType + numCols + 2 * lRows, lRows, 'B');

    //set matrix, row bounds, row senses and lower and upper row indices
    CoinShallowPackedVector origVec;
    CoinPackedVector modifiedVec;
    CoinPackedVector appendedVec;
    CoinPackedMatrix *newMatrix = new CoinPackedMatrix(false, 0, 0);
    double *newRowLB = new double[newNumRows];
    double *newRowUB = new double[newNumRows];
    char *newRowSense = new char[newNumRows];
    int *newLRowInd = new int[newLRows];
    int *newURowInd = NULL;
    if(newURows > 0){
	newURowInd = new int[newURows];
    }
    
    matrix.reverseOrdering();
    newMatrix->setDimensions(0, newNumCols);
    numAddedArtCols = 0;
    numAddedRows = 0;
    numAddedLrows = 0;
    numAdddedURows = 0;
    for(i = 0; i < numRows; i++){
	origVec = matrix.getVector(i);
	if(model_->binarySearch(0, lRows - 1, i, lRowInd) < 0){
	    newMatrix->appendRow(origVec);
	    newRowLB[numAddedRows] = rowLB[i];
	    newRowUB[numAddedRows] = rowUB[i];
	    newRowSense[numAddedRows] = rowSense[i];
	    newURowInd[numAdddedURows] = numAddedRows;
	    numAddedRows++;
	    numAdddedURows++;
	}
	else{
	    if(rowSense[i] == 'L'){
		appendedVec.insert(numCols + numAddedArtCols, -1);
	    }
	    else if(rowSense[i] == 'G'){
		appendedVec.insert(numCols + numAddedArtCols, 1);
	    }
	    else{
		assert(0);
	    }
	    modifiedVec = origVec;
	    modifiedVec.append(appendedVec);
	    newMatrix->appendRow(modifiedVec);
	    newRowLB[numAddedRows] = rowLB[i];
	    newRowUB[numAddedRows] = rowUB[i];
	    newRowSense[numAddedRows] = rowSense[i];
	    newLRowInd[numAddedLrows] = numAddedRows;
	    numAddedRows++;
	    numAddedLrows++;
	    newMatrix->appendRow(modifiedVec);
	    if(rowSense[i] == 'L'){
		newRowLB[numAddedRows] = rowUB[i];
		newRowUB[numAddedRows] = infinity;
		newRowSense[numAddedRows] = 'G';
	    }
	    else{
		newRowUB[numAddedRows] = rowLB[i];
		newRowLB[numAddedRows] = -infinity;
		newRowSense[numAddedRows] = 'L';
	    }
	    newLRowInd[numAddedLrows] = numAddedRows;
	    numAddedRows++;
	    numAddedLrows++;
	    appendedVec.clear();
	    modifiedVec.clear();
	    tmpIndex = numCols + numAddedArtCols;
	    modifiedVec.insert(tmpIndex, -1);
	    modifiedVec.insert(tmpIndex + lRows, 1);
	    modifiedVec.insert(tmpIndex + 2 * lRows, -bigMConst);
	    newMatrix->appendRow(modifiedVec);
	    newRowUB[numAddedRows] = 0.0;
	    newRowLB[numAddedRows] = -infinity;
	    newRowSense[numAddedRows] = 'L';
	    newLRowInd[numAddedLrows] = numAddedRows;
	    numAddedRows++;
	    numAddedLrows++;
	    modifiedVec.clear();
	    modifiedVec.insert(tmpIndex, 1);
	    modifiedVec.insert(tmpIndex + lRows, 1);
	    modifiedVec.insert(tmpIndex + 2 * lRows, bigMConst);
	    newMatrix->appendRow(modifiedVec);
	    newRowUB[numAddedRows] = bigMConst;
	    newRowLB[numAddedRows] = -infinity;
	    newRowSense[numAddedRows] = 'L';
	    newLRowInd[numAddedLrows] = numAddedRows;
	    numAddedRows++;
	    numAddedLrows++;
	    modifiedVec.clear();
	    modifiedVec.insert(tmpIndex + lRows, 1);
	    modifiedVec.insert(tmpIndex + 2 * lRows, bigMConst);
	    newMatrix->appendRow(modifiedVec);
	    newRowUB[numAddedRows] = bigMConst;
	    newRowLB[numAddedRows] = -infinity;
	    newRowSense[numAddedRows] = 'L';
	    newLRowInd[numAddedLrows] = numAddedRows;
	    numAddedRows++;
	    numAddedLrows++;
	    modifiedVec.clear();
	    numAddedArtCols++;
	}
    }

    //set col bounds
    double *newColLB = new double[newNumCols];
    double *newColUB = new double[newNumCols];
    memcpy(newColLB, model_->getOrigColLb(), sizeof(double) * numCols);
    memcpy(newColUB, model_->getOrigColUb(), sizeof(double) * numCols);
    //sahar:modify it
    CoinFillN(newColLB + numCols, lRows, artColLower);
    CoinFillN(newColUB + numCols, lRows, artColUpper);
    CoinFillN(newColLB + numCols + lRows, 2 * lRows, 0.0);
    CoinFillN(newColUB + numCols + lRows, lRows, artColUpper);
    CoinFillN(newColUB + numCols + 2 * lRows, lRows, 1.0);

    //set the upper- and lower-level objectives
    double *newUObjCoeffs = new double[newNumCols];
    double *newLObjCoeffs = new double[newLCols];
    CoinZeroN(newUObjCoeffs, newNumCols);
    CoinZeroN(newLObjCoeffs, newLCols);
    memcpy(newUObjCoeffs, oSolver->getObjCoefficients(), sizeof(double) * numCols);
    CoinFillN(newUObjCoeffs + numCols + lRows, lRows, -1 * bigMObj);
    memcpy(newLObjCoeffs, lObjCoeffs, sizeof(double) * lCols);
    CoinFillN(newLObjCoeffs + lCols + lRows, lRows, bigMObj);

    //set upper and lower column indices
    int *newUColInd = new int[newUCols];
    int *newLColInd = new int[newLCols];
    memcpy(newUColInd, model_->getUpperColInd(), sizeof(int) * newUCols);
    memcpy(newLColInd, lColInd, sizeof(int) * lCols);
    CoinIotaN(newLColInd + lCols, 3 * lRows, numCols);

    //set structure row indices
    int *newStructRowInd = new int[newNumRows];
    CoinIotaN(newStructRowInd, newNumRows, 0);

    double *firstPhaseSol = new double[newNumCols];
    CoinZeroN(firstPhaseSol, newNumCols);

    double objVal(0.0), realLObj(0.0), lObj(0.0);
    int argc = 1;
    char **argv = new char *[1];
    argv[0] = (char *) "mibs";
    
    while(!isFirstPhaseFinished){
	MibSModel *firstPhaseModel = new MibSModel();
	firstPhaseModel->setSolver(&lpSolver);

	CoinPackedMatrix tmpMatrix;
	tmpMatrix = *newMatrix;
	tmpMatrix.reverseOrdering();

	remainingTime = timeLimit - model_->broker_->subTreeTimer().getTime();

	firstPhaseModel->AlpsPar()->setEntry(AlpsParams::timeLimit, remainingTime);
	firstPhaseModel->BlisPar()->setEntry(BlisParams::heurStrategy, PARAM_OFF);
	//firstPhaseModel->BlisPar()->setEntry(BlisParams::heurRound, PARAM_OFF);
	//firstPhaseModel->MibSPar()->setEntry(MibSParams::checkInstanceStructure, false);
	firstPhaseModel->MibSPar()->setEntry(MibSParams::usePreprocessor, PARAM_OFF);
	firstPhaseModel->MibSPar()->setEntry(MibSParams::useLowerObjHeuristic, PARAM_OFF);
	firstPhaseModel->MibSPar()->setEntry(MibSParams::useObjCutHeuristic, PARAM_OFF);
	firstPhaseModel->MibSPar()->setEntry(MibSParams::useWSHeuristic, PARAM_OFF);
	firstPhaseModel->MibSPar()->setEntry(MibSParams::useGreedyHeuristic, PARAM_OFF);
	firstPhaseModel->MibSPar()->setEntry(MibSParams::bilevelCutTypes, 0);
	firstPhaseModel->MibSPar()->setEntry(MibSParams::usePureIntegerCut, PARAM_OFF);
	firstPhaseModel->MibSPar()->setEntry(MibSParams::useNoGoodCut, PARAM_OFF);
	firstPhaseModel->MibSPar()->setEntry(MibSParams::useIncObjCut, PARAM_OFF);
	firstPhaseModel->MibSPar()->setEntry(MibSParams::cutStrategy, 2);
	firstPhaseModel->MibSPar()->setEntry(MibSParams::objBoundStrategy, 0);
	firstPhaseModel->MibSPar()->setEntry(MibSParams::blisCutStrategy, -2);
	firstPhaseModel->MibSPar()->setEntry(MibSParams::branchStrategy,
					     MibSBranchingStrategyLinking);
	firstPhaseModel->MibSPar()->setEntry(MibSParams::useBendersCut, PARAM_OFF);
	firstPhaseModel->MibSPar()->setEntry(MibSParams::useGeneralNoGoodCut, PARAM_OFF);
	firstPhaseModel->MibSPar()->setEntry(MibSParams::useIntersectionCut, PARAM_ON);
	firstPhaseModel->MibSPar()->setEntry(MibSParams::intersectionCutType,
					     MibSIntersectionCutTypeHypercubeIC);
	firstPhaseModel->MibSPar()->setEntry(MibSParams::feasCheckSolver,
					     feasCheckSolver.c_str());
	firstPhaseModel->MibSPar()->setEntry(MibSParams::solveSecondLevelWhenXVarsInt,
					     PARAM_ON);
	firstPhaseModel->MibSPar()->setEntry(MibSParams::solveSecondLevelWhenLVarsFixed,
					     PARAM_ON);
	firstPhaseModel->MibSPar()->setEntry(MibSParams::computeBestUBWhenLVarsFixed,
					     PARAM_ON);
	firstPhaseModel->MibSPar()->setEntry(MibSParams::useLinkingSolutionPool, PARAM_ON);
	firstPhaseModel->MibSPar()->setEntry(MibSParams::useZeroSumAlg, PARAM_OFF);

	firstPhaseModel->loadAuxiliaryData(newLCols, newLRows,
					   newLColInd, newLRowInd, 1,
					   newLObjCoeffs, newUCols, newURows,
					   newUColInd, newURowInd, newNumRows,
					   newStructRowInd, 0, NULL);
	firstPhaseModel->loadProblemData(tmpMatrix, newColLB, newColUB, newUObjCoeffs,
					 newRowLB, newRowUB, newColType, 1, infinity,
					 newRowSense, true);

#ifdef  COIN_HAS_MPI
	AlpsKnowledgeBrokerMPI broker(argc, argv, *firstPhaseModel);
#else
	AlpsKnowledgeBrokerSerial broker(argc, argv, *firstPhaseModel);
#endif
	
	broker.search(firstPhaseModel);
	
	if(firstPhaseModel->getNumSolutions() > 0){
	    BlisSolution * optSol =dynamic_cast<BlisSolution*>
		(broker.getBestKnowledge(AlpsKnowledgeTypeSolution).first);
	    memcpy(firstPhaseSol, optSol->getValues(), sizeof(double) * newNumCols);
	    objVal = broker.getBestQuality();
	}
	else{
	    assert(0);
	}

	broker.printBestSolution();

	/*hasPositiveArtCol = false;
	for(i = 0; i < lRows; i++){
	    if(firstPhaseSol[i + numCols] > etol){
		hasPositiveArtCol = true;
		break;
	    }
	}*/

	isFirstPhaseFinished = true;
	
	if(timeLimit > model_->broker_->subTreeTimer().getTime()){
	    if(objVal < -10000000000){
		isUnbounded_ = true;
		if(!optimalSol_){
		    optimalSol_ = new double[numCols];
		}
		memcpy(optimalSol_, firstPhaseSol, sizeof(double) * numCols);
	    }
	    else{
		foundOptimal_ = true;
		if(!optimalSol_){
		    optimalSol_ = new double[numCols];
		}
		memcpy(optimalSol_, firstPhaseSol, sizeof(double) * numCols);
	    }
	}
	else{
	    returnedUpperBound_ = true;
	}
        
	delete firstPhaseModel;
    }

    //delete matrixA2;
    delete [] newColType;
    delete newMatrix;
    delete [] newRowLB;
    delete [] newRowUB;
    delete [] newRowSense;
    delete [] newLRowInd;
    if(newURowInd){
	delete [] newURowInd;
    }
    delete [] newColLB;
    delete [] newColUB;
    delete [] newLObjCoeffs;
    delete [] newUObjCoeffs;
    delete [] newUColInd;
    delete [] newLColInd;
    delete [] newStructRowInd;
    delete [] firstPhaseSol;
    
}

//#############################################################################
double
MibSZeroSum::findBigMObj()
{
    std::string feasCheckSolver(model_->MibSPar_->
				entry(MibSParams::feasCheckSolver));
    int maxThreadsLL(model_->MibSPar_->entry
		     (MibSParams::maxThreadsLL));
    int whichCutsLL(model_->MibSPar_->entry
		    (MibSParams::whichCutsLL));

    int i;
    int intCnt(0);
    double bigM(0.0), value(0.0);
    double etol(model_->etol_);
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
			 oSolver->getObjCoefficients(), model_->getOrigRowLb(),
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
	    dynamic_cast<OsiCpxSolverInterface*>(nSolver)->getEnvironmentPtr();
	assert(cpxEnv);
	CPXsetintparam(cpxEnv, CPX_PARAM_SCRIND, CPX_OFF);
	CPXsetintparam(cpxEnv, CPX_PARAM_THREADS, maxThreadsLL);
#endif
    }

    nSolver->branchAndBound();

    if(nSolver->isProvenOptimal()){
	value = nSolver->getObjValue();
	if(value < etol){
	    bigM = -1 * 100000 * value + 5;
	}
	else{
	    bigM = 5;
	}
    }
    //this should not happen because we assume that the relaxation
    //is bounded and feasible
    else{
	assert(0);
    }

    delete nSolver;
    delete [] integerVars;

    return bigM;
}

//#############################################################################
double
MibSZeroSum::findBigMConst(double &artColLower, double &artColUpper)
{

    int i , j;
    int numElements(0), rowIndex(0), colIndex(0);
    double mult(0.0), value(0.0), rhs(0.0);
    double tmpMin(0.0), tmpMax(0.0), minValue(0.0), maxValue(0.0);
    double bigM(0.0);
    double etol(model_->etol_);
    double infinity(model_->solver()->getInfinity());
    int lRows(model_->getLowerRowNum());
    int *lRowInd(model_->getLowerRowInd());
    double *origColLb(model_->getOrigColLb());
    double *origColUb(model_->getOrigColUb());
    double *origRowLb(model_->getOrigRowLb());
    double *origRowUb(model_->getOrigRowUb());
    const char *rowSense(model_->getOrigRowSense());
    
    CoinPackedMatrix matrix = *model_->getOrigConstCoefMatrix();
    matrix.reverseOrdering();

    CoinShallowPackedVector origVec;

    minValue = infinity;
    maxValue = -infinity;

    for(i = 0; i < lRows; i++){
	rowIndex = lRowInd[i];
	if(rowSense[rowIndex] == 'L'){
	    mult = 1;
	    rhs = origRowUb[rowIndex];
	}
	else{
	    mult = -1;
	    rhs = -1 * origRowLb[rowIndex];
	}
	origVec = matrix.getVector(rowIndex);
	numElements = origVec.getNumElements();
	const int *indices = origVec.getIndices();
	const double *elements = origVec.getElements();
	tmpMin = -1 * rhs;
	tmpMax = -1 * rhs;
	for(j = 0; j < numElements; j++){
	    value = elements[j] * mult;
	    colIndex = indices[j];
	    if(value > 0){
		tmpMin += value * origColLb[colIndex];
		tmpMax += value * origColUb[colIndex];
	    }
	    else{
		tmpMin += value* origColUb[colIndex];
		tmpMax += value* origColLb[colIndex];
	    }
	}
	if(minValue - tmpMin > etol){
	    minValue = tmpMin;
	}
	if(maxValue - tmpMax < -etol){
	    maxValue = tmpMax;
	}
    }

    artColLower = minValue;
    artColUpper = maxValue;

    value = CoinMax(fabs(minValue), fabs(maxValue));

    bigM = 2 * value + 1;

    return bigM;
}

//#############################################################################   
