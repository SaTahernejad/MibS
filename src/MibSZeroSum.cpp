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
    lSolver_ = 0;
    isAlgStarted_ = false;
    isUnbounded_ = false;
    foundOptimal_ = false;
    returnedNothing_ = false;
    returnedLowerBound_ = false;
    optimalSol_ = NULL;
    matrixA2_ = NULL;
}

//#############################################################################
MibSZeroSum::~MibSZeroSum()
{
    if(lSolver_) delete lSolver_;
    if(optimalSol_) delete [] optimalSol_;
    if(matrixA2_) delete matrixA2_;
}

//#############################################################################
void
MibSZeroSum::solveZeroSum(MibSModel *mibs, double *sol)
{
    
    model_ = mibs;
    isAlgStarted_ = true;

    int numCols(model_->getNumOrigVars());
    int lCols(model_->getLowerDim());
    double initialBoundSecondPh(0.0);

    double *lColLb = new double[lCols];
    double *lColUb = new double[lCols];
    CoinZeroN(lColLb, lCols);
    CoinZeroN(lColUb, lCols);

    double *LLSol = new double[lCols];
    memcpy(LLSol, sol, sizeof(double) * lCols);
    
    // first-phase of the algorithm: finding a bound on lower-level
    //variables so that the lower-level problem will be feasible for
    //all feasible first-level variables
    doFirstPhase(LLSol, lColLb, lColUb, initialBoundSecondPh);

    if((!foundOptimal_) && (!isUnbounded_) && (!returnedNothing_)){
	//start second phase
	doSecondPhase(lColLb, lColUb, initialBoundSecondPh);
    }

    delete [] lColLb;
    delete [] lColUb;
    delete [] LLSol;

}

//#############################################################################
void
    MibSZeroSum::doFirstPhase(double *LLSol, double *lColLb, double *lColUb,
			      double &initialBoundSecondPh)
{

    OsiSolverInterface * oSolver = model_->getSolver();

    std::string feasCheckSolver(model_->MibSPar_->
				entry(MibSParams::feasCheckSolver));

    double timeLimit(model_->AlpsPar()->entry(AlpsParams::timeLimit));

    bool isFirstPhaseFinished(false);
    int i;
    int colIndex(0);
    int numAddedArtCols(0);
    double remainingTime(0.0);
    double infinity(oSolver->getInfinity());
    double etol(model_->etol_);
    int numRows(model_->getNumOrigCons());
    int numCols(model_->getNumOrigVars());
    int lRows(model_->getLowerRowNum());
    int lCols(model_->getLowerDim());
    int newNumCols(numCols + lRows);
    int newNumRows(numRows);
    int newUCols(model_->getUpperDim());
    int newLCols(lCols + lRows);
    int newURows(model_->getOrigUpperRowNum());
    int newLRows(lRows);
    int *lColInd(model_->getLowerColInd());
    int *lRowInd(model_->getLowerRowInd());
    double *lObjCoeffs(model_->getLowerObjCoeffs());
    const char *rowSense(model_->getOrigRowSense());
    CoinPackedMatrix matrix = *model_->getOrigConstCoefMatrix();
    matrixA2_ = new CoinPackedMatrix(false, 0, 0);


    findBoundsOnLLCols(lColLb, lColUb, LLSol, true);

    double bigMObj = findBigMObj();

    double artColBound = findBoundArtCol();

    //Set up lp solver
    /*OsiClpSolverInterface lpSolver;
    lpSolver.getModelPtr()->setDualBound(1.0e10);
    lpSolver.messageHandler()->setLogLevel(0);*/

    OsiCpxSolverInterface lpSolver;
    lpSolver.messageHandler()->setLogLevel(0);
    CPXENVptr cpxEnv = lpSolver.getEnvironmentPtr();
    assert(cpxEnv);
    CPXsetintparam(cpxEnv, CPX_PARAM_SCRIND, CPX_OFF);
    CPXsetintparam(cpxEnv, CPX_PARAM_THREADS, 1);

    //set col type
    char * newColType = new char[newNumCols];
    memcpy(newColType, model_->colType_, numCols);
    CoinFillN(newColType + numCols, lRows, 'C');

    //set matrix
    CoinShallowPackedVector origVec;
    CoinPackedVector modifiedVec;
    CoinPackedVector appendedVec;
    CoinPackedMatrix *newMatrix = new CoinPackedMatrix(false, 0, 0);
    matrix.reverseOrdering();
    newMatrix->setDimensions(0, newNumCols);
    numAddedArtCols = 0;
    for(i = 0; i < numRows; i++){
	origVec = matrix.getVector(i);
	if(model_->binarySearch(0, lRows - 1, i, lRowInd) < 0){
	    newMatrix->appendRow(origVec);
	}
	else{
	    if(rowSense[i] == 'L'){
		appendedVec.insert(numCols + numAddedArtCols, -1);
	    }
	    else{
		appendedVec.insert(numCols + numAddedArtCols, 1);
	    }
	    modifiedVec = origVec;
	    modifiedVec.append(appendedVec);
	    newMatrix->appendRow(modifiedVec);
	    modifiedVec.clear();
	    appendedVec.clear();
	    numAddedArtCols ++;
	}
    }

    double *newRowLB = new double[newNumRows];
    double *newRowUB = new double[newNumRows];
    memcpy(newRowLB, model_->getOrigRowLb(), sizeof(double) * newNumRows);
    memcpy(newRowUB, model_->getOrigRowUb(), sizeof(double) * newNumRows);
    char *newRowSense = new char[newNumRows];
    memcpy(newRowSense, model_->getOrigRowSense(), newNumRows);

    int *newLRowInd = new int[newLRows];
    int *newURowInd = NULL;
    if(newURows > 0){
	newURowInd = new int[newURows];
	memcpy(newURowInd, model_->getUpperRowInd(), sizeof(int) * newURows);
    }
    memcpy(newLRowInd, model_->getLowerRowInd(), sizeof(int) * newLRows);
    

    //set col bounds
    double *newColLB = new double[newNumCols];
    double *newColUB = new double[newNumCols];
    memcpy(newColLB, model_->getOrigColLb(), sizeof(double) * numCols);
    memcpy(newColUB, model_->getOrigColUb(), sizeof(double) * numCols);
    CoinFillN(newColLB + numCols, lRows, 0.0);
    CoinFillN(newColUB + numCols, lRows, artColBound);

    //set the upper- and lower-level objectives
    double *newUObjCoeffs = new double[newNumCols];
    double *newLObjCoeffs = new double[newLCols];
    CoinZeroN(newUObjCoeffs, newNumCols);
    CoinZeroN(newLObjCoeffs, newLCols);
    memcpy(newUObjCoeffs, oSolver->getObjCoefficients(), sizeof(double) * numCols);
    CoinFillN(newUObjCoeffs + numCols, lRows, -1 * bigMObj);
    memcpy(newLObjCoeffs, lObjCoeffs, sizeof(double) * lCols);
    CoinFillN(newLObjCoeffs + lCols, lRows, 0.0);

    //set upper and lower column indices
    int *newUColInd = new int[newUCols];
    int *newLColInd = new int[newLCols];
    memcpy(newUColInd, model_->getUpperColInd(), sizeof(int) * newUCols);
    memcpy(newLColInd, lColInd, sizeof(int) * lCols);
    CoinIotaN(newLColInd + lCols, lRows, numCols);

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
	for(i = 0; i < lCols; i++){
	    colIndex = lColInd[i];
	    newColLB[colIndex] = lColLb[i];
	    newColUB[colIndex] = lColUb[i];
	}

	MibSModel *firstPhaseModel = new MibSModel();
	firstPhaseModel->setSolver(&lpSolver);

	CoinPackedMatrix tmpMatrix;
	tmpMatrix = *newMatrix;
	tmpMatrix.reverseOrdering();

	remainingTime = timeLimit - model_->broker_->subTreeTimer().getTime();

	firstPhaseModel->AlpsPar()->setEntry(AlpsParams::timeLimit, remainingTime);
	//firstPhaseModel->AlpsPar()->setEntry(AlpsParams::msgLevel, 1000);
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

	if(timeLimit > model_->broker_->subTreeTimer().getTime()){

	    if(!lSolver_){
		initialSetUpLowerSolver(matrixA2_);
	    }

	    //correct the lower-level solver
	    modifyLowerSolver(firstPhaseSol, matrixA2_);

	    //solve the lower-level problem
	    solveLowerProblem();

	    if(objVal < -10000000000){
		if(!lSolver_->isProvenOptimal()){
		    isUnbounded_ = true;
		    isFirstPhaseFinished = true;
		    if(!optimalSol_){
			optimalSol_ = new double[numCols];
		    }
		    memcpy(optimalSol_, firstPhaseSol, sizeof(double) * numCols);
		}
		else{
		    memcpy(LLSol, lSolver_->getColSolution(), sizeof(double) * lCols);
		    findBoundsOnLLCols(lColLb, lColUb, LLSol, false);
		}
	    }
	    else{
		isFirstPhaseFinished = true;
		if(lSolver_->isProvenOptimal()){
		    realLObj = lSolver_->getObjValue();
		    lObj = 0.0;
		    for(i = 0; i < lCols; i++){
			colIndex = lColInd[i];
			lObj += lObjCoeffs[i] * firstPhaseSol[lColInd[i]];
		    }
		    if(fabs(realLObj - lObj) <= etol){
			foundOptimal_ = true;
			if(!optimalSol_){
			    optimalSol_ = new double[numCols];
			}
			memcpy(optimalSol_, firstPhaseSol, sizeof(double) * numCols);
		    }
		    else{
			memcpy(LLSol, lSolver_->getColSolution(), sizeof(double) * lCols);
			findBoundsOnLLCols(lColLb, lColUb, LLSol, false);
			initialBoundSecondPh = -1 * lObj;
		    }
		}
		else{
		    assert(0);
		}
	    }
	}
	else{
	    returnedNothing_ = true;
	    isFirstPhaseFinished = true;
	}
	delete firstPhaseModel;
    }

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
void
MibSZeroSum::doSecondPhase(double *lColLb, double *lColUb, double objBound)
{
    OsiSolverInterface * oSolver = model_->getSolver();

    std::string feasCheckSolver(model_->MibSPar_->
				entry(MibSParams::feasCheckSolver));

    double timeLimit(model_->AlpsPar()->entry(AlpsParams::timeLimit));
    
    int i;
    int colIndex(0);
    double value(0.0);
    bool isSecondPhaseFinished(false);
    double remainingTime(0.0);
    double infinity(oSolver->getInfinity());
    double etol(model_->etol_);
    int newNumCols(model_->getNumOrigVars());
    int numRows(model_->getNumOrigCons());
    int newNumRows(numRows + 1);
    int newUCols(model_->getUpperDim());
    int newLCols(model_->getLowerDim());
    int uRows(model_->getOrigUpperRowNum());
    int newURows(uRows + 1);
    int newLRows(model_->getLowerRowNum());
    int *newUColInd(model_->getUpperColInd());
    int *newLColInd(model_->getLowerColInd());
    int *newLRowInd(model_->getLowerRowInd());
    const double *newUObjCoeffs(oSolver->getObjCoefficients());
    double *newLObjCoeffs(model_->getLowerObjCoeffs());
    char *newColType(model_->colType_);

    int *newURowInd = new int[newURows];
    memcpy(newURowInd, model_->getOrigUpperRowInd(), sizeof(double) * uRows);
    newURowInd[uRows] = numRows;

    int *newStructRowInd = new int[newNumRows];
    CoinIotaN(newStructRowInd, newNumRows, 0);

    CoinPackedMatrix newMatrix = *model_->getOrigConstCoefMatrix();
    newMatrix.reverseOrdering();

    CoinPackedVector appendedVec;

    for(i = 0; i < newLCols; i++){
	colIndex = newLColInd[i];
	value = newUObjCoeffs[colIndex];
	if(fabs(value) > etol){
	    appendedVec.insert(colIndex, value);
	}
    }

    newMatrix.appendRow(appendedVec);
    newMatrix.reverseOrdering();

    appendedVec.clear();

    double *newColLB = new double[newNumCols];
    double *newColUB = new double[newNumCols];
    memcpy(newColLB, model_->getOrigColLb(), sizeof(double) * newNumCols);
    memcpy(newColUB, model_->getOrigColUb(), sizeof(double) * newNumCols);

    double *newRowLB = new double[newNumRows];
    double *newRowUB = new double[newNumRows];
    memcpy(newRowLB, model_->getOrigRowLb(), sizeof(double) * numRows);
    memcpy(newRowUB, model_->getOrigRowUb(), sizeof(double) * numRows);
    newRowUB[numRows] = infinity;

    char *newRowSense = new char[newNumRows];
    memcpy(newRowSense, model_->getOrigRowSense(), numRows);
    newRowSense[numRows] = 'G';

    double *secondPhaseSol = new double[newNumCols];
    CoinZeroN(secondPhaseSol, newNumCols);

    double *LLSol = new double[newLCols];
    CoinZeroN(LLSol, newLCols);

    //Set up lp solver
    /*OsiClpSolverInterface lpSolver;
    lpSolver.getModelPtr()->setDualBound(1.0e10);
    lpSolver.messageHandler()->setLogLevel(0);*/

    OsiCpxSolverInterface lpSolver;
    lpSolver.messageHandler()->setLogLevel(0);
    CPXENVptr cpxEnv = lpSolver.getEnvironmentPtr();
    assert(cpxEnv);
    CPXsetintparam(cpxEnv, CPX_PARAM_SCRIND, CPX_OFF);
    CPXsetintparam(cpxEnv, CPX_PARAM_THREADS, 1);

    double objVal(0.0), realLObj(0.0), lObj(0.0);
    int argc = 1;
    char **argv = new char *[1];
    argv[0] = (char *) "mibs";

    while(!isSecondPhaseFinished){
	for(i = 0; i < newLCols; i++){
	    colIndex = newLColInd[i];
	    newColLB[colIndex] = lColLb[i];
	    newColUB[colIndex] = lColUb[i];
	}

	newRowLB[numRows] = objBound;

	MibSModel *secondPhaseModel = new MibSModel();
	secondPhaseModel->setSolver(&lpSolver);

	remainingTime = timeLimit - model_->broker_->subTreeTimer().getTime();
	
	secondPhaseModel->AlpsPar()->setEntry(AlpsParams::timeLimit, remainingTime);
	secondPhaseModel->BlisPar()->setEntry(BlisParams::heurStrategy, PARAM_OFF);
	//firstPhaseModel->BlisPar()->setEntry(BlisParams::heurRound, PARAM_OFF);
	//firstPhaseModel->MibSPar()->setEntry(MibSParams::checkInstanceStructure, false);
	secondPhaseModel->MibSPar()->setEntry(MibSParams::usePreprocessor, PARAM_OFF);
	secondPhaseModel->MibSPar()->setEntry(MibSParams::useLowerObjHeuristic, PARAM_OFF);
	secondPhaseModel->MibSPar()->setEntry(MibSParams::useObjCutHeuristic, PARAM_OFF);
	secondPhaseModel->MibSPar()->setEntry(MibSParams::useWSHeuristic, PARAM_OFF);
	secondPhaseModel->MibSPar()->setEntry(MibSParams::useGreedyHeuristic, PARAM_OFF);
	secondPhaseModel->MibSPar()->setEntry(MibSParams::bilevelCutTypes, 0);
	secondPhaseModel->MibSPar()->setEntry(MibSParams::usePureIntegerCut, PARAM_OFF);
	secondPhaseModel->MibSPar()->setEntry(MibSParams::useNoGoodCut, PARAM_OFF);
	secondPhaseModel->MibSPar()->setEntry(MibSParams::useIncObjCut, PARAM_OFF);
	secondPhaseModel->MibSPar()->setEntry(MibSParams::cutStrategy, 2);
	secondPhaseModel->MibSPar()->setEntry(MibSParams::objBoundStrategy, 0);
	secondPhaseModel->MibSPar()->setEntry(MibSParams::blisCutStrategy, -2);
	secondPhaseModel->MibSPar()->setEntry(MibSParams::branchStrategy,
					     MibSBranchingStrategyLinking);
	secondPhaseModel->MibSPar()->setEntry(MibSParams::useBendersCut, PARAM_OFF);
	secondPhaseModel->MibSPar()->setEntry(MibSParams::useGeneralNoGoodCut, PARAM_OFF);
	secondPhaseModel->MibSPar()->setEntry(MibSParams::useIntersectionCut, PARAM_ON);
	secondPhaseModel->MibSPar()->setEntry(MibSParams::intersectionCutType,
					     MibSIntersectionCutTypeHypercubeIC);
	secondPhaseModel->MibSPar()->setEntry(MibSParams::feasCheckSolver,
					     feasCheckSolver.c_str());
	secondPhaseModel->MibSPar()->setEntry(MibSParams::solveSecondLevelWhenXYVarsInt,
					     PARAM_ON);
	secondPhaseModel->MibSPar()->setEntry(MibSParams::solveSecondLevelWhenLVarsFixed,
					     PARAM_ON);
	secondPhaseModel->MibSPar()->setEntry(MibSParams::computeBestUBWhenLVarsFixed,
					     PARAM_ON);
	secondPhaseModel->MibSPar()->setEntry(MibSParams::useLinkingSolutionPool, PARAM_ON);
       
	secondPhaseModel->loadAuxiliaryData(newLCols, newLRows, newLColInd,
					    newLRowInd, 1, newLObjCoeffs,
					    newUCols, newURows, newUColInd,
					    newURowInd, newNumRows, newStructRowInd,
					    0, NULL);
	
	secondPhaseModel->loadProblemData(newMatrix, newColLB, newColUB, newUObjCoeffs,
					 newRowLB, newRowUB, newColType, 1, infinity,
					 newRowSense, false);

#ifdef  COIN_HAS_MPI
	AlpsKnowledgeBrokerMPI broker(argc, argv, *secondPhaseModel);
#else
	AlpsKnowledgeBrokerSerial broker(argc, argv, *secondPhaseModel);
#endif

	broker.search(secondPhaseModel);

	if(secondPhaseModel->getNumSolutions() > 0){
	    BlisSolution * optSol =dynamic_cast<BlisSolution*>
		(broker.getBestKnowledge(AlpsKnowledgeTypeSolution).first);
	    memcpy(secondPhaseSol, optSol->getValues(), sizeof(double) * newNumCols);
	    objVal = broker.getBestQuality();
	}
	else{
	    assert(0);
	}

	broker.printBestSolution();

	if(timeLimit > model_->broker_->subTreeTimer().getTime()){

	    //correct the lower-level solver
	    modifyLowerSolver(secondPhaseSol, matrixA2_);

	    //solve the lower-level problem
	    solveLowerProblem();

	    if(lSolver_->isProvenOptimal()){
		realLObj = lSolver_->getObjValue();
	        lObj = 0.0;
	        for(i = 0; i < newLCols; i++){
		    lObj += newLObjCoeffs[i] * secondPhaseSol[newLColInd[i]];
		}
	    
	        if(fabs(realLObj - lObj) <= etol){
		    isSecondPhaseFinished = true;
		    foundOptimal_ = true;
		    if(!optimalSol_){
			optimalSol_ = new double[newNumCols];
		    }
		    memcpy(optimalSol_, secondPhaseSol, sizeof(double) * newNumCols);
		}
	        else{
		    memcpy(LLSol, lSolver_->getColSolution(), sizeof(double) * newLCols);
		    findBoundsOnLLCols(lColLb, lColUb, LLSol, false);
		    objBound = -1 * lObj;
		}
	    }
	    else{
		assert(0);
	    }
	}
	else{
	    returnedLowerBound_ = true;
	    isSecondPhaseFinished = true;
	}
	delete secondPhaseModel;
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
MibSZeroSum::findBoundArtCol()
{

    int i, j;
    int numElements(0), rowIndex(0), colIndex(0);
    double mult(0.0), value(0.0), rhs(0.0);
    double tmpMax(0.0), maxValue(0.0);
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
	tmpMax = -1 * rhs;
	for(j = 0; j < numElements; j++){
	    value = elements[j] * mult;
	    colIndex = indices[j];
	    if(value > 0){
		tmpMax += value * origColUb[colIndex];
	    }
	    else{
		tmpMax += value* origColLb[colIndex];
	    }
	}
	if(maxValue - tmpMax < -etol){
	    maxValue = tmpMax;
	}
    }

    if(maxValue <= 0){
	maxValue = 5;
    }
    
    return maxValue;


}

//#############################################################################   
void
MibSZeroSum::findBoundsOnLLCols(double *newLb, double *newUb, double *LLSol,
				bool isFirstTime)
{

    bool areBoundsChanged(false);
    int i;
    double value(0.0);
    int colIndex(0);
    double etol(model_->etol_);
    int lCols(model_->getLowerDim());
    const double *origColLb(model_->getOrigColLb());
    const double *origColUb(model_->getOrigColUb());
    int *lColInd(model_->getLowerColInd());
    double *currentLb = new double[lCols];
    double *currentUb =new double[lCols];

    memcpy(currentLb, newLb, sizeof(double) * lCols);
    memcpy(currentUb, newUb, sizeof(double) * lCols);
   
    for(i = 0; i < lCols; i++){
	colIndex = lColInd[i];
	value = LLSol[i];
	if((isFirstTime) || (value - currentLb[i] < -etol)){
	    areBoundsChanged = true;
	    if(fabs(floor(value + 0.5) - value) <= etol){
		value = floor(value + 0.5);
		if(1 - value + origColLb[colIndex] <= 0){
		    newLb[i] = value - 1;
		}
		else{
		    newLb[i] = origColLb[colIndex];
		}
	    }
	    else{
		newLb[i] = floor(value);
	    }	
	}
	if((isFirstTime) || (value - currentUb[i] > etol)){
	    areBoundsChanged = true;
	    if(fabs(floor(value + 0.5) - value) <= etol){
		value = floor(value + 0.5);
		if(value + 1 - origColUb[colIndex] <= 0){
		    newUb[i] = value + 1;
		}
		else{
		    newUb[i] = origColUb[colIndex];
		}
	    }
	    else{
		newUb[i] = ceil(value);
	    }
	}
    }

    assert(areBoundsChanged);

    delete [] currentLb;
    delete [] currentUb;
    
}

//#############################################################################
void
MibSZeroSum::initialSetUpLowerSolver(CoinPackedMatrix *matrixA2)
{

    std::string feasCheckSolver =
	model_->MibSPar_->entry(MibSParams::feasCheckSolver);
    
    OsiSolverInterface * oSolver = model_->solver();
    int i, j;
    int colIndex(0), rowIndex(0), numElements(0), pos(0), cntInt(0);
    int uCols(model_->getUpperDim());
    int lCols(model_->getLowerDim());
    int lRows(model_->getLowerRowNum());
    int numArtCols(lRows/2);
    int startArtCols(lCols - numArtCols);
    double *origColLb(model_->getOrigColLb());
    double *origColUb(model_->getOrigColUb());
    double *origRowLb(model_->getOrigRowLb());
    double *origRowUb(model_->getOrigRowUb());
    int *uColInd(model_->getUpperColInd());
    int *lColInd(model_->getLowerColInd());
    int *lRowInd(model_->getLowerRowInd());

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
	matrixG2->appendRow(vecG2);
	matrixA2->appendRow(vecA2);
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

    if (feasCheckSolver == "Cbc"){
	lSolver_ = new OsiCbcSolverInterface();
    }else if (feasCheckSolver == "SYMPHONY"){
#ifdef COIN_HAS_SYMPHONY
	lSolver_ = new OsiSymSolverInterface();
#else
	throw CoinError("SYMPHONY chosen as solver, but it has not been enabled",
			"initialSetUpLowerSolver", "MibSZeroSum");
#endif
    }else if (feasCheckSolver == "CPLEX"){
#ifdef COIN_HAS_CPLEX
	lSolver_ = new OsiCpxSolverInterface();
#else
	throw CoinError("CPLEX chosen as solver, but it has not been enabled",
			"initialSetUpLowerSolver", "MibSZeroSum");
#endif
    }else{
	throw CoinError("Unknown solver chosen",
			"initialSetUpLowerSolver", "MibSZeroSum");
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
    int colIndex(0), rowIndex(0);
    double value(0.0);
    int uCols(model_->getUpperDim());
    int lRows(model_->getLowerRowNum());
    int *uColInd(model_->getUpperColInd());
    int *lColInd(model_->getLowerColInd());
    int *lRowInd(model_->getLowerRowInd());
    double *origRowLb(model_->getOrigRowLb());
    double *origRowUb(model_->getOrigRowUb());

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
	rowIndex = lRowInd[i];
	value = mult[i];
	lSolver_->setRowBounds(i, origRowLb[rowIndex] - value,
			       origRowUb[rowIndex] - value);
    }

    delete [] mult;
    delete [] uSol;
}

//#############################################################################
void
MibSZeroSum::solveLowerProblem()
{

    std::string feasCheckSolver(model_->MibSPar_->
				entry(MibSParams::feasCheckSolver));
    int maxThreadsLL(model_->MibSPar_->entry
		     (MibSParams::maxThreadsLL));
    int whichCutsLL(model_->MibSPar_->entry
		    (MibSParams::whichCutsLL));

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
    lSolver_->branchAndBound();
}

