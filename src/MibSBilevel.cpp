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

#include "CoinPackedVector.hpp"
#include "OsiCbcSolverInterface.hpp"

#include "MibSBilevel.hpp"
#include "MibSModel.hpp"
#include "MibSTreeNode.hpp"
#include "MibSSolution.hpp"
#include "MibSHeuristic.hpp"
#include "MibSConfig.hpp"

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
MibSSolType
MibSBilevel::createBilevel(CoinPackedVector* sol, 
			   MibSModel *mibs)
{
  
  /** Splits sol into two parts, upper- and lower-level components **/

  if(!mibs) return MibSNoSol;
  
  model_ = mibs;

  //heuristic_ = new MibSHeuristic(mibs);

  int i(0),j(0);
  int N(model_->numVars_);
  int lN(model_->lowerDim_); // lower-level dimension
  int stocLN(model_->stocLowerDim_);
  int uN(model_->upperDim_); // upper-level dimension
  double etol(model_->BlisPar()->entry(BlisParams::integerTol));
  
  MibSBranchingStrategy branchPar = static_cast<MibSBranchingStrategy>
      (model_->MibSPar_->entry(MibSParams::branchStrategy));

  int solveSecondLevelWhenXYVarsInt(model_->MibSPar_->entry
			   (MibSParams::solveSecondLevelWhenXYVarsInt));
  int solveSecondLevelWhenXVarsInt(model_->MibSPar_->entry
			  (MibSParams::solveSecondLevelWhenXVarsInt));
  int solveSecondLevelWhenLVarsInt(model_->MibSPar_->entry
			  (MibSParams::solveSecondLevelWhenLVarsInt));
  int solveSecondLevelWhenLVarsFixed(model_->MibSPar_->entry
			    (MibSParams::solveSecondLevelWhenLVarsFixed));
  int cutStrategy(model_->MibSPar_->entry
		  (MibSParams::cutStrategy));

  int useLinkingSolutionPool(model_->MibSPar_->entry
			   (MibSParams::useLinkingSolutionPool));

  //saharSto: Delete this
  useLinkingSolutionPool = 0;
  solveSecondLevelWhenXVarsInt = 1;

  MibSSolType storeSol(MibSNoSol);
  
  assert(N == model_->solver()->getNumCols());

  model_->setNumRows(model_->solver()->getNumRows());
  model_->setUpperRowNum(model_->solver()->getNumRows() - model_->getStocLowerRowNum());
  model_->setUpperRowData();
  
  int *indices = sol->getIndices();
  double *values = sol->getElements();
  int numElements(sol->getNumElements()); // number of nonzero elements
  int * fixedInd = model_->fixedInd_;

  
  //saharStoc: In the stochastic version, we do not have
  //upper(lower)SolutionOrd_ and optUpper(Lower)SolutionOrd_
  //will not be NULL only if we return a solution (heuristic or
  //verify that the optimal solution of relaxation is bilevel feasible).
  //if(!upperSolutionOrd_)
  //  upperSolutionOrd_ = new double[uN];
  //if(!lowerSolutionOrd_)
  //  lowerSolutionOrd_ = new double[stocLN];
  if(!optUpperSolutionOrd_)
      optUpperSolutionOrd_ = new double[uN];
  if(!optLowerSolutionOrd_)
      optLowerSolutionOrd_ = new double[stocLN];
  
  //CoinZeroN(upperSolutionOrd_, uN);
  //CoinZeroN(lowerSolutionOrd_, stocLN);
  CoinZeroN(optUpperSolutionOrd_, uN);
  CoinZeroN(optLowerSolutionOrd_, stocLN);
  
  isIntegral_ = true;
  isUpperIntegral_ = true;
  isLinkVarsIntegral_ = true;
  LPSolStatus_ = MibSLPSolStatusUnknown;
  isLinkVarsFixed_ = true;
  shouldPrune_ = false;
  isLowerSolved_ = false;
  isUBSolved_ = false;
  isContainedInLinkingPool_ = false;
  useBilevelBranching_ = true;
  isProvenOptimal_ = true;
  tagInSeenLinkingPool_ = MibSLinkingPoolTagIsNotSet;

  model_->countIteration_ ++;
  /*std::cout << "countIteration = " << model_->countIteration_ << std::endl;
  if(model_->countIteration_ == 821){
      std::cout << "Stop here!" << std::endl;
      }*/
  
  //int * lowerColInd = mibs->getLowerColInd();
  //int * upperColInd = mibs->getUpperColInd();
  int index(0); //uCount(0), lCount(0);

  const double *lower = model_->solver()->getColLower();
  const double *upper = model_->solver()->getColUpper();
  double value;
  for(i = 0; i < numElements; i++){
    index = indices[i];
    value = CoinMax(values[i], lower[index]);
    value = CoinMin(value, upper[index]);
    //if(binarySearch(0, uN - 1, index, upperColInd) >= 0){
    if(index <= uN){
       if(fabs(floor(value + 0.5) - value) > etol){
#if 1
	   if(fixedInd[index] == 1){
	       isLinkVarsIntegral_ = false;
	   }
	   if(mibs->solver()->isInteger(index)){
	       isUpperIntegral_ = false;
	       isIntegral_ = false;
	       LPSolStatus_ = MibSLPSolStatusInfeasible;
	   }
#endif
       }  
    }
    else{
       if(fabs(floor(value + 0.5) - value) > etol){
#if 1
	  //This check is failing when Blis has already declared the solution integral
	  //It's not really needed
	  if(mibs->solver()->isInteger(index)){
	     isIntegral_ = false;
	     LPSolStatus_ = MibSLPSolStatusInfeasible;
	  }
#endif
       }    
    }
  }

  //for(i = 0; i < N; i ++){
  for(i = 0; i < uN; i++){
      //if(binarySearch(0, uN - 1, i, upperColInd) >= 0){
      if((fixedInd[i] == 1) && (fabs(upper[i] - lower[i]) > etol)){
	  isLinkVarsFixed_ = false;
	  break;
      }
  }
  
  /* put the solution in order by integers first */

  int pos(0);

  for (i = 0; i < numElements; i++){
      index = indices[i];
      //pos = binarySearch(0, lN - 1, index, lowerColInd);
      //if (pos < 0){
      if(index <= uN){
	  //pos = binarySearch(0, uN - 1, index, upperColInd);
	  pos = index;
	  if ((mibs->solver()->isInteger(index)) &&
	      (((values[i] - floor(values[i])) < etol) ||
	       ((ceil(values[i]) - values[i]) < etol))){
	      //upperSolutionOrd_[pos] = (double) floor(values[i] + 0.5);
	      optUpperSolutionOrd_[pos] = (double) floor(values[i] + 0.5);
	  }else{
	      //upperSolutionOrd_[pos] = values[i];
	      optUpperSolutionOrd_[pos] = values[i]; 
	  }
	  //optUpperSolutionOrd_[pos] = upperSolutionOrd_[pos];
      }else{
	  pos = index - uN;
	  if ((mibs->solver()->isInteger(index)) &&
	      (((values[i] - floor(values[i])) < etol) ||
	       ((ceil(values[i]) - values[i]) < etol))){
	      //lowerSolutionOrd_[pos] = (double) floor(values[i] + 0.5);
	      optLowerSolutionOrd_[pos] = (double) floor(values[i] + 0.5);
	  }else{
	      //lowerSolutionOrd_[pos] = values[i];
	      optLowerSolutionOrd_[pos] = values[i];
	  }
	  //optLowerSolutionOrd_[pos] = lowerSolutionOrd_[pos];	
      }
  }

  int solType(0), numSolvedScenarios(0);
  std::vector<double> linkSol;
  //saharSto: check the below part for linking pool
  if(isLinkVarsIntegral_){
      //std::vector<double> linkSol;
      for(i = 0; i < uN; i++){
	  //index = upperColInd[i];
	  if(fixedInd[i] == 1){
	      linkSol.push_back(optUpperSolutionOrd_[i]);
	  }
      }
      if(model_->seenLinkingSolutions.find(linkSol) !=
	 model_->seenLinkingSolutions.end()){
	  isContainedInLinkingPool_ = true;
	  solType = model_->seenLinkingSolutions.find(linkSol)->second.tag;
	  numSolvedScenarios =
	      model_->seenLinkingSolutions.find(linkSol)->second.numSolvedScenarios;
      }
  }
	 
  if(isContainedInLinkingPool_){
      tagInSeenLinkingPool_ = static_cast<MibSLinkingPoolTag>(solType);
  }
  if(tagInSeenLinkingPool_ == MibSLinkingPoolTagLowerIsInfeasible){
      LPSolStatus_ = MibSLPSolStatusInfeasible;
  }

  //steps 5-6
  if((isLinkVarsFixed_) && ((tagInSeenLinkingPool_ ==
			     MibSLinkingPoolTagLowerIsInfeasible) ||
			 (tagInSeenLinkingPool_ == MibSLinkingPoolTagUBIsSolved))){
      useBilevelBranching_ = false;
      shouldPrune_ = true;
  }

  //step 7
  if(!shouldPrune_){
      if(((tagInSeenLinkingPool_ == MibSLinkingPoolTagLowerIsFeasible)
	  || (tagInSeenLinkingPool_ == MibSLinkingPoolTagUBIsSolved)) ||
	 (((!isContainedInLinkingPool_) ||
	   (tagInSeenLinkingPool_ == MibSLinkingPoolTagLowerIsIncomplete)) &&
	  (((branchPar == MibSBranchingStrategyLinking) &&
	    (isIntegral_) && (isLinkVarsFixed_)) ||
	   ((branchPar == MibSBranchingStrategyFractional)
	    && (isIntegral_)) ||
	   ((solveSecondLevelWhenXYVarsInt == PARAM_ON) && (isIntegral_)) ||
	   ((solveSecondLevelWhenXVarsInt == PARAM_ON) && (isUpperIntegral_)) ||
	   ((solveSecondLevelWhenLVarsInt == PARAM_ON) && (isLinkVarsIntegral_)) ||
	   ((solveSecondLevelWhenLVarsFixed == PARAM_ON) && (isLinkVarsFixed_ ))))){
	  storeSol = checkBilevelFeasiblity(mibs->isRoot_, numSolvedScenarios, linkSol);
      }
  }
  if(cutStrategy == 1){
      useBilevelBranching_ = false;
  }
  
  //heuristic_->findHeuristicSolutions();

  //delete heuristic_;

  if(storeSol == MibSNoSol){
      delete [] optUpperSolutionOrd_;
      delete [] optLowerSolutionOrd_;
  }

  return storeSol;
}


//#############################################################################
MibSSolType
MibSBilevel::checkBilevelFeasiblity(bool isRoot, int numSolvedScenarios,
				    std::vector<double> linkSol)
{
    //bool warmStartLL(model_->MibSPar_->entry
    //		     (MibSParams::warmStartLL));
    int maxThreadsLL(model_->MibSPar_->entry
		     (MibSParams::maxThreadsLL));
    int whichCutsLL(model_->MibSPar_->entry
		    (MibSParams::whichCutsLL));
    //int probType(model_->MibSPar_->entry
    //		 (MibSParams::bilevelProblemType));
    std::string feasCheckSolver(model_->MibSPar_->entry
				(MibSParams::feasCheckSolver));
    MibSBranchingStrategy branchPar = static_cast<MibSBranchingStrategy>
	(model_->MibSPar_->entry(MibSParams::branchStrategy));
    int computeBestUBWhenXVarsInt(model_->MibSPar_->entry
			      (MibSParams::computeBestUBWhenXVarsInt));
    int computeBestUBWhenLVarsInt(model_->MibSPar_->entry
				     (MibSParams::computeBestUBWhenLVarsInt));
    int computeBestUBWhenLVarsFixed(model_->MibSPar_->entry
				    (MibSParams::computeBestUBWhenLVarsFixed));
    int useLinkingSolutionPool(model_->MibSPar_->entry
			    (MibSParams::useLinkingSolutionPool));

    //saharSto: Delete this
    computeBestUBWhenLVarsInt = 1;
    useLinkingSolutionPool = 0;
    
    double timeLimit(model_->AlpsPar()->entry(AlpsParams::timeLimit));
    double remainingTime(0.0), startTimeVF(0.0), startTimeUB(0.0);
    int beg(0);
    //saharSto: define a parameter for this 
    bool shouldSolveAllScenarios(true);
    bool foundInfeasibility(false);
    double *tmpOptLowerSol;
    double infinity(model_->solver()->getInfinity());
    MibSSolType storeSol(MibSNoSol);
    int lN(model_->lowerDim_); 
    int stocLN(model_->stocLowerDim_);
    int uN(model_->upperDim_); 
    int i(0), index(0), length(0), pos(0);
    int sizeFixedInd(model_->sizeFixedInd_);
    double etol(model_->etol_), objVal(0.0), lowerObj(0.0);
    int numScenarios(model_->getNumScenarios());
    int * fixedInd = model_->fixedInd_;
    int * lowerColInd = model_->getLowerColInd();
    int * upperColInd = model_->getUpperColInd();
    //double *lowerSol = new double[stocLN];
    //CoinFillN(lowerSol, stocLN, 0.0);

    std::vector<double> shouldStoreLowerObjs;
    //std::vector<double> shouldStoreValuesUBSol(stocLN + uN);

    const double * sol = model_->solver()->getColSolution();
    
    /*std::vector<double> linkSol;
    for(i = 0; i < uN; i++){
	index = upperColInd[i];
	if(fixedInd[index] == 1){
	    linkSol.push_back(optUpperSolutionOrd_[i]);
	}
    }*/

    //saharSto: checkThis
    if(((branchPar == MibSBranchingStrategyLinking) &&
       (isIntegral_ == true) && (isLinkVarsFixed_ == true)) ||
       ((branchPar == MibSBranchingStrategyFractional) &&
	(isIntegral_ == true))){
	shouldSolveAllScenarios = true;
    }

    if((shouldSolveAllScenarios == false) && (isIntegral_ == false)){
	shouldSolveAllScenarios = true;
    }

    isProvenOptimal_ = true; 

	//isProvenOptimal_ = true;
    
	/*if (warmStartLL && (feasCheckSolver == "SYMPHONY") && solver_){
	    solver_ = setUpModel(model_->getSolver(), false);
	}else{
	    //if (solver_){
	//	delete solver_;
	   // }
	    solver_ = setUpModel(model_->getSolver(), true);
	}*/

    if(isContainedInLinkingPool_ == true){
	shouldStoreLowerObjs = model_->seenLinkingSolutions.find(linkSol)->second.lowerObjValues;
    }
	    
    assert(numSolvedScenarios <= numScenarios);

    if(numSolvedScenarios <= 0){
	tmpOptLowerSol = new double[stocLN];
	CoinZeroN(tmpOptLowerSol, stocLN); 
    }

    if(isIntegral_ == false){
	foundInfeasibility = true;
    }
    
    for(i = 0; i < numScenarios; i++){
	if(i >= numSolvedScenarios){
	    if(lSolver_){
		lSolver_ = setUpModel(model_->getSolver(), false, i);
	    }
	    else{
		lSolver_ = setUpModel(model_->getSolver(), true, i);
	    }
	    //}

	    OsiSolverInterface *lSolver = lSolver_;

	    remainingTime = timeLimit - model_->broker_->subTreeTimer().getTime();
	    if(remainingTime <= etol){
		shouldPrune_ = true;
	        storeSol = MibSNoSol;
		goto TERM_CHECKBILEVELFEAS;
	    }
	
	    remainingTime = CoinMax(remainingTime, 0.00);
	
	    if (feasCheckSolver == "Cbc"){
		dynamic_cast<OsiCbcSolverInterface *>
		    (lSolver)->getModelPtr()->messageHandler()->setLogLevel(0);
	    }else if (feasCheckSolver == "SYMPHONY"){
#if COIN_HAS_SYMPHONY
	    //dynamic_cast<OsiSymSolverInterface *> 
	    // (lSolver)->setSymParam("prep_level", -1);

		sym_environment *env = dynamic_cast<OsiSymSolverInterface *>
		    (lSolver)->getSymphonyEnvironment();
	    
		/*if (warmStartLL){
		sym_set_int_param(env, "keep_warm_start", TRUE);
		if (probType == 1){ //Interdiction
		    sym_set_int_param(env, "should_use_rel_br", FALSE);
		    sym_set_int_param(env, "use_hot_starts", FALSE);
		    sym_set_int_param(env, "should_warmstart_node", TRUE);
		    sym_set_int_param(env, "sensitivity_analysis", TRUE);
		    sym_set_int_param(env, "sensitivity_bounds", TRUE);
		    sym_set_int_param(env, "set_obj_upper_lim", FALSE);
		}
	    }*/
	    //Always uncomment for debugging!!
		sym_set_dbl_param(env, "time_limit", remainingTime);
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
		lSolver->setHintParam(OsiDoReducePrint);
	        lSolver->messageHandler()->setLogLevel(0);
		CPXENVptr cpxEnv =
		    dynamic_cast<OsiCpxSolverInterface*>(lSolver)->getEnvironmentPtr();
		assert(cpxEnv);
	        CPXsetintparam(cpxEnv, CPX_PARAM_SCRIND, CPX_OFF);
	        CPXsetintparam(cpxEnv, CPX_PARAM_THREADS, maxThreadsLL);
#endif
	    }

	    //step 8
	    /*if (warmStartLL && feasCheckSolver == "SYMPHONY"){
	    lSolver->resolve();
	    setWarmStart(lSolver->getWarmStart());
	}else{*/
	    startTimeVF = model_->broker_->subTreeTimer().getTime();
	    lSolver->branchAndBound();
	    model_->timerVF_ += model_->broker_->subTreeTimer().getTime() - startTimeVF;
	  //}

	    model_->counterVF_ ++;
	    //saharSto: check this
	    isLowerSolved_ = true;

	    if((feasCheckSolver == "SYMPHONY") && (sym_is_time_limit_reached
						   (dynamic_cast<OsiSymSolverInterface *>
						    (lSolver)->getSymphonyEnvironment()))){
		shouldPrune_ = true;
	        storeSol = MibSNoSol;
	        goto TERM_CHECKBILEVELFEAS;
	    }
	    else if(!lSolver->isProvenOptimal()){
		LPSolStatus_ = MibSLPSolStatusInfeasible;
		isProvenOptimal_ = false;
		if(useLinkingSolutionPool){
	    //step 10
	    //Adding x_L to set E
		    addSolutionToSeenLinkingSolutionPool
			(MibSLinkingPoolTagLowerIsInfeasible, 0,
			 shouldStoreLowerObjs, 0.0, linkSol);
		}
		if(isLinkVarsFixed_){
		    useBilevelBranching_ = false;
	            shouldPrune_ = true;
		}
		break;
	    }
	    else{
	    //const double * sol = model_->solver()->getColSolution();
		objVal = lSolver->getObjValue() * model_->getLowerObjSense();

	         objVal_ = objVal;

	         const double * values = lSolver->getColSolution();
	    
	         if(useLinkingSolutionPool){
		//std::copy(values, values + lN, shouldStoreValuesLowerSol.begin());
		
		//step 12
		//Adding x_L to set E  
		//addSolutionToSeenLinkingSolutionPool
		//(MibSLinkingPoolTagLowerIsFeasible, shouldStoreValuesLowerSol, objVal_);
		//shouldStoreValuesLowerSol.clear();
		     shouldStoreLowerObjs.push_back(objVal_);
		 }
		 if(tmpOptLowerSol != NULL){
		     CoinDisjointCopyN(values, lN, tmpOptLowerSol + i * lN);
		 }
	    }
	}
	else{
	    objVal = shouldStoreLowerObjs[i];
	    objVal_ = objVal;
	}

	MibSTreeNode * node = static_cast<MibSTreeNode *>(model_->activeNode_);
	MibSTreeNode * parent =
	    static_cast<MibSTreeNode *>(model_->activeNode_->getParent());
	    
	if((!node->isBoundSet())
	   && (node->getIndex() != 0)){
	    double parentBound = parent->getLowerUB();
	    node->setLowerUB(parentBound);
	    node->setIsBoundSet(true);
	}

	if(objVal > node->getLowerUB()){
	    node->setLowerUB(objVal);
	    node->setIsBoundSet(true);
	}
    
	/*if (!warmStartLL){
	    delete solver_;
	}*/


    //step 13
    /*if(((!useLinkingSolutionPool) && (isProvenOptimal_)) ||
       ((tagInSeenLinkingPool_ == MibSLinkingPoolTagLowerIsFeasible) ||
	(tagInSeenLinkingPool_ == MibSLinkingPoolTagUBIsSolved))){*/

	//OsiSolverInterface *UBSolver;
	
	//double *lowerSol = new double[lN];
	//CoinFillN(lowerSol, lN, 0.0);

	/*if(useLinkingSolutionPool){
	    //get optimal value  of (VF) from solution pool
	    //model_->it = seenLinkingSolutions.find(linkSol);
	    //objVal = model_->it->second.lowerObjVal1;
	    objVal = model_->seenLinkingSolutions[linkSol].lowerObjValue; 
	    //objVal = seenLinkingSolutions.find(linkSol).
	    objVal_ = objVal;
	    std::copy(model_->seenLinkingSolutions[linkSol].lowerSolution.begin(),
		      model_->seenLinkingSolutions[linkSol].lowerSolution.end(), lowerSol);
		      }*/
  
	if(foundInfeasibility == false){
	    //lowerObj = getLowerObj(sol, model_->getLowerObjSense());
	    lowerObj = 0.0;
	    beg = i * lN;
	    for(j = beg; j < beg + lN; j++){
		lowerObj += sol[j + uN] * lObjCoeffs[j - beg];
	    }
	    assert((objVal - lowerObj) <= etol);
	}
       	
	LPSolStatus_ = MibSLPSolStatusInfeasible;

	
	if((foundInfeasibility == false) && (fabs(objVal - lowerObj) > etol)){
	    foundInfeasibility = true;
	    if((shouldSolveAllScenarios == false) && (i < numScenarios)){
		delete [] tmpOptLowerSol;
		addSolutionToSeenLinkingSolutionPool
		    (MibSLinkingPoolTagLowerIsIncomplete, i + 1,
		     shouldStoreLowerObjs, 0.0, linkSol);
		break;
	    }
	}
    }

    if(i == numScenarios){
	if(foundInfeasibility == false){
	    LPSolStatus_ = MibSLPSolStatusFeasible;
	    useBilevelBranching_ = false;
	    shouldPrune_ = true;
	    storeSol = MibSRelaxationSol;
	    }
	if(numSolvedScenarios < numScenarios){
	    addSolutionToSeenLinkingSolutionPool
		(MibSLinkingPoolTagLowerIsFeasible,
		 numScenarios, shouldStoreLowerObjs, 0.0, linkSol);
	}
    }

    bool seenAllScenarios(false);

    if((i == numScenarios) || (numSolvedScenarios == numScenarios)){
	seenAllScenarios = true;
    }
    
	//step 15
	/** Current solution is bilevel feasible **/
	/*if((fabs(objVal - lowerObj) < etol) && (isIntegral_)){
	    LPSolStatus_ = MibSLPSolStatusFeasible;
	    useBilevelBranching_ = false;
	    shouldPrune_ = true;
	    storeSol = MibSRelaxationSol;
	}
	else{
	    memcpy(optLowerSolutionOrd_, lowerSol, sizeof(double) * lN);
	}*/

    if(!shouldPrune_){
	//step 18
	if(((seenAllScenarios == true) && (tagInSeenLinkingPool_ != MibSLinkingPoolTagUBIsSolved)) &&
	   (((branchPar == MibSBranchingStrategyLinking) &&
	     (isIntegral_) && (isLinkVarsFixed_)) ||
	    ((computeBestUBWhenXVarsInt == PARAM_ON) && (isUpperIntegral_)) ||
	    ((computeBestUBWhenLVarsInt == PARAM_ON)) ||
	    ((computeBestUBWhenLVarsFixed == PARAM_ON) && (isLinkVarsFixed_)))){
	    OsiSolverInterface *UBSolver;
	    if(UBSolver_){
		UBSolver_ = setUpUBModel(model_->getSolver(), shouldStoreLowerObjs,
					 false);
	    }
	    else{
		UBSolver_ =setUpUBModel(model_->getSolver(), shouldStoreLowerObjs,
					true);
	    }

	    UBSolver = UBSolver_;

	    remainingTime = timeLimit - model_->broker_->subTreeTimer().getTime();

	    if(remainingTime <= etol){
		shouldPrune_ = true;
		storeSol = MibSNoSol;
		goto TERM_CHECKBILEVELFEAS;
	    }
		
	    remainingTime = CoinMax(remainingTime, 0.00);

	    if (feasCheckSolver == "Cbc"){
		dynamic_cast<OsiCbcSolverInterface *>
		    (UBSolver)->getModelPtr()->messageHandler()->setLogLevel(0);
	    }else if (feasCheckSolver == "SYMPHONY"){
#if COIN_HAS_SYMPHONY
		    //dynamic_cast<OsiSymSolverInterface *>
		    // (lSolver)->setSymParam("prep_level", -1);
		sym_environment *env = dynamic_cast<OsiSymSolverInterface *>
		    (UBSolver)->getSymphonyEnvironment();
		    //Always uncomment for debugging!!
		sym_set_dbl_param(env, "time_limit", remainingTime);
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
		UBSolver->setHintParam(OsiDoReducePrint);
		UBSolver->messageHandler()->setLogLevel(0);
		CPXENVptr cpxEnv =
		    dynamic_cast<OsiCpxSolverInterface*>(UBSolver)->getEnvironmentPtr();
		assert(cpxEnv);
		CPXsetintparam(cpxEnv, CPX_PARAM_SCRIND, CPX_OFF);
		CPXsetintparam(cpxEnv, CPX_PARAM_THREADS, maxThreadsLL);
#endif
	    }

		//step 19
	    startTimeUB = model_->broker_->subTreeTimer().getTime(); 
	    UBSolver->branchAndBound();
	    model_->timerUB_ += model_->broker_->subTreeTimer().getTime() - startTimeUB;
	    model_->counterUB_ ++;
	    isUBSolved_ = true;
	    if((feasCheckSolver == "SYMPHONY") && (sym_is_time_limit_reached
						   (dynamic_cast<OsiSymSolverInterface *>
						    (UBSolver)->getSymphonyEnvironment()))){
		shouldPrune_ = true;
		storeSol = MibSNoSol;
		goto TERM_CHECKBILEVELFEAS;
	    }
	    else if (UBSolver->isProvenOptimal()){
		isProvenOptimal_ = true;
		const double * valuesUB = UBSolver->getColSolution();
		    //std::copy(valuesUB, valuesUB + uN + lN, shouldStoreValuesUBSol.begin());
		 for (i = 0; i < uN + lN; i++){
			//pos = binarySearch(0, uN - 1, i, upperColInd);
		     if (i < uN){
			 if ((UBSolver->isInteger(i)) &&
			     (((valuesUB[i] - floor(valuesUB[i])) < etol) ||
			      ((ceil(valuesUB[i]) - valuesUB[i]) < etol))){
			     optUpperSolutionOrd_[i] = (double) floor(valuesUB[i] + 0.5);
			 }else{
			     optUpperSolutionOrd_[i] = (double) valuesUB[i];
			 }
		     }else{
			    //pos = binarySearch(0, lN - 1, i, lowerColInd);
			 if ((UBSolver->isInteger(i)) &&
			     (((valuesUB[i] - floor(valuesUB[i])) < etol) ||
			      ((ceil(valuesUB[i]) - valuesUB[i]) < etol))){
			     optLowerSolutionOrd_[i - uN] = (double) floor(valuesUB[i] + 0.5);
			 }else{
			     optLowerSolutionOrd_[i - uN] = (double) valuesUB[i];
			 }
		     }
		 }
		 objVal = UBSolver->getObjValue() * model_->solver()->getObjSense();
		 storeSol = MibSHeurSol;
		}else{
		    isProvenOptimal_ = false;
		    objVal = infinity;
		}
		//step 22
		//Adding x_L to set E
		if(useLinkingSolutionPool){
		    addSolutionToSeenLinkingSolutionPool
		  (MibSLinkingPoolTagUBIsSolved, numScenarios, shouldStoreLowerObjs,
		   objVal, linkSol);
		}
		//shouldStoreValuesUBSol.clear();
	    
		//step 23
		if(isLinkVarsFixed_){
		    useBilevelBranching_ = false;
		    //isProvenOptimal_ = false;
		    shouldPrune_ = true;
		}	
	    }
	else if ((tmpOptLowerSol != NULL) &&
		 (tagInSeenLinkingPool_ != MibSLinkingPoolTagUBIsSolved)){
	    for (i = 0; i < lN; i++){
		index = uN + i;
		value = tmpOptLowerSol[i];
		if ((model_->solver()->isInteger(index)) &&
		    (((value - floor(tmpOptLowerSol[i])) < etol) ||
		     ((ceil(value) - value) < etol))){
		    optLowerSolutionOrd_[i] = (double) floor(value + 0.5);
		}else{
		    optLowerSolutionOrd_[i] = (double) value;
		}
	    }
	    if(isUpperIntegral_ == true){
		storeSol = MibSHeurSol;
	    }
	}
    }
	//}

 TERM_CHECKBILEVELFEAS:
    if(tmpOptLowerSol != NULL){
	delete [] tmpOptLowerSol;
    }
    //delete [] lowerSol;

    return storeSol;
}
	    
//#############################################################################
void 
MibSBilevel::gutsOfDestructor()
{

  if(optUpperSolutionOrd_) delete [] optUpperSolutionOrd_;
  if(optLowerSolutionOrd_) delete [] optLowerSolutionOrd_;
  if(upperSolutionOrd_) delete [] upperSolutionOrd_;
  if(lowerSolutionOrd_) delete [] lowerSolutionOrd_;
  if(lSolver_) delete lSolver_;
  if(UBSolver_) delete UBSolver_;
  //delete heuristic_;
}


//#############################################################################
OsiSolverInterface *
    MibSBilevel::setUpUBModel(OsiSolverInterface * oSolver,
			      std::vector<double> stocObjValLL,
			      bool newOsi, const double *lpSol)
{
    
    std::string feasCheckSolver =
	model_->MibSPar_->entry(MibSParams::feasCheckSolver);

    OsiSolverInterface * nSolver;

    if (!lpSol){
	lpSol = oSolver->getColSolution();
    }

    int * fixedInd = model_->fixedInd_;

    int i(0), j(0), index1(0), beg(0);
    double value(0.0);
    double minusInf(-1 * oSolver->getInfinity());

    int numScenarios(model_->getNumScenarios());
    int uCols(model_->getUpperDim());
    //int uRows(model_->getOrigUpperRowNum());
    //int lRows(model_->getLowerRowNum());
    int lCols(model_->getLowerDim());
    //int rowNum(model_->getNumOrigCons() + 1);
    int rowNum(model_->getNumOrigCons() + numScenarios);
    int colNum(model_->getNumOrigVars());
    //int * uColIndices(model_->getUpperColInd());


    if(newOsi){
	double objSense(model_->getLowerObjSense());
	double uObjSense(1);
        //int * lColIndices(model_->getLowerColInd());
	//int * uRowIndices(model_->getOrigUpperRowInd());
	//int * lRowIndices(model_->getLowerRowInd());
	const double * origColLb(model_->getOrigColLb());
	const double * origColUb(model_->getOrigColUb());
	const double * origRowLb(model_->getOrigRowLb());
	const double * origRowUb(model_->getOrigRowUb());
	double * lObjCoeffs(model_->getLowerObjCoeffs());
	const double * uObjCoeffs(oSolver->getObjCoefficients());
	int tmpRowNum(rowNum - numScenarios);

	CoinPackedMatrix matrix = *model_->origConstCoefMatrix_;
	matrix.reverseOrdering();

	double * rowUb = new double[rowNum];
        double * rowLb = new double[rowNum];
        double * colUb = new double[colNum];
        double * colLb = new double[colNum];

        CoinZeroN(colUb, colNum);
        CoinZeroN(colLb, colNum);

	/** Set the row bounds **/
	memcpy(rowLb, origRowLb, sizeof(double) * tmpRowNum);
	memcpy(rowUb, origRowUb, sizeof(double) * tmpRowNum);
	CoinFillN(rowLb + tmpRowNum, rowNum - tmpRowNum, minusInf);

	//Set the col bounds
	memcpy(colLb, origColLb, sizeof(double) * colNum);
	memcpy(colUb, origColUb, sizeof(double) * colNum);
	for(i = 1; i < numScenarios; i++){
	    beg = uCols + i * lCols;
	    CoinDisjointCopyN(origColLb + uCols, lCols, colLb + beg);
	    CoinDisjointCopyN(origColUb + uCols, lCols, colUb + beg);
	}

	for(i = 0; i < uCols; i++){
	    index1 = i;
	    if(fixedInd[index1] == 1){
		colLb[index1] = floor(lpSol[index1] + 0.5);
	        colUb[index1] = colLb[index1];
	    }
	}
    
        if (feasCheckSolver == "Cbc"){
	    nSolver = new OsiCbcSolverInterface();
	}else if (feasCheckSolver == "SYMPHONY"){
#ifdef COIN_HAS_SYMPHONY
	    nSolver = new OsiSymSolverInterface();
#else
	    throw CoinError("SYMPHONY chosen as solver, but it has not been enabled",
			    "setUpUBModel", "MibsBilevel");
#endif
	}else if (feasCheckSolver == "CPLEX"){
#ifdef COIN_HAS_CPLEX
	    nSolver = new OsiCpxSolverInterface();
#else
	    throw CoinError("CPLEX chosen as solver, but it has not been enabled",
			    "setUpUBModel", "MibsBilevel");
#endif
	}else{
	    throw CoinError("Unknown solver chosen",
			    "setUpUBModel", "MibsBilevel");
	}
    
        int * integerVars = new int[colNum];
        double * objCoeffs = new double[colNum];
        //double * newRow = new double[lCols];
        CoinFillN(integerVars, colNum, 0);
        CoinFillN(objCoeffs, colNum, 0.0);
        //CoinFillN(newRow, lCols, 0.0);

        int intCnt(0);

	/** Fill in array of integer variables **/
	for(i = 0; i < colNum; i++){
	    if(oSolver->isInteger(i)){
		integerVars[intCnt] = i;
	        intCnt ++;
	    }
	}

        CoinDisjointCopyN(uObjCoeffs, colNum, objCoeffs);
        //CoinDisjointCopyN(lObjCoeffs, lCols, newRow);

        //for(i = 0; i < lCols; i++){
	//  newRow[i] = newRow[i] * objSense;
	//}

        //rowUb[rowNum-1] = objValLL;
        //rowLb[rowNum-1] = -1 * (oSolver->getInfinity());

        CoinPackedMatrix * newMat = new CoinPackedMatrix(false, 0, 0);
        newMat->setDimensions(0, colNum);

        for(i = 0; i < rowNum - 1; i++){
	    newMat->appendRow(matrix.getVector(i));
	}

        CoinPackedVector row;
	for(i = 0; i < numScenarios; i++){
	    index1 = uCols + i * lCols;
	    for(j = 0; j < lCols; j++){
		index1 ++;
		row.insert(index1, lObjCoeffs[j] * objSense);
	    }
	    newMat->appendRow(row);
	    row.clear();
	    rowUb[tmpRowNum + i] = stocObjValLL[i];
	}
	    
	/*for(i = 0; i < lCols; i++){
	    index1 = lColIndices[i];
	    row.insert(index1, newRow[i]);
	}
        newMat->appendRow(row);*/

        nSolver->loadProblem(*newMat, colLb, colUb,
			     objCoeffs, rowLb, rowUb);

        for(i = 0; i < intCnt; i++){
	    nSolver->setInteger(integerVars[i]);
	}

        nSolver->setObjSense(uObjSense); //1 min; -1 max

        nSolver->setHintParam(OsiDoReducePrint, true, OsiHintDo);
    
        delete [] integerVars;
        delete [] rowUb;
        delete [] rowLb;
        delete [] colUb;
        delete [] colLb;
        delete [] objCoeffs;
        //delete [] newRow;
        delete newMat;
    }
    else{
	nSolver = UBSolver_;
	//nSolver->setRowUpper(rowNum-1, objValLL);
	for(i = 0; i < numScenarios; i++){
	    nSolver->setRowUpper(tmpRowNum + i,
				 stocObjValLL[i]);
	}
	for(i = 0; i < uCols; i++){
	    index1 = i;
	    if(fixedInd[index1] == 1){
		value = floor(lpSol[index1] + 0.5);
		nSolver->setColLower(index1, value);
		nSolver->setColUpper(index1, value);
	    }
	}
    }

    return nSolver;

}

//#############################################################################
OsiSolverInterface *
MibSBilevel::setUpModel(OsiSolverInterface * oSolver, bool newOsi, int numScenario,
			const double *lpSol)
{

  /** Create lower-level model with fixed upper-level vars **/

  //int probType =
    //model_->MibSPar_->entry(MibSParams::bilevelProblemType);

    //bool warmStartLL =
    //model_->MibSPar_->entry(MibSParams::warmStartLL);

    //bool doDualFixing =
    //model_->MibSPar_->entry(MibSParams::doDualFixing);

  std::string feasCheckSolver =
    model_->MibSPar_->entry(MibSParams::feasCheckSolver);

  OsiSolverInterface * nSolver;

  bool signsChanged(false);
  double etol(model_->etol_);
  int i(0), j(0), index1(0), index2(0);
  double mult(0.0);
  int uCols(model_->getUpperDim());
  int origURows(model_->getOrigUpperRowNum());
  int lRows(model_->getLowerRowNum());
  int lCols(model_->getLowerDim());
  int * uColIndices = model_->getUpperColInd();
  int * lColIndices = model_->getLowerColInd();
  int * lRowIndices = model_->getLowerRowInd();

  double objSense(model_->getLowerObjSense());  
  double * lObjCoeffs = model_->getLowerObjCoeffs();
     
  const CoinPackedMatrix * matrix = oSolver->getMatrixByRow();
  double coeff(0.0);
  if (!lpSol){
     lpSol = oSolver->getColSolution();   
  }
  const double * origRowLb = model_->getOrigRowLb();
  const double * origRowUb = model_->getOrigRowUb();
  const double * origColLb = model_->getOrigColLb();
  const double * origColUb = model_->getOrigColUb();
  const double * lColLbInLProb = model_->getLColLbInLProb();
  const double * lColUbInLProb = model_->getLColUbInLProb();
  const char *origRowSense = model_->getOrigRowSense();
  double * rowUb = new double[lRows];
  double * rowLb = new double[lRows];
  //CoinZeroN(rowUb, lRows);
  //CoinZeroN(rowLb, lRows);
  //CoinZeroN(colUb, lCols);
  //CoinZeroN(colLb, lCols);

  
  if (newOsi){
      /** Set the row bounds **/
      for(i = 0; i < lRows; i++){
	  //index1 = lRowIndices[i];
	  index1 = origURows + numScenario * lRows + i;
	  rowLb[i] = floor(origRowLb[index1] + 0.5);
	  rowUb[i] = floor(origRowUb[index1] + 0.5);
      }
      
     if (feasCheckSolver == "Cbc"){
	nSolver = new OsiCbcSolverInterface();
     }else if (feasCheckSolver == "SYMPHONY"){
#ifdef COIN_HAS_SYMPHONY
	nSolver = new OsiSymSolverInterface();
#else
	throw CoinError("SYMPHONY chosen as solver, but it has not been enabled",
			"setUpModel", "MibsBilevel");
#endif
     }else if (feasCheckSolver == "CPLEX"){
#ifdef COIN_HAS_CPLEX
	nSolver = new OsiCpxSolverInterface();
#else
	throw CoinError("CPLEX chosen as solver, but it has not been enabled",
			"setUpModel", "MibsBilevel");
#endif
     }else{
	throw CoinError("Unknown solver chosen",
			"setUpModel", "MibsBilevel");
     }

     double * colUb = new double[lCols];
     double * colLb = new double[lCols];
     int * integerVars = new int[lCols];
     double * objCoeffs = new double[lCols];
     
     CoinFillN(integerVars, lCols, 0);
     //CoinZeroN(objCoeffs, lCols);

     // set the col bounds
     if(lColLbInLProb == NULL){
	 for(i = 0; i < lCols; i++){
	     colLb[i] = origColLb[lColIndices[i]];
	     //colUb[i] = origColUb[lColIndices[i]];
	 }
     }
     else{
	 memcpy(colLb, lColLbInLProb, sizeof(double) * lCols);
     }

     if(lColUbInLProb == NULL){
	 for(i = 0; i < lCols; i++){
	     colUb[i] = origColUb[lColIndices[i]];
	 }
     }
     else{
	 memcpy(colUb, lColUbInLProb, sizeof(double) * lCols);
     }
     
     int intCnt(0);
  
     /** Fill in array of lower-level integer variables **/
     for(i = 0; i < lCols; i++){
	index1 = lColIndices[i];
	if(oSolver->isInteger(index1)){
	   integerVars[intCnt] = i;
	   intCnt++;
	}
     }
     
     CoinDisjointCopyN(lObjCoeffs, lCols, objCoeffs);
     /*for(i = 0; i < lCols; i++){
	 objCoeffs[i] = objCoeffs[i] * objSense;
	 }*/
     
     CoinPackedMatrix * newMat = new CoinPackedMatrix(false, 0, 0);
     newMat->setDimensions(0, lCols);
     CoinShallowPackedVector origRow;
     CoinPackedVector row;
     CoinPackedVector rowG2;
     int numElements(0), pos(0);

     for(i = 0; i < lRows; i++){
	 index1 = lRowIndices[i];
	 origRow = matrix->getVector(index1);
	 row = origRow;
	 numElements = row.getNumElements();
	 const int *indices = row.getIndices();
	 const double *elements = row.getElements();
	 for(j = 0; j < numElements; j++){
	     //pos = binarySearch(0, lCols -1, indices[j], lColIndices);
	     if(indices[j] >= uCols){
		 pos = indices[j];
		 //if(pos >= 0){
		 rowG2.insert(pos, elements[j]);
	     }
	 }
	 newMat->appendRow(rowG2);
	 rowG2.clear();
	 row.clear();
     }

     /*
       nSolver->assignProblem(newMat, colLb, colUb,
       objCoeffs, rowLb, rowUb);
     */
     
     nSolver->loadProblem(*newMat, colLb, colUb,
			  objCoeffs, rowLb, rowUb);
     
     for(i = 0; i < intCnt; i++){
	nSolver->setInteger(integerVars[i]);
     }
     //nSolver->setInteger(integerVars, intCnt);
     
     nSolver->setObjSense(objSense); //1 min; -1 max
     
     nSolver->setHintParam(OsiDoReducePrint, true, OsiHintDo);

     /*#if 0
     if(0){
	dynamic_cast<OsiCbcSolverInterface *> 
	   (nSolver)->getModelPtr()->messageHandler()->setLogLevel(0);
     }
     else{
#ifdef COIN_HAS_SYMPHONY
	dynamic_cast<OsiSymSolverInterface *> 
	   (nSolver)->setSymParam("prep_level", -1);
	
	dynamic_cast<OsiSymSolverInterface *> 
	   (nSolver)->setSymParam("verbosity", -2);
	
	dynamic_cast<OsiSymSolverInterface *> 
	   (nSolver)->setSymParam("max_active_nodes", 1);
#endif
     }
#endif*/
     delete [] colLb;
     delete [] colUb;
     delete [] integerVars;
     delete [] objCoeffs;
     delete newMat;

  }else{
     nSolver = lSolver_;
     
     signsChanged = false;
     if(feasCheckSolver == "SYMPHONY"){
	 for(i = 0; i < lRows; i++){
	     index1 = lRowIndices[i];
	     if(origRowSense[index1] == 'G'){
		     signsChanged = true;
		     break;
	     }
	 }
     }

     if(!signsChanged){
	 for(i = 0; i < lRows; i++){
	     //index1 = lRowIndices[i];
	     index1 = origURows + numScenario * lRows + i;
	     rowLb[i] = origRowLb[index1];
	     rowUb[i] = origRowUb[index1];
	 }
     }
     else{
	 for(i = 0; i < lRows; i++){
	     //index1 = lRowIndices[i];
	     index1 = origURows + numScenario * lRows + i;
	     if(origRowSense[index1] == 'G'){
		 rowLb[i] = -origRowUb[index1];
		 rowUb[i] = -origRowLb[index1];
	     }
	     else{
		 rowLb[i] = origRowLb[index1];
		 rowUb[i] = origRowUb[index1];
	     }
	 }
     }
     
  }
  
  /*#define SYM_VERSION_IS_WS strcmp(SYMPHONY_VERSION, "WS")  

#if SYMPHONY_VERSION_IS_WS
  if (feasCheckSolver == "SYMPHONY" && probType == 1 && warmStartLL &&
      !newOsi && doDualFixing){ //Interdiction

     // Get upper bound from best known (feasible) lower level solution and 
    //	 try to fix additional variables by sensitivity analysis

     std::vector<std::pair<AlpsKnowledge*, double> > solutionPool;
     model_->getKnowledgeBroker()->
	getAllKnowledges(AlpsKnowledgeTypeSolution, solutionPool);

     const double * sol; 
     double objval, Ub(objSense*nSolver->getInfinity());
     BlisSolution* blisSol;
     std::vector<std::pair<AlpsKnowledge*, double> >::const_iterator si;
     for (si = solutionPool.begin(); si != solutionPool.end(); ++si){
	blisSol = dynamic_cast<BlisSolution*>(si->first);
	sol = blisSol->getValues();
	for (i = 0; i < uCols; i++){
	   if (lpSol[uColIndices[i]] > 1 - etol &&
	       sol[lColIndices[i]] > 1-etol){
	      break;
	   }
	}
	if (i == uCols && -objSense*blisSol->getQuality() < Ub){
	   Ub = -objSense*blisSol->getQuality();
	}
     }

     // Figure out which variables get fixed by the upper level solution

     int *newUbInd = new int[uCols];
     int *newLbInd = new int[uCols];
     double *newUbVal = new double[uCols];
     double *newLbVal = new double[uCols];
     double newLb;
     for (i = 0; i < uCols; i++){
	newUbInd[i] = uColIndices[i];
	newLbInd[i] = uColIndices[i];
	newLbVal[i] = 0;
	if (lpSol[uColIndices[i]] > 1 - etol){
	   newUbVal[i] = 0;
	}else{
	   newUbVal[i] = 1;
	}
     }
	
     // If we found an upper bound, then do the dual fixing

     if (Ub < objSense*nSolver->getInfinity()){
	sym_environment *env = 
	   dynamic_cast<OsiSymSolverInterface *>(nSolver)->getSymphonyEnvironment();

	for (i = 0; i < uCols; i++){
	   if (newUbVal[i] == 1){
	      // Try fixing it to zero
	      newUbVal[i] = 0; 
	      sym_get_lb_for_new_rhs(env, 0, NULL, NULL, uCols, newLbInd,
				     newLbVal, uCols, newUbInd, newUbVal,
				     &newLb);
	      if (objSense*newLb > Ub + etol){
		 //Victory! This variable can be fixed to 1 permanently
		 newLbVal[i] = 1;
	      }
	      //Set upper bound back to 1
	      newUbVal[i] = 1;
	      
	      if (newLbVal[i] == 0){
		 // Try fixing it to one
		 newLbVal[i] = 1;
		 sym_get_lb_for_new_rhs(env, 0, NULL, NULL, uCols, newLbInd,
					newLbVal, uCols, newUbInd, newUbVal,
					&newLb);
		 if (objSense*newLb > Ub + etol){
		    //Victory! This variable can be fixed to 0 permanently
		    newUbVal[i] = 0;
		 }
		 newLbVal[i] = 0;
	      }
	   }
	   
	}
     }
     
     // Now set the row bounds to account for fixings 
     
     // This is probably very frgaile. Assuming interdiction 
     //	 rows come last. It would be better to set variable 
    //	 bounds directly, but this doesn't seem to work right now.
     int iRowStart = lRows-uCols; 

#if 1
     for(i = iRowStart; i < lRows; i++){
	nSolver->setRowLower(i, newLbVal[i-iRowStart]);
	nSolver->setRowUpper(i, newUbVal[i-iRowStart]);
     }
#else
     for(i = 0; i < uCols; i++){
	nSolver->setColLower(i, newLbVal[i]);
	nSolver->setColUpper(i, newUbVal[i]);
     }
#endif     
     delete[] newUbInd;
     delete[] newLbInd;
     delete[] newUbVal;
     delete[] newLbVal;
  }else{
#endif*/
     //FIXME: NEED TO GET ROW SENSE HERE
     
     /** Get contribution of upper-level columns **/
      double * upComp = new double[lRows];
      CoinFillN(upComp, lRows, 0.0);

      for(i = 0; i < lRows; i++){
          index1 = origURows + numScenario * lRows + i;
	  mult = 1;
	  if((signsChanged) && (origRowSense[index1] == 'G')){
	      mult = -1;
	  }
	  for(j = 0; j < uCols; j++){
	      index2 = uColIndices[j];
	      coeff = matrix->getCoefficient(index1, index2);
	      if (coeff != 0){
		  upComp[i] += mult * coeff * lpSol[index2];
	      }
	  }
      }
     
     /** Correct the row bounds to account for fixed upper-level vars **/
     
     for(i = 0; i < lRows; i++){
	 nSolver->setRowLower(i, floor(rowLb[i] - upComp[i] + 0.5));
	 nSolver->setRowUpper(i, floor(rowUb[i] - upComp[i] + 0.5));
     }
     
     delete [] upComp;

     /*#if SYMPHONY_VERSION_IS_WS
  }
  #endif*/

  //I don't think this is needed
  //if(!getWarmStart())
  //  setWarmStart(nSolver->getWarmStart());

  delete [] rowLb;
  delete [] rowUb;
  
  return nSolver;

}

//#############################################################################
int
MibSBilevel::findIndex(int index, int size, int * indices)
{

  int i(0);
  int pos(-1);

  for(i = 0; i < size; i++){
    if(indices[i] == index)
       pos = i;
  }

  return pos;

}

//#############################################################################
int 
MibSBilevel::binarySearch(int start, int stop, int index, int * indexArray)
{
   int i(0);
   int pos(-1);

  //=========================================================================
  // If the list is short enough, finish with sequential search
  // Otherwise, call binary search recursively
  //=========================================================================

   //FIXME: CHANGED THIS 7/15
   //NOW USES ONLY SEQUENTIAL SEARCH
   //BINARY REQUIRES ORDERED ARRAY
   //OUR ARRAY IS ORDERED BY INTEGERS FIRST, NOT INDICES


   //if((stop - start) < 4){
   if(1){
      for(i = start; i < stop + 1; ++i){
	 if(indexArray[i] == index){
	    pos = i;
	    break;
	 }
      }
      
      return pos;
   }
   else{ 
      
      int midpoint((start + stop)/2);
      int val(indexArray[midpoint]);
      
      if(val == index){
	 pos = midpoint;
	 return pos;
      }
      else if(val > index){
	 pos = binarySearch(start, midpoint - 1, index, indexArray);
      }
      else{
	 pos = binarySearch(midpoint + 1, stop, index, indexArray);
      }
   }
   return pos;
}



//#############################################################################
double
MibSBilevel::getLowerObj(const double * sol, double objSense)
{

   int lCols(model_->getLowerDim());
   int * lColIndices = model_->getLowerColInd();
   double * lObjCoeffs = model_->getLowerObjCoeffs();

   int i(0), index(0);
   double objVal(0.0);


   for(i = 0; i < lCols; i++){
      index = lColIndices[i];
      if(0){
	std::cout << "sol[" << index << "]: " << sol[index] << std::endl;
	std::cout << "lObjCoeff[" << i << "]: " << lObjCoeffs[i] << std::endl;
      }      
      objVal += lObjCoeffs[i] * sol[index];
   }

   return objVal * objSense;

}
//#############################################################################
void
    MibSBilevel::addSolutionToSeenLinkingSolutionPool(MibSLinkingPoolTag solTag,
						      int numSolvedScenarios, std::vector<double>
						      &shouldStoreLowerObjs, double ubObjValue,
						      std::vector<double> linkSol)
{
    int i(0),index(0);
    int uN(model_->upperDim_);
    int lN(model_->lowerDim_);
    int * upperColInd = model_->getUpperColInd();
    int * fixedInd = model_->fixedInd_;
    int solType = static_cast<int>(solTag);

    LINKING_SOLUTION linkingSolution;
    
    /*std::vector<double> linkSol;
    for(i = 0; i < uN; i++){
	index = upperColInd[i];
	if(fixedInd[index] == 1){
	    linkSol.push_back(upperSolutionOrd_[i]);
	}
	}*/

    tagInSeenLinkingPool_ = solTag;
    linkingSolution.lowerObjValues(0);
    linkingSolution.lowerObjValues.clear();
    //linkingSolution.lowerSolution.push_back(0);
    //linkingSolution.UBSolution.push_back(0);
    //linkingSolution.lowerSolution.clear();
    //linkingSolution.UBSolution.clear();
    
    switch(solTag){
    case MibSLinkingPoolTagLowerIsInfeasible:
	{
	    if(isContainedInLinkingPool_ == false){
		linkingSolution.tag = solType;
		linkingSolution.numSolvedScenarios = 0;
		linkingSolution.UBObjValue = 0.0;
		linkingSolution.lowerObjValues.push_back(0);
	        model_->seenLinkingSolutions[linkSol] = linkingSolution;
	    }
	    else{
		model_->seenLinkingSolutions[linkSol].tag = MibSLinkingPoolTagLowerIsInfeasible;
		model_->seenLinkingSolutions[linkSol].numSolvedScenarios = 0;
		model_->seenLinkingSolutions[linkSol].lowerObjValues.clear();
		model_->seenLinkingSolutions[linkSol].lowerObjValues.push_back(0);
		model_->seenLinkingSolutions[linkSol].UBObjValue = 0.0;
	    }
	    break;
	}
    case MibSLinkingPoolTagLowerIsIncomplete:
	{
	    if(isContainedInLinkingPool_ == false){
		linkingSolution.tag = solType;
		linkingSolution.numSolvedScenarios = numSolvedScenarios;
		linkingSolution.UBObjValue = 0.0;
		linkingSolution.lowerObjValues = shouldStoreLowerObjs;
		model_->seenLinkingSolutions[linkSol] = linkingSolution;
	    }
	    else{
		model_->seenLinkingSolutions[linkSol].numSolvedScenarios = numSolvedScenarios;
		model_->seenLinkingSolutions[linkSol].UBObjValue = 0.0;
		model_->seenLinkingSolutions[linkSol].lowerObjValues.clear();
		model_->seenLinkingSolutions[linkSol].lowerObjValues.push_back(0);
		model_->seenLinkingSolutions[linkSol].lowerObjValues = shouldStoreLowerObjs;
	    }
	    break;
	}
    case MibSLinkingPoolTagLowerIsFeasible:
	{
	    if(isContainedInLinkingPool_ == false){
		linkingSolution.tag = solType;
		linkingSolution.numSolvedScenarios = numSolvedScenarios;
		linkingSolution.UBObjValue = 0.0;
		linkingSolution.lowerObjValues = shouldStoreLowerObjs;
		model_->seenLinkingSolutions[linkSol] = linkingSolution;
	    }
	    else{
		model_->seenLinkingSolutions[linkSol].tag = MibSLinkingPoolTagLowerIsFeasible;
		model_->seenLinkingSolutions[linkSol].numSolvedScenarios = numSolvedScenarios;
		model_->seenLinkingSolutions[linkSol].UBObjValue = 0.0;
		model_->seenLinkingSolutions[linkSol].lowerObjValues.clear();
		model_->seenLinkingSolutions[linkSol].lowerObjValues = shouldStoreLowerObjs;
	    }
	    break;
	}
    case MibSLinkingPoolTagUBIsSolved:
	{
	    model_->seenLinkingSolutions[linkSol].tag = MibSLinkingPoolTagUBIsSolved;
	    model_->seenLinkingSolutions[linkSol].numSolvedScenarios = numSolvedScenarios;
	    //if UB is infeasible, ubObjValue is infinity.
	    model_->seenLinkingSolutions[linkSol].UBObjValue = ubObjValue;
	    break;
	}
    }
}
	    
	
	
	    






    
 
