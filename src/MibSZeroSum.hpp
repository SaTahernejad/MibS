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

#ifndef MibSZeroSum_h_
#define MibSZeroSum_h_

#include "MibSModel.hpp"

//#############################################################################

class MibSZeroSum {

    friend class MibSModel;

private:

    MibSModel *model_;
    OsiSolverInterface *lSolver_;
    bool isAlgStarted_;
    bool isFirstPhaseFinished_;

public:

    MibSZeroSum();

    ~MibSZeroSum();

    /** Solve the zero-sum problem by using the available feasible solution **/
    void solveZeroSum(MibSModel *mibs, double *sol);

    /** First phase of the algorithm **/
    void doFirstPhase(double *LLSol, double *lColLb, double *lColUb,
		      double *optimalSol, bool &isUnbounded,
		      bool &foundOptimal);

    /** Find the value of big M **/
    double findBigM();

    /** Find the bounds on the lower-level variables **/
    double *findBoundsOnLLCols(double *newLb, double* newUb, double *LLSol);

    /** Set the initial part of lower-level solver **/
    void initialSetUpLowerSolver(CoinPackedMatrix *matrixA2);

    /** Modify the initial lower-level solver **/
    void modifyLowerSolver(double *sol);

    /** Finds the lower part of the solution **/
    void findLowerSol(double *sol, double *LLSol);

};

#endif
