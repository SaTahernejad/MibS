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

class MibSModel;

//#############################################################################

class MibSZeroSum {

    friend class MibSModel;
    friend class MibSBilevel;

private:

    MibSModel *model_;
    bool isAlgStarted_;
    bool isUnbounded_;
    bool foundOptimal_;
    bool returnedUpperBound_;
    double *optimalSol_;

public:

    MibSZeroSum();

    ~MibSZeroSum();

    /** Solve the zero-sum problem by using the available feasible solution **/
    void solveZeroSum(MibSModel *mibs, double *sol);

    /** First phase of the algorithm **/
    void doFirstPhase();

    /** Find the value of big M objective **/
    double findBigMObj();

    /** Find bound on artificial variable **/
    double findBoundArtCol();

};

#endif
