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

#ifndef MibSTreeNode_h_
#define MibSTreeNode_h_

#include "Alps.h"
#include "AlpsEncoded.h"

#include "BlisTreeNode.h"

//#############################################################################

class MibSTreeNode : public BlisTreeNode {

 private:
   
  double lowerUpperBound_;
  bool boundSet_;
  BlisLpStatus lpStatus_;
  bool useUBObj_;
  double *dual_;
  double *dj_;
  //FIXME: remove following hard coding of bounds after fixing leafBranchPath in WS
  double *lb_;
  double *ub_;
  double boundCutRhs_;
     
 public:
  
  MibSTreeNode();
  MibSTreeNode(AlpsNodeDesc *&desc);
  MibSTreeNode(BlisModel *m);

  ~MibSTreeNode();
  
  void setIsBoundSet(bool val) {boundSet_ = val;}
  void setLowerUB(double bound) {lowerUpperBound_ = bound;}
  void setLpStatus(BlisLpStatus lpStatus) {lpStatus_ = lpStatus;}
  void setDuals(double *dual) {dual_ = dual;}
  void setDjs(double *dj) {dj_ = dj;}
  void setLb(double *lb) {lb_ = lb;}
  void setUb(double *ub) {ub_ = ub;}
  void setBoundCutRhs(double rhs) {boundCutRhs_ = rhs;}
  inline bool isBoundSet() {return boundSet_;}
  inline double getLowerUB() {return lowerUpperBound_;}
  AlpsTreeNode* createNewTreeNode(AlpsNodeDesc *&desc) const;
  int process(bool isRoot, bool rampUp);

    
  BlisLpStatus getLpStatus() {return lpStatus_;}
  bool getUseUBObj() {return useUBObj_;}
  double * getDuals() {return dual_;}
  double * getDjs() {return dj_;}
  double * getLb() {return lb_;}
  double * getUb() {return ub_;}
  double getBoundCutRhs() {return boundCutRhs_;}

private:

    void setDualDj(BlisModel *model, double tol);
};

#endif
