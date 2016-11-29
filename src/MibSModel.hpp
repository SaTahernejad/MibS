/*===========================================================================*/
/* This file is part of a Mixed Integer Bilevel Solver                       */
/* developed using the Discrete Conic Optimization (DisCO).                  */
/*                                                                           */
/* Authors: Scott DeNegre, Lehigh University                                 */
/*          Ted Ralphs, Lehigh University                                    */
/*                                                                           */
/* Copyright (C) 2007-2015 Lehigh University, Scott DeNegre, and Ted Ralphs. */
/* All Rights Reserved.                                                      */
/*                                                                           */
/* This software is licensed under the Eclipse Public License. Please see    */
/* accompanying file for terms.                                              */
/*===========================================================================*/

#ifndef MibSModel_h_
#define MibSModel_h_

#include "OsiClpSolverInterface.hpp"

#include "DcoModel.hpp"
#include "DcoSolution.hpp"
#include "MibSBilevel.hpp"
#include "MibSParams.hpp"

class MibSBilevel;
class MibSCutGenerator;

//#############################################################################
class MibSModel : public DcoModel {

    friend class MibSCutGenerator;
    friend class MibSBilevel;
    friend class MibSBranchStrategyMaxInf;
    friend class MibSTreeNode;
    friend class MibSHeuristic;

private:

    /** Data file specifying LL problem **/
    std::string lowerDataFile_;

    /** Data file specifying Omega and UL objective **/
    std::string upperDataFile_;

    /** AMPL model file specifying Omega and UL objective **/
    /** (may or may not have data) **/
    std::string upperAmplModelFile_;

    /** AMPL data file specifying instance **/
    /** (may or may not have data) **/
    std::string upperAmplDataFile_;

    /** Original number of total variables **/
    int numOrigVars_;

    /** Original number of total constraints **/
    int numOrigCons_;

    /** Current number of total variables **/
    int numVars_;

    /** Current number of total constraints **/
    int numCons_;

    //sahar:I removed objSense_ because
    //it is a member of DcoModel class
    /** Objective sense of UL problem **/
    //double objSense_;

    /** Objective sense of LL problem **/
    double lowerObjSense_;

    /** Number of UL variables **/
    int upperColNum_;

    /** Number of UL constraints **/
    int upperRowNum_;

    /** Number of LL variables **/
    int lowerColNum_;

    /** Number of LL constraints **/
    int lowerRowNum_;

    /** Number of structural constraints **/
    int structRowNum_;

    /** Determines type of problem(general or interdiction) **/
    bool isInterdict_;

    /** Determines if problem is pure integer or not **/
    bool isPureInteger_;

    /** Determines if all UL coefficients are integer or not**/
    bool isUpperCoeffInt_;

    /** Determines if all LL coefficients are integer or not**/
    bool isLowerCoeffInt_;

    /** Determines if all variables of UL problem are binary or not **/
    bool allUpperBin_;

    /** Determines if all variables of LL problem are binary or not **/
    bool allLowerBin_;

    /** Determines if matrix A1 is positive or not **/
    bool positiveA1_;

    /** Determines if matrix A2 is positive or not **/
    bool positiveA2_;

    /** Determines if matrix G1 is positive or not **/
    bool positiveG1_;

    /** Determines if matrix G2 is positive or not **/
    bool positiveG2_;

    /** the left (negative) slope of the LL value function **/
    double leftSlope_;

    /** the right (positive) slope of the LL value function **/
    double rightSlope_;

    //sahar:I added this memebr and member function setColType
    //here because BlisModel has this member, but DcoModel
    //does not have it.
    /** Type of variables: continuous, integer, boolean **/
    char *colType_;

    /** Indices of UL variables **/
    int * upperColInd_;

    /** Indices of UL rows **/
    int * upperRowInd_;

    /** Indices of LL variables **/
    int * lowerColInd_;

    /** Indices of LL rows **/
    int * lowerRowInd_;

    /** Indices of structural (non-vub) rows **/
    int * structRowInd_;

    /** Indices of UL variables in LL constraints **/
    int * fixedInd_;

    /** LL objective coefficients **/
    double * lowerObjCoeffs_;

    /** Interdiction coefficients **/
    double * interdictCost_;

    /** Interdiction budget **/
    double interdictBudget_;

    /** Original column lower bounds from Omega **/
    double * origColLb_;

    /** Original column upper bounds from Omega **/
    double * origColUb_;

    /** Original row lower bounds from Omega **/
    double * origRowLb_;

    /** Original row upper bounds from Omega **/
    double * origRowUb_;

    /** MibSBilevel object **/
    MibSBilevel *bS_;

    /** Method for determining binding cons **/
    std::string bindingMethod_;

    /** MibS Parameters **/
    MibSParams *MibSPar_;

    /** Indicator telling whether solution has been updated **/
    bool solIsUpdated_;

public:

    MibSModel();
    ~MibSModel();

    /** Read in the problem data **/
    void readInstance(const char * dataFile);

    /** Set the UL file **/
    inline void setUpperFile(std::string infile)
    {
	upperDataFile_ = infile;
    }

    /** Set the UL AMPL model file **/
    inline void setUpperAmplModelFile(std::string infile)
    {
	upperAmplModelFile_ = infile;
    }

    /** Set the UL AMPL data file **/
    inline void setUpperAmplDataFile(std::string infile)
    {
	upperAmplDataFile_ = infile;
    }

    /** Set the LL file **/
    inline void setLowerFile(std::string infile)
    {
	lowerDataFile_ = infile;
    }

    /** Set the MibsBilevel pointer **/
    inline void setMibSBilevel(MibSBilevel *bs) {bS_ = bs;}

    /** Set the number of rows **/
    inline void setNumRows(int val) {numRows_ = val;}

    /** Set the number of columns **/
    inline void setNumCols(int val) {numCols_ = val;}

    /** Set the LL dimension **/
    inline void setLowerDim(int val) {lowerColNum_ = val;}

    /** Set the UL dimension **/
    inline void setUpperDim(int val) {upperColNum_ = val;}

    /** Set the UL row number **/
    inline void setUpperRowNum(int val) {upperRowNum_ = val;}

    /** Set the LL row number **/
    inline void setLowerRowNum(int val) {lowerRowNum_ = val;}

    /** Set the number of structural rows **/
    inline void setStructRowNum(int val) {structRowNum_ = val;}

    /** Pass variable types */
    void setColType(char *colType){
	colType_ = colType;
    }

    /** Set the interdiction cost **/
    inline void setInterdictCost(double *ptr) {interdictCost_ = ptr;}

    /** Set the interdiction budget **/
    inline void setInterdictBudget(double val) {interdictBudget_ = val;}

    /** Set UL column indices **/
    void setUpperColInd(int *ptr) {upperColInd_ = ptr;}

    /** Set UL column data **/
    void setUpperColData();

    /** Set UL row indices **/
    void setUpperRowInd(int *ptr) {upperRowInd_ = ptr;}

    /** Set UL row indices **/
    void setUpperRowData();

    /** Set pointer to array of LL column indices **/
    void setLowerColInd(int *ptr) {lowerColInd_ = ptr;}

    /** Set pointer to array of LL row indices **/
    void setLowerRowInd(int *ptr) {lowerRowInd_ = ptr;}

    /** Set pointer to array of structural row indices **/
    void setStructRowInd(int *ptr) {structRowInd_ = ptr;}

    /** Set pointer to array of LL objective coefficients **/
    void setLowerObjCoeffs(double *ptr) {lowerObjCoeffs_ = ptr;}

    /** Set objective sense of LL problem **/
    void setLowerObjSense(double os) {lowerObjSense_ = os;}

    /** Set the number of original variables **/
    void setNumOrigCols(int num) {numOrigCols_ = num;}

    /** Set the number of original constraints **/
    void setNumOrigRows(int num) {numOrigRows_ = num;}

    /** set the slopes of the LL value function **/
    void setValFuncSlopes();

    /** Get the UL file **/
    std::string getUpperFile() {return upperDataFile_;}

    /** Get the UL AMPL model file **/
    std::string getUpperAmplModelFile() {return upperAmplModelFile_;}

    /** Get the UL AMPL data file **/
    std::string getUpperAmplDataFile() {return upperAmplDataFile_;}

    /** Get the LL file **/
    std::string getLowerFile()
    {
	return MibSPar_->entry(MibSParams::auxiliaryInfoFile);
    }

    /** Get the LL dimension **/
    int getLowerDim() {return lowerColNum_;}

    /** Get the UL dimension **/
    int getUpperDim() {return upperColNum_;}

    /** Get the number of original variables **/
    int getNumOrigCols() {return numOrigCols_;}

    /** Get the number of original constraints **/
    int getNumOrigRows() {return numOrigRows_;}

    /** Get the UL row number **/
    int getUpperRowNum() {return upperRowNum_;}

    /** Get the LL row number **/
    int getLowerRowNum() {return lowerRowNum_;}

    /** Get objective sense of LL problem **/
    double getLowerObjSense() {return lowerObjSense_;}

    /** Get the tolerance **/
    double getTolerance() {return etol_;}

    /** Get pointer to the UL column index array **/
    int * getUpperColInd() {return upperColInd_;}

    /** Get pointer to the UL row index array **/
    int * getUpperRowInd() {return upperRowInd_;}

    /** Get pointer to the LL column index array **/
    int * getLowerColInd() {return lowerColInd_;}

    /** Get pointer to the LL row index array **/
    int * getLowerRowInd() {return lowerRowInd_;}

    /** Get pointer to the UL columns in LL problem array **/
    int * getFixedInd() {return fixedInd_;}

    /** Get pointer to the array of original column lower bounds **/
    double * getOrigColLb() const {return origColLb_;}

    /** Get pointer to the array of original column upper bounds **/
    double * getOrigColUb() const {return origColUb_;}

    /** Get pointer to the array of original row lower bounds **/
    double * getOrigRowLb() const {return origRowLb_;}

    /** Get pointer to the array of original row upper bounds **/
    double * getOrigRowUb() const {return origRowUb_;}

    /** Get pointer to the LL objective coefficient array **/
    double * getLowerObjCoeffs() {return lowerObjCoeffs_;}

    /** Get pointer to the interdiction coefficient array **/
    double * getInterdictCost() {return interdictCost_;}

    /** Get the interdiction budget **/
    double getInterdictBudget() {return interdictBudget_;}

    /** Get the pointer to MibsBilevel **/
    inline MibSBilevel *getMibSBilevel() {return bS_;}

    /** Get the parameters **/
    MibSParams *MibSPar() {return MibSPar_;}

    /** Set the DisCO parameters **/
    void setDisCOParameters();

    /** Read auxiliary data file **/
    void readAuxiliaryData();

    /** Set auxiliary data directly when using MibS as a library **/
    void loadAuxiliaryData(int lowerColNum, int lowerRowNum,
			   const int *lowerColInd,
			   const int *lowerRowInd,
			   double lowerObjSense,
			   const double *lowerObjCoef,
			   int upperColNum, int upperRowNum,
			   const int *upperColInd,
			   const int *upperRowInd,
			   int structRowNum,
			   const int *structRowInd,
			   double interdictBudget,
			   const double *interdictCost);

    /** Read problem description file **/
    void readProblemData();

    /** Set problem data directly when using MibS as a library **/
    void loadProblemData(const CoinPackedMatrix& matrix,
			 const double* colLB, const double* colUB,
			 const double* obj, const double* rowLB,
			 const double* rowUB, const char *types,
			 double objSense, double infinity,  const char *rowSense);

    /** Check for solution feasiblity **/
    DcoSolution * userFeasibleSolution(const double * solution,
					bool &userFeasible);

    /** Check if a LL solution satisfies UL constraints **/
    bool checkUpperFeasibility(double * solution);

    /** Get solution information **/
    CoinPackedVector * getSolution();

    /** Calls MibSBilevel::createBilevel(CoinPackedVector *vec) **/
    void createBilevel(CoinPackedVector *vec);

    /** Print current solution **/
    void printCurSol();

    /** Read in Alps, DisCO, MibS parameters. */
    virtual void readParameters(const int argnum, const char * const *arglist);

    /** Pack MibS portion of the model into an encoded object. */
    AlpsReturnStatus encodeMibS(AlpsEncoded *encoded) const;

    /** Unpack MibS portion of the model from an encoded object. */
    AlpsReturnStatus decodeMibS(AlpsEncoded &encoded);

    /** The method that encodes the model into an encoded object. */
    virtual AlpsEncoded* encode() const;

    /** The method that decodes the model from an encoded object. */
    virtual void decodeToSelf(AlpsEncoded&);

    /** Determine the list of first-stage variables participate in second-stage constraints */
    void setRequiredFixedList(const CoinPackedMatrix *newMatrix);

    /** Determines the properties of instance. */
    void instanceStructure(const CoinPackedMatrix *newMatrix, const double* rowLB,
			   const double* rowUB, const char *rowSense);

    AlpsTreeNode * createRoot();

    virtual bool setupSelf();

    void setBounds();

    void runPreprocessor();

    void runPreprocessor1();

    double getObjectiveBound();

    double lowerObjectiveBound();

    double interdictionBound();

private:

    /** Initialize the object data **/
    void initialize();

    bool findIndex(int index, int size, int * indices);

    int binarySearch(int start, int stop, int index, int * indexArray);

    OsiSolverInterface * setUpModel(int * fixed);
};

#endif

    

    
    

    

    
    

    
