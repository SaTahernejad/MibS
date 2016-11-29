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

#include "DcoModel.hpp"
#include "DcoLinearConstraint.hpp"
#include "DcoVariable.hpp"
#include "BlisNodeDesc.h"
#include "BlisConfig.h"

#include "BcpsConfig.h"

//#############################################################################
MibSModel::MibSModel()
{
    initialize();
}

//#############################################################################
MibSModel::~MibSModel()
{
    if(colType_) delete [] colType_;
    if(upperColInd_) delete [] upperColInd_;
    if(lowerColInd_) delete [] lowerColInd_;
    if(upperRowInd_) delete [] upperRowInd_;
    if(lowerRowInd_) delete [] lowerRowInd_;
    if(structRowInd_) delete [] structRowInd_;
    if(fixedInd_) delete [] fixedInd_;
    if(interdictCost_) delete [] interdictCost_;
    if(origColLb_) delete [] origColLb_;
    if(origColUb_) delete [] origColUb_;
    if(origRowLb_) delete [] origRowLb_;
    if(origRowUb_) delete [] origRowUb_;
    if(lowerObjCoeffs_) delete [] lowerObjCoeffs_;
    if(MibSPar_) delete MibSPar_;
    if(bS_) delete bS_;
}

//#############################################################################   
void
MibSModel::initialize()
{
    lowerDataFile_ = "";
    upperDataFile_ = "";
    upperAmplModelFile_ = "";
    upperAmplDataFile_ = "";
    etol_ = 1e-5;
    numOrigVars_ = 0;
    numOrigCons_ = 0;
    numVars_ = 0;
    numCons_ = 0;
    objSense_ = 0.0;
    lowerObjSense_ = 0.0;
    upperColNum_ = 0;
    upperRowNum_ = 0;
    lowerColNum_ = 0;
    lowerRowNum_ = 0;
    structRowNum_ = 0;
    leftSlope_ = 0.0;
    rightSlope_ = 0.0;
    isInterdict_ = false;
    isPureInteger_ = true;
    isUpperCoeffInt_ = true;
    isLowerCoeffInt_ = true;
    allUpperBin_ = true;
    allLowerBin_ = true;
    positiveA1_ = true;
    positiveA2_ = true;
    positiveG1_ = true;
    positiveG2_ = true;
    colType_ = NULL;
    upperColInd_ = NULL;
    lowerColInd_ = NULL;
    upperRowInd_ = NULL;
    lowerRowInd_ = NULL;
    structRowInd_ = NULL;
    fixedInd_ = NULL;
    origColLb_ = NULL;
    origColUb_ = NULL;
    origRowLb_ = NULL;
    origRowUb_ = NULL;
    lowerObjCoeffs_ = NULL;
    interdictCost_ = NULL;
    //To Do:sahar: check the below lines
    bS_ = new MibSBilevel();
    MibSPar_ = new MibSParams;
    solIsUpdated_ = false;
    MibSCutGenerator *cg = new MibSCutGenerator(this);
    cg->setStrategy(DcoCutStrategyPeriodic);
    cg->setCutGenerationFreq(1);  // Generate cuts at every node
    addCutGenerator(cg);
    setDcoParameters();
}

//#############################################################################
void
MibSModel::setDcoParameters()
{
    int dcoCuts(MibSPar_->entry(MibSParams::dcoCutStrategy));
    int dcoBranch(MibSPar_->entry(MibSParams::dcoBranchStrategy));

    /* Set Blis Parameters to keep cutting until no cut is found */
    dcoPar()->setEntry(DcoParams::cutFactor, ALPS_DBL_MAX);
    dcoPar()->setEntry(DcoParams::cutPass, ALPS_INT_MAX);
    dcoPar()->setEntry(DcoParams::tailOff, -10000);
    dcoPar()->setEntry(DcoParams::denseConFactor, ALPS_DBL_MAX);

    //sahar:I removed the line for setting cutgenerationfrequency to 1
    //because its default value in Dco is 1.

    /* Set Dco cut strategy using MibS parameters dcoCutStrategy */
    dcoPar()->setEntry(DcoParams::cutStrategy, dcoCuts);
    /* Set Blis branch strategy using MibS parameters blisBranchStrategy */
    dcoPar()->setEntry(DcoParams::branchStrategy, dcoBranch);
}

//#############################################################################  
void
MibSModel::readMibSInstance(const char* dataFile)
{
    readInstance(dataFile);
    setUpperFile(dataFile);
    readAuxiliaryData(); // reads in lower-level vars, rows, obj coeffs
    readProblemData(); // reads in max c^1x + d^1y s.t. (x,y) in Omega^I
}

//#############################################################################
void
MibSModel::readProblemData()
{
    CoinMpsIO *mps = new CoinMpsIO;

    int i(0);

    // mps file reader
    mps->readMps(getUpperFile().c_str(), "");

    mps->messageHandler()->setLogLevel(msgLevel);

    CoinPackedMatrix matrix = *(mps->getMatrixByCol());

    //Set colType_
    char *colType (NULL);
    colType = new char [numCols_];
    
    for (i = 0; i < numCols_; ++i) {
	if (mps->isContinuous(i)) {
	    colType[i] = 'C';
	}
	else{
	    if (colLB_[i] == 0 && colUB_[i] == 1.0) {
		colType[i] = 'B';
	    }
	    else {
		colType[i] = 'I';
	    }
	}
    }

    setColType(colType);

    CoinPackedMatrix matrix = *(mps->getMatrixByCol());

    const char* rowSense = mps->getRowSense();

    loadProblemData(matrix, colLB_, colUB_, objCoef_, rowLB_, rowUB_,
		    colType, objSense_, mps->getInfinity()); 

    delete [] colType;
    delete mps;
}

//############################################################################# 
void
MibSModel::loadProblemData(const CoinPackedMatrix& matrix,
			   const double* colLB, const double* colUB,
			   const double* obj, const double* rowLB,
			   const double* rowUB, const char *types,
			   double objSense, double infinity, const char *rowSense)
{
    int i(0), j(0);
    
    int problemType(MibSPar_->entry(MibSParams::bilevelProblemType));

    if(isInterdict_ == true){
	if(problemType == PARAM_NOTSET){
	    MibSPar()->setEntry(MibSParams::bilevelProblemType, INTERDICT);
	}
	else if(problemType == GENERAL){
	    std::cout << "Wrong value for MibSProblemType. Automatically modifying its value." << std::endl;
	    MibSPar()->setEntry(MibSParams::bilevelProblemType, INTERDICT);
	}
    }
    else{
	if(problemType == PARAM_NOTSET){
	    MibSPar()->setEntry(MibSParams::bilevelProblemType, GENERAL);
	}
	else if(problemType == INTERDICT){
	    std::cout << "Wrong value for MibSProblemType. Automatically modifying its value." << std::endl;
	    MibSPar()->setEntry(MibSParams::bilevelProblemType, GENERAL);
	}
    }

    problemType = MibSPar_->entry(MibSParams::bilevelProblemType);

    int rowNum(0), colNum(0);

    double *varLB(NULL);
    double *varUB(NULL);
    double *conLB(NULL);
    double *conUB(NULL);
    double *objCoef(NULL);
    char   *colType(NULL);

    CoinPackedMatrix *newMatrix = NULL;

    switch (problemType) {

    case GENERAL:
	
	rowNum = matrix.getNumRows();
	colNum = matrix.getNumCols();
	
	//set structural rows
	rowNum = matrix.getNumRows();
	colNum = matrix.getNumCols();
	if (!structRowInd_) {
	    structRowInd_ = new int[rowNum];
	    CoinIotaN(structRowInd_, rowNum, 0);
	    structRowNum_ = rowNum;
	}

	//Make copies of the data
	newMatrix = new CoinPackedMatrix();
	*newMatrix = matrix;

	varLB = new double [colNum];
	varUB = new double [colNum];
	conLB = new double [rowNum];
	conUB = new double [rowNum];
	objCoef = new double [colNum];
	colType = new char [colNum];

	CoinDisjointCopyN(colLB, colNum, varLB);
	CoinDisjointCopyN(colUB, colNum, varUB);
	CoinDisjointCopyN(rowLB, rowNum, conLB);
	CoinDisjointCopyN(rowUB, rowNum, conUB);
	CoinDisjointCopyN(obj, colNum, objCoef);
	memcpy(colType, types, colNum);

	break;

    case INTERDICT:
	//sahar: in interdiction problems, we have
	//three sets of constraints
	//1. UL rows -> mainUpperRow
	//2. LL rows which consists
	//2.1. main LL rows -> mainLowerRow
	//2.2. interdiction LL rows -> interdictRow
	//we have three sets of variables
	//1. ULCols
	//2. LLCOls (the number of UL and LL cols is the same).
	//2.auxCols

	//sahar:The objective function in MPS file shows
	//the LL objective.
	
	//sahar:FIXME: Allow more than one  UL Row
	//Add interdict vars and aux UL rows
	int upperColNum(matrix.getNumCols());
	int lowerColNum(matrix.getNumCols());
	double *intCosts = getInterdictCost();
	int numInterdictNZ(0);
	int mainUpperRowNum(1);
	int mainLowerRowNum(matrix.getNumRows());
	int interdictRowNum(upperColNum);
	int auxRowNum(mainUpperRowNum + interdictRowNum);
	//sahar:why do we add aux col?
	int auxColNum(1);

	//set total number of rows and variables
	rowNum = mainUpperRowNum + mainLowerRowNum + interdictRowNum;
	colNum = upperColNum + lowerColNum + auxColNum;

	//set the number of non-zero coeffs in UL row
	for (i = 0; i < upperColNum; i++) {
	    if (fabs(intCosts[i]) > etol_) {
		numInterdictNZ++;
	    }
	}

	//set structural rows: it consists only main UL and LL rows
	int structRowNum(mainUpperRowNum + mainLowerRowNum);
	structRowInd_ = new int[structRowNum];
	CoinIotaN(structRowInd_, structRowNum, 0);
	structRowNum_ = structRowNum;

	varLB = new double [colNum];
	varUB = new double [colNum];

	conLB = new double [rowNum];
	conUB = new double [rowNum];

	CoinDisjointCopyN(colLB, upperColNum, varLB + upperColNum);
	CoinDisjointCopyN(colUB, upperColNum, varUB + upperColNum);

	CoinFillN(varLB, upperColNum, 0.0);
	CoinFillN(varUB, upperColNum, 1.0);

	CoinFillN(varLB + (upperColNum + lowerColNum), auxColNum, 0.0);
	CoinFillN(varUB + (upperColNum + lowerColNum), auxColNum, 1.0);

	CoinDisjointCopyN(rowLB, mainLowerRowNum, conLB + mainUpperRowNum);
	CoinDisjointCopyN(rowUB, mainLowerRowNum, conUB + mainUpperRowNum);

	//add UL row (interdiction budget row) 
	CoinFillN(conLB, mainUpperRowNum, - 1 * infinity);
	CoinFillN(conUB, mainUpperRowNum, getInterdictBudget());

	//add interdiction rows
	CoinFillN(conLB + (mainUpperRowNum + mainLowerRowNum),
		  interdictRowNum, - 1 * infinity);
	CoinDisjointCopyN(colUB, interdictRowNum, conUB + (mainUpperRowNum + mainLowerRowNum));

	//set UL objective
	objCoef = new double [colNum];
	
	CoinZeroN(objCoef, colNum);
	
	for (j = 0; j < colNum; j++) {
	    objCoef[j + colNum] = -1 * obj[j];
	}

	//Set colType_   
	colType = new char [colNum];

	for (j = 0; j < upperColNum; ++j) {
	    colType[j] = 'B';
	}

	CoinDisjointCopyN(types, lowerColNum, colType + upperColNum);

	//auxilliary indicator columns, used later
	for(j = 0; j < auxColNum; j++) {
	    colType[j + upperColNum + lowerColNum] = 'B';
	}

	//set the matrix
	CoinPackedMatrix rowMatrix;
	rowMatrix = matrix;
	rowMatrix.reverseOrdering();
	const double * matElements = rowMatrix.getElements();
	const int * matIndices = rowMatrix.getIndices();
	const int * matStarts = rowMatrix.getVectorStarts();

	newMatrix = new CoinPackedMatrix(false, 0, 0);
	newMatrix->setDimensions(0, colNum);
	int start(0), end(0), tmp(0), index(0);

	CoinPackedVector row;
	//add UL row to matrix
	for (i = 0; i < mainUpperRowNum; i++) {
	    for (j = 0; j < upperColNum; j++) {
		row.insert(j, intCosts[j]);
	    }
	    newMatrix->appendRow(row);
	    row.clear();
	}

	//add main LL rows to matrix
	for (i = 0; i < mainLowerRowNum; i++) {
	    start = matStarts[i];
	    end = start + rowMatrix.getVectorSize(i);
	    for (j = start; j < end; j++) {
		index = matIndices[j] + upperColNum;
		row.insert(index, matElements[j]);
	    }
	    newMatrix->appendRow(row);
	    row.clear();
	}

	//add interdiction rows
	for (i = 0; i < interdictRowNum; i++) {
	    row.insert(i, colUB[i]);
	    row.insert(i + upperColNum, 1.0);
	    newMatrix->appendRow(row);
	    row.clear();
	}

	newMatrix->reverseOrdering();

	//sahar: in previous version of MibS, the information for
	//UL problem were set here, but I removed them because
	//at the end of this function, we call setUpperColData()
	//and setUpperRowData().
	
	break;
    }

    setColMatrix(newMatrix);

    
    //sahar: do we need numORigVars_
    // and numOrigCons_? I think
    //we can use numCoreVariables_
    //and numCoreConstraints_ instead.
    //even if they are required, I think
    //the way that they are defined in
    //previous version is not correct.
    //in previous version, these variables
    //are defined in memver functions
    //setUpperColData() and setUpperRowData();
    numCoreConstraints_ = numCols_ = numCons_
	= newMatrix->getMinorDim();
    
    numCoreVariables_   = numRows_ = numLinearRows_ =
	numVars_ = newMatrix->getMajorDim();

    rowLB_ = conLB;
    rowUB_ = conUB;
   
    colLB_ = varLB;
    colUB_ = varUB;

    objCoef_ = objCoef;

    setColType(colType);

    dcoPar_->setEntry(DcoParams::objSense, objSense);

    //create variables and constraints
    //sahar:In previous version, the codes for creating
    //variables and constraints in BlisModel were rewritten
    //here, I used the member functions of DcoModel
    setupAddVariables();
    setupAddLinearConstraints();

    //set the information of UL problem
    setUpperColData();
    setUpperRowData();

    //sahar:I think we can remove member function setBounds()
    setBounds();

    //sahar: I removed the member function setProblemType()

    setRequiredFixedList(newMatrix);
    instanceStructure(newMatrix, rowLB, rowUB, rowSense);
}


//#############################################################################  
void
MibSModel::setUpperColData()
{
    int lowerColNum(lowerColNum_);
    

    

    

    


	
	
	

	
	
	
     
	
	
	
    



    
    
    

    
    
    
    
    
    
    
