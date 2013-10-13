/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2013   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include "staticfracture.h"
#include "fracturemanager.h"
#include "dofmanager.h"
#include "enrichmentdomain.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"
#include "isolinearelasticmaterial.h"
#include "exportmodulemanager.h"

#include "oofemtxtdatareader.h" // for reading .in files
#include "util.h"               // for creating eng models
#include <vector>

namespace oofem {

REGISTER_EngngModel( StaticFracture );

StaticFracture :: StaticFracture(int i, EngngModel *_master) : NonLinearStatic(i, _master)
{
    ndomains = 1;
}


IRResultType
StaticFracture :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, this->numberOfSteps, _IFT_OptimizationProblem_nsteps);
    if ( this->numberOfSteps <= 0 ) {
        _error("instanciateFrom: nsteps not specified, bad format");
    }

    //if ( ir->hasField(_IFT_StaggeredProblem_deltat) ) {
    //    EngngModel :: initializeFrom(ir);
    //    IR_GIVE_FIELD(ir, deltaT, _IFT_StaggeredProblem_deltat);
    //    dtTimeFunction = 0;
    //} else if ( ir->hasField(_IFT_StaggeredProblem_prescribedtimes) ) {
    //    EngngModel :: initializeFrom(ir);
    //    IR_GIVE_FIELD(ir, discreteTimes, _IFT_StaggeredProblem_prescribedtimes);
    //    dtTimeFunction = 0;
    //} else {
    //    IR_GIVE_FIELD(ir, timeDefinedByProb, _IFT_StaggeredProblem_timeDefinedByProb);
    //}
    //if ( dtTimeFunction < 1 ) {
    //    ndomains = 0;
    //}
    //IR_GIVE_OPTIONAL_FIELD(ir, dtTimeFunction, _IFT_StaggeredProblem_dtf);
    //IR_GIVE_OPTIONAL_FIELD(ir, stepMultiplier, _IFT_StaggeredProblem_stepmultiplier);
    //if ( stepMultiplier < 0 ) {
    //    _error("stepMultiplier must be > 0")
    //}
  
    
    IR_GIVE_FIELD(ir, this->numObjFunc, _IFT__IFT_OptimizationProblem_NumObjFunc);
    emodelList = new AList< EngngModel >(1);
    inputStreamNames = new std :: string [ 1 ];
    //for (int i = 1; i <= this->numObjFunc; i++) {
        //IR_GIVE_FIELD(ir, inputStreamNames [ i-1 ], _IFT__IFT_OptimizationProblem_Name_prob1);
    //}
    
    IR_GIVE_FIELD(ir, inputStreamNames [ 0 ], _IFT__IFT_OptimizationProblem_Name_prob1);
//    IR_GIVE_FIELD(ir, inputStreamNames [ 1 ], "top.in");


    this->objFuncList.resize(1);
    this->objFuncList[0] = new MinCompliance();



    return IRRT_OK;
}

///////////
int
StaticFracture :: instanciateYourself(DataReader *dr, InputRecord *ir, const char *dataOutputFileName, const char *desc)
{

    this->Instanciate_init(dataOutputFileName, this->ndomains);



    fprintf(outputStream, "%s", PRG_HEADER);
    this->startTime = time(NULL);
    ////this->startClock= this-> getClock();
    fprintf( outputStream, "\nStarting analysis on: %s\n", ctime(& this->startTime) );
    fprintf(outputStream, "%s\n", desc);


    // instanciate receiver
    this->initializeFrom(ir);
    this->instanciateDefaultMetaStep(ir); // after numSteps has been set

    //exportModuleManager->initializeFrom(ir);
    //initModuleManager->initializeFrom(ir);

    int result;
    //result = EngngModel :: instanciateYourself(dr, ir, dataOutputFileName, desc);
    ir->finish();
    // instanciate slave problems
    result &= this->instanciateSlaveProblems();



    // should loop over slave problems and obj. fncs.
    EngngModel *sp = emodelList->at(1);
    Domain *d = sp->giveDomain(1);   
    int numMat = d->giveNumberOfMaterialModels();

    MinCompliance *objFunc = dynamic_cast< MinCompliance * >( this->objFuncList[0] ); 


    // Assumed that there is only one design parameter per element
    objFunc->designVarList.resize(numMat);
    objFunc->sensitivityList.resize(numMat);
    for (int i = 1; i <= numMat; i++) {
        objFunc->designVarList.at(i) = objFunc->volFrac;
        objFunc->sensitivityList.at(i) = 0.0;
    }



    return result;
}


//int EngngModel :: instanciateYourself(DataReader *dr, InputRecord *ir, const char *dataOutputFileName, const char *desc)
//// simple input - only number of steps variable is read
//{
//    bool inputReaderFinish = true;
//
//    this->Instanciate_init(dataOutputFileName, this->ndomains);
//
//    fprintf(outputStream, "%s", PRG_HEADER);
//    this->startTime = time(NULL);
//    //this->startClock= this-> getClock();
//    fprintf( outputStream, "\nStarting analysis on: %s\n", ctime(& this->startTime) );
//
//    fprintf(outputStream, "%s\n", desc);
//
//#  ifdef VERBOSE
//    OOFEM_LOG_DEBUG( "Reading all data from input file %s\n", dr->giveDataSourceName() );
//#  endif
//#ifdef __PARALLEL_MODE
//    if ( this->isParallel() ) {
//        fprintf(outputStream, "Problem rank is %d/%d on %s\n\n", this->rank, this->numProcs, this->processor_name);
//    }
//
//#endif
//
//    // instanciate receiver
//    this->initializeFrom(ir);
//    exportModuleManager->initializeFrom(ir);
//    initModuleManager->initializeFrom(ir);
//
//    if ( this->nMetaSteps == 0 ) {
//        inputReaderFinish = false;
//        this->instanciateDefaultMetaStep(ir);
//    } else {
//        this->instanciateMetaSteps(dr);
//    }
//
//    // instanciate initialization module manager
//    initModuleManager->instanciateYourself(dr, ir);
//    // instanciate export module manager
//    exportModuleManager->instanciateYourself(dr, ir);
//    this->instanciateDomains(dr);
//
//    exportModuleManager->initialize();
//
//    // Milan ??????????????????
//    //GPImportModule* gim = new GPImportModule(this);
//    //gim -> getInput();
//    // Milan ??????????????????
//
//#endif
//    // check emodel input record if no default metastep, since all has been read
//    if ( inputReaderFinish ) {
//        ir->finish();
//    }
//
//    return 1;
//}

int
StaticFracture :: instanciateSlaveProblems()
{
    EngngModel *timeDefProb = NULL, *slaveProb;

    int i=1;
    this->nModels = 1;
    OOFEMTXTDataReader dr( inputStreamNames [ i - 1 ].c_str() );
    //the slave problem dictating time needs to have attribute master=NULL, other problems point to the dictating slave
    //slaveProb = oofem :: InstanciateProblem(& dr, this->pMode, this->contextOutputMode, this);
    slaveProb = oofem :: InstanciateProblem(& dr, this->pMode, this->contextOutputMode, NULL);
    emodelList->put(i, slaveProb);
   
    return 1;
}

void
StaticFracture :: solveYourself()
{
    MetaStep *activeMStep;
    FILE *out = this->giveOutputStream();
    this->timer.startTimer(EngngModelTimer :: EMTT_AnalysisTimer);
    this->giveNumberOfSlaveProblems();
    
    this->timer.startTimer(EngngModelTimer :: EMTT_SolutionStepTimer);
    this->timer.initTimer(EngngModelTimer :: EMTT_NetComputationalStepTimer);


    int numMetaSteps = this->giveNumberOfMetaSteps();
    for (int imstep = 1; imstep <= numMetaSteps; imstep++) {
        activeMStep = this->giveMetaStep(imstep);

        //this->initMetaStepAttributes(activeMStep); // give numerical method and read other input

        int nTimeSteps = activeMStep->giveNumberOfSteps();
        for ( int tStepNum = 1; tStepNum <= nTimeSteps; tStepNum++ ) { //loop over time steps in opt analysis

            EngngModel *sp = this->giveSlaveProblem(1); 
           
            //sp->initMetaStepAttributes( sp->giveMetaStep(imstep) );
            
            // Resetting the time step number for each sp after each optimization time step

            
            sp->solveYourself();

            this->updateYourself( this->giveCurrentStep());
            

            //this->terminate( this->giveCurrentStep() );
        
            //sp->giveExportModuleManager()->doOutput( this->giveCurrentStep() );
            
            // optimization
            this->optimize( this->giveCurrentStep() );    


            

            if ( tStepNum > 0) {
                TimeStep *tStep =  sp->giveCurrentStep();

                tStep->setNumber(tStepNum); 
                sp->giveExportModuleManager()->doOutput(tStep);


                tStep->setNumber(0); 
            }

        }
    }
    

}


void
StaticFracture :: solveYourselfAt(TimeStep *tStep)
{
    for ( int subProb = 1; subProb <= this->giveNumberOfSlaveProblems(); subProb++ ) {
        EngngModel *sp = this->giveSlaveProblem(subProb);

        //sp->giveNextStep();
        sp->solveYourself();
        //sp->solveYourselfAt( sp->giveCurrentStep() );
        
    }   
}


EngngModel *
StaticFracture :: giveSlaveProblem(int i)
{
    if ( ( i > 0 ) && ( i <= this->nModels ) ) {
        return this->emodelList->at(i);
    } else {
        _error("giveSlaveProblem: Undefined problem");
    }

    return NULL;
}

void
StaticFracture :: optimize(TimeStep *tStep)
{
    // Main optimization loop

    MinCompliance *objFunc = dynamic_cast< MinCompliance * >( this->objFuncList[0] ); 

    for ( int subProb = 1; subProb <= this->giveNumberOfSlaveProblems(); subProb++ ) {
        
        EngngModel *sp = this->giveSlaveProblem(subProb);
        Domain *d = sp->giveDomain(1);   
        int numMat = d->giveNumberOfMaterialModels();


        double cost = 0.0;
        int numEl = d->giveNumberOfElements(); 

        double dce = 0.0;
        for (int i = 1; i <= numEl; i++) {
            
            Element *el = d->giveElement(i);
            cost += objFunc->evaluateYourself(el, dce, sp->giveCurrentStep() ); // add cost for each element
            objFunc->sensitivityList.at(i) = dce; // save derivative of cost function

        }
        
        // Update design variables based on some method. For now use the 'standard' optimality
        // criteria
        this->optimalityCriteria(objFunc);

        
        objFunc->designVarList.printYourself();

        // Update material parameters
        for (int i = 1; i <= numMat; i++) {
            DynamicInputRecord ir;
            Material *mat = d->giveMaterial(i);
            mat->giveInputRecord(ir);
            double E0 = 1.0;
            double fac = pow( objFunc->designVarList.at(i), objFunc->penalty);
            ir.setField( E0 * fac, _IFT_IsotropicLinearElasticMaterial_e);
            mat->initializeFrom(&ir);
            if (i==1)
            {
                mat->printYourself();
            }
        }

        //objFunc->designVarList.printYourself();
    
        printf("\n costfunction %e and sum design %e \n \n", cost, objFunc->designVarList.sum());

    }
}


void 
StaticFracture :: updateYourself(TimeStep *tStep)
{
    for ( int subProb = 1; subProb <= this->giveNumberOfSlaveProblems(); subProb++ ) {    
        EngngModel *sp = this->giveSlaveProblem(subProb);
        //sp->giveCurrentStep()
//        sp->updateYourself(sp->giveCurrentStep());
    }
}

double
MinCompliance :: evaluateYourself(Element *el, double &dc, TimeStep *tStep)
{
    FloatArray help, ae;
    FloatMatrix Ke;
    el->giveCharacteristicMatrix(Ke, ElasticStiffnessMatrix, tStep);
    double fac = pow( this->designVarList.at(el->giveNumber()), this->penalty);
    el->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, ae);
    help.beProductOf(Ke,ae);
    
    double cost = ae.dotProduct(help);
    dc = -this->penalty * pow( designVarList.at(el->giveNumber()), this->penalty-1.0) * cost/fac;
    return cost;
}





void 
StaticFracture :: optimalityCriteria(ObjectiveFunction *objFunc)
{
    FloatArray &x = objFunc->designVarList; 
    FloatArray &dc = objFunc->sensitivityList;

    EngngModel *sp = this->giveSlaveProblem(1);
    Domain *d = sp->giveDomain(1);   
    int numEl = d->giveNumberOfElements(); 
    
   
    // Solver paramters
    double l1 = 0; 
    double l2 = 100000; 
    double move = 0.2;

    FloatArray xOld = x;
    while (l2-l1 > 1.0e-4) {
        
        double lmid = 0.5*(l2+l1);
        double totalVolume = 0.0;
        double help = 0.0;
        //for (int i = 1; i < x.giveSize(); i++) {
        for (int i = 1; i <= numEl; i++) { // should only be el that contribute
            Element *el = d->giveElement(i);
            // xnew = max(0.001,max(x-move,min(1.,min(x+move,x.*sqrt(-dc./lmid)))));
            double temp1 = max(0.001, max(xOld.at(i)-move, min(1., min(xOld.at(i) + move, xOld.at(i) * sqrt(-dc.at(i)/lmid) ) )));
            //double temp1 = x.at(i) + move;
        //double temp2 = x.at(i) * sqrt(-dc.at(i)/lmid);
        //double temp3 = min(temp1, temp2);
        //temp2 = min(1.0, temp3);
        //temp1 = x.at(i) - move;
        //temp3 = max(temp1, temp2);
        //temp1 = max(0.001,temp3);
        ////designVarList.at(i) = max(0.001, max(x-move, min(1., min(x+move,x.*sqrt(-dCostFunction.at(i)/lmid)))));
            
            x.at(i) = temp1;
            double dV = el->computeVolumeAreaOrLength();
            totalVolume += dV;
            help += temp1 * dV; 
        }
        double averageDensity = help / totalVolume; 
        
        if ( averageDensity - objFunc->constraint > 0 ) {
            l1 = lmid;
        } else {
            l2 = lmid;
        }
    }

}

void
StaticFracture :: terminate(TimeStep *tStep)
{    
    for ( int subProb = 1; subProb <= this->giveNumberOfSlaveProblems(); subProb++ ) {
        EngngModel *sp = this->giveSlaveProblem(subProb);
       // sp->terminate(tStep);
    }

}




//remove
void
StaticFracture :: updateLoadVectors(TimeStep *tStep)
{
 NonLinearStatic ::updateLoadVectors(tStep);
}







} // end namespace oofem
