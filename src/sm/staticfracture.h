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

#ifndef staticfracture_h
#define staticfracture_h
#include "nlinearstatic.h"
#include "metastep.h"
#include "xfemmanager.h"
#include "fracturemanager.h"
#include "Element.h"

#include "engngm.h"

///@name Input fields for StaggeredProblem
//@{
#define _IFT_OptimizationProblem_nsteps "nsteps"
#define _IFT_StaticFracture_Name "optimizationproblem"
#define _IFT__IFT_OptimizationProblem_NumObjFunc "numobjfunc"
//#define _IFT__IFT_OptimizationProblem_Name_dtf "dtf"
//#define _IFT__IFT_OptimizationProblem_Name_timeDefinedByProb "timedefinedbyprob"
//#define _IFT__IFT_OptimizationProblem_Name_stepmultiplier "stepmultiplier"
//#define _IFT__IFT_OptimizationProblem_Name_prescribedtimes "prescribedtimes"
#define _IFT__IFT_OptimizationProblem_Name_prob1 "prob1"
#define _IFT__IFT_OptimizationProblem_Name_prob2 "prob2"
//@}
namespace oofem {
/**
 * This class implements a nonlinear static fracture problem.
 * It provides (or will) functionality for updating the model after each time step according to some propagation model
 */
class StaticFracture : public NonLinearStatic
{
protected:

    // from nlinearstatic
    virtual void solveYourselfAt(TimeStep *tStep);
    virtual void terminate(TimeStep *tStep);
    virtual void updateLoadVectors(TimeStep *tStep);
    virtual void updateYourself(TimeStep *tStep);
    void optimize(TimeStep *tStep);

    /// number of objective functions to minimize
    int numObjFunc;

    /// Number of engineering models to run.
    int nModels;
    
    /// List of slave models to which this model is coupled    
    IntArray coupledModels;

    /// List of engineering models to solve sequentially.
    AList< EngngModel > *emodelList;
    
    std :: string *inputStreamNames;

    int instanciateSlaveProblems();

public:
    StaticFracture(int i, EngngModel *_master = NULL);
    virtual ~StaticFracture(){};

    // new topology opt

    void optimalityCriteria(int numElX, int numElY, FloatArray &x, FloatArray &dc);
    FloatArray designVarList;
    double penalty;
    double volFrac;
    double min(double a, double b) { return a < b ? a : b; };
    double max(double a, double b) { return a > b ? a : b; };
    void costFunctionAndDerivative(Element *el, double &costFunction, double &dCostFunction, TimeStep *tStep);


    virtual void solveYourself();
    //virtual void initializeYourself(TimeStep *tStep) { }
    
    //virtual void doStepOutput(TimeStep *tStep);

    virtual int instanciateYourself(DataReader *dr, InputRecord *ir, const char *outFileName, const char *desc);

    virtual IRResultType initializeFrom(InputRecord *ir);
    //virtual void updateAttributes(MetaStep *mStep);
    virtual int giveNumberOfSlaveProblems() { return nModels; }
    virtual EngngModel *giveSlaveProblem(int i);

};

} // end namespace oofem
#endif // nlinearstatic_h
