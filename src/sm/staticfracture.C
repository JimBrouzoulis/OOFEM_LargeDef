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
#include <vector>

namespace oofem {

REGISTER_EngngModel( StaticFracture );

StaticFracture :: StaticFracture(int i, EngngModel *_master) : NonLinearStatic(i, _master)
{
    updateStructureFlag = false; // if true, then the internal structure needs to be updated
}


void
StaticFracture :: solveYourselfAt(TimeStep *tStep)
{

    NonLinearStatic :: solveYourselfAt(tStep);

}

#include "isolinearelasticmaterial.h"
void 
StaticFracture :: updateYourself(TimeStep *tStep)
{

    NonLinearStatic :: updateYourself(tStep);

    Domain *d = this->giveDomain(1);   
    int numMat = d->giveNumberOfMaterialModels();

    double volFrac = 0.3;
    double penalty = 3.0;
    double prop = 30000;

    if (tStep->isTheFirstStep() ) {
        for (int i = 1; i < numMat; i++) {
            designVarList.resize(numMat);
            designVarList.at(i) = volFrac;
        }
    }

    double costFunction = 0.0;
    
    int numEl = d->giveNumberOfElements();
    FloatMatrix Ke; 
    FloatArray ae, help, dCostFunction(numMat);
    for (int i = 1; i < numMat; i++)
    {
        Element *el = d->giveElement(i);
        el->giveCharacteristicMatrix(Ke, SecantStiffnessMatrix, tStep);
        el->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, ae);
        help.beProductOf(Ke,ae);
        double temp =  ae.dotProduct(help);
        costFunction += pow( designVarList.at(i), penalty) * temp;
        dCostFunction.at(i) = -penalty * pow( designVarList.at(i), penalty-1.0) * temp;

    }

    this->optimalityCriteria(20, 20, designVarList, volFrac, dCostFunction);
    
    for (int i = 1; i < numMat; i++) {
        DynamicInputRecord ir;
        Material *mat = d->giveMaterial(i);
        mat->giveInputRecord(ir);
        //ir.setField(prop * pow( designVarList.at(i), penalty), _IFT_IsotropicLinearElasticMaterial_e);
        ir.setField(prop * designVarList.at(i), _IFT_IsotropicLinearElasticMaterial_e);
        mat->initializeFrom(&ir);
    }

    //designVarList.printYourself();
    
    printf("costfunction %e and sum design %e\n", costFunction, designVarList.sum());
}

void 
StaticFracture :: optimalityCriteria(int numElX, int numElY, FloatArray &x, double volFrac, FloatArray dCostFunction)
{
 double l1 = 0; 
 double l2 = 100000; 
 double move = 0.15;
 while (l2-l1 > 1.0e-4) {
    double lmid = 0.5*(l2+l1);

    for (int i = 1; i < x.giveSize(); i++)
    {
        //xnew = max(0.001,max(x-move,min(1.,min(x+move,x.*sqrt(-dc./lmid)))));
        double temp1 = x.at(i) + move;
        double temp2 = x.at(i) * sqrt(-dCostFunction.at(i)/lmid);
        double temp3 = temp1 < temp2 ? temp1 : temp2;
        temp2 = 1.0 < temp3 ? 1.0 : temp3;
        temp1 = x.at(i) - move;
        temp3 = temp1 > temp2 ? temp1 : temp2;
        temp1 = 0.001 > temp3 ? 0.001 : temp3;
        //designVarList.at(i) = max(0.001, max(x-move, min(1., min(x+move,x.*sqrt(-dCostFunction.at(i)/lmid)))));
        x.at(i) = temp1;
        
    }
    double density = x.sum();
    if ( density - volFrac*numElX*numElY > 0 ) {
        l1 = lmid;
    } else {
        l2 = lmid;
    }
 }

 //printf("l1 and l2 %e %e \n", l1,l2);
}

// remove
void
StaticFracture :: terminate(TimeStep *tStep)
{
    NonLinearStatic :: terminate(tStep);
}



void
StaticFracture :: updateLoadVectors(TimeStep *tStep)
{
 NonLinearStatic ::updateLoadVectors(tStep);
}



// Updating dofs and evaluating unknowns
#if 1

double 
StaticFracture ::  giveUnknownComponent(ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof)
{
    // Returns the unknown quantity corresponding to the dof
    if ( this->requiresUnknownsDictionaryUpdate() ) {
        int hash = this->giveUnknownDictHashIndx(mode, tStep);
        if ( dof->giveUnknowns()->includes(hash) ) {
            return dof->giveUnknowns()->at(hash);
        } else { // Value is not initiated in UnknownsDictionary
            return 0.0; ///@todo: how should one treat newly created dofs?
            //OOFEM_ERROR2( "giveUnknown:  Dof unknowns dictionary does not contain unknown of value mode (%s)", __ValueModeTypeToString(mode) );
        }
    } else {
        return NonLinearStatic ::  giveUnknownComponent(mode, tStep, d, dof);
    }
    
}


void
StaticFracture :: initializeDofUnknownsDictionary(TimeStep *tStep) 
{
    // Initializes all dof values to zero
    
    Domain *domain;
    Dof *iDof;
    DofManager *node;
    
    int nDofs;
    for ( int idomain = 1; idomain <= this->giveNumberOfDomains(); idomain++ ) {
        domain = this->giveDomain(idomain);
        int nnodes = domain->giveNumberOfDofManagers();
        for ( int inode = 1; inode <= nnodes; inode++ ) {
            node = domain->giveDofManager(inode);
            nDofs = node->giveNumberOfDofs();
            for ( int i = 1; i <= nDofs; i++ ) {
                iDof = node->giveDof(i);
                iDof->updateUnknownsDictionary(tStep->givePreviousStep(), VM_Total, 0.0);
            }
        }
    }
}



void
StaticFracture :: updateDofUnknownsDictionary(DofManager *inode, TimeStep *tStep)
{

    // update DoF unknowns dictionary. 
    Dof *iDof;
    double val;
    for ( int i = 1; i <= inode->giveNumberOfDofs(); i++ ) {
        iDof = inode->giveDof(i);
        int eqNum = iDof->__giveEquationNumber();
        if ( iDof->hasBc(tStep) ) { 
            val = iDof->giveBcValue(VM_Total, tStep);
        } else {
            if ( eqNum > 0 ) {
                val = totalDisplacement.at(eqNum);
            } else { // new eq number
                val = 0.0;
            }
        }

        iDof->updateUnknownsDictionary(tStep, VM_Total, val);
    }
}


void
StaticFracture :: setTotalDisplacementFromUnknownsInDictionary(EquationID type, ValueModeType mode, TimeStep *tStep) 
{
    // Sets the values in the displacement vector based on stored values in the unknowns dictionaries.
    // Used in the beginning of each time step.
    Domain *domain;
    DofManager *inode;
    Dof *iDof;
    for ( int idomain = 1; idomain <= this->giveNumberOfDomains(); idomain++ ) {
        domain = this->giveDomain(idomain);
        for ( int j = 1; j <= domain->giveNumberOfDofManagers(); j++ ) {
            inode = domain->giveDofManager(j);
            int eqNum;
            for ( int i = 1; i <= inode->giveNumberOfDofs(); i++ ) {
                iDof = inode->giveDof(i);
                eqNum = iDof->giveEqn();
                if ( eqNum > 0 ) {
                    double val = iDof->giveUnknown(mode, tStep);
                    totalDisplacement.at(eqNum) = val;
                }
            }
        }   
    }

}

#endif




} // end namespace oofem
