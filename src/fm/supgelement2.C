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

#include "supgelement2.h"
#include "domain.h"
#include "timestep.h"
#include "load.h"
#include "intarray.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "fluidmodel.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "connectivitytable.h"
#endif

namespace oofem {
SUPGElement2 :: SUPGElement2(int n, Domain *aDomain) :
    SUPGElement(n, aDomain)
    // Constructor. Creates an element with number n, belonging to aDomain.
{
}


SUPGElement2 :: ~SUPGElement2()
// Destructor.
{ }

IRResultType
SUPGElement2 :: initializeFrom(InputRecord *ir)
{
    //const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    //IRResultType result;                               // Required by IR_GIVE_FIELD macro

    return SUPGElement :: initializeFrom(ir);
}



void
SUPGElement2 :: giveCharacteristicMatrix(FloatMatrix &answer,
                                          CharType mtrx, TimeStep *tStep)
//
// returns characteristics matrix of receiver according to mtrx
//
{
#if 0
    if ( mtrx == StiffnessMatrix ) {
        // support for stokes solver
        IntArray vloc, ploc;
        FloatMatrix h;
        int size = this->computeNumberOfDofs(EID_MomentumBalance_ConservationEquation);
        this->giveLocalVelocityDofMap(vloc);
        this->giveLocalPressureDofMap(ploc);
        answer.resize(size, size);
        answer.zero();
        this->computeAdvectionDerivativeTerm_MB(h, tStep);
        answer.assemble(h, vloc);
        this->computeDiffusionDerivativeTerm_MB(h, TangentStiffness, tStep);
        answer.assemble(h, vloc);
        this->computePressureTerm_MB(h, tStep);
        answer.assemble(h, vloc, ploc);
        this->computeLinearAdvectionTerm_MC(h, tStep);
        answer.assemble(h, ploc, vloc);
        this->computeAdvectionDerivativeTerm_MC(h, tStep);
        answer.assemble(h, ploc, vloc);
        this->computeDiffusionDerivativeTerm_MC(h, tStep);
        answer.assemble(h, ploc, vloc);
        this->computePressureTerm_MC(h, tStep);
        answer.assemble(h, ploc);
        this->computeLSICStabilizationTerm_MB(h, tStep);
        answer.assemble(h, vloc);
    } else
#endif
    {
        _error("giveCharacteristicMatrix: Unknown Type of characteristic mtrx.");
    }
}


void
SUPGElement2 :: giveCharacteristicVector(FloatArray &answer, CharType mtrx, ValueModeType mode,
                                          TimeStep *tStep)
//
// returns characteristics vector of receiver according to requested type
//
{
    if ( mtrx == ExternalForcesVector ) {
        // stokes flow
        IntArray vloc, ploc;
        FloatArray h;
        int size = this->computeNumberOfDofs(EID_MomentumBalance_ConservationEquation);
        this->giveLocalVelocityDofMap(vloc);
        this->giveLocalPressureDofMap(ploc);
        answer.resize(size);
        answer.zero();
        this->computeBCRhsTerm_MB(h, tStep);
        answer.assemble(h, vloc);
        this->computeBCRhsTerm_MC(h, tStep);
        answer.assemble(h, ploc);
    } else
#if 0
        if ( mtrx == InternalForcesVector ) {
        // stokes flow
        IntArray vloc, ploc;
        FloatArray h;
        int size = this->computeNumberOfDofs(EID_MomentumBalance_ConservationEquation);
        this->giveLocalVelocityDofMap(vloc);
        this->giveLocalPressureDofMap(ploc);
        answer.resize(size);
        answer.zero();
        this->computeAdvectionTerm_MB(h, tStep);
        answer.assemble(h, vloc);
        this->computeAdvectionTerm_MC(h, tStep);
        answer.assemble(h, ploc);
        this->computeDiffusionTerm_MB(h, tStep);
        answer.assemble(h, vloc);
        this->computeDiffusionTerm_MC(h, tStep);
        answer.assemble(h, ploc);

        FloatMatrix m1;
        FloatArray v, p;
        // add lsic stabilization term
        this->giveCharacteristicMatrix(m1, LSICStabilizationTerm_MB, tStep);
        //m1.times( lscale / ( dscale * uscale * uscale ) );
        this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, v);
        h.beProductOf(m1, v);
        answer.assemble(h, vloc);
        this->giveCharacteristicMatrix(m1, LinearAdvectionTerm_MC, tStep);
        //m1.times( 1. / ( dscale * uscale ) );
        h.beProductOf(m1, v);
        answer.assemble(h, ploc);


        // add pressure term
        this->giveCharacteristicMatrix(m1, PressureTerm_MB, tStep);
        this->computeVectorOf(EID_ConservationEquation, VM_Total, tStep, v);
        h.beProductOf(m1, v);
        answer.assemble(h, vloc);

        // pressure term
        this->giveCharacteristicMatrix(m1, PressureTerm_MC, tStep);
        this->computeVectorOf(EID_ConservationEquation, VM_Total, tStep, v);
        h.beProductOf(m1, v);
        answer.assemble(h, ploc);
    } else
#endif
    {
        _error("giveCharacteristicVector: Unknown Type of characteristic mtrx.");
    }
}


double
SUPGElement2 :: giveCharacteristicValue(CharType mtrx, TimeStep *tStep)
{
    if ( mtrx == CriticalTimeStep ) {
        return this->computeCriticalTimeStep(tStep);
    } else {
        _error("giveCharacteristicValue: Unknown Type of characteristic mtrx.");
    }

    return 0.0;
}


int
SUPGElement2 :: checkConsistency()
//
// check internal consistency
// mainly tests, whether material and crossSection data
// are safe for conversion to "Structural" versions
//
{
    int result = 1;
    /*
     * if (!this->giveMaterial()->testMaterialExtension(Material_TransportCapability)) {
     * _warning("checkConsistency : material without support for transport problems");
     * result =0;
     * }
     */
    return result;
}


void
SUPGElement2 :: updateInternalState(TimeStep *stepN)
{
    FloatArray stress;

    // force updating strains & stresses
    for ( int i = 0; i < numberOfIntegrationRules; i++ ) {
        IntegrationRule *iRule = integrationRulesArray [ i ];
        for ( int j = 0; j < iRule->giveNumberOfIntegrationPoints(); j++ ) {
            computeDeviatoricStress(stress, iRule->getIntegrationPoint(j), stepN);
        }
    }
}

void
SUPGElement2 :: printOutputAt(FILE *file, TimeStep *stepN)
// Performs end-of-step operations.
{
    int i;

#ifdef __PARALLEL_MODE
    fprintf( file, "element %d [%8d] :\n", this->giveNumber(), this->giveGlobalNumber() );
#else
    fprintf(file, "element %d :\n", number);
#endif

    for ( i = 0; i < numberOfIntegrationRules; i++ ) {
        integrationRulesArray [ i ]->printOutputAt(file, stepN);
    }
}


#ifdef __OOFEG
int
SUPGElement2 :: giveInternalStateAtNode(FloatArray &answer, InternalStateType type, InternalStateMode mode,
                                        int node, TimeStep *atTime)
{
    int indx = 1;
    Node *n = this->giveNode(node);

    if ( type == IST_Velocity ) {
        answer.resize( this->giveSpatialDimension() );
        int dofindx;
        if ( ( dofindx = n->findDofWithDofId(V_u) ) ) {
            answer.at(indx++) = n->giveDof(dofindx)->giveUnknown(VM_Total, atTime);
        }

        if ( ( dofindx = n->findDofWithDofId(V_v) ) ) {
            answer.at(indx++) = n->giveDof(dofindx)->giveUnknown(VM_Total, atTime);
        }

        if ( ( dofindx = n->findDofWithDofId(V_w) ) ) {
            answer.at(indx++) = n->giveDof(dofindx)->giveUnknown(VM_Total, atTime);
        }

        return 1;
    } else if ( type == IST_Pressure ) {
        int dofindx;
        if ( ( dofindx = n->findDofWithDofId(P_f) ) ) {
            answer.resize(1);
            answer.at(1) = n->giveDof(dofindx)->giveUnknown(VM_Total, atTime);
            return 1;
        } else {
            return 0;
        }
    } else {
        return Element :: giveInternalStateAtNode(answer, type, mode, node, atTime);
    }
}

#endif

void
SUPGElement2 :: computeAccelerationTerm_MB(FloatMatrix &answer, TimeStep *atTime)
{
    FloatMatrix n, b;
    double dV, rho;
    int k, undofs = this->computeNumberOfDofs(EID_MomentumBalance);
    GaussPoint *gp;

    answer.resize(undofs, undofs);
    answer.zero();

    int rule = 2;
    IntegrationRule *iRule = this->integrationRulesArray [ rule ];
    /* consistent part + supg stabilization term */
    for ( k = 0; k < iRule->giveNumberOfIntegrationPoints(); k++ ) {
        gp = iRule->getIntegrationPoint(k);
        this->computeNuMatrix(n, gp);
        this->computeUDotGradUMatrix( b, gp, atTime->givePreviousStep() );
        dV  = this->computeVolumeAround(gp);
        rho = this->giveMaterial()->give('d', gp);
        /* consistent part */
        answer.plusProductUnsym(n, n, rho * dV);
        /* supg stabilization */
        answer.plusProductUnsym(b, n, rho * t_supg * dV);
    }
}

void
SUPGElement2 :: computeAdvectionTerm_MB(FloatArray &answer, TimeStep *atTime)
{
    FloatMatrix n, b, bn;
    FloatArray u, v(3);
    double dV, rho, coeff, sum;
    int i, k, undofs = this->computeNumberOfDofs(EID_MomentumBalance);
    int isd, nsd = this->giveNumberOfSpatialDimensions();
    GaussPoint *gp;

    answer.resize(undofs);
    answer.zero();

    this->computeVectorOf(EID_MomentumBalance, VM_Total, atTime, u);

    int rule = 2;
    IntegrationRule *iRule = this->integrationRulesArray [ rule ];
    /* consistent part + supg stabilization term */
    for ( k = 0; k < iRule->giveNumberOfIntegrationPoints(); k++ ) {
        gp = iRule->getIntegrationPoint(k);
        this->computeNuMatrix(n, gp);
        this->computeUDotGradUMatrix( bn, gp, atTime->givePreviousStep() );
        this->computeUDotGradUMatrix(b, gp, atTime);
        v.beProductOf(b, u);
        dV  = this->computeVolumeAround(gp);
        rho = this->giveMaterial()->give('d', gp);
        /* consistent part */
        coeff = rho * dV;
        for ( i = 1; i <= undofs; i++ ) {
            for ( sum = 0.0, isd = 1; isd <= nsd; isd++ ) {
                sum += n.at(isd, i) * v.at(isd);
            }

            answer.at(i) += coeff * sum;
        }

        /* supg stabilization */
        coeff = t_supg * rho * dV;
        for ( i = 1; i <= undofs; i++ ) {
            for ( sum = 0, isd = 1; isd <= nsd; isd++ ) {
                sum += bn.at(isd, i) * v.at(isd);
            }

            answer.at(i) += coeff * sum;
        }
    }
}

void
SUPGElement2 :: computeAdvectionDerivativeTerm_MB(FloatMatrix &answer, TimeStep *atTime)
{
    FloatMatrix n, b, bn, grad_u, grad_uN, N;
    double dV, rho;
    int k, undofs = this->computeNumberOfDofs(EID_MomentumBalance);
    GaussPoint *gp;

    answer.resize(undofs, undofs);
    answer.zero();

    int rule = 2;
    IntegrationRule *iRule = this->integrationRulesArray [ rule ];
    /* consistent part + supg stabilization term */
    for ( k = 0; k < iRule->giveNumberOfIntegrationPoints(); k++ ) {
        gp = iRule->getIntegrationPoint(k);
        this->computeNuMatrix(n, gp);
        this->computeUDotGradUMatrix( bn, gp, atTime->givePreviousStep() );
        this->computeUDotGradUMatrix(b, gp, atTime);
        dV  = this->computeVolumeAround(gp);
        rho = this->giveMaterial()->give('d', gp);

        this->computeGradUMatrix(grad_u, gp, atTime);


        /* consistent part */
        answer.plusProductUnsym(n, b, rho * dV);

        grad_uN.beProductOf(grad_u, n);
        answer.plusProductUnsym(n, grad_uN, rho * dV);
        /* supg stabilization */
        answer.plusProductUnsym(bn, b, t_supg * rho * dV);
        answer.plusProductUnsym(bn, grad_uN, t_supg * rho * dV);
    }
}

void
SUPGElement2 :: computeDiffusionTerm_MB(FloatArray &answer, TimeStep *atTime)
{
    int undofs = this->computeNumberOfDofs(EID_MomentumBalance);
    FloatArray u, eps, stress, bs, dDB_u;
    FloatMatrix b, un_gu, dDB;
    GaussPoint *gp;
    double coeff, sum, dV, Re = static_cast<FluidModel*>(domain->giveEngngModel())->giveReynoldsNumber();
    int isd, nsd;

    answer.resize(undofs);
    answer.zero();

    nsd = this->giveNumberOfSpatialDimensions();
    this->computeVectorOf(EID_MomentumBalance, VM_Total, atTime, u);

    int rule = 1;
    IntegrationRule *iRule = this->integrationRulesArray [ rule ];
    for ( int k = 0; k < iRule->giveNumberOfIntegrationPoints(); k++ ) {
        gp = iRule->getIntegrationPoint(k);
        dV  = this->computeVolumeAround(gp);
        this->computeBMatrix(b, gp);
        this->computeDivTauMatrix(dDB, gp, atTime);
        this->computeUDotGradUMatrix( un_gu, gp, atTime->givePreviousStep() );
        eps.beProductOf(b, u);
        static_cast< FluidDynamicMaterial * >( this->giveMaterial() )->computeDeviatoricStressVector(stress, gp, eps, atTime);
        dDB_u.beProductOf(dDB, u);
        /* consistent part */
        answer.plusProduct(b, stress, dV / Re);

        /* SUPG term */
        //answer.plusProductUnsym(un_gu,dDB_u, t_supg * dV * (-1.0) * (1./Re));

        coeff = ( -1.0 ) * t_supg * dV;
        for ( int i = 1; i <= undofs; i++ ) {
            for ( sum = 0, isd = 1; isd <= nsd; isd++ ) {
                sum += un_gu.at(isd, i) * dDB_u.at(isd);
            }

            answer.at(i) += coeff * sum;
        }
    }
}

void
SUPGElement2 :: computeDiffusionDerivativeTerm_MB(FloatMatrix &answer, MatResponseMode mode, TimeStep *atTime)
{
    int undofs = this->computeNumberOfDofs(EID_MomentumBalance);

    answer.resize(undofs, undofs);
    answer.zero();
    FloatMatrix _db, _d, _b, dDB, un_gu;
    double dV, Re = static_cast<FluidModel*>(domain->giveEngngModel())->giveReynoldsNumber();
    GaussPoint *gp;
    FloatArray dDB_u;

    int rule = 1;
    IntegrationRule *iRule = this->integrationRulesArray [ rule ];
    for ( int k = 0; k < iRule->giveNumberOfIntegrationPoints(); k++ ) {
        gp = iRule->getIntegrationPoint(k);
        dV  = this->computeVolumeAround(gp);
        this->computeBMatrix(_b, gp);
        static_cast< FluidDynamicMaterial * >( this->giveMaterial() )->giveDeviatoricStiffnessMatrix(_d, mode, gp, atTime);

        this->computeDivTauMatrix(dDB, gp, atTime);
        this->computeUDotGradUMatrix( un_gu, gp, atTime->givePreviousStep() );
        /* standard term */
        _db.beProductOf(_d, _b);
        answer.plusProductUnsym(_b, _db, dV / Re); //answer.plusProduct (_b,_db,area);

        /* SUPG term */

        answer.plusProductUnsym( un_gu, dDB, t_supg * dV * ( -1.0 ) * ( 1. / Re ) );
    }
}


void
SUPGElement2 :: computePressureTerm_MB(FloatMatrix &answer, TimeStep *atTime)
{
    double dV;
    GaussPoint *gp;
    FloatMatrix gu, np, b;
    int undofs = this->computeNumberOfDofs(EID_MomentumBalance);
    int pndofs = this->computeNumberOfDofs(EID_ConservationEquation);

    answer.resize(undofs, pndofs);
    answer.zero();

    int rule = 1;
    IntegrationRule *iRule = this->integrationRulesArray [ rule ];
    for ( int k = 0; k < iRule->giveNumberOfIntegrationPoints(); k++ ) {
        gp = iRule->getIntegrationPoint(k);
        dV  = this->computeVolumeAround(gp);
        this->computeDivUMatrix(gu, gp);
        this->computeNpMatrix(np, gp);

        /*alternative computing*/
        //this->computeNuMatrix(gu, gp);
        //this->computeGradPMatrix(np, gp);

        /* standard term */


        answer.plusProductUnsym(gu, np, ( -1.0 ) * dV);
    }

    iRule = this->integrationRulesArray [ 1 ];
    for ( int k = 0; k < iRule->giveNumberOfIntegrationPoints(); k++ ) {
        gp = iRule->getIntegrationPoint(k);
        dV  = this->computeVolumeAround(gp);
        this->computeUDotGradUMatrix( b, gp, atTime->givePreviousStep() );
        this->computeGradPMatrix(np, gp);

        // supg term
        answer.plusProductUnsym(b, np, t_supg * dV);
    }
}
void
SUPGElement2 :: computeLSICStabilizationTerm_MB(FloatMatrix &answer, TimeStep *atTime)
{
    int undofs = this->computeNumberOfDofs(EID_MomentumBalance);
    double dV, rho;
    GaussPoint *gp;
    FloatMatrix b;

    answer.resize(undofs, undofs);
    answer.zero();

    int rule = 0;
    IntegrationRule *iRule = this->integrationRulesArray [ rule ];
    for ( int k = 0; k < iRule->giveNumberOfIntegrationPoints(); k++ ) {
        gp = iRule->getIntegrationPoint(k);
        dV  = this->computeVolumeAround(gp);
        rho = this->giveMaterial()->give('d', gp);
        this->computeDivUMatrix(b, gp);

        answer.plusProductSymmUpper(b, b, dV * rho * t_lsic);
    }

    answer.symmetrized();
}

void
SUPGElement2 :: computeLinearAdvectionTerm_MC(FloatMatrix &answer, TimeStep *atTime)
{
    double dV;
    GaussPoint *gp;
    FloatMatrix gu, np;
    int undofs = this->computeNumberOfDofs(EID_MomentumBalance);
    int pndofs = this->computeNumberOfDofs(EID_ConservationEquation);

    answer.resize(pndofs, undofs);
    answer.zero();

    int rule = 0;
    IntegrationRule *iRule = this->integrationRulesArray [ rule ];
    for ( int k = 0; k < iRule->giveNumberOfIntegrationPoints(); k++ ) {
        gp = iRule->getIntegrationPoint(k);
        dV  = this->computeVolumeAround(gp);
        this->computeDivUMatrix(gu, gp);
        this->computeNpMatrix(np, gp);

        /* standard term */

        answer.plusProductUnsym(np, gu, dV);
    }
}

void
SUPGElement2 :: computeAdvectionTerm_MC(FloatArray &answer, TimeStep *atTime)
{
    // N_epsilon (due to PSPG stabilization)
    FloatMatrix g, b;
    FloatArray u, v;
    GaussPoint *gp;
    double dV, coeff, sum;
    int pndofs = this->computeNumberOfDofs(EID_ConservationEquation);
    int isd, nsd = this->giveNumberOfSpatialDimensions();

    answer.resize(pndofs);
    answer.zero();

    int rule = 1;
    IntegrationRule *iRule = this->integrationRulesArray [ rule ];
    this->computeVectorOf(EID_MomentumBalance, VM_Total, atTime, u);

    /* pspg stabilization term */
    for ( int k = 0; k < iRule->giveNumberOfIntegrationPoints(); k++ ) {
        gp = iRule->getIntegrationPoint(k);
        this->computeGradPMatrix(g, gp);
        this->computeUDotGradUMatrix(b, gp, atTime);
        v.beProductOf(b, u);
        dV  = this->computeVolumeAround(gp);
        coeff = dV * t_pspg;
        for ( int i = 1; i <= pndofs; i++ ) {
            for ( sum = 0.0, isd = 1; isd <= nsd; isd++ ) {
                sum += g.at(isd, i) * v.at(isd);
            }

            answer.at(i) += coeff * sum;
        }
    }
}


void
SUPGElement2 :: computeAdvectionDerivativeTerm_MC(FloatMatrix &answer, TimeStep *atTime)
{
    int pndofs = this->computeNumberOfDofs(EID_ConservationEquation);
    int undofs = this->computeNumberOfDofs(EID_MomentumBalance);
    GaussPoint *gp;
    FloatMatrix g, b;
    double dV, coeff;

    answer.resize(pndofs, undofs);
    answer.zero();

    int rule = 1;
    IntegrationRule *iRule = this->integrationRulesArray [ rule ];
    /* pspg stabilization term */
    for ( int k = 0; k < iRule->giveNumberOfIntegrationPoints(); k++ ) {
        gp = iRule->getIntegrationPoint(k);
        this->computeGradPMatrix(g, gp);
        this->computeUDotGradUMatrix(b, gp, atTime);
        dV  = this->computeVolumeAround(gp);
        coeff = dV * t_pspg;
        answer.plusProductUnsym(g, b, coeff);
    }
}

void
SUPGElement2 :: computeDiffusionDerivativeTerm_MC(FloatMatrix &answer, TimeStep *atTime)
{
    int pndofs = this->computeNumberOfDofs(EID_ConservationEquation);
    int undofs = this->computeNumberOfDofs(EID_MomentumBalance);

    answer.resize(pndofs, undofs);
    FloatMatrix dDB, _d, g;
    double dV, coeff, rho;
    GaussPoint *gp;

    int rule = 1;
    IntegrationRule *iRule = this->integrationRulesArray [ rule ];

    for ( int k = 0; k < iRule->giveNumberOfIntegrationPoints(); k++ ) {
        gp = iRule->getIntegrationPoint(k);
        dV  = this->computeVolumeAround(gp);
        rho = this->giveMaterial()->give('d', gp);

        coeff = ( -1.0 ) * dV * t_pspg / rho;
        //( ( FluidDynamicMaterial * ) this->giveMaterial() )->giveDeviatoricStiffnessMatrix(_d, TangentStiffness,gp, atTime);

        this->computeDivTauMatrix(dDB, gp, atTime);
        this->computeGradPMatrix(g, gp);

        answer.plusProductUnsym(g, dDB, coeff);
    }


    answer.zero();
}

void
SUPGElement2 :: computeDiffusionTerm_MC(FloatArray &answer, TimeStep *atTime)
{
    int pndofs = this->computeNumberOfDofs(EID_ConservationEquation);
    answer.resize(pndofs);

    /*
     * FloatMatrix dDB, _d, g;
     * double sum, dV, coeff, rho, isd, nsd = this->giveNumberOfSpatialDimensions();
     * GaussPoint *gp;
     * IntegrationRule *iRule = this->integrationRulesArray [ 1 ];
     * int i,k;
     * FloatArray u, dDB_u;
     *
     * this->computeVectorOf(EID_MomentumBalance, VM_Total, atTime, u);
     *
     * for ( k = 0; k < iRule->giveNumberOfIntegrationPoints(); k++ ) {
     *  gp = iRule->getIntegrationPoint(k);
     *  dV  = this->computeVolumeAround(gp);
     *  rho = this->giveMaterial()->give('d', gp);
     *
     *  coeff = (-1.0) * dV * t_pspg / rho;
     *
     *  this->computeGradPMatrix(g, gp);
     *  this->computeDivTauMatrix(dDB, gp, atTime);
     *
     *  dDB_u.beProductOf(dDB, u);
     *
     *  for ( i = 1; i <= pndofs; i++ ) {
     *      for ( sum = 0, isd = 1; isd <= nsd; isd++ ) {
     *          sum += g.at(isd, i) * dDB_u.at(isd);
     *      }
     *
     *      answer.at(i) += coeff * sum;
     *  }
     *
     *
     *
     * }
     */

    answer.zero();
}


void
SUPGElement2 :: computeAccelerationTerm_MC(FloatMatrix &answer, TimeStep *atTime)
{
    int pndofs = this->computeNumberOfDofs(EID_ConservationEquation);
    int undofs = this->computeNumberOfDofs(EID_MomentumBalance);
    double dV, coeff;
    FloatMatrix g, n;
    GaussPoint *gp;

    answer.resize(pndofs, undofs);
    answer.zero();
    // pspg stabilization term: M_\epsilon term


    int rule = 0;
    IntegrationRule *iRule = this->integrationRulesArray [ rule ];

    for ( int k = 0; k < iRule->giveNumberOfIntegrationPoints(); k++ ) {
        gp = iRule->getIntegrationPoint(k);
        this->computeGradPMatrix(g, gp);
        this->computeNuMatrix(n, gp);
        dV  = this->computeVolumeAround(gp);
        coeff = dV * t_pspg;
        answer.plusProductUnsym(g, n, coeff);
    }
}

void
SUPGElement2 :: computePressureTerm_MC(FloatMatrix &answer, TimeStep *atTime)
{
    int pndofs = this->computeNumberOfDofs(EID_ConservationEquation);
    double dV, rho, coeff;
    GaussPoint *gp;
    FloatMatrix g;

    answer.resize(pndofs, pndofs);
    answer.zero();
    int rule = 0;
    IntegrationRule *iRule = this->integrationRulesArray [ rule ];

    for ( int k = 0; k < iRule->giveNumberOfIntegrationPoints(); k++ ) {
        gp = iRule->getIntegrationPoint(k);
        this->computeGradPMatrix(g, gp);
        dV  = this->computeVolumeAround(gp);
        rho = this->giveMaterial()->give('d', gp);
        coeff = dV * t_pspg / rho;
        answer.plusProductSymmUpper(g, g, coeff);
    }

    answer.symmetrized();
}


void
SUPGElement2 :: computeBCRhsTerm_MB(FloatArray &answer, TimeStep *atTime)
{
    int undofs = this->computeNumberOfDofs(EID_MomentumBalance);

    answer.resize(undofs);
    answer.zero();

    int id, n, nLoads;
    double dV, rho;
    Load *load;
    bcGeomType ltype;

    int rule = 0;
    IntegrationRule *iRule = this->integrationRulesArray [ rule ];
    GaussPoint *gp;
    FloatArray un, gVector, s, helpLoadVector;
    FloatMatrix b, nu;

    // add body load (gravity) termms
    nLoads = this->giveBodyLoadArray()->giveSize();
    for ( int i = 1; i <= nLoads; i++ ) {
        load  = domain->giveLoad( bodyLoadArray.at(i) );
        ltype = load->giveBCGeoType();
        if ( ( ltype == BodyLoadBGT ) && ( load->giveBCValType() == ForceLoadBVT ) ) {
            load->computeComponentArrayAt(gVector, atTime, VM_Total);
            if ( gVector.giveSize() ) {
                for ( int k = 0; k < iRule->giveNumberOfIntegrationPoints(); k++ ) {
                    gp = iRule->getIntegrationPoint(k);
                    this->computeUDotGradUMatrix( b, gp, atTime->givePreviousStep() );
                    this->computeNuMatrix(nu, gp);
                    dV  = this->computeVolumeAround(gp);
                    rho = this->giveMaterial()->give('d', gp);
                    answer.plusProduct(b, gVector, t_supg * rho * dV);
                    answer.plusProduct(nu, gVector, rho * dV);
                }
            }
        }
    }

    // integrate tractions
    // if no traction bc applied but side marked as with traction load
    // then zero traction is assumed !!!

    // loop over boundary load array
    nLoads = this->giveBoundaryLoadArray()->giveSize() / 2;
    for ( int i = 1; i <= nLoads; i++ ) {
        n = boundaryLoadArray.at(1 + ( i - 1 ) * 2);
        id = boundaryLoadArray.at(i * 2);
        load  = domain->giveLoad(n);
        ltype = load->giveBCGeoType();
        if ( ltype == EdgeLoadBGT ) {
            this->computeEdgeLoadVector_MB(helpLoadVector, load, id, atTime);
            if ( helpLoadVector.giveSize() ) {
                answer.add(helpLoadVector);
            }
        } else if ( ltype == SurfaceLoadBGT ) {
            this->computeSurfaceLoadVector_MB(helpLoadVector, load, id, atTime);
            if ( helpLoadVector.giveSize() ) {
                answer.add(helpLoadVector);
            }
        } else {
            _error("computeForceLoadVector : unsupported load type class");
        }
    }
}



void
SUPGElement2 :: computeBCRhsTerm_MC(FloatArray &answer, TimeStep *atTime)
{
    int nLoads;
    double dV;
    Load *load;
    bcGeomType ltype;
    FloatArray s, gVector, helpLoadVector;
    FloatMatrix g;
    int pndofs = this->computeNumberOfDofs(EID_ConservationEquation);

    int rule = 1;
    IntegrationRule *iRule = this->integrationRulesArray [ rule ];
    GaussPoint *gp;

    answer.resize(pndofs);
    answer.zero();
    nLoads = this->giveBodyLoadArray()->giveSize();
    for ( int i = 1; i <= nLoads; i++ ) {
        load  = domain->giveLoad( bodyLoadArray.at(i) );
        ltype = load->giveBCGeoType();
        if ( ( ltype == BodyLoadBGT ) && ( load->giveBCValType() == ForceLoadBVT ) ) {
            load->computeComponentArrayAt(gVector, atTime, VM_Total);
            if ( gVector.giveSize() ) {
                for ( int k = 0; k < iRule->giveNumberOfIntegrationPoints(); k++ ) {
                    gp = iRule->getIntegrationPoint(k);
                    this->computeGradPMatrix(g, gp);
                    dV  = this->computeVolumeAround(gp);
                    answer.plusProduct(g, gVector, t_pspg * dV);
                }
            }
        }
    }

    // integrate tractions
    // if no traction bc applied but side marked as with traction load
    // then zero traction is assumed !!!

    // loop over boundary load array
    nLoads = this->giveBoundaryLoadArray()->giveSize() / 2;
    for ( int i = 1; i <= nLoads; i++ ) {
        int n = boundaryLoadArray.at(1 + ( i - 1 ) * 2);
        int id = boundaryLoadArray.at(i * 2);
        load = domain->giveLoad(n);
        ltype = load->giveBCGeoType();
        if ( ltype == EdgeLoadBGT ) {
            this->computeEdgeLoadVector_MC(helpLoadVector, load, id, atTime);
            if ( helpLoadVector.giveSize() ) {
                answer.add(helpLoadVector);
            }
        } else if ( ltype == SurfaceLoadBGT ) {
            this->computeSurfaceLoadVector_MC(helpLoadVector, load, id, atTime);
            if ( helpLoadVector.giveSize() ) {
                answer.add(helpLoadVector);
            }
        } else {
            _error("computeForceLoadVector : unsupported load type class");
        }
    }
}

void
SUPGElement2 :: computeEdgeLoadVector_MB(FloatArray &answer, Load *load, int id, TimeStep *stepN)
{
    _error("computeEdgeLoadVectorAt_MB: not implemented");
}

void
SUPGElement2 :: computeSurfaceLoadVector_MB(FloatArray &answer, Load *load, int id, TimeStep *stepN)
{
    _error("computeSurfaceLoadVectorAt_MB: not implemented");
}

void
SUPGElement2 :: computeEdgeLoadVector_MC(FloatArray &answer, Load *load, int id, TimeStep *stepN)
{
    _error("computeEdgeLoadVectorAt_MC: not implemented");
}

void
SUPGElement2 :: computeSurfaceLoadVector_MC(FloatArray &answer, Load *load, int id, TimeStep *stepN)
{
    _error("computeEdgeLoadVectorAt_MC: not implemented");
}


void
SUPGElement2 :: computeDeviatoricStrain(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
{
    FloatArray u;
    FloatMatrix b;
    this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, u);

    this->computeBMatrix(b, gp);
    answer.beProductOf(b, u);
}

void
SUPGElement2 :: computeDeviatoricStress(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
{
    FloatArray eps;

    // compute deviatoric strain
    this->computeDeviatoricStrain(eps, gp, tStep);
    // call material to compute stress
    static_cast< FluidDynamicMaterial * >( this->giveMaterial() )->computeDeviatoricStressVector(answer, gp, eps, tStep);
}

} // end namespace oofem
