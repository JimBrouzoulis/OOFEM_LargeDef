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

#include "quad1mindlin.h"
#include "node.h"
#include "material.h"
#include "crosssection.h"
#include "structuralms.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "load.h"
#include "structuralcrosssection.h"
#include "mathfem.h"
#include "classfactory.h"

namespace oofem {

REGISTER_Element( Quad1Mindlin );

FEI2dQuadLin Quad1Mindlin :: interp_lin(1, 2);

Quad1Mindlin :: Quad1Mindlin(int n, Domain *aDomain) :
    NLStructuralElement(n, aDomain)
{
	numberOfGaussPoints = 4;
    numberOfDofMans = 4;
}


FEInterpolation *
Quad1Mindlin :: giveInterpolation(DofIDItem id) const
{
    return & interp_lin;
}


void
Quad1Mindlin :: computeGaussPoints()
// Sets up the array containing the four Gauss points of the receiver.
{
    if ( !integrationRulesArray ) {
        numberOfIntegrationRules = 1;
        integrationRulesArray = new IntegrationRule * [ numberOfIntegrationRules ];
        integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 5);
        this->giveCrossSection()->setupIntegrationPoints( *integrationRulesArray[0], numberOfGaussPoints, this );
    }
}


void
Quad1Mindlin :: computeBodyLoadVectorAt(FloatArray &answer, Load *forLoad, TimeStep *stepN, ValueModeType mode)
{
    // Only gravity load
    double dV, load;
    GaussPoint *gp;
    FloatArray force, gravity, n;

    if ( ( forLoad->giveBCGeoType() != BodyLoadBGT ) || ( forLoad->giveBCValType() != ForceLoadBVT ) ) {
        _error("computeBodyLoadVectorAt: unknown load type");
    }

    // note: force is assumed to be in global coordinate system.
    forLoad->computeComponentArrayAt(gravity, stepN, mode);

    force.resize(0);
    if ( gravity.giveSize() ) {
        IntegrationRule *ir = integrationRulesArray [ 0 ]; ///@todo Other/higher integration for lumped mass matrices perhaps?
        for ( int i = 0; i < ir->giveNumberOfIntegrationPoints(); ++i) {
            gp = ir->getIntegrationPoint(i);

            this->interp_lin.evalN(n, *gp->giveCoordinates(), FEIElementGeometryWrapper(this));
            dV = this->computeVolumeAround(gp) * this->giveCrossSection()->give(CS_Thickness);
            load = this->giveMaterial()->give('d', gp) * gravity.at(3) * dV;

            force.add(load, n);
        }

        answer.resize(12);
        answer.zero();

        answer.at(1)  = force.at(1);
        answer.at(4)  = force.at(2);
        answer.at(7)  = force.at(3);
        answer.at(10) = force.at(4);

    } else {
        answer.resize(0);
    }
}


void
Quad1Mindlin :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
// Returns the [5x9] strain-displacement matrix {B} of the receiver,
// evaluated at gp.
{
    FloatArray n;
    FloatMatrix dn;

    this->interp_lin.evaldNdx(dn, *gp->giveCoordinates(),  FEIElementGeometryWrapper(this));
    this->interp_lin.evalN(n, *gp->giveCoordinates(),  FEIElementGeometryWrapper(this));

    answer.resize(5, 12);
    answer.zero();

    ///@todo Check sign here
    for (int i = 0; i < 4; ++i) {
        answer(0, 1 + i*3) = dn(i,0);
        answer(1, 2 + i*3) = dn(i,1);
        answer(2, 1 + i*3) = dn(i,1);
        answer(2, 2 + i*3) = dn(i,0);

        answer(3, 0 + i*3) = -dn(i,1);
        answer(3, 2 + i*3) = n(i);
        answer(4, 0 + i*3) = -dn(i,0);
        answer(4, 1 + i*3) = n(i);
    }
}


void
Quad1Mindlin :: computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer)
// Returns the [3x9] displacement interpolation matrix {N} of the receiver,
// evaluated at gp.
{
    FloatArray n;
    this->interp_lin.evalN(n, *gp->giveCoordinates(), FEIElementGeometryWrapper(this));
    answer.beNMatrixOf(n, 3);
}


IRResultType
Quad1Mindlin :: initializeFrom(InputRecord *ir)
{
    this->numberOfGaussPoints = 4;
    return this->NLStructuralElement :: initializeFrom(ir);
}


void
Quad1Mindlin :: giveDofManDofIDMask(int inode, EquationID, IntArray &answer) const
{
    answer.setValues(3, D_w, R_u, R_v);
}


void
Quad1Mindlin :: computeMidPlaneNormal(FloatArray &answer, const GaussPoint *gp)
{
    FloatArray u, v;
    u.beDifferenceOf( * this->giveNode(2)->giveCoordinates(), * this->giveNode(1)->giveCoordinates() );
    v.beDifferenceOf( * this->giveNode(3)->giveCoordinates(), * this->giveNode(1)->giveCoordinates() );

    answer.beVectorProductOf(u, v);
    answer.normalize();
}


double
Quad1Mindlin :: giveCharacteristicLenght(GaussPoint *gp, const FloatArray &normalToCrackPlane)
{
    return this->giveLenghtInDir(normalToCrackPlane);
}


double
Quad1Mindlin :: computeVolumeAround(GaussPoint *gp)
{
    double detJ, weight;

    weight = gp->giveWeight();
    detJ = fabs( this->interp_lin.giveTransformationJacobian( * gp->giveCoordinates(), FEIElementGeometryWrapper(this) ) );
    return detJ * weight;
}


void
Quad1Mindlin :: computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep)
// Returns the lumped mass matrix of the receiver.
{
    GaussPoint *gp;
    double dV, mass = 0.;

    IntegrationRule *ir = integrationRulesArray [ 0 ]; ///@todo Other/higher integration for lumped mass matrices perhaps?
    for ( int i = 0; i < ir->giveNumberOfIntegrationPoints(); ++i) {
        gp = ir->getIntegrationPoint(i);

        dV = this->computeVolumeAround(gp);
        mass += dV * this->giveMaterial()->give('d', gp);
    }

    answer.resize(12, 12);
    answer.zero();
    answer.at(1, 1) = mass*0.25;
    answer.at(4, 4) = mass*0.25;
    answer.at(7, 7) = mass*0.25;
    answer.at(10, 10) = mass*0.25;
}


int
Quad1Mindlin :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *atTime)
{
    if ( type == IST_ShellForceMomentumTensor ) {
        answer = static_cast< StructuralMaterialStatus * >( this->giveMaterial()->giveStatus(gp) )->giveStressVector();
        return 1;
    } else if ( type == IST_ShellStrainCurvatureTensor ) {
        answer = static_cast< StructuralMaterialStatus * >( this->giveMaterial()->giveStatus(gp) )->giveStrainVector();
        return 1;
    } else {
        return NLStructuralElement::giveIPValue(answer, gp, type, atTime);
    }
}


void
Quad1Mindlin :: computeEgdeNMatrixAt(FloatMatrix &answer, int iedge, GaussPoint *gp)
{
    IntArray edgeNodes;
    FloatArray n;

    this->interp_lin.edgeEvalN( n, iedge, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );
    this->interp_lin.computeLocalEdgeMapping(edgeNodes, iedge);

    answer.beNMatrixOf(n, 3);
}


void
Quad1Mindlin :: giveEdgeDofMapping(IntArray &answer, int iEdge) const
{
    if ( iEdge == 1 ) { // edge between nodes 1 2
        answer.setValues(6, 1,2,3,4,5,6);
    } else if ( iEdge == 2 ) { // edge between nodes 2 3
        answer.setValues(6, 4,5,6,7,8,9);
    } else if ( iEdge == 3 ) { // edge between nodes 3 4
        answer.setValues(6, 7,8,9,10,11,12);
    } else if ( iEdge == 4 ) { // edge between nodes 4 1
        answer.setValues(6, 10,11,12,1,2,3);
    } else {
        _error("giveEdgeDofMapping: wrong edge number");
    }
}


double
Quad1Mindlin :: computeEdgeVolumeAround(GaussPoint *gp, int iEdge)
{
    double detJ = this->interp_lin.edgeGiveTransformationJacobian( iEdge, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );
    return detJ *gp->giveWeight();
}


void
Quad1Mindlin :: computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge)
{
    this->interp_lin.edgeLocal2global( answer, iEdge, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );
}


int
Quad1Mindlin :: computeLoadLEToLRotationMatrix(FloatMatrix &answer, int iEdge, GaussPoint *gp)
{
    double dx, dy, length;
    IntArray edgeNodes;
    Node *nodeA, *nodeB;

    answer.resize(3, 3);
    answer.zero();

    this->interp_lin.computeLocalEdgeMapping(edgeNodes, iEdge);

    nodeA = this->giveNode(edgeNodes.at(1));
    nodeB = this->giveNode(edgeNodes.at(2));

    dx = nodeB->giveCoordinate(1) - nodeA->giveCoordinate(1);
    dy = nodeB->giveCoordinate(2) - nodeA->giveCoordinate(2);
    length = sqrt(dx * dx + dy * dy);

    answer.at(1, 1) = 1.0;
    answer.at(2, 2) = dx / length;
    answer.at(2, 3) = -dy / length;
    answer.at(3, 2) = -answer.at(2, 3);
    answer.at(3, 3) = answer.at(2, 2);

    return 1;
}

} // end namespace oofem
