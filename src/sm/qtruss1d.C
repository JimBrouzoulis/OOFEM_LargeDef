/* $Header: /home/cvs/bp/oofem/sm/src/truss1d.C,v 1.6 2003/04/06 14:08:32 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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

//   file Truss1d.C

#include "qtruss1d.h"
#include "fei1dquad.h"
#include "domain.h"
#include "node.h"
#include "material.h"
#include "crosssection.h"
#include "gausspnt.h"
#include "gaussintegrationrule.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "intarray.h"

#include "engngm.h"
#ifndef __MAKEDEPEND
 #include <stdlib.h>
 #include <math.h>
#endif

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {
  FEI1dQuad QTruss1d :: interpolation(1);

QTruss1d :: QTruss1d(int n, Domain *aDomain) : StructuralElement(n, aDomain)
    // Constructor.
{
    numberOfDofMans     = 3;
    length              = 0.;
}


IRResultType
QTruss1d :: initializeFrom(InputRecord *ir)
{
  const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
  IRResultType result;                 // Required by IR_GIVE_FIELD macro
  this->StructuralElement :: initializeFrom(ir);
  IR_GIVE_OPTIONAL_FIELD(ir, numberOfGaussPoints, IFT_QTruss1d_nip, "nip"); // Macro
  
  this->computeGaussPoints();
  return IRRT_OK;
}

void
QTruss1d ::   giveDofManDofIDMask(int inode, EquationID, IntArray &answer) const
// returns DofId mask array for inode element node.
// DofId mask array determines the dof ordering requsted from node.
// DofId mask array contains the DofID constants (defined in cltypes.h)
// describing physical meaning of particular DOFs.
//IntArray* answer = new IntArray (2);
{
  answer.resize(1);
  
  answer.at(1) = D_u;

  return;
}


double
QTruss1d :: computeVolumeAround(GaussPoint *aGaussPoint)
// Returns the length of the receiver. This method is valid only if 1
// Gauss point is used.
{
  double J = this->interpolation.giveTransformationJacobian(* aGaussPoint->giveCoordinates(), FEIElementGeometryWrapper(this), 0.0);
  double weight  = aGaussPoint->giveWeight();
  return  J * weight * this->giveCrossSection()->give(CS_Area);
}

int
QTruss1d :: computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords)
{

  this->interpolation.local2global(answer,lcoords,FEIElementGeometryWrapper(this),0.0);
  return 1;
}


void QTruss1d :: computeGaussPoints()
// Sets up the array of Gauss Points of the receiver.XF
{
    if ( !integrationRulesArray ) 
      {
        numberOfIntegrationRules = 1;
        integrationRulesArray = new IntegrationRule * [ 1 ];
        integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 3);
        integrationRulesArray [ 0 ]->setUpIntegrationPoints(_Line, numberOfGaussPoints, _1dMat);
      }
}

void
QTruss1d :: computeNmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer)
// Returns the displacement interpolation matrix {N} of the receiver, eva-
// luated at aGaussPoint.
{
  FloatArray n;
  answer.resize(1,3);
  answer.zero();

  this->interpolation.evalN(n, * aGaussPoint->giveCoordinates(), FEIElementGeometryWrapper(this), 0.0);
  answer.at(1,1) = n.at(1);
  answer.at(1,2) = n.at(2);
  answer.at(1,3) = n.at(3);
  
  return;
}

void
QTruss1d :: computeBmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer, int li, int ui)
//
// Returns linear part of geometrical equations of the receiver at gp.
// Returns the linear part of the B matrix
//
{
  answer.resize(1,3);
  answer.zero();
  
  this->interpolation.evaldNdx(answer, * aGaussPoint->giveCoordinates(), FEIElementGeometryWrapper(this), 0.0);
  
  return;
}




  double
QTruss1d :: giveLength()
  // Returns the length of the receiver.
{
  length = interpolation.computeLength(FEIElementGeometryWrapper(this));
  return length;
  }




int
QTruss1d :: computeLocalCoordinates(FloatArray &answer, const FloatArray &gcoords)
{  
  int ok;
   ok = interpolation.global2local(answer,gcoords, FEIElementGeometryWrapper(this),0);
  return ok;
}

} // end namespace oofem