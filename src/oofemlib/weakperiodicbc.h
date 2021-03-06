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

#ifndef WEAKPERIODICBC_H_
#define WEAKPERIODICBC_H_

#include <vector>
#include <iostream>

#include "activebc.h"
#include "inputrecord.h"
#include "gaussintegrationrule.h"

///@name Input fields for WeakPeriodicBoundaryCondition
//@{
#define _IFT_WeakPeriodicBoundaryCondition_Name "weakperiodicbc"
#define _IFT_WeakPeriodicBoundaryCondition_order "order"
#define _IFT_WeakPeriodicBoundaryCondition_descritizationType "descritizationtype"
#define _IFT_WeakPeriodicBoundaryCondition_dofid "dofid"
#define _IFT_WeakPeriodicBoundaryCondition_ngp "ngp"
#define _IFT_WeakPeriodicBoundaryCondition_elementSidesPositive "elementsidespositive"
#define _IFT_WeakPeriodicBoundaryCondition_elementSidesNegative "elementsidesnegative"
#define _IFT_WeakPeriodicBoundaryCondition_elementSidesPositiveSet "elementsidespositiveset"
#define _IFT_WeakPeriodicBoundaryCondition_elementSidesNegativeSet "elementsidesnegativeset"
//@}

namespace oofem {
enum basisType { monomial=0, trigonometric=1, legendre=2 };
/**
 * Imposes weak periodicity on the doftype of choice. 2D. It is required that the two edges are parallel and either horizontal or vertical.
 *
 * @author Carl Sandström
 */
class WeakPeriodicBoundaryCondition : public ActiveBoundaryCondition
{
private:

    basisType useBasisType;
    int bcID;
    int orderOfPolygon;

    /** Direction of normal. 1 if normal in x, 2 if y and 3 if z. */
    int direction;

    /** Keeps info on which coordinates varies over the surface. Depends on number of spatial dimensions and normal direction */
    IntArray surfaceIndexes;

    FloatArray smax, smin;
    bool doUpdateSminmax;

    /** Number of Gausspoints used when integrating along the element edges */
    int ngp;

    /** Number of degrees of freedom */
    int ndof;

    /** Set containing positive side */
    int posSet;

    /** Set containing negative side */
    int negSet;

    /** ID of dofs on which weak periodicity is imposed */
    int dofid;

    /** sideSign is the sign of the normal for each side */
    signed int sideSign [ 2 ];

    /** side[] keeps track of which side of the triangle is located along the boundary. element[] keeps track of what element is located along the boundary */
    std :: vector< int >side [ 2 ], element [ 2 ];

    /** Keeps track of which coordinate(s) are changing on the surface/edge */
    std :: vector< double> directions;

    /** Type of surface/edge */
    integrationDomain sideGeom;

    void giveEdgeNormal(FloatArray &answer, int element, int side);

    void updateSminmax();

    void updateDirection();

    double computeBaseFunctionValue(int baseID, double coordinate);

    Node *gammaDman;

    double factorial(int n);

    double binomial(double n , int k);

    void getExponents(int n, int &i, int &j);
public:
    WeakPeriodicBoundaryCondition(int n, Domain *d);
    virtual ~WeakPeriodicBoundaryCondition() { };

    virtual IRResultType initializeFrom(InputRecord *ir);

    basisType giveBasisType() {return useBasisType; };

    virtual void assemble(SparseMtrx *answer, TimeStep *tStep, EquationID eid, CharType type, 
                          const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s);

    virtual void assembleVector(FloatArray &answer, TimeStep *tStep, EquationID eid,
                                  CharType type, ValueModeType mode,
                                  const UnknownNumberingScheme &s, FloatArray *eNorm = NULL);

    virtual int giveNumberOfInternalDofManagers();

    virtual DofManager *giveInternalDofManager(int i);

    virtual void addElementSide(int elem, int side);

    virtual const char *giveClassName() const { return "WeakPeriodicBoundaryCondition"; }
    virtual const char *giveInputRecordName() const { return _IFT_WeakPeriodicBoundaryCondition_Name; }

protected:
    void computeElementTangent(FloatMatrix &answer, Element *e, int boundary);
};
}
#endif /* WEAKPERIODICBC_H_ */
