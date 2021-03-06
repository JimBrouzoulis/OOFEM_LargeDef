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

#include "gravitypressure.h"
#include "load.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "classfactory.h"

namespace oofem {

REGISTER_BoundaryCondition( GravityPressure );

IRResultType
GravityPressure :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    Load :: initializeFrom(ir);

    normalVector.resize(3);
    normalVector.zero();
    normalVector.at(3) = -1.;
    IR_GIVE_OPTIONAL_FIELD(ir, normalVector, _IFT_GravityPressure_normal);

    zeroLevel = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, zeroLevel, _IFT_GravityPressure_zerolevel);
    return IRRT_OK;
}

void
GravityPressure :: computeValueAt(FloatArray &answer, TimeStep *atTime, FloatArray &coords, ValueModeType mode)
{
    //Need to include the information on the fluid
    //This assumes that the z-direction represents gravity.

    //Use normal input normal vector to set up the local coordinate system
    FloatArray s(3), t(3);

    if ( this->normalVector.at(1) == 0 ) {
        s.at(1) = 0.;
        s.at(2) = this->normalVector.at(3);
        s.at(3) = -this->normalVector.at(2);
    } else if ( this->normalVector.at(2) == 0 ) {
        s.at(1) = this->normalVector.at(3);
        s.at(2) = 0.;
        s.at(3) = -this->normalVector.at(1);
    } else {
        s.at(1) = this->normalVector.at(2);
        s.at(2) = -this->normalVector.at(1);
        s.at(3) = 0.;
    }

    s.normalize();

    t.beVectorProductOf(normalVector, s);
    t.normalize();

    //Locate coordinate matrix
    FloatMatrix lcs(3, 3);
    for ( int i = 1; i <= 3; i++ ) {
        lcs.at(1, i) = normalVector.at(i);
        lcs.at(2, i) = s.at(i);
        lcs.at(3, i) = t.at(i);
    }

    //Express integration point in local coordinate system
    FloatArray coordsLocal(3);
    coordsLocal.beProductOf(lcs, coords);

    //Subtract zero level. Assume that pressure is positive.
    double pressureHead = coordsLocal.at(1);
    computeComponentArrayAt(answer, atTime, mode);
    answer.times(pressureHead);
}
} // end namespace oofem
