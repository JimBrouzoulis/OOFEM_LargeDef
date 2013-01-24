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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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

#ifndef Tr2Shell7_h
#define Tr2Shell7_h


#include "eleminterpmapperinterface.h"
#include "nodalaveragingrecoverymodel.h"
#include "layeredcrosssection.h"

#include "nlstructuralelement.h"
#include "shell7base.h"

namespace oofem {

class FEI3dTrQuad;
class BoundaryLoad;

/**
 * This class represent a 7 parameter shell element. 
 * Each node has 7 degrees of freedom (displ. vec., director vec., inhomogeneous thickness strain ).
 * Add ref. to paper!
 * @author Jim Brouzoulis
 * @date 2012-11-01
 */

class Tr2Shell7 : public Shell7Base
{
protected:
    int numberOfGaussPoints;	
    static FEI3dTrQuad interpolation;
    static bool __initialized;
	static IntArray ordering_phibar;
	static IntArray ordering_m;
	static IntArray ordering_gam;
    static IntArray ordering_all;

    static bool initOrdering() {
        ordering_phibar.setValues(18, 1, 2, 3, 8, 9, 10, 15, 16, 17, 22, 23, 24, 29, 30 ,31, 36, 37, 38);
		ordering_m.setValues(18, 4, 5, 6, 11, 12, 13, 18, 19, 20, 25, 26, 27, 32, 33 ,34, 39, 40, 41);
		ordering_gam.setValues(6, 7, 14, 21, 28, 35, 42);
        ordering_all.setValues(42, 1, 2, 3, 8, 9, 10, 15, 16, 17, 22, 23, 24, 29, 30 ,31, 36, 37, 38,
                                   4, 5, 6, 11, 12, 13, 18, 19, 20, 25, 26, 27, 32, 33 ,34, 39, 40, 41,
                                   7, 14, 21, 28, 35, 42);
        return true;
    }

	//specific
    void giveSurfaceDofMapping(IntArray &answer, int iSurf) const;
	void giveEdgeDofMapping(IntArray &answer, int iEdge) const;

    virtual double computeVolumeAround(GaussPoint *gp);
    virtual double computeVolumeAroundLayer(GaussPoint *mastergp, int layer);
    virtual double computeAreaAround(GaussPoint *gp);

    virtual void computeGaussPoints();
    virtual void giveLocalNodeCoords(FloatArray &nodeLocalXiCoords, FloatArray &nodeLocalEtaCoords);

	//only used for debuging 
    void compareMatrices(const FloatMatrix &matrix1, const FloatMatrix &matrix2, FloatMatrix &answer);

    virtual FEInterpolation *giveInterpolation();

    
public:
    Tr2Shell7(int n, Domain *d);	// constructor
    virtual ~Tr2Shell7() { }		// destructor -> declaring as virtual will make each subclass call their respective destr.
    // definition & identification
    virtual int giveNumberOfDofs()           { return 42; }
    virtual int giveNumberOfEdgeDofs()       { return 21; }
    virtual int giveNumberOfEdgeDofManagers(){ return 3;  }
    virtual const char *giveClassName()                const { return "Tr2Shell7"; }
    virtual classType giveClassID()                    const { return Tr2Shell7Class; }
    virtual Element_Geometry_Type giveGeometryType()   const { return EGT_triangle_2; }
	virtual integrationDomain  giveIntegrationDomain() const { return _Triangle; } // write new wedge-like type 'layeredWedge'


	virtual void printOutputAt(FILE *file, TimeStep *tStep);

};



} // end namespace oofem
#endif 
