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

#ifndef qspace_h
#define qspace_h

#include "nlstructuralelement.h"
#include "fei3dhexaquad.h"
#include "zznodalrecoverymodel.h"
#include "nodalaveragingrecoverymodel.h"
#include "eleminterpmapperinterface.h"
#include "huertaerrorestimator.h"
#include "sprnodalrecoverymodel.h"

#define _IFT_QSpace_Name "qspace"

namespace oofem {
/**
 * This class implements an Quadratic 3d  20 - node element. Each node has 3 degrees of freedom.
 *
 * One single additional attribute is needed for Gauss integration purpose :
 * 'jacobianMatrix'. This 3x3 matrix contains polynomials.
 * Tasks:
 * - calculating its Gauss points
 * - calculating its B,D,N matrices and dV
 *
 * @author Ladislav Svoboda
 */
class QSpace : public NLStructuralElement, public SPRNodalRecoveryModelInterface, public ZZNodalRecoveryModelInterface, public NodalAveragingRecoveryModelInterface
{
protected:
    static FEI3dHexaQuad interpolation;

public:
    QSpace(int n, Domain *d);
    virtual ~QSpace() {}

    virtual FEInterpolation *giveInterpolation() const { return &interpolation; }

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const;
    virtual double computeVolumeAround(GaussPoint *);

    virtual double giveCharacteristicLenght(GaussPoint *gp, const FloatArray &normalToCrackPlane);

    virtual Interface *giveInterface(InterfaceType);
    virtual int testElementExtension(ElementExtension ext) { return ( ( ext == Element_SurfaceLoadSupport ) ? 1 : 0 ); }

    virtual Element *ZZNodalRecoveryMI_giveElement() { return this; }

    virtual void SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap);
    virtual void SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap);
    virtual int SPRNodalRecoveryMI_giveNumberOfIP();
    virtual SPRPatchType SPRNodalRecoveryMI_givePatchType();

    virtual void NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node, InternalStateType type, TimeStep *tStep);
    virtual void NodalAveragingRecoveryMI_computeSideValue(FloatArray &answer, int side, InternalStateType type, TimeStep *tStep);

    // definition & identification
    virtual const char *giveInputRecordName() const { return _IFT_QSpace_Name; }
    virtual const char *giveClassName() const { return "QSpace"; }
    virtual classType giveClassID() const { return QSpaceClass; }
    virtual int computeNumberOfDofs(EquationID ut) { return 60; }
    virtual MaterialMode giveMaterialMode();

protected:
    virtual void computeGaussPoints();
    virtual void computeNmatrixAt(GaussPoint *, FloatMatrix &);
    
    virtual void computeBmatrixAt(GaussPoint *, FloatMatrix &, int = 1, int = ALL_STRAINS);
    //void computeBFmatrixAt(GaussPoint *, FloatMatrix &);
    void computeBHmatrixAt(GaussPoint *, FloatMatrix &);

    virtual int giveApproxOrder() { return 2; }
    virtual int giveNumberOfIPForMassMtrxIntegration() { return 27; }

    /**
     * @name Surface load support
     */
    //@{
    virtual IntegrationRule *GetSurfaceIntegrationRule(int);
    virtual void computeSurfaceNMatrixAt(FloatMatrix &answer, int iSurf, GaussPoint *gp);
    virtual void giveSurfaceDofMapping(IntArray &answer, int) const;
    virtual double computeSurfaceVolumeAround(GaussPoint *gp, int);
    virtual void computeSurfIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int);
    virtual int computeLoadLSToLRotationMatrix(FloatMatrix &answer, int, GaussPoint *gp);
    //@}
};
} // end namespace oofem
#endif
