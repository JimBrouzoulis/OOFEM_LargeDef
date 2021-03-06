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

#ifndef surfacetensionbc_h
#define surfacetensionbc_h

#include "activebc.h"

#include <utility>
#include <list>

///@name Input fields for surface tension boundary condition
//@{
#define _IFT_SurfaceTensionBoundaryCondition_Name "surfacetension"
#define _IFT_SurfaceTensionBoundaryCondition_gamma "gamma"
#define _IFT_SurfaceTensionBoundaryCondition_useTangent "usetangent"
//@}

namespace oofem {
class Element;
class ElementSide;

/**
 * Computes the load (and possibly tangent) for surface tension. This boundary condition applicable to both solid and flow problems.
 * Computing this load and tangent is involves some fairly complicated expression with lengthy derivations, especially for 3D surfaces.
 */
class SurfaceTensionBoundaryCondition : public ActiveBoundaryCondition
{
    double gamma; ///< Surface tension.
    bool useTangent; ///< Determines if tangent should be used.

public:
    /**
     * Constructor. Creates boundary an active condition with given number, belonging to given domain.
     * @param n Boundary condition number.
     * @param d Domain to which new object will belongs.
     */
    SurfaceTensionBoundaryCondition(int n, Domain *d) : ActiveBoundaryCondition(n, d) { }
    /// Destructor.
    virtual ~SurfaceTensionBoundaryCondition() { }

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual void assemble(SparseMtrx *answer, TimeStep *tStep, EquationID eid,
                          CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s);

    virtual void assembleVector(FloatArray &answer, TimeStep *tStep, EquationID eid,
                                  CharType type, ValueModeType mode,
                                  const UnknownNumberingScheme &s, FloatArray *eNorms = NULL);

    virtual void giveLocationArrays(std::vector<IntArray> &rows, std::vector<IntArray> &cols, EquationID eid, CharType type,
                                    const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s);

    virtual classType giveClassID() const { return SurfaceTensionBoundaryConditionClass; }
    virtual const char *giveClassName() const { return "SurfaceTensionBoundaryCondition"; }
    virtual const char *giveInputRecordName() const { return _IFT_SurfaceTensionBoundaryCondition_Name; }

protected:
    /**
     * Helper function for computing the contributions to the load vector.
     */
    void computeLoadVectorFromElement(FloatArray &answer, Element *e, int side, TimeStep *tStep);
    /**
     * Helper function for computing the tangent (@f$ K = \frac{\mathrm{d}F}{\mathrm{d}u} @f$)
     */
    void computeTangentFromElement(FloatMatrix &answer, Element *e, int side, TimeStep *tStep);
};
} // end namespace oofem
#endif // surfacetensionbc_h
