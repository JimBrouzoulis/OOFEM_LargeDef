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

#ifndef enrichmentfunction_h
#define enrichmentfunction_h

#include "intarray.h"
#include "classfactory.h"
#include "femcmpnn.h"

#define _IFT_DiscontinuousFunction_Name "discontinuousfunction"
#define _IFT_RampFunction_Name "rampfunction"

namespace oofem {
class EnrichmentItem;
class EnrichmentDomain;
class BasicGeometry;
class GaussPoint;

/**
 * Abstract class representing global shape function
 * Base class declares abstract interface common to all enrichment functions.
 * Particularly, evaluateEnrFuncAt() and evaluateEnrFuncDerivAt()
 * evaluate the value and spatial derivatives at a given point.
 * @author chamrova
 * @author Jim Brouzoulis
 * @author Erik Svenning
 */
class EnrichmentFunction : public FEMComponent
{
public:
    /**
     * Constructor.
     * @param n Number associated with receiver.
     * @param aDomain Reference to domain.
     */
    EnrichmentFunction(int n, Domain *aDomain) : FEMComponent(n, aDomain) { }
    /// Destructor
    virtual ~EnrichmentFunction() { };

    // New interface
    virtual void evaluateEnrFuncAt(double &oEnrFunc, const FloatArray &iPos, const double &iLevelSet, const EnrichmentDomain *ipEnrDom) const = 0;
    virtual void evaluateEnrFuncDerivAt(FloatArray &oEnrFuncDeriv, const FloatArray &iPos, const double &iLevelSet, const FloatArray &iGradLevelSet, const EnrichmentDomain *ipEnrDom) const = 0;

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);
    virtual const char *giveClassName() const { return "EnrichmentFunction"; }
    /// Accessor.
    int giveNumberOfDofs() const { return numberOfDofs; }

protected:

    int numberOfDofs;
};

/** Class representing Heaviside EnrichmentFunction. */
class DiscontinuousFunction : public EnrichmentFunction
{
public:
    DiscontinuousFunction(int n, Domain *aDomain) : EnrichmentFunction(n, aDomain) {
        this->numberOfDofs = 1;
    }

    virtual void evaluateEnrFuncAt(double &oEnrFunc, const FloatArray &iPos, const double &iLevelSet, const EnrichmentDomain *ipEnrDom) const;
    virtual void evaluateEnrFuncDerivAt(FloatArray &oEnrFuncDeriv, const FloatArray &iPos, const double &iLevelSet, const FloatArray &iGradLevelSet, const EnrichmentDomain *ipEnrDom) const;

    virtual const char *giveClassName() const { return "DiscontinuousFunction"; }
    virtual const char *giveInputRecordName() const { return _IFT_DiscontinuousFunction_Name; }
};

/** Class representing the four classical linear elastic branch functions. */
class LinElBranchFunction
{
public:
    LinElBranchFunction() {}
    virtual ~LinElBranchFunction() {}

    virtual void evaluateEnrFuncAt(std :: vector< double > &oEnrFunc, const double &iR, const double &iTheta) const;
    virtual void evaluateEnrFuncDerivAt(std :: vector< FloatArray > &oEnrFuncDeriv, const double &iR, const double &iTheta) const;
};

/** Class representing bimaterial interface. */
class RampFunction : public EnrichmentFunction
{
public:

    RampFunction(int n, Domain *aDomain) : EnrichmentFunction(n, aDomain) {
        this->numberOfDofs = 1;
    }

    virtual void evaluateEnrFuncAt(double &oEnrFunc, const FloatArray &iPos, const double &iLevelSet, const EnrichmentDomain *ipEnrDom) const;
    virtual void evaluateEnrFuncDerivAt(FloatArray &oEnrFuncDeriv, const FloatArray &iPos, const double &iLevelSet, const FloatArray &iGradLevelSet, const EnrichmentDomain *ipEnrDom) const;

    virtual const char *giveClassName() const { return "RampFunction"; }
    virtual const char *giveInputRecordName() const { return _IFT_RampFunction_Name; }
};
} // end namespace oofem
#endif  // enrichmentfunction_h
