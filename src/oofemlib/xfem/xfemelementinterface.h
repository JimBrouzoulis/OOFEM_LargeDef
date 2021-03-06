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

#ifndef xfemelementinterface_h
#define xfemelementinterface_h

#include "interface.h"
#include "alist.h"
#include "xfemmanager.h"

namespace oofem {
class FloatArray;
class FloatMatrix;
class Triangle;
class Element;
class GaussPoint;
class Element;
class XfemManager;

//#define XFEM_DEBUG_VTK 1

/**
 * Provides Xfem interface for an element.
 * @author Erik Svenning
 */
class XfemElementInterface : public Interface
{
protected:
    Element *element;

public:
    /// Constructor.
    XfemElementInterface(Element *e) : Interface() { this->element = e; }

    virtual ~XfemElementInterface() {}

    /// Creates enriched B-matrix.
    void XfemElementInterface_createEnrBmatrixAt(FloatMatrix &oAnswer, GaussPoint &iGP, Element &iEl);

    /// Creates enriched N-matrix.
    void XfemElementInterface_createEnrNmatrixAt(FloatMatrix &oAnswer, const FloatArray &iLocCoord, Element &iEl);

    /// Partitions the element into patches by a triangulation.
    virtual void XfemElementInterface_partitionElement(std::vector< Triangle > &oTriangles, const std :: vector< FloatArray > &iPoints);
    /// Updates integration rule based on the triangulation.
    virtual void XfemElementInterface_updateIntegrationRule();

    /// Helpful routine to put the nodes for triangulation together, should be in protected members probably.
    /// Returns an array of array of points. Each array of points defines the points of a subregion of the element.
    virtual void XfemElementInterface_prepareNodesForDelaunay(std :: vector< std :: vector< FloatArray > > &oPointPartitions);

    /**
     * If the enrichment evolves in time, the element subdivision
     * need to be updated. That is done by recomputeGaussPoints.
     */
    virtual void recomputeGaussPoints();
};
} // end namespace oofem
#endif // xfemelementinterface_h
