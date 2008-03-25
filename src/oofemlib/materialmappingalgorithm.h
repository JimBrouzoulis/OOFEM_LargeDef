/* $Header: /home/cvs/bp/oofem/oofemlib/src/materialmappingalgorithm.h,v 1.10 2003/05/19 13:03:57 bp Exp $ */
/*

                   *****    *****   ******  ******  ***   ***                            
                 **   **  **   **  **      **      ** *** **                             
                **   **  **   **  ****    ****    **  *  **                              
               **   **  **   **  **      **      **     **                               
              **   **  **   **  **      **      **     **                                
              *****    *****   **      ******  **     **         
            
                                                                   
               OOFEM : Object Oriented Finite Element Code                 
                    
                 Copyright (C) 1993 - 2000   Borek Patzak                                       



         Czech Technical University, Faculty of Civil Engineering,
     Department of Structural Mechanics, 166 29 Prague, Czech Republic
                                                                               
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                                                                              
*/

//   **********************************************
//   *** CLASS MATERIAL MODEL MAPPING ALGORITHM ***
//   **********************************************

#ifndef materilmappingalgorithm_h 

#include "compiler.h"
#include "cltypes.h"
#include "interface.h"
#include "intarray.h"
#include "inputrecord.h"

class Domain;
class Element;
class TimeStep;
class FloatArray;
class GaussPoint;

/**
 The class representing the general material model mapping algorithm.
 The basic task is to map the selected internal variable from old mesh to specific integration point
 on new mesh. 
 Unfortunately, there are two principal ways how to implement adaptive mapping:
  1) The material service is called for each IP to perform all mappings. The mapper object is used (created locally 
 or class object) to map all necessary variables.  The problem is that some mappers can be used efficiently if
 they are initialized for specific position and then they can be reused for all internal variables. The others,
 hovewer, are efficient if if they are initialized for certain variable (of smoothing type) and then are used 
 for all possible points. This is the key problem, whether to initialize mapper for position or for variable first. 
 This also prevents to create some general efficient interface.
 2) The other possibility is to make mapper object responsible for mapping variables in all IPs of material model.
 Mapper is responsible for mapping required variables, which are specified. The advantage is, that mapper can decide,
 whether to initialize yourself for certain variable and map all points or initialize for given point and map all 
 necessary variables. Such approach allows efficient handling of the problem, but some additional support is 
 required: there must exists standard way, how to access and update required history variable, elements must provide 
 access to all integration points, and there is another tradeoff - mapper should perform loop only over
 elements belonging to corresponding material model and there is no efficient support for this at this time.

 The idea used here is to use option (1):  material model is responsible. The mappers are all assumed to 
 reguire first initialization for all variables in given point. So the mappers requiring the initialization 
 based on position can ignore given variable list and can setup for given position. The mappers requiring
 variable init must in general support initialization for multiple variables, but even if this operation is global, i.e., 
 recovered vars are computed in whole domain, they should be kept, since other gp's will have typically the same requests.
 This also assomes that all mappers are are stored as class variables of material model in order to reuse 
 their initialization.

*/
class MaterialMappingAlgorithm  {

protected:

public:
 /// Constructor
 MaterialMappingAlgorithm () {}
 /// Destructor
 virtual ~MaterialMappingAlgorithm() {}


 // add general virtual services here!
 /**
  Initializes the receiver state before mapping. The idea is to place some 
  comon global oprerations before mapping particular IP's if necessary.
  Stores Times stamp of last initialization, so multiple calls for same step
  do not initialize receiver again.
  @param dold old domain
  @param varTypes array of InternalStateType values, identifiing all vars to be mapped
  @param tStep time step
  */
 void init (Domain* dold, IntArray& varTypes, GaussPoint* gp, TimeStep* tStep);
 /**
  Initializes the receiver state before mapping. The idea is to place some 
  comon global oprerations before mapping particular IP's if necessary.
  Stores Times stamp of last initialization, so multiple calls for same step
  do not initialize receiver again.
  @param dold old domain
  @param varTypes array of InternalStateType values, identifiing all vars to be mapped
  @param coords coordinates of the receiver point
  @param region if > 0 region id of receiver point,, if < 0 ignore regions.
  @param tStep time step
  */
 virtual void __init (Domain* dold, IntArray& varTypes, FloatArray& coords, int region, TimeStep* tStep) = 0;
 /**
  Finishes the mapping for given time step. Used to perform cleanup. 
  Typically some mappers reguire to compute some global mesh data related to
  current step, which are valid for all IPs - so they are computed only once for
  all IPs. 
  */
 virtual void finish (TimeStep* tStep) = 0;
 /** Maps and update the unknown of given type from 
  old mesh oldd to new mesh to which gp belongs to. The result is stored in answer array.
  @param answer contains result
  @param type determines the type of internal variable
  @param gp Integration point belonging to new domain to which mapping occur
  @param tStep time step
  @return nonzero if o.k.
  */
 virtual int mapVariable (FloatArray& answer, GaussPoint* gp, InternalStateType type, TimeStep* tStep) ;
 /** Maps and update the unknown of given type from 
  old mesh oldd to new mesh to which gp belongs to. The result is stored in answer array.
  @param answer contains result
  @param type determines the type of internal variable
  @param coords coordinates of receiver point to which mapping occur
  @param tStep time step
  @return nonzero if o.k.
  */
 virtual int __mapVariable (FloatArray& answer, FloatArray& coords, InternalStateType type, TimeStep* tStep) = 0;
 /** Initializes receiver acording to object description stored in input record.
  InitString can be imagined as data record in component database
  belonging to receiver. Receiver may use value-name extracting functions 
  to extract particular field from record. 
  @see readInteger, readDouble and similar functions */
 virtual IRResultType initializeFrom (InputRecord* ir) {return IRRT_OK;}
 /// Returns class name of the receiver.
 virtual const char* giveClassName () const  = 0;
 

protected:
};

#define materilmappingalgorithm_h
#endif





