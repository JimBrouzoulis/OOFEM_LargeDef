/* $Header: /home/cvs/bp/oofem/oofemlib/src/combuff.h,v 1.5 2003/04/06 14:08:23 bp Exp $ */
/*

                   *****    *****   ******  ******  ***   ***                            
                 **   **  **   **  **      **      ** *** **                             
                **   **  **   **  ****    ****    **  *  **                              
               **   **  **   **  **      **      **     **                               
              **   **  **   **  **      **      **     **                                
              *****    *****   **      ******  **     **         
            
                                                                   
               OOFEM : Object Oriented Finite Element Code                 
                    
                 Copyright (C) 1993 - 2006   Borek Patzak                                       



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
#ifndef loadballancer_h
#ifdef __PARALLEL_MODE
#include "inputrecord.h"
#include "interface.h"
#include "clock.h"

class Domain;
class EngngModel;
class ProcessCommunicator;
#define MIGRATE_LOAD_TAG       9998
#define LOADBALLANCER_END_DATA 9999

/**
   Abstract base class representing general load ballancer monitor. The task of the monitor is to
   detect the imbalance and to make the decision, whether to redistribute the work or to continue 
   with existing partitioning.
 */
class LoadBallancerMonitor
{
 protected:
  EngngModel* emodel;
 public:
  enum LoadBallancerDecisionType {LBD_CONTINUE, LBD_RECOVER};
  
  LoadBallancerMonitor (EngngModel* em) {emodel=em;}
  virtual ~LoadBallancerMonitor() {}

  /**@name Load evaluation and imbalance detection methods*/
  //@{
  virtual LoadBallancerDecisionType decide () = 0;
  //@}
  
};

/**
   Implementation of simple wall-clock based monitor. 
   It detect imbalance based on waal clock diference required for slotion step
   on particular nodes. When difference in wall clock solution times is greater
   than a treshold value, the load migration is performed.
*/
class WallClockLoadBallancerMonitor : public LoadBallancerMonitor
{
 public:
  WallClockLoadBallancerMonitor (EngngModel* em): LoadBallancerMonitor(em) {}
  LoadBallancerDecisionType decide ();
};


/**
   Abstract base class representing general load ballancer. The task of load ballancer is to 
   recover load ballance when running in parallel. This is achieved by moving work from busy 
   nodes to other nodes to achieve an equal distribution of work. 
   In general load ballancer should repartition the problem domain, taking into account several 
   criteria:
   - It should take into account different computational requirement of different elements
   - The new partitioning should minimize the cut (to minimize the communication)
   - The new partitioning should minimize data movement (the cost of repartitioning) by 
   preserving the locality as much as possible. In other words the new and existing partitioning
   should be "similar".
 */
class LoadBallancer
{
 public:
  /**
     Describes the state of dofmanager after load ballancing 
     on the local partition:
     DM_NULL  - undefined (undetermined) state, if assigned means internal error
     DM_Local - local dofman that remains local
     DM_Remote- local dofman that becames remote (becames local on remote partition)
     //DM_SharedNew - local shared that became shared
     DM_Shared- shared dofman that remains shared.
     //DM_SharedExlude - shared dofman that remains shared, 
                       possibly with changed partitions, but local partition
		       is no more in shared list (should be exluded on remote partitions).
     //DM_SharedUpdate - Shared dofman that remains shared, the partition list may changed.
  */
  enum DofManMode { DM_NULL, DM_Local, DM_Shared, DM_Remote};
 protected:
  Domain* domain;

 public:

  LoadBallancer (Domain* d);
  virtual ~LoadBallancer () {}
  
  

 /**@name Work transfer calculation methods  */
 //@{
  virtual void calculateLoadTransfer () = 0;
 //@}
  
 /**@name Work migration methods  */
 //@{
  void migrateLoad (Domain* d);
 //@}
  
 /**@name Query methods after work transfer calculation */
 //@{
  /// Returns the label of dofmanager after load ballancing
  virtual DofManMode giveDofManState (int idofman) = 0;
  
  /// Returns the partition list of given dofmanager after load balancing
  virtual IntArray* giveDofManPartitions (int idofman) = 0;
  
  /// Returns the new partition number assigned to local element after LB
  virtual int giveElementPartition(int ielem) = 0;

 //@}
  ///Initializes receiver acording to object description stored in input record.
  virtual IRResultType initializeFrom (InputRecord* ir) {return IRRT_OK;}

 protected:

  int packMigratingData (Domain*, ProcessCommunicator& pc) ;
  int unpackMigratingData (Domain*, ProcessCommunicator& pc) ;
  void deleteRemoteDofManagers (Domain* );
  void deleteRemoteElements (Domain* );
  

};


class LoadBallancerElementInterface : public Interface
{
 public:
  LoadBallancerElementInterface () {}
  
  virtual double predictRelativeComputationalCost ();
};

#endif
#define loadballancer_h
#endif