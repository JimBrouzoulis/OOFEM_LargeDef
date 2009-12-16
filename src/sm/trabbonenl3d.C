/* $Header: /home/cvs/bp/oofem/sm/src/idmnl1.C,v 1.7 2003/04/06 14:08:30 bp Exp $ */
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


#include "trabbonenl3d.h"
#include "gausspnt.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "structuralcrosssection.h"
#include "mathfem.h"

#include "sparsemtrx.h"
#include "isolinearelasticmaterial.h"
#include "dynalist.h"
#include "error.h"
#include "nonlocalmaterialext.h"
#ifndef __MAKEDEPEND
#include <math.h>
#endif

#ifdef __PARALLEL_MODE
#include "idmnl1.h"
#include "combuff.h"
#endif

#ifdef __OOFEG
#include "oofeggraphiccontext.h"
#include "conTable.h"
#endif


/////////////////////////////////////////////////////////////////
////////////TRABECULAR BONE NONLOCAL MATERIAL////////////////////
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: CONSTRUCTOR
//

TrabBoneNL3D::TrabBoneNL3D(int n, Domain *d):TrabBone3D(n,d), StructuralNonlocalMaterialExtensionInterface(d), NonlocalMaterialStiffnessInterface()
{
  R = 0.;
}

//
// END: CONSTRUCTOR
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: DESTRUCTOR
//

TrabBoneNL3D::~TrabBoneNL3D()
{
}

//
// END: DESTRUCTOR
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: SUBROUTINE FOR UPDATE BEFORE NON LOCAL AVERAGE
// update local values of accumulated pastic strain

void 
TrabBoneNL3D::updateBeforeNonlocAverage(const FloatArray& strainVector, GaussPoint* gp, TimeStep* atTime)
{
  FloatArray SDstrainVector, fullSDStrainVector;
  double cumPlastStrain;
  TrabBoneNL3DStatus *nlStatus = (TrabBoneNL3DStatus*) this -> giveStatus (gp);
 
  this->initTempStatus(gp);
  this->initGpForNewStep(gp);
  this->giveStressDependentPartOfStrainVector(SDstrainVector, gp, strainVector,atTime, VM_Total);

  nlStatus->letTempStrainVectorBe(strainVector); 

  StrainVector strain(strainVector, gp->giveMaterialMode());
 
  this->performPlasticityReturn(gp, strain);
  this->computeLocalCumPlastStrain (cumPlastStrain, strain, gp, atTime);
  nlStatus->setLocalCumPlastStrainForAverage (cumPlastStrain);
}

//
// END: SUBROUTINE FOR UPDATE BEFORE NON LOCAL AVERAGE
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: EVALUATION OF LOCAL STIFFNESS MATRIX
//

void
TrabBoneNL3D::give3dMaterialStiffnessMatrix (FloatMatrix& answer,
                                   MatResponseForm form,MatResponseMode mode,GaussPoint* gp,
                                   TimeStep* atTime)
{
  TrabBoneNL3DStatus *nlStatus = (TrabBoneNL3DStatus*) this -> giveStatus (gp);

  double tempDam, beta, deltaKappa, nlKappa;
  FloatArray tempEffectiveStress, tempTensor2, prodTensor, plasFlowDirec;
  FloatMatrix elasticity, compliance, SSaTensor, secondTerm, thirdTerm, tangentMatrix;

  if (mode == ElasticStiffness)
  {
    this->constructAnisoComplTensor(compliance, m1, m2, (3.0-(m1+m2)), rho, eps0, nu0, mu0, expk, expl);
    elasticity.beInverseOf(compliance);

    answer = elasticity;
  }
  else if (mode == SecantStiffness)
  {
    this->constructAnisoComplTensor(compliance, m1, m2, (3.0-(m1+m2)), rho, eps0, nu0, mu0, expk, expl);
    elasticity.beInverseOf(compliance);
    tempDam = nlStatus -> giveTempDam();

    answer = elasticity;
    answer.times(1.0-tempDam);
  }
  else if (mode == TangentStiffness)
  {
    deltaKappa = nlStatus -> giveDeltaKappa();
    if (deltaKappa > 0.0)
    {
// Imports
      tempEffectiveStress = *nlStatus -> giveTempEffectiveStress();
      this->computeCumPlastStrain (nlKappa, gp, atTime);
      tempDam = nlStatus -> giveTempDam();
      plasFlowDirec = *nlStatus -> givePlasFlowDirec();
      SSaTensor = *nlStatus -> giveSSaTensor();
      beta = nlStatus -> giveBeta();
// Construction of the dyadic product tensor
      prodTensor.beTProductOf(SSaTensor,plasFlowDirec);
// Construction of the tangent stiffness third term
      thirdTerm.beDyadicProductOf(tempEffectiveStress,prodTensor);
      thirdTerm.times(-expDam*critDam*exp(-expDam*nlKappa)*(1.0-mParam)/beta);
// Construction of the tangent stiffness second term
      tempTensor2.beProductOf(SSaTensor,plasFlowDirec);
      secondTerm.beDyadicProductOf(tempTensor2,prodTensor);
      secondTerm.times(-(1.0-tempDam)/beta);
// Construction of the tangent stiffness
      tangentMatrix = SSaTensor;
      tangentMatrix.times(1.0-tempDam);
      tangentMatrix.plus(secondTerm);
      tangentMatrix.plus(thirdTerm);

      answer = tangentMatrix;

    }
    else
    {
//  printf("\n Element %d, GP %d", gp->giveElement()->giveNumber(), gp->giveNumber());
// Import of state variables
      tempDam = nlStatus -> giveTempDam();
// Construction of the tangent stiffness
      this->constructAnisoComplTensor(compliance, m1, m2, (3.0-(m1+m2)), rho, eps0, nu0, mu0, expk, expl);
      elasticity.beInverseOf(compliance);

      answer = elasticity;
      answer.times(1.0-tempDam);
    }
  }

  nlStatus -> setSmtrx(answer);

}

//
// END: EVALUATION OF LOCAL STIFFNESS MATRIX
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: EVALUATION OF NONLOCAL CONTRIBUTION TO STIFFNESS MATRIX
//

void 
TrabBoneNL3D::NonlocalMaterialStiffnessInterface_addIPContribution(SparseMtrx& dest, GaussPoint* gp, TimeStep* atTime)
{
  TrabBoneNL3DStatus *nlStatus = (TrabBoneNL3DStatus*) this -> giveStatus (gp);
  dynaList<localIntegrationRecord> *list = nlStatus->giveIntegrationDomainList();
  dynaList<localIntegrationRecord>::iterator pos;
  TrabBoneNL3D* rmat;

  double coeff;
  FloatArray rcontrib, lcontrib;
  IntArray   loc, rloc;

  FloatMatrix contrib;

  if (this->giveLocalNonlocalStiffnessContribution (gp, loc, lcontrib, atTime) == 0) return;

  for (pos = list->begin(); pos!= list->end(); ++pos)
  {
    rmat = (TrabBoneNL3D*) ((*pos).nearGp)->giveMaterial();
    if (rmat->giveClassID() == this->giveClassID())
    {
      rmat->giveRemoteNonlocalStiffnessContribution ((*pos).nearGp, rloc, rcontrib, atTime);
      coeff = gp->giveElement()->computeVolumeAround(gp) * (*pos).weight / nlStatus->giveIntegrationScale();

////       printf ("\nLocal Vector Contribution Element %d (GP %d )", gp->giveElement()->giveNumber(), gp->giveNumber());
////       lcontrib.printYourself();
////       printf ("Remote Vector Contribution Element %d (GP %d ), weight factor %f \n", (*pos).nearGp->giveElement()->giveNumber(), (*pos).nearGp->giveNumber(), coeff);
////       rcontrib.printYourself();

      int i,j,dim1 = loc.giveSize(), dim2 = rloc.giveSize();
      contrib.resize (dim1, dim2);
      for (i=1; i<=dim1; i++)
      {
        for (j=1; j<=dim2; j++)
        {
          contrib.at(i,j) = -1.0*lcontrib.at(i)*rcontrib.at(j)*coeff;
        }
      }

      dest.assemble (loc, rloc, contrib);

////       printf ("\nLocal Element %d (GP %d )", gp->giveElement()->giveNumber(), gp->giveNumber());
////       printf (", Remote Element %d (GP %d ), weight factor %f \n", (*pos).nearGp->giveElement()->giveNumber(), (*pos).nearGp->giveNumber(), coeff);
////       printf ("Remote Contribution Output\n");
////       contrib.printYourself();

    }
  }
}

dynaList<localIntegrationRecord>* 
TrabBoneNL3D::NonlocalMaterialStiffnessInterface_giveIntegrationDomainList(GaussPoint* gp)
{
  TrabBoneNL3DStatus *nlStatus = (TrabBoneNL3DStatus*) this -> giveStatus (gp);
  this->buildNonlocalPointTable (gp);
  return nlStatus->giveIntegrationDomainList();
}

int
TrabBoneNL3D::giveLocalNonlocalStiffnessContribution (GaussPoint* gp, IntArray& loc, FloatArray& lcontrib, TimeStep* atTime)
{
  TrabBoneNL3DStatus *nlStatus = (TrabBoneNL3DStatus*) this -> giveStatus (gp);
  StructuralElement* elem = (StructuralElement*)(gp->giveElement());

  int nrows, nsize, i, j;
  double sum, nlKappa, dDamFunc, dam, tempDam;
  FloatArray localNu;
  FloatMatrix b;

  this->computeCumPlastStrain (nlKappa, gp, atTime);
  dam = nlStatus -> giveDam();
  tempDam = nlStatus -> giveTempDam();

//  if ((gp->giveElement()->giveNumber()==1) && (gp->giveNumber()==1))
//  {
//    printf("DeltaDam: %10.16e , deltaKappa: %10.16e \n", tempDam - dam, deltaKappa);
//  }

  if ((tempDam - dam) > 0.0)
  {
    localNu = *nlStatus -> giveTempEffectiveStress();

    elem -> giveLocationArray (loc, EID_MomentumBalance, EModelDefaultEquationNumbering());
    elem -> computeBmatrixAt (gp, b);
    dDamFunc = expDam*critDam*exp(-expDam*nlKappa);

    nrows = b.giveNumberOfColumns();
    nsize = localNu.giveSize();
    lcontrib.resize(nrows);
    for (i=1; i<= nrows; i++)
    {
      sum = 0.0;
      for (j=1; j<=nsize; j++) sum+= b.at(j,i)*localNu.at(j);
      lcontrib.at(i) = dDamFunc*mParam*sum;
    }
  }
  else
  {
    loc.resize(0);
    return 0;
  }
 return 1;
}

void
TrabBoneNL3D::giveRemoteNonlocalStiffnessContribution (GaussPoint* gp, IntArray& rloc, FloatArray& rcontrib, TimeStep* atTime)
{
  TrabBoneNL3DStatus *nlStatus = (TrabBoneNL3DStatus*) this -> giveStatus (gp);
  StructuralElement* elem = (StructuralElement*) (gp->giveElement());

  int ncols, nsize, i,j;
  double sum, deltaKappa, beta;
  FloatArray remoteNu, plasFlowDirec, prodTensor;
  FloatMatrix b, SSaTensor;

  deltaKappa = nlStatus -> giveDeltaKappa();

  elem -> giveLocationArray(rloc, EID_MomentumBalance, EModelDefaultEquationNumbering());
  elem -> computeBmatrixAt(gp, b);

  ncols = b.giveNumberOfColumns();
  rcontrib.resize(ncols);

  if (deltaKappa > 0.0)
  {
    plasFlowDirec = *nlStatus -> givePlasFlowDirec();
    SSaTensor = *nlStatus -> giveSSaTensor();
    beta = nlStatus -> giveBeta();

    prodTensor.beTProductOf(SSaTensor,plasFlowDirec);
    remoteNu = 1/beta*prodTensor;
    nsize = remoteNu.giveSize();

    for (i=1; i<=ncols; i++)
    {
      sum = 0.0;
      for (j=1; j<=nsize; j++) sum += remoteNu.at(j)*b.at(j, i);
      rcontrib.at(i) = sum;
    }
  }
  else
  {
    for (i=1; i<=ncols; i++) rcontrib.at(i) = 0.;
  }

}

//
// END: EVALUATION OF NONLOCAL CONTRIBUTION TO STIFFNESS MATRIX
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: SUBROUTINE FOR EVALUATION OF TOTAL STRESS
//

void  
TrabBoneNL3D::giveRealStressVector(FloatArray& answer, MatResponseForm form, GaussPoint* gp,
                                  const FloatArray& totalStrain, TimeStep* atTime)
{
  TrabBoneNL3DStatus *nlStatus = (TrabBoneNL3DStatus*) this -> giveStatus (gp);
  this->initGpForNewStep(gp);

  int i;
  double tempDam;
  FloatArray effStress, totalStress, densStress;

  performPlasticityReturn(gp, totalStrain);
  tempDam = computeDamage(gp, atTime);
  effStress = *nlStatus -> giveTempEffectiveStress();

  totalStress = (1-tempDam)*effStress;

  for(i=1;i<=6;i++){
     if (sqrt(totalStress.at(i)*totalStress.at(i))< pow(10,-8))
	totalStress.at(i)=0.;
    }

  computePlasStrainEnerDensity (gp, totalStrain, totalStress);

  if (JCrit != 0.) computeDensificationStress(densStress, gp, totalStrain, atTime);
  else densStress.resize(6);

  answer = totalStress + densStress;

  nlStatus -> setTempDam(tempDam);
  nlStatus -> letTempStrainVectorBe(totalStrain);
  nlStatus -> letTempStressVectorBe(answer);

////   double tempKappa = nlStatus -> giveTempKappa(), nlKappa;
////   this->computeCumPlastStrain (nlKappa, gp, atTime);
////   double deltaKappa = nlStatus -> giveDeltaKappa();
////   printf("\nInternal Variable, Element %d (GP %d )\n", gp->giveElement()->giveNumber(), gp->giveNumber());
////   printf("deltaKappa %10.16e , kappa %10.16e , damage %10.16e , nlKappa %10.16e \n", deltaKappa, tempKappa, tempDam, nlKappa);
////   printf("Plastic Stress");
////   effStress.printYourself();

  return ;
}

//
// END: SUBROUTINE FOR EVALUATION OF TOTAL STRESS
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: SUBROUTINE OF NONLOCAL ALPHA EVALUATION
//

void 
TrabBoneNL3D::computeCumPlastStrain (double& kappa, GaussPoint* gp, TimeStep* atTime)
{
  double nonlocalContribution, nonlocalCumPlastStrain = 0.0, coeff =0.0;
  TrabBoneNL3DStatus *nonlocStatus, *nlStatus = (TrabBoneNL3DStatus*) this -> giveStatus (gp);

  this->buildNonlocalPointTable(gp);
  this->updateDomainBeforeNonlocAverage(atTime);
  
  dynaList<localIntegrationRecord>* list = nlStatus->giveIntegrationDomainList();
  dynaList<localIntegrationRecord>::iterator pos;

  for (pos = list->begin(); pos!= list->end(); ++pos)
  {

////     printf ("\nStart Element %d (GP %d )", gp->giveElement()->giveNumber(), gp->giveNumber());
    coeff = gp->giveElement()->computeVolumeAround(gp) * (*pos).weight / nlStatus->giveIntegrationScale();
////     printf (", Remote Element %d (GP %d ), weight factor %f \n", (*pos).nearGp->giveElement()->giveNumber(), (*pos).nearGp->giveNumber(), coeff);

    nonlocStatus = (TrabBoneNL3DStatus*) this -> giveStatus ((*pos).nearGp);
    nonlocalContribution = nonlocStatus->giveLocalCumPlastStrainForAverage();
    nonlocalContribution *= (*pos).weight;
    nonlocalCumPlastStrain += nonlocalContribution;
  }

  nonlocalCumPlastStrain *= 1./nlStatus->giveIntegrationScale();

  double localCumPlastStrain = nlStatus->giveLocalCumPlastStrainForAverage();
  kappa = mParam*nonlocalCumPlastStrain +(1-mParam)*localCumPlastStrain ;
}

//
// END: SUBROUTINE OF NONLOCAL ALPHA EVALUATION
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: INTERFACE ????
//

Interface* 
TrabBoneNL3D::giveInterface (InterfaceType type)
{
 if (type == NonlocalMaterialExtensionInterfaceType) return (StructuralNonlocalMaterialExtensionInterface*) this;
 else if (type == NonlocalMaterialStiffnessInterfaceType) return (NonlocalMaterialStiffnessInterface*) this;
 //
 else return NULL;
}

//
// END: INTERFACE ????
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: PARAMETERS OF INPUT FILE
//

IRResultType
TrabBoneNL3D::initializeFrom (InputRecord* ir)
{

  const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
  IRResultType result;                               // Required by IR_GIVE_FIELD macro
 
  TrabBone3D::initializeFrom (ir);
  StructuralNonlocalMaterialExtensionInterface::initializeFrom (ir);

  IR_GIVE_FIELD (ir, R, IFT_TrabBoneNL3D_r, "r"); // Macro
  if (R < 0.0) R = 0.0;

  mParam = 2.;
  IR_GIVE_OPTIONAL_FIELD (ir, mParam, IFT_TrabBoneNL3D_m, "m"); // Macro
 
  return IRRT_OK;
}

//
// END: PARAMETERS OF INPUT FILE
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: ??????????????
//

int
TrabBoneNL3D::giveInputRecordString(std::string &str, bool keyword)
{
  char buff[1024];

  TrabBone3D::giveInputRecordString(str, keyword);
  StructuralNonlocalMaterialExtensionInterface::giveInputRecordString(str, false);
  sprintf(buff, " r %e", this -> R);
  str += buff;

  return 1;
}

//
// END: ????????????????
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: WEIGHT FUNCTION FOR NONLOCAL CONTRIBUTION 
//

double 
TrabBoneNL3D::computeWeightFunction (const FloatArray& src, const FloatArray& coord)
{

  double dist = src.distance (coord);

  if ((dist >= 0.) && (dist <= this-> R))
  {
    double help = (1.-dist*dist/(R*R));
    return help*help;
  }
  return 0.0;
}

//
// END: WEIGHT FUNCTION FOR NONLOCAL CONTRIBUTION
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
/////////TRABECULAR BONE NONLOCAL MATERIAL STATUS////////////////
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: CONSTRUCTOR
// init state variables

TrabBoneNL3DStatus::TrabBoneNL3DStatus (int n, Domain*d, GaussPoint *g) 
  : TrabBone3DStatus(n,d,g), StructuralNonlocalMaterialStatusExtensionInterface()
{
 localCumPlastStrainForAverage = 0.0;
}

//
// END: CONSTRUCTOR
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: DESTRUCTOR
//

TrabBoneNL3DStatus::~TrabBoneNL3DStatus ()
{
}

//
// END: DESTRUCTOR
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: PRINTOUT
//

void 
TrabBoneNL3DStatus::printOutputAt  (FILE *file, TimeStep* tStep)
{
 StructuralMaterialStatus :: printOutputAt (file, tStep);
  fprintf (file,"status {");
  fprintf (file," plastrains ") ;
  int n = plasDef.giveSize() ;
  for (int i=1 ; i<=n ; i++) 
    fprintf (file," % .4e", plasDef.at(i)) ;
  fprintf (file, " ,");
  fprintf ( file," kappa % .4e ,", kappa ) ;
  fprintf ( file," dam  % .4e ,", tempDam ) ;
  fprintf ( file," esed  % .4e ,", this->tempTSED - this->tempPSED ) ;
  fprintf ( file," psed  % .4e ,", this->tempPSED ) ;
  fprintf ( file," tsed  % .4e", this->tempTSED ) ;
  fprintf (file,"}\n");
}

//
// END: PRINTOUT
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: INITIALIZE TEMP VARIABLE (UPDATED DURING ITERATIONS)
// initialize temporary state variables according to equilibriated state vars

void 
TrabBoneNL3DStatus::initTempStatus ()
{
 TrabBone3DStatus::initTempStatus();
}

//
// END: INITIALIZE TEMP VARIABLE (UPDATED DURING ITERATIONS)
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: SETS VARIABLE EQUAL TO TEMP VARIABLE AT THE END OF THE STEP
// Called when equlibrium reached, set equilibriated vars according to temporary (working) ones.

void 
TrabBoneNL3DStatus::updateYourself(TimeStep* atTime)
{
 TrabBone3DStatus::updateYourself(atTime);
}

//
// END: SETS VARIABLE EQUAL TO TEMP VARIABLE AT THE END OF THE STE
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: INTERFACE
//

Interface* 
TrabBoneNL3DStatus::giveInterface (InterfaceType type)
{
  if (type == NonlocalMaterialStatusExtensionInterfaceType) return this;
  else return NULL;
}

//
// END: INTERFACE ???????????????
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: PARALLEL MODE OPTION
//

#ifdef __PARALLEL_MODE
int
TrabBoneNL3D::packUnknowns (CommunicationBuffer& buff, TimeStep* stepN, GaussPoint* ip)
{
	abort();
#if 0
  IDNLMaterialStatus *nlStatus = (IDNLMaterialStatus*) this -> giveStatus (ip);

  this->buildNonlocalPointTable(ip);
  this->updateDomainBeforeNonlocAverage(stepN);

  return buff.packDouble (nlStatus->giveLocalEquivalentStrainForAverage());
#endif
}

int 
TrabBoneNL3D::unpackAndUpdateUnknowns (CommunicationBuffer& buff, TimeStep* stepN, GaussPoint* ip)
{
	abort();
#if 0
  int result ;
  IDNLMaterialStatus *nlStatus = (IDNLMaterialStatus*) this -> giveStatus (ip);
  double localEquivalentStrainForAverage;

  result = buff.unpackDouble (localEquivalentStrainForAverage);
  nlStatus->setLocalEquivalentStrainForAverage (localEquivalentStrainForAverage);
  return result;
#endif
}

int 
TrabBoneNL3D::estimatePackSize (CommunicationBuffer& buff, GaussPoint* ip)
{
	abort();
#if 0
  // Note: nlStatus localStrainVectorForAverage memeber must be properly sized!
  // IDNLMaterialStatus *nlStatus = (IDNLMaterialStatus*) this -> giveStatus (ip);
  return buff.givePackSize(MPI_DOUBLE,1);
#endif
}
#endif

//
// END: PARALLEL MODE OPTION
/////////////////////////////////////////////////////////////////
