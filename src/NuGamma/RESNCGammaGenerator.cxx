//____________________________________________________________________________
/*
 Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author:  Pierre Lasorak <p.j.j.lasorak \at qmul.ac.uk>
          Queen Mary University, London

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include <TMath.h>
#include <cstdlib>

#include "Algorithm/AlgConfigPool.h"

#include "NuGamma/RESNCGammaGenerator.h"

#include "EVGCore/EVGThreadException.h"
#include "EVGCore/EventGeneratorI.h"
#include "EVGCore/RunningThreadInfo.h"

#include "Conventions/GBuild.h"
#include "Conventions/KineVar.h"
#include "Conventions/Constants.h"
#include "Conventions/Controls.h"
#include "Conventions/KinePhaseSpace.h"

#include "GHEP/GHepStatus.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepRecord.h"
#include "GHEP/GHepFlags.h"

#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"

#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGLibrary.h"

#include "Utils/Range1.h"
#include "Utils/KineUtils.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::controls;

//___________________________________________________________________________
RESNCGammaGenerator::RESNCGammaGenerator() :
KineGeneratorWithCache("genie::RESNCGammaGenerator")
{
  LOG("RESNCgKinematic", pINFO) << "RESNCGammaGenerator::RESNCGammaGenerator()";
}

//___________________________________________________________________________
RESNCGammaGenerator::RESNCGammaGenerator(string config) :
KineGeneratorWithCache("genie::RESNCGammaGenerator", config)
{
  LOG("RESNCgKinematic", pINFO) << "RESNCGammaGenerator::()RESNCGammaGenerator(string config)";
}

//___________________________________________________________________________
RESNCGammaGenerator::~RESNCGammaGenerator()
{
  LOG("RESNCgKinematic", pINFO) << "RESNCGammaGenerator::~RESNCGammaGenerator()";
}

//___________________________________________________________________________
void RESNCGammaGenerator::ProcessEventRecord(GHepRecord * evrec) const
{

  LOG("RESNCgKinematic", pINFO) << "-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/";
  LOG("RESNCgKinematic", pINFO) << "-/-RESNCGammaGenerator::ProcessEventRecord(GHepRecord * evrec) const-/-/";
  LOG("RESNCgKinematic", pINFO) << "-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/";
  //exit(1);
  // this->AddPhoton(evrec);
  // this->AddFinalStateNeutrino(evrec);
  // this->AddRecoilNucleon(evrec);
  // this->AddTargetRemnant(evrec);
  // for
  // if(fGenerateUniformly) {
  //   LOG("RESNCgKinematic", pNOTICE)
  //         << "Generating kinematics uniformly over the allowed phase space";
  // }

  //-- Get the random number generators
  RandomGen * rnd = RandomGen::Instance();

  //-- Access cross section algorithm for running thread
  RunningThreadInfo * rtinfo = RunningThreadInfo::Instance();
  const EventGeneratorI * evg = rtinfo->RunningThread();
  fXSecModel = evg->CrossSectionAlg();

  //-- Get the interaction from the GHEP record
  Interaction * interaction = evrec->Summary();
  interaction->SetBit(kISkipProcessChk);

  //-- Check for EM or CC process
  bool is_not_valid = interaction->ProcInfo().IsEM() || interaction->ProcInfo().IsWeakCC();
  if(is_not_valid){
    LOG("RESNCgKinematic", pFATAL) << "The interaction is EM or CC";
    //evrec->EventFlags()->SetBitNumber(kKineGenErr, true);
    genie::exceptions::EVGThreadException exception;
    exception.SetReason("Bad configuration of interaction");
    exception.SwitchOnFastForward();
    throw exception;
  }
  
  //-- Compute the phase space limits
  const KPhaseSpace & kps = interaction->PhaseSpace();
  Range1D_t range_W  = kps.Limits(kKVW);
  Range1D_t range_Q2 = kps.Limits(kKVQ2);

  // Luis's model can go till W<~2GeV
  // Check is the limit are reasonable
  if(range_W.max <= 0 || range_W.min >= range_W.max)
    {
      LOG("RESNCgKinematic", pWARN) << "No available phase space";
      LOG("RESNCgKinematic", pWARN) << "range_W.min " << range_W.min;
      LOG("RESNCgKinematic", pWARN) << "range_W.max " << range_W.max;
      evrec->EventFlags()->SetBitNumber(kKineGenErr, true);
      genie::exceptions::EVGThreadException exception;
      exception.SetReason("No available phase space");
      exception.SwitchOnFastForward();
      throw exception;
    }

  // Set initial state
  const InitialState & init_state = interaction -> InitState();
  double E = init_state.ProbeE(kRfHitNucRest);
  //  double M = init_state.Tgt().HitNucP4().M();
  //  double ml  = interaction->FSPrimLepton()->Mass();

  //-- For the subsequent kinematic selection with the rejection method:
  //   Calculate the max differential cross section or retrieve it from the
  //   cache. Throw an exception and quit the evg thread if a non-positive
  //   value is found.
  //   If the kinematics are generated uniformly over the allowed phase
  //   space the max xsec is irrelevant
  double xsec_max = this->MaxXSec(evrec);

  //-- Try to select a valid W, Q2 pair using the rejection method
  double dW   = range_W.max - range_W.min;
  
  double xsec = -1;

  unsigned int iter = 0;
  bool accept = false;

  // Variables that we are throwing in the rejection method
  double gW        = 0; // hadronic invariant mass
  double gQ2       = 0; // momentum transfer
  double gEGamma   = 0; // energy of the photon
  double gPhiGamma = 0; // cosine of the angle between the scattering plane and the photon emission

  while(1) {
    iter++;
    if(iter > kRjMaxIterations) {
      LOG("RESNCgKinematic", pWARN)
	<< "*** Could not select a valid (W,Q2, EGamma, cosPhi) pair after "
	<< iter << " iterations";
      evrec->EventFlags()->SetBitNumber(kKineGenErr, true);
      genie::exceptions::EVGThreadException exception;
      exception.SetReason("Couldn't select kinematics");
      exception.SwitchOnFastForward();
      throw exception;
    }
    
    // We do firstly a dumb rejection method (yes, it will be slow)
    gW        = rnd->RndKine().Rndm() * xsec_max;
    gQ2	      = rnd->RndKine().Rndm() * xsec_max;
    gEGamma   = rnd->RndKine().Rndm() * xsec_max;
    gPhiGamma = rnd->RndKine().Rndm() * xsec_max;

    LOG("RESNCgKinematic", pINFO) << "Trying: W        = " << gW;	   
    LOG("RESNCgKinematic", pINFO) << "Trying: Q2       = " << gQ2;	   
    LOG("RESNCgKinematic", pINFO) << "Trying: EGamma   = " << gEGamma; 
    LOG("RESNCgKinematic", pINFO) << "Trying: PhiGamma = " << gPhiGamma;

    //-- Set kinematics for current trial
    interaction->KinePtr()->SetW(gW);		     
    interaction->KinePtr()->SetQ2(gQ2);	     
    interaction->KinePtr()->SetKV(kKVEGamma,   gEGamma);    
    interaction->KinePtr()->SetKV(kKVPhiGamma, gPhiGamma);
  
    //-- Computing cross section for the current kinematics
    xsec = fXSecModel->XSec(interaction);

    double t =  xsec_max * rnd->RndKine().Rndm();
    //This is probably something I will need to implement?
    //this->AssertXSecLimits(interaction, xsec, xsec_max);
    LOG("RESNCgKinematic", pINFO) << "t = " << t;
    LOG("RESNCgKinematic", pINFO) << "xsec = " << xsec;
    
    accept = (t < xsec);

    //-- If the generated kinematics are accepted, finish-up module's job
    if(accept) {        
      LOG("RESNCgKinematic", pINFO) << "Selected: W        = " << gW;
      LOG("RESNCgKinematic", pINFO) << "          Q2       = " << gQ2;
      LOG("RESNCgKinematic", pINFO) << "          EGamma   = " << gEGamma;
      LOG("RESNCgKinematic", pINFO) << "          PhiGamma = " << gPhiGamma;

      // set the cross section for the selected kinematics
      //evrec->SetDiffXSec(xsec);

      // lock selected kinematics & clear running values
      interaction->KinePtr()->SetW(gW);
      interaction->KinePtr()->SetQ2(gQ2);
      interaction->KinePtr()->SetKV(kKVEGamma,   gEGamma);
      interaction->KinePtr()->SetKV(kKVPhiGamma, gPhiGamma);

      interaction->KinePtr()->ClearRunningValues();

      return;
    }
  } 
}


//___________________________________________________________________________
void RESNCGammaGenerator::AddPhoton(GHepRecord * evrec) const
{
  std::cout << "void RESNCGammaGenerator::AddPhoton(GHepRecord * evrec) const" << std::endl;
  // Adding the final state photon
  //
  LOG("RESNCgKinematic", pINFO) << "Adding final state photon";

}

//___________________________________________________________________________
void RESNCGammaGenerator::AddFinalStateNeutrino(GHepRecord * evrec) const
{
  std::cout << "void RESNCGammaGenerator::AddFinalStateNeutrino(GHepRecord * evrec) const" << std::endl;
  // Adding the final state neutrino
  // Just use 4-momentum conservation (init_neutrino = photon + final_neutrino)
  
  LOG("RESNCgKinematic", pINFO) << "Adding final state neutrino";

 
}

//___________________________________________________________________________
void RESNCGammaGenerator::AddRecoilNucleon(GHepRecord * evrec) const
{
  // Adding the recoil nucleon.
  LOG("RESNCgKinematic", pINFO) << "Adding recoil nucleon";

}

//___________________________________________________________________________
void RESNCGammaGenerator::AddTargetRemnant(GHepRecord * evrec) const
{
  // Add the remnant nuclear target at the GHEP record
  std::cout << "void RESNCGammaGenerator::AddTargetRemnant(GHepRecord * evrec) const" << std::endl;
  LOG("RESNCgKinematic", pINFO) << "Adding final state nucleus";


}
//___________________________________________________________________________
double RESNCGammaGenerator::Energy(const Interaction * interaction) const
{
  // Override the base class Energy() method to cache the max xsec for the
  // neutrino energy in the LAB rather than in the hit nucleon rest frame.
  
  const InitialState & init_state = interaction->InitState();
  double E = init_state.ProbeE(kRfLab);
  return E;
}
//___________________________________________________________________________
double RESNCGammaGenerator::ComputeMaxXSec(GHepRecord * evrec) const
{
  std::cout << "double RESNCGammaGenerator::ComputeMaxXSec(GHepRecord * evrec) const" << std::endl;
  return 2.;
}
//___________________________________________________________________________
//___________________________________________________________________________
void RESNCGammaGenerator::Configure(const Registry & config)
{
  LOG("RESNCgKinematic", pINFO) << "registry " << config;
  Algorithm::Configure(config);
  LOG("RESNCgKinematic", pINFO) << "registry " << config;
  this->LoadConfig();
}
//____________________________________________________________________________
void RESNCGammaGenerator::Configure(string config)
{
  LOG("RESNCgKinematic", pINFO) << "registry " << config;
  Algorithm::Configure(config);
  LOG("RESNCgKinematic", pINFO) << "registry " << config;

  this->LoadConfig();
}
//____________________________________________________________________________
void RESNCGammaGenerator::LoadConfig(void)
{
  fGenerateUniformly = fConfig->GetBoolDef("UniformOverPhaseSpace", false);
}
//____________________________________________________________________________

