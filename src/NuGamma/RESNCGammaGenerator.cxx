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
using namespace TMath;

//___________________________________________________________________________
RESNCGammaGenerator::RESNCGammaGenerator() :
KineGeneratorWithCache("genie::RESNCGammaGenerator")
{
  Gamma             = new TLorentzVector();
  OutgoingNucleon   = new TLorentzVector();
  OutgoingNeutrino  = new TLorentzVector();
  Resonance         = new TLorentzVector();



  
}

//___________________________________________________________________________
RESNCGammaGenerator::RESNCGammaGenerator(string config) :
KineGeneratorWithCache("genie::RESNCGammaGenerator", config)
{
  LOG("RESNCgKinematic", pINFO) << "RESNCGammaGenerator::RESNCGammaGenerator(string config)";
}

//___________________________________________________________________________
RESNCGammaGenerator::~RESNCGammaGenerator()
{
  LOG("RESNCgKinematic", pINFO) << "RESNCGammaGenerator::~RESNCGammaGenerator()";
  delete Gamma;
  delete OutgoingNucleon;
  delete OutgoingNeutrino;
  delete Resonance;

}

//___________________________________________________________________________
void RESNCGammaGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
  ThrowKinematics(evrec);
  AddPhoton             (evrec);
  AddFinalStateNeutrino (evrec);
  AddTargetRemnant      (evrec);
  AddRecoilNucleon      (evrec);
  
}

void RESNCGammaGenerator::ThrowKinematics(GHepRecord * evrec) const{
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
  
  //-- Compute the phase space limits only for W
  const KPhaseSpace & kps = interaction->PhaseSpace();
  Range1D_t range_W  = kps.Limits(kKVW);
  // Range1D_t range_Q2 = kps.Limits(kKVQ2);

  // Luis's model can go till W<~2GeV
  // Check is the limit are reasonable
  if(range_W.max <= 0 || range_W.min >= range_W.max)
    {
      LOG("RESNCgKinematic", pWARN) << "No available W phase space";
      LOG("RESNCgKinematic", pWARN) << "range_W.min " << range_W.min;
      LOG("RESNCgKinematic", pWARN) << "range_W.max " << range_W.max;
      evrec->EventFlags()->SetBitNumber(kKineGenErr, true);
      genie::exceptions::EVGThreadException exception;
      exception.SetReason("No available W phase space");
      exception.SwitchOnFastForward();
      throw exception;
    }

  // Set initial state
  const InitialState & init_state = interaction -> InitState();
  //  double E = init_state.ProbeE(kRfHitNucRest);
  //  double M = init_state.Tgt().HitNucP4().M();
  //  double ml  = interaction->FSPrimLepton()->Mass();

  //-- For the subsequent kinematic selection with the rejection method:
  //   Calculate the max differential cross section or retrieve it from the
  //   cache. Throw an exception and quit the evg thread if a non-positive
  //   value is found.
  //   If the kinematics are generated uniformly over the allowed phase
  //   space the max xsec is irrelevant
  double xsec_max = this->MaxXSec(evrec);

  //-- Try to select a valid W using the rejection method
  double dW   = range_W.max - range_W.min;
  
  double xsec = -1;

  unsigned int iter = 0;
  bool acceptW = false;
  bool acceptQ2 = false;
  //bool acceptEGamma = false;
  //bool acceptPhiGamma = false;

  // Variables that we are throwing in the rejection method
  double gW        = 0; // hadronic invariant mass
  double gQ2       = 0; // momentum transfer
  double gEGamma   = 0; // energy of the photon in LAB FRAME
  double gPhiGamma = 0; // cosine of the angle between the scattering plane and the photon emission in LAB FRAME

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
    
    // We firstly use a (simple) rejection method to choose a value for W

    gW        = range_W.min + dW * rnd->RndKine().Rndm();

    LOG("RESNCgKinematic", pINFO) << "Trying: W        = " << gW;	   

    //-- Set kinematics for current trial
    interaction->KinePtr()->SetW(gW);		     
    //interaction->KinePtr()->SetQ2(gQ2);	     
    //interaction->KinePtr()->SetKV(kKVEGamma,   gEGamma);    
    //interaction->KinePtr()->SetKV(kKVPhiGamma, gPhiGamma);
  
    //-- Computing cross section for the current kinematics
    xsec = fXSecModel->XSec(interaction);

    double tW =  xsec_max * rnd->RndKine().Rndm();
    //This is probably something I will need to implement?
    //this->AssertXSecLimits(interaction, xsec, xsec_max);
    LOG("RESNCgKinematic", pINFO) << "tW = " << tW;
    LOG("RESNCgKinematic", pINFO) << "xsec = " << xsec;
    
    acceptW = (tW < xsec);

    if(acceptW) {
	// We now do the same for Q2, given our fixed value for W

	while(acceptQ2 == false){
	Range1D_t range_Q2 = kps.Q2Lim_W();
	double dQ2   = range_Q2.max - range_Q2.min;
	gQ2 = range_Q2.min + dQ2 * rnd->RndKine().Rndm();

	interaction->KinePtr()->SetQ2(gQ2);
	xsec = fXSecModel->XSec(interaction);
	
	double tQ2 = xsec_max * rnd->RndKine().Rndm();
	LOG("RESNCgKinematic", pINFO) << "tQ2 = " << tQ2;
    	LOG("RESNCgKinematic", pINFO) << "xsec = " << xsec;

	acceptQ2 = (tQ2 < xsec);
	}
	
     }

    //-- If the generated kinematics are accepted, we randomly generate the photon
    //-- kinematics, and finish
    if(acceptQ2) {

	// We calculate the photon 4-momentum in the lab frame. This involves
	// working in the hadronic centre of mass frame and disregarding any
	// polarisation effects, thus allowing the photons to be generated
	// isotropically. We then boost back to the lab frame

	double TgtNucM = -1;

	if(init_state.IsNuN()) {
		TgtNucM = 0.939565346;
	} else if(init_state.IsNuP()) {
		TgtNucM = 0.938272046;
	} else {
		LOG("RESNCgKinematic", pWARN) << "*** Hit nucleon not a nucleon???";
	}

	// A quick check

	std::cout << "Lower limit of range of W: " << gW << std::endl;
	std::cout << "Hit nucleon rest mass: " << TgtNucM << std::endl;
	
	// With nucleon rest mass determined, easy now to calculate hadr frame 4-momentum
	// In the hadr frame (i.e. rest frame of resonance):

	double EGamma_Hadr = ( Power(gW,2) - Power(TgtNucM,2) ) / (2 * gW);
	double CosThetaGamma_Hadr = ( 2 * rnd->RndKine().Rndm() ) -1;
	double PhiGamma_Hadr = 2*TMath::Pi()* rnd->RndKine().Rndm();
	Gamma->SetPxPyPzE( EGamma_Hadr * Sqrt(1 - Power(CosThetaGamma_Hadr,2)) * Cos(PhiGamma_Hadr) ,
                           EGamma_Hadr * Sqrt(1 - Power(CosThetaGamma_Hadr,2)) * Sin(PhiGamma_Hadr) ,
                           EGamma_Hadr * CosThetaGamma_Hadr ,
                           EGamma_Hadr );

	// We now need to determine the 4-momentum of the resonance, in the lab frame

	// First, we get the probe 4-momentum in the hit nucleon rest frame. We first do this assuming
	// the probe travels along the z axis (which it doesn't necessaerily in this frame)

	TLorentzVector * ProbeP4 = init_state.GetProbeP4(kRfHitNucRest);
	TLorentzVector HitNucP4 = init_state.Tgt().HitNucP4();

	// In the hit nucleon rest frame, we use W and Q2 to calculate the resonance
	// 4-momentum, randomising the angle about which the system is symmetric
	// We first do this assuming the probe travels along the z axis (which it 
	// doesn't necessarily in this frame)

	double HitNucM = init_state.Tgt().HitNucP4().M();
	double OutNeuE = ProbeP4->E() + (1/(2*HitNucM))*( Power(HitNucM,2) - Power(gW,2) - gQ2 );
	double OutNeuCosTheta = 1 - (1/(2*ProbeP4->E()*OutNeuE))*gQ2;
	double OutNeuSinTheta = Sqrt(1-Power(OutNeuCosTheta,2));
	double OutNeuPhi = 2*TMath::Pi() * rnd->RndKine().Rndm();
	
	Resonance->SetPxPyPzE( -OutNeuE * OutNeuSinTheta * Cos(OutNeuPhi) ,
                               -OutNeuE * OutNeuSinTheta * Sin(OutNeuPhi) ,
                               ProbeP4->E() - OutNeuE * OutNeuCosTheta ,
                               ProbeP4->E() - OutNeuE + HitNucM );

	// We must now rotate to align with the probe angle of approach

	double ProbeTheta = ACos( ProbeP4->Z() / ProbeP4->E() );
	double ProbePhi = ATan2( ProbeP4->Y() , ProbeP4->X() );

	Resonance->RotateY(ProbeTheta);
	Resonance->RotateZ(ProbePhi);

	// Now, we boost the resonance 4-momementum by the hit nucleon velocity,
	// which gives us the resonance 4-momentum in the lab frame			

	Resonance->Boost(HitNucP4.BoostVector());

	// Finally, we can boost the photon 4-momentum by the resonance velocity,
	// giving us the photon 4-momentum in the lab frame

	Gamma->Boost(Resonance->BoostVector());

	// From this 4-momentum, we calculate the energy and phi angle of the
	// photon, in the lab frame

	gEGamma = Gamma->E();
	gPhiGamma = ATan2(Gamma->Y(),Gamma->X());

	if(gPhiGamma < 0) {
		gPhiGamma = gPhiGamma + 2*Pi();
	}
        
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
  Gamma;
}

//___________________________________________________________________________
void RESNCGammaGenerator::AddFinalStateNeutrino(GHepRecord * evrec) const
{
  std::cout << "void RESNCGammaGenerator::AddFinalStateNeutrino(GHepRecord * evrec) const" << std::endl;
  // Adding the final state neutrino
  // Just use 4-momentum conservation (init_neutrino = photon + final_neutrino)
  
  LOG("RESNCgKinematic", pINFO) << "Adding final state neutrino";
  (void)evrec;
 
}

//___________________________________________________________________________
double RESNCGammaGenerator::ComputeMaxXSec (const Interaction * in) const{
  (void)in;
  return 2.;
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
