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

  //-- Get the random number generators
  RandomGen * rnd = RandomGen::Instance();

  //-- Access cross section algorithm for running thread
  RunningThreadInfo * rtinfo = RunningThreadInfo::Instance();
  const EventGeneratorI * evg = rtinfo->RunningThread();
  fXSecModel = evg->CrossSectionAlg();

  //-- Get the interaction from the GHEP record
  Interaction * interaction = evrec->Summary();
  interaction->SetBit(kISkipProcessChk);

  // Set initial state
  const InitialState & init_state = interaction -> InitState();
  
  //-- Compute the phase space limits for W, given input parameters, hit nucleon momentum
  const KPhaseSpace & kps = interaction->PhaseSpace();
  Range1D_t range_W  = kps.Limits(kKVW);

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

  // We now check these limits against what we'd expect

//  TLorentzVector * ProbeP4p = init_state.GetProbeP4(kRfLab);
//  TLorentzVector ProbeP4 = *ProbeP4p;
//  double W2UpperLimCheck = Power( init_state.Tgt().HitNucP4().M() , 2 ) + 2 * ProbeP4.Dot(init_state.Tgt().HitNucP4());
//  double W2LowerLimCheck =  ProbeP4.E() + init_state.Tgt().HitNucP4().E() + Sqrt( 2 * ProbeP4.Dot(init_state.Tgt().HitNucP4()) + Power( init_state.Tgt().HitNucP4().M() , 2 ) - Power( ProbeP4.E() + init_state.Tgt().HitNucP4().E() , 2 ) ) ;

//  std::cout << std::endl;
//  std::cout << "Automatically calculated limits for W: " << std::endl;
//  std::cout << "Lower limit: " << range_W.min << std::endl;
//  std::cout << "Upper limit: " << range_W.max << std::endl;
//  std::cout << std::endl;
//  std::cout << "Check limits for W: " << std::endl;
//  std::cout << "Lower limit: " << Sqrt(W2LowerLimCheck) << std::endl;
//  std::cout << "Upper limit: " << Sqrt(W2UpperLimCheck) << std::endl;

  //-- For the subsequent kinematic selection with the rejection method:
  //   Calculate the max differential cross section or retrieve it from the
  //   cache. Throw an exception and quit the evg thread if a non-positive
  //   value is found.
  //   If the kinematics are generated uniformly over the allowed phase
  //   space the max xsec is irrelevant
  double xsec_max = this->MaxXSec(evrec);

  // Set W range length
  double dW   = range_W.max - range_W.min;
  
  double xsec = -1;

  unsigned int iter = 0;
  bool accept = false;

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

    gW = range_W.min + dW * rnd->RndKine().Rndm();

    LOG("RESNCgKinematic", pINFO) << "Trying: W        = " << gW;	   

    // Set kinematics for current trial
    interaction->KinePtr()->SetW(gW);

    // Given chosen gW, find the range of Q2		
    
    Range1D_t range_Q2 = kps.Q2Lim_W();
    double dQ2   = range_Q2.max - range_Q2.min;
    gQ2 = range_Q2.min + dQ2 * rnd->RndKine().Rndm();


    // Set kinematics
    interaction->KinePtr()->SetQ2(gQ2);
  
    // Computing cross section for the current kinematics
    xsec = fXSecModel->XSec(interaction);


    // Randomly pick point for rejection method
    double tQ2W =  xsec_max * rnd->RndKine().Rndm();
    //This is probably something I will need to implement?
    //this->AssertXSecLimits(interaction, xsec, xsec_max);
    LOG("RESNCgKinematic", pINFO) << "tQ2W = " << tQ2W;
    LOG("RESNCgKinematic", pINFO) << "xsec = " << xsec;
    
    accept = (tQ2W < xsec);

    // If the generated kinematics are accepted, we randomly generate the photon
    // kinematics, and finish
    if(accept) {

	// We calculate the photon 4-momentum in the lab frame. This involves
	// working in the hadronic centre of mass frame and disregarding any
	// polarisation effects, thus allowing the photons to be generated
	// isotropically. We then boost back to the lab frame

	double TgtNucM = -1;

	if(init_state.IsNuN()) {
		TgtNucM = PDGLibrary::Instance()->Find(2112)->Mass();
	} else if(init_state.IsNuP()) {
		TgtNucM = PDGLibrary::Instance()->Find(2212)->Mass();
	} else {
		LOG("RESNCgKinematic", pWARN) << "*** Hit nucleon not a nucleon???";
	}

	// A quick check

	//std::cout << "Lower limit of range of W: " << gW << std::endl;
	//std::cout << "Hit nucleon rest mass: " << TgtNucM << std::endl;
	
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
	std::cout << "Check - Hit nucleon 4-momentum: ";
	HitNucP4.Print();

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

	// We also calculate the outgoing nucleon 4-momentum

	OutgoingNucleon->SetPxPyPzE( Resonance->Px() - Gamma->Px() , Resonance->Py() - Gamma->Py() , Resonance->Pz() - Gamma->Pz() , Resonance->E() - Gamma->E() );

	// and also the outgoing neutrino

	TLorentzVector * ProbeP4Lab = init_state.GetProbeP4(kRfLab);
	OutgoingNeutrino->SetPxPyPzE( ProbeP4Lab->Px() + HitNucP4.Px() - Resonance->Px() , ProbeP4Lab->Py() + HitNucP4.Py() - Resonance->Py() , ProbeP4Lab->Pz() + HitNucP4.Pz() - Resonance->Pz() , ProbeP4Lab->E() + HitNucP4.E() - Resonance->E() );
        
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

  TLorentzVector x4(0,0,0,0);
  evrec->AddParticle(22, kIStStableFinalState, 2,-1,-1,-1, *Gamma, x4);
  
}

//___________________________________________________________________________
void RESNCGammaGenerator::AddFinalStateNeutrino(GHepRecord * evrec) const
{
  std::cout << "void RESNCGammaGenerator::AddFinalStateNeutrino(GHepRecord * evrec) const" << std::endl;
  // Adding the final state neutrino
  
  LOG("RESNCgKinematic", pINFO) << "Adding final state neutrino";

  TLorentzVector x4(0,0,0,0);

  // Get neutrino flavour
 
  Interaction * interaction = evrec->Summary();
  const InitialState & init_state = interaction->InitState();
  int ProbePdg = init_state.ProbePdg();

  evrec->AddParticle(ProbePdg, kIStStableFinalState, 0,-1,-1,-1, *OutgoingNeutrino, x4);

}

//___________________________________________________________________________
double RESNCGammaGenerator::ComputeMaxXSec (const Interaction * in) const{
  (void)in;
  return 3.;
}

//___________________________________________________________________________
void RESNCGammaGenerator::AddRecoilNucleon(GHepRecord * evrec) const
{
  // Adding the recoil nucleon.
  LOG("RESNCgKinematic", pINFO) << "Adding recoil nucleon";

  TLorentzVector x4(0,0,0,0);

  // Get nucleon pdg
 
  Interaction * interaction = evrec->Summary();
  const InitialState & init_state = interaction->InitState();
  int HitNucPdg = init_state.Tgt().HitNucPdg();

  evrec->AddParticle(HitNucPdg, kIStStableFinalState, 2,-1,-1,-1, *OutgoingNucleon, x4);

}

//___________________________________________________________________________
void RESNCGammaGenerator::AddTargetRemnant(GHepRecord * evrec) const
{
  // Add the remnant nuclear target at the GHEP record
  std::cout << "void RESNCGammaGenerator::AddTargetRemnant(GHepRecord * evrec) const" << std::endl;
  LOG("RESNCgKinematic", pINFO) << "Adding final state nucleus";


}
