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
#include <string>

#include "Algorithm/AlgConfigPool.h"

#include "NuGamma/RESNCGammaGenerator.h"
#include "NuGamma/LARNCGammaXSec.h"

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

#include <iostream>
#include <fstream>
#include <sstream>

using namespace genie;
using namespace genie::constants;
using namespace genie::controls;
using namespace TMath;
using namespace std;
using std::string;

string ToString (double k);


//TFile* RESNCGammaGenerator::FileOfEvents = new TFile("file.root", "RECREATE");
//TTree* RESNCGammaGenerator::TreeOfEvents = new TTree("events");
  
//___________________________________________________________________________
RESNCGammaGenerator::RESNCGammaGenerator() :
KineGeneratorWithCache("genie::RESNCGammaGenerator")
{
  Gamma             = new TLorentzVector(0,0,0,0);
  OutgoingNucleon   = new TLorentzVector(0,0,0,0);
  OutgoingNeutrino  = new TLorentzVector(0,0,0,0);
  Resonance	    = new TLorentzVector(0,0,0,0);
//  FileOfEvents->cd();
//  TreeOfEvents->SetBranchAddress(Q);
//  TreeOfEvents->SetBranchAddress(Q);
//  TreeOfEvents->SetBranchAddress(Q);

  
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

  // Set initial state
  const InitialState & init_state = interaction -> InitState();

  //-- For the subsequent kinematic selection with the rejection method:
  //   Calculate the max differential cross section or retrieve it from the
  //   cache. Throw an exception and quit the evg thread if a non-positive
  //   value is found.
  //   If the kinematics are generated uniformly over the allowed phase
  //   space the max xsec is irrelevant  
  double xsec_max = this->ComputeMaxXSec(interaction);
  
  
  double xsec = -1;

  unsigned int iter = 0;
  bool env_accept = false;

  // Variables that we are throwing in the rejection method
  double gW        = -900; // hadronic invariant mass
  double gQ2       = -900; // momentum transfer
  double gEGamma   = -900; // energy of the photon in LAB FRAME
  double gPhiGamma = -900; // cosine of the angle between the scattering plane and the photon emission in LAB FRAME
  double gThetaGamma = -900;
//  double separationW = 1.4;
//  Range1D_t range_Wsep[2] = {range_W,range_W};
//  range_Wsep[0].max = separationW;
//  range_Wsep[1].min = separationW;
  // Set W range length
//  double dW[2]   = {range_Wsep[0].max - range_Wsep[0].min,
//                    range_Wsep[1].max - range_Wsep[1].min};
    
//  double area[2] = {(separationW - range_W.min) * 200.,
//                    (range_W.max - separationW) * 30.,};
//  double totalarea= area[0]+area[1];
//  double trialarea = rnd->RndKine().Rndm() * totalarea;
//  int narea=-1;
//  if(trialarea< area[0]) narea = 0;
//  else                   narea = 1;

double dW = range_W.max - range_W.min;
  
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
//	gW = 1.232;

    LOG("RESNCgKinematic", pINFO) << "Trying: W        = " << gW;	   

    // Set kinematics for current trial
    interaction->KinePtr()->SetW(gW);

    // Given chosen gW, find the range of Q2		
    
    Range1D_t range_Q2 = kps.Q2Lim_W();
    LOG("RESNCgKinematic", pINFO) << "Q2 min = " << range_Q2.min;
    double dQ2   = range_Q2.max - range_Q2.min;

//    if( range_Q2.min > 0.25 || range_Q2.max < 0.25){
//	break;
//    }
//  FileOfEvents->cd();
    gQ2 = range_Q2.min + dQ2 * rnd->RndKine().Rndm();
//    gQ2 = 0.25;

    // Set kinematics
    interaction->KinePtr()->SetQ2(gQ2);
  
    // Computing cross section for the current kinematics
//    xsec = fXSecModel->XSec(interaction);

//    TreeOfEvent->Fill();
    // Randomly pick point for rejection method
//    double tQ2W =  xsec_max * rnd->RndKine().Rndm();
    //This is probably something I will need to implement?
    //this->AssertXSecLimits(interaction, xsec, xsec_max);
//    LOG("RESNCgKinematic", pINFO) << "tQ2W = " << tQ2W;
//    LOG("RESNCgKinematic", pINFO) << "xsec = " << xsec;
    
//    if(xsec > xsec_max) {
//	LOG("RESNCgKinematic", pWARN) << "@~@~@~@~@~@ WARNING: MAX XSEC SET TOO LOW! @~@~@~@~@~@~@";
//   }

    // Check that we lie below (naive) envelope
    double xsec_test = xsec_max * rnd->RndKine().Rndm();
    double env_height = 0;
    double xsec_max_W = 0;
    
    if (gW > 1.208){
	xsec_max_W = xsec_max * ( 1 - 0.5 * (1/( range_W.max - 1.208 )) * ( gW - 1.208 ) );
    } else {
	xsec_max_W = xsec_max * ( 1 + 0.5 * (1/1.208) * ( gW - 1.208 ) );
    }

    if (gQ2 < 1){
	env_height = xsec_max_W;
    } else {
	env_height = xsec_max_W * ( 1 - 0.4 * (1/(range_Q2.max-1)) * (gQ2-1) );
    }

    env_height = xsec_max;
    env_accept = ( xsec_test < env_height );
    
//    env_accept = true;

    // If the generated kinematics are accepted, we generate the photon
    // kinematics based on Luis code
    if(env_accept) {
   
	gThetaGamma = acos( 2 * rnd->RndKine().Rndm() -1 );
 	gPhiGamma = ( 2 * rnd->RndKine().Rndm() -1 ) * TMath::Pi();
	((LARNCGammaXSec*)fXSecModel)->SetThetaPhoton(gThetaGamma);
	((LARNCGammaXSec*)fXSecModel)->SetPhiPhoton(gPhiGamma);
	xsec = fXSecModel->XSec(interaction);
	if ( xsec > env_height ) {
		LOG("RESNCgKinematic", pWARN) << "@~@~@~@~@~@~@~@~@~@ CROSS SECTION ABOVE ENVELOPE!! @~@~@~@~@~@~@~@~@~@";
		ofstream failfile;
		failfile.open ("failfile.txt");
		string failstring = "Envelope was lower than xsec. Q2 = " + ToString(gQ2) + " , W = " + ToString(gW) + " , Theta = " + ToString(gThetaGamma) + " , Phi = " + ToString(gPhiGamma) + ". xsec = " + ToString(xsec) + ". env_height = " + ToString(env_height); 
		failfile << failstring.c_str();
		failfile.close();
	}
	LOG("RESNCgKinematic", pINFO) << "env_height = " << env_height;
	LOG("RESNCgKinematic", pINFO) << "xsec = " << xsec;
	LOG("RESNCgKinematic", pINFO) << "xsec_test = " << xsec_test;
	
	if ( xsec_test < xsec){
		
	
        
	double TgtNucM = -1;

	if(init_state.IsNuN()) {
		TgtNucM = PDGLibrary::Instance()->Find(2112)->Mass();
	} else if(init_state.IsNuP()) {
		TgtNucM = PDGLibrary::Instance()->Find(2212)->Mass();
	} else {
		LOG("RESNCgKinematic", pWARN) << "*** Hit nucleon not a nucleon???";
	}

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

	// Now with the resonance 4-mom, we can calculate the photon 4-momentum

	gEGamma = 0.5 * ( Resonance->Dot(*Resonance) - Power(TgtNucM,2) ) * (1 / ( Resonance->E() - Resonance->Px() * Sin(gThetaGamma) * Cos(gPhiGamma) - Resonance->Py() * Sin(gThetaGamma) * Sin(gPhiGamma) - Resonance->Pz() * Cos(gThetaGamma) ) );

	Gamma->SetPxPyPzE( gEGamma * Sin(gThetaGamma) * Cos(gPhiGamma) , gEGamma * Sin(gThetaGamma) * Sin(gPhiGamma) , gEGamma * Cos(gThetaGamma) , gEGamma );

	OutgoingNucleon->SetPxPyPzE( Resonance->Px() - Gamma->Px() , Resonance->Py() - Gamma->Py() , Resonance->Pz() - Gamma->Pz() , Resonance->E() - Gamma->E() );

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
double RESNCGammaGenerator::ComputeMaxXSec (const Interaction * interaction) const{
//  (void)interaction;
//  return 60.;

// Naively find the max xsec, assuming it lies near Q2=0 and W=1.232

  LOG("RESNCgKinematic", pINFO) << "Calculating max xsec";

double max_xsec = 0;
int iter = 0;
double maxQ2 = 0;
double maxW = 0;
double maxTheta = 0;
double maxPhi = 0;

	for (int i=0; i<100; i++){
	  for (int j=0; j<1; j++){
	    for (int k=0; k<1; k++){
	      for (int l=0; l<1; l++){
		interaction->KinePtr()->SetQ2(0.001 + i*0.05);
		interaction->KinePtr()->SetW(1.232 + 0.0001*j);
		((LARNCGammaXSec*)fXSecModel)->SetThetaPhoton(0.7);
		((LARNCGammaXSec*)fXSecModel)->SetPhiPhoton( 1.4 );
		double xsec = fXSecModel->XSec(interaction);
		if ( xsec > max_xsec ){
			max_xsec = xsec;
			maxQ2 = 0.001 + i*0.01;
			maxW = 1.805 + 0.0001*j;
			maxTheta = k * (TMath::Pi()/10);
			maxPhi =  (l-5) * (TMath::Pi()/5) ;
		}
		interaction->KinePtr()->ClearRunningValues();
		iter++;
		LOG("RESNCgKinematic", pINFO) << "Iteration " << iter << " of " << 1*401*11*11;
		LOG("RESNCgKinematic", pINFO) << "At ( Q2 , W , Theta , Phi ) = ( " << 0.001 + i*0.01 << " , " << 1.805 + 0.0001*j << " , " << k * (TMath::Pi()/10) << " , " <<  (l-5) * (TMath::Pi()/5)  << " )";
		LOG("RESNCgKinematic", pINFO) << "Xsec so far:  " << xsec;
		LOG("RESNCgKinematic", pINFO) << "At ( Q2 , W , Theta , Phi ) = ( " << maxQ2 << " , " << maxW << " , " << maxTheta << " , " << maxPhi << " )";
	      }
	    }
	  }
	}

//max_xsec = 1.2 * max_xsec;

LOG("RESNCgKinematic", pFATAL) << "max xsec: " << max_xsec;
exit(0);
return max_xsec;
//return 28.9127; // Correct for 1 GeV neutrino
//return 50;
//return 5;
//return 24.3454;
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
//_________________________________________________________________________________
string ToString(double k)
{
  stringstream stream;
  stream << k;
  return stream.str();
}
