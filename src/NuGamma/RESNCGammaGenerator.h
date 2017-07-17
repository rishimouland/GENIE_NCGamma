//____________________________________________________________________________
/*!

\class    genie::RESNCGammaGenerator

\brief    

\author   Pierre Lasorak <p.j.j.lasorak \at qmul.ac.uk>
          Queen Mary University, London

\created  November 24, 2015

\cpright  Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _RESNCGAMMA_GENERATOR_H_
#define _RESNCGAMMA_GENERATOR_H_

#include "EVGModules/KineGeneratorWithCache.h"

namespace genie {

class RESNCGammaGenerator : public KineGeneratorWithCache {

public :
  RESNCGammaGenerator();
  RESNCGammaGenerator(string config);
  ~RESNCGammaGenerator();
  
  void Configure(const Registry & config);
  void Configure(string config);
  
  // methods to load sub-algorithms and config data from the Registry
  void LoadConfig (void);
  
  //-- implement the EventRecordVisitorI interface
  void ProcessEventRecord (GHepRecord * event_rec) const;

  // Overload KineGeneratorWithCache functions
  double ComputeMaxXSec      (GHepRecord * event_rec) const;
  double Energy              (const Interaction * in) const;


private:
  void AddPhoton             (GHepRecord * event_rec) const;
  void AddFinalStateNeutrino (GHepRecord * event_rec) const;
  void AddTargetRemnant      (GHepRecord * event_rec) const;
  void AddRecoilNucleon      (GHepRecord * event_rec) const;
  bool fGenerateUniformly;
};

}      // genie namespace
#endif // _RESNCGAMMA_GENERATOR_H_
