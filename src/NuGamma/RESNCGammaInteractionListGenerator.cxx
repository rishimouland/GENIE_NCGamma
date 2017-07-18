//____________________________________________________________________________
/*
 Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Pierre Lasorak <p.j.j.lasorak \at qmul.ac.uk>
         Queen Mary University, London

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include "EVGCore/InteractionList.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "NuGamma/RESNCGammaInteractionListGenerator.h"
#include "PDG/PDGCodes.h"

using namespace genie;

//___________________________________________________________________________
RESNCGammaInteractionListGenerator::RESNCGammaInteractionListGenerator() :
InteractionListGeneratorI("genie::RESNCGammaInteractionListGenerator")
{
  std::cout << "RESNCGammaInteractionListGenerator" << std::endl;
}
//___________________________________________________________________________
RESNCGammaInteractionListGenerator::RESNCGammaInteractionListGenerator(string config):
  InteractionListGeneratorI("genie::RESNCGammaInteractionListGenerator", config)
{

}
//___________________________________________________________________________
RESNCGammaInteractionListGenerator::~RESNCGammaInteractionListGenerator()
{

}
//___________________________________________________________________________
InteractionList * RESNCGammaInteractionListGenerator::CreateInteractionList(
                                       const InitialState & init_state) const
{
  LOG("IntLst", pINFO)
     << "InitialState = " << init_state.AsString();
  std::cout << "READING RECNCGAMMAINTERACTION LIST GENERATOR"  << std::endl;

  int nupdg   = init_state.ProbePdg();
  int tgtpdg  = init_state.Tgt().Pdg();
  
  InteractionList * intlist = new InteractionList;
  
  const Target & target = init_state.Tgt();
  
  if(target.Z()>0) {
    Interaction * interaction = Interaction::RESNCGamma(tgtpdg,kPdgProton,nupdg,0);
    intlist->push_back(interaction);
  }
  if(target.N()>0) {
    Interaction * interaction = Interaction::RESNCGamma(tgtpdg,kPdgNeutron,nupdg,0);
    intlist->push_back(interaction);
  }

  return intlist;
}
