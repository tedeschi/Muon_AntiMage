//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//

#include "StackingAction.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

StackingAction::StackingAction()
{}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
StackingAction::~StackingAction()
{}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4ClassificationOfNewTrack
StackingAction::ClassifyNewTrack(const G4Track * aTrack)
{
  G4ClassificationOfNewTrack classification = fWaiting;


  switch(stage)
  {
  case 0: // Stage 0 : Primary muons only
    {
    if(aTrack->GetParentID()==0) {
      G4ParticleDefinition * particleType = aTrack->GetDefinition();
      if(particleType==G4MuonMinus::MuonMinusDefinition()){
		    classification = fUrgent;
	      }
      }
    break;
    }
  default: //Stage 1: all other particles, dont start the ones below a certain depth (-1.1m is lower veto panel)
    {
   	classification = fUrgent;
    if(aTrack->GetPosition().getZ() < -2*CLHEP::m){
	    classification = fKill;
	    break;
      } 
    }
  }
  return classification;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void StackingAction::NewStage()
{
  stage++;
  //G4cout<< "new stage " << stage << G4endl;
  if (stage == 1){
	if(!(reachedveto)) {
      stackManager->clear();
      //G4cout << "++++++++ event aborted" << G4endl;
      return;
    }
    stackManager->ReClassify();
    return;
  }
  else
  {
  // Other stage change : just re-classify
    stackManager->ReClassify();
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void StackingAction::PrepareNewEvent()
{
  stage = 0;
  reachedveto = FALSE;
}
