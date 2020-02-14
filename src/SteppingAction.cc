//
// ********************************************************************
// * License and Disclaimer                                           *
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class

#include "SteppingAction.hh"

#include "DetectorConstruction.hh"
#include "EventAction.hh"

#include "G4Step.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4ThreeVector.hh"
#include "G4TrackStatus.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction* det, EventAction* evt, StackingAction* sta)
:detector(det), eventaction(evt), stackingaction(sta)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{

  // get volume of the current step
  G4VPhysicalVolume* volume
    = step->GetPreStepPoint()->GetTouchableHandle()
      ->GetVolume();

  G4String volumename = volume->GetName();
  G4String particleType = step->GetTrack()->GetDefinition()->GetParticleType();
  G4double edep = step->GetTotalEnergyDeposit();
  G4ThreeVector position = step->GetTrack()->GetPosition();
  
  // check if we deposit energy in a scoring volume
  if (edep>0){                                   //doesnt need a unit, as long as it not zero
	  if ((volumename(0,3) == "Det")||(volumename(0,4) == "pane")){
      G4int det_no = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber();
      eventaction->FillEnergy(edep,det_no);
	    //G4cout << "-- Hit in "<< volume->GetName() << " no "
      //       << step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber() << G4endl;

  	}
  }
  if((particleType="mu-")&&(position.mag() < 5*CLHEP::m)){
   stackingaction->Veto();  //only follow up with full simulation if incident muon is within 5 m of the origin (between M1 and M2)
  }
  //kill everyting below -5m to speed up things
  if(position.getZ() < -5*CLHEP::m){
    step->GetTrack()->SetTrackStatus(fStopAndKill);
  }
}


