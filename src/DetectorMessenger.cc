#include "DetectorMessenger.hh"

#include "DetectorConstruction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

#include "G4UIcmdWithABool.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction* Det)
:Detector(Det)
{
  detDir = new G4UIdirectory("/exp/");
  detDir->SetGuidance("experiment control");

  UpdateCmd = new G4UIcmdWithoutParameter("/exp/update",this);
  UpdateCmd->SetGuidance("Update geometry.");
  UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateCmd->SetGuidance("if you changed geometrical value(s).");
  UpdateCmd->AvailableForStates(G4State_Idle);

  innerCmd = new G4UIcmdWithABool("/exp/setInnerShieldOn",this);
  innerCmd->SetGuidance("Set inner shield");

  outerCmd = new G4UIcmdWithABool("/exp/setOuterShieldOn",this);
  outerCmd->SetGuidance("Set outer shield");

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete UpdateCmd;
  delete detDir;
  //delete Dir;
  delete innerCmd;
  delete outerCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  if( command == UpdateCmd )
   { Detector->UpdateGeometry(); }

  if( command == innerCmd )
   { Detector->SetInnerShield(innerCmd->GetNewBoolValue(newValue));}

  if( command == outerCmd )
   { Detector->SetOuterShield(outerCmd->GetNewBoolValue(newValue));}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......