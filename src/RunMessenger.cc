//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunMessenger.hh"

#include <sstream>

#include "RunAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunMessenger::RunMessenger(RunAction* manager)
:outputAction (manager)
{
  fileDir = new G4UIdirectory("/file/");
  fileDir->SetGuidance("output control");

  factoryCmd = new G4UIcmdWithAString("/file/setFileName",this);
  factoryCmd->SetGuidance("set name for the root file");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunMessenger::~RunMessenger()
{
  delete fileDir;
  delete factoryCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RunMessenger::SetNewValue(G4UIcommand* command,
                                          G4String newValue)
{ 
  if(command == factoryCmd) outputAction->filename=newValue;
              
   
}
