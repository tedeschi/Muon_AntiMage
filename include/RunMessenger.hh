//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef RunMessenger_h
#define RunMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class RunAction;
class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWithAString;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class RunMessenger: public G4UImessenger
{
  public:

   RunMessenger(RunAction* );
  ~RunMessenger();

   void SetNewValue(G4UIcommand*, G4String);
   
  private:

   RunAction*          outputAction;
   
   G4UIdirectory*          fileDir;  
   G4UIcmdWithAString*     factoryCmd;
   

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

