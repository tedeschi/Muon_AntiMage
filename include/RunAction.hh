/// \file RunAction.hh
/// \brief Definition of the RunAction class

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"


//root includes
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include <vector>
#include <sys/time.h>


class G4Run;
class DetectorConstruction;
class PrimaryGeneratorAction;
class RunMessenger;


class RunAction : public G4UserRunAction
{
  public:
    RunAction(DetectorConstruction*, PrimaryGeneratorAction*);
    virtual ~RunAction();

    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);

	double energy[301];
	void Fill();
		double GetTimeDiff();

    G4String filename;
  private:

	TH1F*  energyHistoGe;
	TH1F*  energyHistoGeWithVeto;
	TH1F*  energyHistoVeto;
	TH1F*  MultiplicityGe;
	TH1F*  MultiplicityGeWithVeto;
	TH1F*  MultiplicityVeto;

	TFile* fRootFout;
	TTree* fRootTree;

	DetectorConstruction   *detector;
	PrimaryGeneratorAction *primary;
    RunMessenger           *runMessenger;


    timeval start, end, now;

	std::vector<double> det_energy;
	std::vector<int> det_channel;
	double start_energy;
	double start_phi;
	double start_theta;
  double start_x;
  double start_y;
  double start_z;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
