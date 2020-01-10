#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "RunMessenger.hh"
//root includes

#include <iostream>
using namespace std;
#include <time.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(DetectorConstruction* det,PrimaryGeneratorAction* prim)
:G4UserRunAction(),detector(det),primary(prim)
{
 filename = "Muon_AntiMaGe";
 runMessenger = new RunMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
  delete runMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{

  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  //inform the runManager to save random number seed
  //G4RunManager::GetRunManager()->SetRandomNumberStore(false);

  G4String filename2 = filename + "_"+std::to_string(time(0))+".root";
  fRootFout = new TFile(filename2,"RECREATE",filename) ;
  fRootTree = new TTree("tree",filename2);
  fRootTree->SetMaxTreeSize(1073741824);


  fRootTree->Branch("start_energy",&start_energy,"start_energy/D");
  fRootTree->Branch("start_phi",&start_phi,"start_phi/D");
  fRootTree->Branch("start_theta",&start_theta,"start_theta/D");
  fRootTree->Branch("start_x",&start_x,"start_x/D");
  fRootTree->Branch("start_y",&start_y,"start_y/D");
  fRootTree->Branch("start_z",&start_z,"start_z/D");


  fRootTree->Branch("det_energy",&det_energy);
  fRootTree->Branch("det_channel",&det_channel);

  energyHistoGe = new TH1F("energyHistoGe","energyHistoGe",20000,0,20000);
  energyHistoGeWithVeto = new TH1F("energyHistoGeWithVeto","energyHistoGeWithVeto",20000,0,20000);
  energyHistoVeto = new TH1F("energyHistoVeto","energyHistoVeto",20000,0,20000);
  MultiplicityGe = new TH1F("MultiplicityGe","MultiplicityGe",50,0,50);
  MultiplicityGeWithVeto = new TH1F("MultiplicityGeWithVeto","MultiplicityGeWithVeto",50,0,50);
  MultiplicityVeto = new TH1F("MultiplicityVeto","MultiplicityVeto",50,0,50);


  gettimeofday(&start, 0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* aRun)
{

  G4int nofEvents = aRun->GetNumberOfEvent();
  gettimeofday(&end, 0);

  if (nofEvents == 0) return;

  fRootFout->cd();
  fRootTree->Print();
  fRootTree->Write();

  energyHistoGe->Write();
  energyHistoGeWithVeto->Write();
  energyHistoVeto->Write();
  MultiplicityGe->Write();
  MultiplicityGeWithVeto->Write();
  MultiplicityVeto->Write();

  fRootFout->Close();

  G4cout <<nofEvents<<" Events in "<< end.tv_sec-start.tv_sec<<"s"<< G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::Fill()
{
  Int_t m_veto = 0;
  Int_t m_ge = 0;
  Int_t m_geWithVeto = 0;
  det_channel.clear();
  det_energy.clear();
  start_energy = (primary->fParticleGun->GetParticleEnergy())/keV;    //convert from MeV to keV
  start_phi = primary->fParticleGun->GetParticleMomentumDirection().getPhi();
  start_theta = primary->fParticleGun->GetParticleMomentumDirection().getTheta();
  start_x = primary->fParticleGun->GetParticlePosition().x();
  start_y = primary->fParticleGun->GetParticlePosition().y();
  start_z = primary->fParticleGun->GetParticlePosition().z();


  //multiplicity
  //Ge
  for (G4int k=100;k<300;k++){
	if (energy[k]>0){
	  m_ge++;
	  det_energy.push_back((double)(energy[k]));
	  det_channel.push_back((int)k);
	}
  }
  //Veto
  for (G4int k=0;k<100;k++){
	if (energy[k]>0){
	  m_veto++;
    det_energy.push_back((double)(energy[k]));
	  det_channel.push_back((int)k);
	}
  }
  //together
  if (m_veto>0) m_geWithVeto = m_ge;


  if (m_ge>0) MultiplicityGe->Fill(m_ge);
  if (m_geWithVeto>0) MultiplicityGeWithVeto->Fill(m_geWithVeto);
  if (m_veto>0) MultiplicityVeto->Fill(m_veto);


  //energy
  for (G4int k=100;k<300;k++){
  	if (energy[k]>0) energyHistoGe->Fill(energy[k]);
  	if ((energy[k]>0)&&(m_veto>0)) energyHistoGeWithVeto->Fill(energy[k]);
  }
  for (G4int k=0;k<100;k++){
	  if (energy[k]>0) energyHistoVeto->Fill(energy[k]);
  }

  //if (det_energy.size()>0) fRootTree->Fill();
  fRootTree->Fill();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
double RunAction::GetTimeDiff(){
	gettimeofday(&now, 0);
	double timediff;
	timediff = now.tv_sec-start.tv_sec;
	return timediff;
}

