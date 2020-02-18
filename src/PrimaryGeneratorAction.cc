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
// $Id$
//
/// \file PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class

#include "PrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction* PrimaryGeneratorAction::fgInstance = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const PrimaryGeneratorAction* PrimaryGeneratorAction::Instance()
{
// Static acces function via G4RunManager

  return fgInstance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(0)
{
//  fParticleGun = new G4GeneralParticleSource;
G4int n_particle = 1;
fParticleGun  = new G4ParticleGun(n_particle);

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
  fgInstance = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  G4ParticleDefinition* particle = G4ParticleTable::GetParticleTable()->FindParticle("mu-");
  fParticleGun->SetParticleDefinition(particle);
  double MuonEnergy = 1.0*GeV+199999.0*G4UniformRand()*GeV;  //1GeV-200TeV (1x10^9-200x10^12)

  fParticleGun->SetParticleEnergy(MuonEnergy);

//  fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., 20*m));
//  G4ThreeVector dir(0,0,-1.);
//  fParticleGun->SetParticleMomentumDirection(dir);

  // particle must go through the floor of the hemisphere
  double radius = 10.0*m;
  double theLowerX = radius-2*radius*G4UniformRand();
  double theLowerY = radius-2*radius*G4UniformRand();
  double theLowerZ = 0.0*m;
  //G4cout << "LowerX= " << theLowerX << " " << "LowerY= " << theLowerY << " "<< "LowerZ= " << theLowerZ << " " << G4endl;
  double outsideBound = std::pow(theLowerX,2) + std::pow(theLowerY,2); //eq for the circle inside box
  while(outsideBound > std::pow(radius,2)) { //if the point falls outside the circle, select a new one
    theLowerX = radius-2*radius*G4UniformRand();
    theLowerY = radius-2*radius*G4UniformRand();
    outsideBound = std::pow(theLowerX,2) + std::pow(theLowerY,2);
  }

  // create particle on surface of hemisphere of size radius, that travels through the bottom surface
  // pick bottom point inside a circle
  ////starting position on the hemisphere
  double theUpperPhi = 3.14159-2*3.1415*G4UniformRand();
  double theUpperTheta = acos(G4UniformRand());
  double theUpperX = radius * cos(theUpperPhi) * sin(theUpperTheta);
  double theUpperY = radius * sin(theUpperPhi) * sin(theUpperTheta);
  double theUpperZ = radius * cos(theUpperTheta)+theLowerZ;
  //G4cout << "UpperX= " << theUpperX << " " << "UpperY= " << theUpperY << " "<< "UpperZ= " << theUpperZ << " " << G4endl;
  fParticleGun->SetParticlePosition(G4ThreeVector(theUpperX, theUpperY, theUpperZ));
  //fParticleGun->SetParticlePosition(G4ThreeVector(206.799*CLHEP::mm,0,-262.114*CLHEP::mm));  //on C2P1D3 for tests


  //compute track direction (note left handed coord system)
  double Xtrk = theUpperX-theLowerX;
  double Ytrk = theUpperY-theLowerY;
  double Ztrk = theUpperZ-theLowerZ;

  double Rtrk = std::sqrt(Xtrk*Xtrk + Ytrk*Ytrk + Ztrk*Ztrk);
  double ThetaTrk=acos(Ztrk/Rtrk);
  double PhiTrk=0.;
  if (Xtrk == 0){
    if (Ytrk > 0) PhiTrk = 1.5708;    //90 degrees
    if (Ytrk < 0) PhiTrk = 4.7124;    //270 degrees
  }
  else{
    PhiTrk=atan(std::abs(Ytrk)/std::abs(Xtrk));
    if ((Xtrk < 0) && (Ytrk > 0)) PhiTrk = 3.14159 - PhiTrk;
    if ((Xtrk < 0) && (Ytrk <= 0)) PhiTrk = 3.14159 + PhiTrk;
    if ((Xtrk > 0) && (Ytrk < 0)) PhiTrk = 6.2832 - PhiTrk;
  }
  //G4cout << "ThetaTrk= " << ThetaTrk << " PhiTrk= " << PhiTrk << G4endl << G4endl;
  //set direction of particle (inward)
  G4ThreeVector dir(-std::sin(ThetaTrk)*std::cos(PhiTrk),-std::sin(ThetaTrk)*std::sin(PhiTrk),-std::cos(ThetaTrk));
  fParticleGun->SetParticleMomentumDirection(dir);

  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
