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
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
   // use of stepping action to set the accounting volume

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Polyhedra.hh"
#include "G4Polycone.hh"
#include "G4SubtractionSolid.hh"

#include "G4Torus.hh"
#include "G4UnionSolid.hh"

#include "G4String.hh"
#include "math.h"

#include "G4VisAttributes.hh"
#include "G4Color.hh"


#ifndef M_PI
#define M_PI    3.14159265358979323846f
#endif


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction()
{
  detectorMessenger = new DetectorMessenger(this);
  ShieldInner = true;
  ShieldOuter = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
 delete detectorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
  G4Material* cupper = nist->FindOrBuildMaterial("G4_Cu");
  G4Material* lead = nist->FindOrBuildMaterial("G4_Pb");
  G4Material* iron = nist->FindOrBuildMaterial("G4_Fe");
  G4Material* chromium = nist->FindOrBuildMaterial("G4_Cr");
  G4Material* nickel = nist->FindOrBuildMaterial("G4_Ni");
  G4Material* Det_mat = nist->FindOrBuildMaterial("G4_Ge");

  G4Material* carbon = nist->FindOrBuildMaterial("G4_C");
  G4Material* hydrogen = nist->FindOrBuildMaterial("G4_H");
  G4Material* mat_Silicon = nist->FindOrBuildMaterial("G4_Si");
  G4Material* mat_Oxygen = nist->FindOrBuildMaterial("G4_O");
  G4Material* mat_Uranium = nist->FindOrBuildMaterial("G4_U");
  G4Material* mat_Thorium = nist->FindOrBuildMaterial("G4_Th");
  G4Material* mat_Potassium = nist->FindOrBuildMaterial("G4_K");


  G4String symbol;             //a=mass of a mole;
  G4double a, z, density;      //z=mean number of protons
  G4Element* H  = new G4Element("Hydrogen",symbol="H",  z= 1, a=   1.00794*g/mole);
  G4Element* C  = new G4Element("Carbon",  symbol="C",  z= 6, a=  12.0107*g/mole);

  G4Material* scint = new G4Material("scint", 1.03*g/cm3, 2, kStateSolid);
  scint-> AddMaterial(carbon,12/13.);
  scint-> AddMaterial(hydrogen,1/13.);
  //scint->AddElement(G4Element::GetElement("Carbon"), 1);
  //scint->AddElement(G4Element::GetElement("Hydrogen"), 1);

  G4Material* mat_Rock= new G4Material("Rock", 2.86 * g / cm3, 5);
  G4double sum_minors = 0.22/1000000. + 0.33/1000000. + 0.96/100.;
  mat_Rock->AddMaterial(mat_Silicon, 1./3.*(1-sum_minors));
  mat_Rock->AddMaterial(mat_Oxygen, 2./3.*(1-sum_minors));
  mat_Rock->AddMaterial(mat_Uranium, 0.22/1000000.);
  mat_Rock->AddMaterial(mat_Thorium, 0.33/1000000.);
  mat_Rock->AddMaterial(mat_Potassium, 0.96/100.);


  G4cout << "RadL " << scint->GetRadlen()/cm
         << " Nucl.Int.Length "<< scint->GetNuclearInterLength()/cm
         << " Imean "<< (scint->GetIonisation()->GetMeanExcitationEnergy())/eV << G4endl;


  //
  // World
  //
  G4Box* solidWorld = new G4Box("World",50*m,50*m,50*m);
  G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld,world_mat,"World");
  G4VPhysicalVolume* physWorld = new G4PVPlacement(0,G4ThreeVector(),logicWorld,"World",0,false,0,checkOverlaps);
  //
  // build lab space
  //
	G4Box* solid_Rock = new G4Box("sol_Rock",50*m,50*m,50*m);
	G4Box* solid_Lab = new G4Box("sol_Lab",5*m,5*m,5*m);
  G4SubtractionSolid *solid_Rock2 = new G4SubtractionSolid("sol_Rock2", solid_Rock, solid_Lab);
  G4LogicalVolume* logical_Rock = new G4LogicalVolume(solid_Rock2,mat_Rock,"log_Rock");
	logical_Rock->SetVisAttributes ( new G4VisAttributes(G4Colour(0.7, 0.7, 0.7, 0.5) )); //grey 50% transparent
  G4VPhysicalVolume* physical_Rock = new G4PVPlacement(0,G4ThreeVector(),logical_Rock,"phy_Rock",logicWorld,false,0,checkOverlaps);

  
  //read detector file
  using namespace std;
  vector<string> Detectors;
  string input;
  string s;

  G4int Detector_number;
  G4String Detector_name;
  G4String Detector_name_sol;
  G4String Detector_name_log;
  G4ThreeVector Detector_position;
  G4int Detector_slices;
  vector<G4double> Detector_r;
  vector<G4double> Detector_z;
  G4RotationMatrix* Detrotation = new G4RotationMatrix();
  Detrotation->rotateX(180*deg);
  
  int j =0;

  ifstream in ("Detectorposition.txt");
  while(getline(in,input)){
    Detectors.push_back(input);
  }
  in.close();
  G4cout << "-------" << Detectors.size() <<G4endl;
  
  for (int i =0;i< Detectors.size();i++){
    G4cout << i << " " << Detectors[i] << endl;
    j = 0;
    Detector_number = 0;
    Detector_r.clear();
    Detector_z.clear();
    stringstream ss(Detectors[i]);
    while (ss >> s) {
     j++;
     if (j==1) Detector_number += atoi(s.c_str())*100;
     else if (j==2) Detector_number += atoi(s.c_str())*10;
     else if (j==3) Detector_number += atoi(s.c_str())*1;
     else if (j==4) Detector_position.setX(atof(s.c_str())*mm);
     else if (j==5) Detector_position.setY(atof(s.c_str())*mm);
     else if (j==6) Detector_position.setZ(atof(s.c_str())*mm);
     else if (j==7) Detector_name = "Det_" + s;
     else if (j==8) Detector_slices = atoi(s.c_str());
     else if (j==9) Detector_r.push_back((G4double) atof(s.c_str())*mm);
     else if (j==10) Detector_z.push_back((G4double) atof(s.c_str())*mm);
     else if (j==11) Detector_r.push_back((G4double) atof(s.c_str())*mm);
     else if (j==12) Detector_z.push_back((G4double) atof(s.c_str())*mm);
     else if (j==13) Detector_r.push_back((G4double) atof(s.c_str())*mm);
     else if (j==14) Detector_z.push_back((G4double) atof(s.c_str())*mm);

    }
    
    if(Detector_slices==2){
      Detector_r.push_back(Detector_r[1]+0.00001*mm);
      Detector_z.push_back(Detector_z[1]+0.00001*mm);
    }

    Detector_position.setZ(Detector_position.z() + Detector_z[2] - 20*cm);

    G4cout << j << " --- " << Detector_name << " " << Detector_number << " " << Detector_position <<  G4endl;         
    
    Detector_name_sol = Detector_name + "_sol";
    Detector_name_log = Detector_name + "_log";
    const G4double r_i[] = {0,0,0};
    const G4double r[] = {Detector_r[0],Detector_r[1],Detector_r[2]};
    const G4double z[] = {Detector_z[0],Detector_z[1],Detector_z[2]};
    
    G4Polycone* Det_solid = new G4Polycone(Detector_name_sol,0,2*M_PI,3,z,r_i,r);
    G4LogicalVolume* Det_logical = new G4LogicalVolume(Det_solid,Det_mat,Detector_name_log);
    new G4PVPlacement (Detrotation,Detector_position,Det_logical,Detector_name,logicWorld,false,Detector_number,0);
  }



  //top
  G4double zConePlanes1[] = {0, 0.07*25.4, 0.07*25.4, 0.18*25.4, 0.18*25.4,
    						0.319*25.4, 0.319*25.4, 0.511*25.4, 0.511*25.4, 0.9263*25.4};
                            //0.32, 0.51
  G4double rConeInner1[] = {6.6582*25.4, 6.6459*25.4, 6.555*25.4, 6.555*25.4,
   							6.5*25.4, 6.5*25.4, 6.5*25.4, 6.5*25.4,
    							6.5*25.4, 6.5*25.4};
  G4double rConeOuter1[] = {6.759*25.4, 6.759*25.4, 6.759*25.4, 6.759*25.4, 6.759*25.4,
    						 6.759*25.4, 6.695*25.4, 6.695*25.4, 6.76*25.4, 6.76*25.4};//6.76
  G4Polycone* barrel1 = new G4Polycone("Barrel", 0, 2*M_PI, 10, zConePlanes1, rConeInner1, rConeOuter1);

  G4Sphere* top = new G4Sphere("top", 53.0125*25.4, 53.2725*25.4, 0, 2*M_PI, 0, 6.69*deg);
  G4UnionSolid* body1 = new G4UnionSolid("body1", barrel1, top, 0,
    										G4ThreeVector(0, 0, -51.3625*25.4));
  G4LogicalVolume* toplid_logical = new G4LogicalVolume(body1, cupper, "toplid_log");
  new G4PVPlacement (0,G4ThreeVector(-20.6*cm,0,0*cm),toplid_logical,"toplidM1_phys",logicWorld,false,0,checkOverlaps);
  new G4PVPlacement (0,G4ThreeVector(20.6*cm,0,0*cm),toplid_logical,"toplidM2_phys",logicWorld,false,0,checkOverlaps);



 //middle
 G4Tubs* cryo2_solid = new G4Tubs("cryo2_sol",16.5*cm,17.4625*cm,7.62*cm,0,2*M_PI);
 G4LogicalVolume* cryo2_logical = new G4LogicalVolume(cryo2_solid,cupper,"cryo2_log");
 new G4PVPlacement (0,G4ThreeVector(-20.6*cm,0,-7.62*cm),cryo2_logical,"cryo2M1_phys",logicWorld,false,0,checkOverlaps);
 new G4PVPlacement (0,G4ThreeVector(20.6*cm,0,-7.62*cm),cryo2_logical,"cryo2M2_phys",logicWorld,false,0,checkOverlaps);

//bottom
   G4double zConePlanes[] = {0*25.4, -0.1*25.4, -0.1*25.4, -0.18*25.4, -0.18*25.4,
  							-0.499*25.4, -0.499*25.4, -0.691*25.4, -0.691*25.4, -8.8183*25.4};//-0.5,-0.69
  G4double rConeInner[] = {6.5*25.4, 6.5*25.4, 6.5*25.4, 6.5*25.4, 6.5*25.4,
    							6.5*25.4, 6.5*25.4, 6.5*25.4, 6.5*25.4, 6.5*25.4};
  G4double rConeOuter[] = {6.55*25.4, 6.55*25.4, 6.6459*25.4, 6.66*25.4, 6.759*25.4, //6.76
   							6.759*25.4, 6.695*25.4, 6.695*25.4, 6.76*25.4, 6.76*25.4};
  G4Polycone* barrel = new G4Polycone("Barrel", 0, 2*M_PI, 10, zConePlanes, rConeInner, rConeOuter);
  G4Sphere* bottom = new G4Sphere("bottom", 53.0125*25.4, 53.2725*25.4, 0, 2*M_PI, 0, 6.69*deg);
  G4RotationMatrix* lidrotation = new G4RotationMatrix();
  lidrotation->rotateX(180*deg);
  G4UnionSolid* body2 = new G4UnionSolid("body2", barrel, bottom, lidrotation,
    										G4ThreeVector(0, 0, 43.4705*25.4));

  G4LogicalVolume* bottomlid_logical = new G4LogicalVolume(body2, cupper, "bottomlid_log");
  new G4PVPlacement (0,G4ThreeVector(-20.6*cm,0,-15.28*cm),bottomlid_logical,"bottomlidM1_phys",logicWorld,false,0,checkOverlaps);
  new G4PVPlacement (0,G4ThreeVector(20.6*cm,0,-15.28*cm),bottomlid_logical,"bottomlidM2_phys",logicWorld,false,0,checkOverlaps);

  //coldplate
  G4Tubs* coldplate_solid = new G4Tubs("coldplate_sol",0*cm,15.5*cm,1*cm,0,2*M_PI);
  G4LogicalVolume* coldplate_logical = new G4LogicalVolume(coldplate_solid,cupper,"coldplate_log");
  new G4PVPlacement (0,G4ThreeVector(20.6*cm,0,-10*cm),coldplate_logical,"coldplate_M1phys",logicWorld,false,0,0);
  new G4PVPlacement (0,G4ThreeVector(-20.6*cm,0,-10*cm),coldplate_logical,"coldplate_M2phys",logicWorld,false,0,0);
  

  if (ShieldOuter){
    //outer copper
    G4Box* cbody = new G4Box("cbody", 20*25.4, 12.0*25.4, 14*25.4);
    G4Box* pocket = new G4Box("pocket", 18.031*25.4, 10.031*25.4, 12.031*25.4);
    G4SubtractionSolid* cbody2 = new G4SubtractionSolid("cbody2", cbody, pocket);

    G4Tubs* tubeCut = new G4Tubs("tubeCut", 0, 1.751*25.4, 20*25.4, 0, 2*M_PI);//1.75
    G4RotationMatrix* rotation = new G4RotationMatrix();
    rotation->rotateY(M_PI/2);
    G4SubtractionSolid* cbody3 = new G4SubtractionSolid("cbody3", cbody2, tubeCut, rotation,
    									G4ThreeVector(30*25.4, 0*25.4, 4.4265*25.4));
    G4RotationMatrix* rotation2 = new G4RotationMatrix();
    rotation2->rotateX(M_PI/2);
    G4SubtractionSolid* shield = new G4SubtractionSolid("shield", cbody3, tubeCut, rotation2,
    						G4ThreeVector(-8.1417*25.4, -20*25.4, 4.4265*25.4));

    G4LogicalVolume* outercopper_logical = new G4LogicalVolume(shield, cupper, "outercopper_log");
    new G4PVPlacement (0,G4ThreeVector(0,0,-15.18*cm),outercopper_logical,"outercopper_phys",logicWorld,false,0,checkOverlaps);
  }

  if (ShieldInner){
    //inner copper
    G4Box* cibody = new G4Box("cibody", 18*25.4, 9.875*25.4, 11.923*25.4);
    G4Box* cipocket = new G4Box("cipocket", 16.031*25.4, 7.906*25.4, 9.954*25.4);
    G4SubtractionSolid* cibody2 = new G4SubtractionSolid("cibody2", cibody, cipocket);

    G4Tubs* citubeCut = new G4Tubs("citubeCut", 0, 1.751*25.4, 20*25.4, 0, 2*M_PI);//1.75
	G4RotationMatrix* cirotation = new G4RotationMatrix();
    cirotation->rotateY(M_PI/2);
    G4SubtractionSolid* cibody3 = new G4SubtractionSolid("cibody3", cibody2, citubeCut, cirotation,
    									G4ThreeVector(30*25.4, 0*25.4, 4.4265*25.4));
    G4RotationMatrix* cirotation2 = new G4RotationMatrix();
    cirotation2->rotateX(M_PI/2);
    G4SubtractionSolid* cishield = new G4SubtractionSolid("cishield", cibody3, citubeCut, cirotation2,
    						G4ThreeVector(-8.1417*25.4, -20*25.4, 4.4265*25.4));
	G4LogicalVolume* innercopper_logical = new G4LogicalVolume(cishield, cupper, "innercopper_log");
    new G4PVPlacement (0,G4ThreeVector(0,0,-15.18*cm),innercopper_logical,"innercopper_phys",logicWorld,false,0,checkOverlaps);
  }


    G4Box* pbbody = new G4Box("pbbody", 38.0*25.4, 30.0*25.4, 32.0*25.4);

    G4Box* pbpocket = new G4Box("pbpocket", 20.001*25.4, 12.001*25.4, 14.001*25.4);
    G4SubtractionSolid* pbbody2 = new G4SubtractionSolid("pbbody2", pbbody, pbpocket);
    G4Tubs* pbtubeCut = new G4Tubs("pbtubeCut", 0, 1.751*25.4, 20*25.4, 0, 2*M_PI);//1.75
    G4RotationMatrix* pbrotation = new G4RotationMatrix();
    pbrotation->rotateY(M_PI/2);
    G4SubtractionSolid* pbbody7 = new G4SubtractionSolid("pbbody7", pbbody2, pbtubeCut,pbrotation,
    									G4ThreeVector(30*25.4, 0*25.4, 4.4265*25.4));
    G4RotationMatrix* pbrotation2 = new G4RotationMatrix();
    pbrotation2->rotateX(M_PI/2);
    G4SubtractionSolid* pbleadShield = new G4SubtractionSolid("pbleadShield", pbbody7, pbtubeCut, pbrotation2,
    						G4ThreeVector(-8.1417*25.4, -20*25.4, 4.4265*25.4));
	G4LogicalVolume* pbleadshield_logical = new G4LogicalVolume(pbleadShield, lead, "pbleadshield_log");
    new G4PVPlacement (0,G4ThreeVector(0,0,-15.18*cm),pbleadshield_logical,"pbleadshield_phys",logicWorld,false,0,checkOverlaps);


  for (G4int k=1;k<33;k++){
	//G4cout << k << G4endl;
    G4double panelsize[3];
	panelsize[0]=0;
    if (k<7){
      panelsize[0]= 0.5*97*25.4;
      panelsize[1]= 0.5*13.7*25.4;
      panelsize[2]= 0.5*1*25.4;
    }
    if ((k>6)&&(k<13)){
      panelsize[0]= 0.5*13.7*25.4;
      panelsize[1]= 0.5*97*25.4;
      panelsize[2]= 0.5*1*25.4;
    }
    if ((k==13)||(k==14)||(k==15)||(k==23)){
	  panelsize[0]= 0.5*1*25.4;
      panelsize[1]= 0.5*34.5*25.4;
      panelsize[2]= 0.5*70.25*25.4;
	}
    if ((k==17)||(k==24)||(k==26)||(k==27)){
	  panelsize[0]= 0.5*33.5*25.4;
      panelsize[1]= 0.5*1*25.4;
      panelsize[2]= 0.5*70.25*25.4;
	}
    if ((k==16)||(k==20)||(k==25)||(k==28)){
	  panelsize[0]= 0.5*47.5*25.4;
      panelsize[1]= 0.5*1*25.4;
      panelsize[2]= 0.5*70.25*25.4;
	}
	if ((k==30)||(k==31)){
	  panelsize[0]= 0.5*1*25.4;
      panelsize[1]= 0.5*29.5*25.4;
      panelsize[2]= 0.5*70.25*25.4;
	}
    if ((k==29)||(k==32)){
	  panelsize[0]= 0.5*1*25.4;
      panelsize[1]= 0.5*37.5*25.4;
      panelsize[2]= 0.5*70.25*25.4;
	}
	if ((k==22)||(k==21)){
	  panelsize[0]= 0.5*85*25.4;
      panelsize[1]= 0.5*34*25.4;
      panelsize[2]= 0.5*1*25.4;
	}
	if ((k==19)||(k==18)){
	  panelsize[0]= 0.5*41*25.4;
      panelsize[1]= 0.5*69*25.4;
      panelsize[2]= 0.5*1*25.4;
	}

	if (!(panelsize[0]>0)) continue;
    G4Box* panel_solid = new G4Box("panel_sol", panelsize[0], panelsize[1], panelsize[2]);
    G4LogicalVolume* panel_logical = new G4LogicalVolume(panel_solid,scint,"panel_log");

	G4ThreeVector panelposition;
//floor
	if (k<7)           panelposition = G4ThreeVector(0,25.4*(0.15+6.85+(k-4)*(13.7+0.3)),-111.55*cm);
	if ((k>6)&&(k<13)) panelposition = G4ThreeVector(25.4*(0.15+6.85+(k-10)*(13.7+0.3)),0,-108*cm);
//sides
	if (k==13) 		   panelposition = G4ThreeVector(-40.6*25.4,+panelsize[1],-15.18*cm);
	if (k==14) 		   panelposition = G4ThreeVector(-40.6*25.4,-panelsize[1],-15.18*cm);
	if (k==15) 		   panelposition = G4ThreeVector(-41.6*25.4,+panelsize[1],-15.18*cm);
	if (k==23) 		   panelposition = G4ThreeVector(-41.6*25.4,-panelsize[1],-15.18*cm);
	if (k==24) 		   panelposition = G4ThreeVector(-7.5*25.4-panelsize[0],+35*25.4,-15.18*cm);
	if (k==17) 		   panelposition = G4ThreeVector(-7.5*25.4-panelsize[0],+36*25.4,-15.18*cm);
	if (k==16) 		   panelposition = G4ThreeVector(-7.5*25.4+panelsize[0],+36*25.4,-15.18*cm);
	if (k==20) 		   panelposition = G4ThreeVector(-7.5*25.4+panelsize[0],+35*25.4,-15.18*cm);
	if (k==25) 		   panelposition = G4ThreeVector(7*25.4-panelsize[0],-35*25.4,-15.18*cm);
	if (k==27) 		   panelposition = G4ThreeVector(7*25.4+panelsize[0],-35*25.4,-15.18*cm);
	if (k==26) 		   panelposition = G4ThreeVector(-7.5*25.4-panelsize[0],-36*25.4,-15.18*cm);
	if (k==28) 		   panelposition = G4ThreeVector(-7.5*25.4+panelsize[0],-36*25.4,-15.18*cm);
	if (k==31) 		   panelposition = G4ThreeVector(40.6*25.4,4*25.4+panelsize[1],-15.18*cm);
	if (k==29) 		   panelposition = G4ThreeVector(40.6*25.4,4*25.4-panelsize[1],-15.18*cm);
	if (k==32) 		   panelposition = G4ThreeVector(41.6*25.4,-4*25.4+panelsize[1],-15.18*cm);
	if (k==30) 		   panelposition = G4ThreeVector(41.6*25.4,-4*25.4-panelsize[1],-15.18*cm);
//top
  if (k==22) 		   panelposition = G4ThreeVector(0,+panelsize[1],75.4*cm);
	if (k==21) 		   panelposition = G4ThreeVector(0,-panelsize[1],75.4*cm);
	if (k==19) 		   panelposition = G4ThreeVector(panelsize[0],0,77.94*cm);
	if (k==18) 		   panelposition = G4ThreeVector(-panelsize[0],0,77.94*cm);


	//G4cout << k << " " << panelsize[0] << " " <<panelposition<< G4endl;
	new G4PVPlacement (0,panelposition,panel_logical,"panel_phys",logicWorld,false,k,0);
  }

  G4Sphere* solidlab = new G4Sphere("lab",2.99*m,3*m,0,2*M_PI,0,0.5*M_PI);
  G4LogicalVolume* logiclab = new G4LogicalVolume(solidlab,world_mat,"lab");
  G4VPhysicalVolume* physlab = new G4PVPlacement(0,G4ThreeVector(0,0,-1*m),logiclab,"lab_phys",logicWorld,false,0,checkOverlaps);

  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#include "G4RunManager.hh"

void DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetInnerShield(G4bool a){
  ShieldInner = a;
  UpdateGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetOuterShield(G4bool a){
  ShieldOuter = a;
  UpdateGeometry();
}
