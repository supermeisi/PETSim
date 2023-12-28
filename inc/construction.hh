#ifndef CONSTRUCTION_HH
#define CONSTRUCTION_HH

#include "G4SystemOfUnits.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4PVPlacement.hh"

class MyDetectorConstruction : public G4VUserDetectorConstruction
{
public:
    MyDetectorConstruction();
    ~MyDetectorConstruction();

    G4LogicalVolume *GetScoringVolume() const { return fScoringVolume; }

    virtual G4VPhysicalVolume *Construct();

private:
    G4Box *solidWorld;
    G4LogicalVolume *logicWorld;
    G4VPhysicalVolume *physWorld;
    
    G4LogicalVolume *fScoringVolume;

    G4double xWorld;
    G4double yWorld;
    G4double zWorld;

    G4Material *worldMat;
    
    void DefineMaterials();
};

#endif