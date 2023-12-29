#include "construction.hh"

MyDetectorConstruction::MyDetectorConstruction()
{
    xWorld = 50. * cm;
    yWorld = 50. * cm;
    zWorld = 50. * cm;

    length = 20. * cm;

    radius = 40 * cm;

    xCryst = 6 * mm;
    yCryst = 6 * mm;
    zCryst = 6 * mm;

    xDet = 1. * mm;
    yDet = yCryst;
    zDet = zCryst;

    gap = 1. * mm;

    DefineMaterials();
}

MyDetectorConstruction::~MyDetectorConstruction()
{
}

void MyDetectorConstruction::DefineMaterials()
{
    G4NistManager *nist = G4NistManager::Instance();

    worldMat = nist->FindOrBuildMaterial("G4_AIR");

    G4double energy[2] = {1.239841939 * eV / 0.9, 1.239841939 * eV / 0.2};

    G4double reflectivity[2] = {1.0, 1.0};

    G4double rindexNaI[2] = {1.78, 1.78};
    G4double rindexWorld[2] = {1.0, 1.0};

    G4MaterialPropertiesTable *mptWorld = new G4MaterialPropertiesTable();
    mptWorld->AddProperty("RINDEX", energy, rindexWorld, 2);
    worldMat->SetMaterialPropertiesTable(mptWorld);

    Na = nist->FindOrBuildElement("Na");
    I = nist->FindOrBuildElement("I");
    NaI = new G4Material("NaI", 3.67 * g / cm3, 2);
    NaI->AddElement(Na, 1);
    NaI->AddElement(I, 1);

    G4double fraction[2] = {1.0, 1.0};

    G4MaterialPropertiesTable *mptNaI = new G4MaterialPropertiesTable();
    mptNaI->AddProperty("RINDEX", energy, rindexNaI, 2);
    mptNaI->AddConstProperty("SCINTILLATIONYIELD", 55. / keV);
    mptNaI->AddConstProperty("RESOLUTIONSCALE", 1.0);
    NaI->SetMaterialPropertiesTable(mptNaI);

    mirrorSurface = new G4OpticalSurface("mirrorSurface");

    mirrorSurface->SetType(dielectric_metal);
    mirrorSurface->SetFinish(ground);
    mirrorSurface->SetModel(unified);

    G4MaterialPropertiesTable *mptMirror = new G4MaterialPropertiesTable();
    mptMirror->AddProperty("REFLECTIVITY", energy, reflectivity, 2);

    mirrorSurface->SetMaterialPropertiesTable(mptMirror);
}

G4VPhysicalVolume *MyDetectorConstruction::Construct()
{
    const G4double pi = 3.14159265358979323846;

    const G4bool checkOverlaps = true;

    solidWorld = new G4Box("solidWorld", xWorld, yWorld, zWorld);
    logicWorld = new G4LogicalVolume(solidWorld, worldMat, "logicWorld");
    physWorld = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), logicWorld, "physWorld", 0, false, 0, true);

    solidScintillator = new G4Box("solidScintillator", xCryst, yCryst, zCryst);
    logicScintillator = new G4LogicalVolume(solidScintillator, NaI, "logicalScintillator");

    G4LogicalSkinSurface *skin = new G4LogicalSkinSurface("skin", logicWorld, mirrorSurface);

    G4int nCrystR = (pi * radius) / (yCryst + gap);
    G4double dAngle = 360. / nCrystR;

    G4int nCrystL = length / (zCryst + gap);

    G4cout << "Number of crystals in radius: " << nCrystR << G4endl;
    G4cout << "Angle distance: " << dAngle << G4endl;

    solidDetector = new G4Box("solidDetector", xDet, yDet, zDet);
    logicDetector = new G4LogicalVolume(solidDetector, worldMat, "logicDetector");
    for (G4int i = 0; i < nCrystL; i++)
    {
        for (G4int j = 0; j < nCrystR; j++)
        {
            G4ThreeVector trans = G4ThreeVector(radius, 0., 2 * (-nCrystL / 2 + i) * (zCryst + gap));

            G4Rotate3D rotZ(j * dAngle * deg, G4ThreeVector(0, 0, 1));
            G4Translate3D transScint(trans);
            G4Transform3D transformScint = (rotZ) * (transScint);

            physScintillator = new G4PVPlacement(transformScint, logicScintillator, "physScintillator", logicWorld, false, i * 16 + j, checkOverlaps);

            G4Translate3D transDet(trans + G4ThreeVector(xCryst + xDet, 0., 0.));
            G4Transform3D transformDet = (rotZ) * (transDet);

            physDetector = new G4PVPlacement(transformDet, logicDetector, "physDetector", logicWorld, false, i * 16 + j, checkOverlaps);
        }
    }

    return physWorld;
}

void MyDetectorConstruction::ConstructSDandField()
{
    MySensitiveDetector *sensDet = new MySensitiveDetector("SensitiveDetector");

    logicDetector->SetSensitiveDetector(sensDet);
}
