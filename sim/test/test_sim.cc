// test_sim.cc
#include <gtest/gtest.h>

#include "G4LogicalVolumeStore.hh"
#include "G4MTRunManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4RunManager.hh"
#include "G4SolidStore.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "simActionInitialization.hh"
#include "simDetectorConstruction.hh"
#include "simPhysicsList.hh"

// Function to check for overlaps
bool CheckOverlaps()
{
    G4PhysicalVolumeStore *pvolStore = G4PhysicalVolumeStore::GetInstance();
    for (auto pvol : *pvolStore)
    {
        if (pvol->CheckOverlaps())
        {
            return true; // Overlaps found
        }
    }
    return false; // No overlaps found
}

// Test fixture class for the Geant4 simulation
class Geant4SimulationTest : public ::testing::Test
{
  protected:
    void SetUp() override
    {
        G4cout << "Initializing G4RunManager..." <<  G4endl;

#ifdef G4MULTITHREADED
        runManager = new G4MTRunManager;
#else
        runManager = new G4RunManager;
#endif

        G4cout << "G4RunManager initialized." <<  G4endl;

        // Set up mandatory initialization classes

        G4cout << "Setting up DetectorConstruction..." <<  G4endl;
        runManager->SetUserInitialization(new simDetectorConstruction());

        G4cout << "Setting up PhysicsList..." <<  G4endl;
        runManager->SetUserInitialization(new simPhysicsList());

        G4cout << "Setting up ActionInitialization..." <<  G4endl;
        runManager->SetUserInitialization(new simActionInitialization());

        // Initialize visualization manager
        G4cout << "Initializing G4VisExecutive..." <<  G4endl;
        visManager = new G4VisExecutive();
        visManager->Initialize();
        G4cout << "G4VisExecutive initialized." <<  G4endl;

        // Initialize UI manager
        uiManager = G4UImanager::GetUIpointer();

        // Run the simulation initialization
        uiManager->ApplyCommand("/run/initialize");
    }

    void TearDown() override
    {
        G4cout << "Cleaning up..." << G4endl;
        delete visManager;
        delete runManager;
        G4cout << "Cleaning up done." << G4endl;
    }

#ifdef G4MULTITHREADED
    G4MTRunManager *runManager;
#else
    G4RunManager *runManager;
#endif
    G4VisExecutive *visManager;
    G4UImanager *uiManager;
};

// Additional geometry tests
TEST_F(Geant4SimulationTest,
       MaterialAssignmentTest_ShouldAssignCorrectMaterials)
{
    auto logicalVolumes = G4LogicalVolumeStore::GetInstance();
    ASSERT_GT(logicalVolumes->size(), 0)
        << "Logical volume store should not be empty.";

    for (const auto &volume : *logicalVolumes)
    {
        ASSERT_NE(volume->GetMaterial(), nullptr)
            << "Material should be assigned to volume: " << volume->GetName();
    }
}

// Test case to check for overlaps in the geometry
TEST_F(Geant4SimulationTest, GeometryOverlapTest_ShouldNotHaveOverlaps)
{
    // Check for overlaps
    bool overlaps = CheckOverlaps();
    EXPECT_FALSE(overlaps) << "No overlaps should be present in the geometry.";
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
