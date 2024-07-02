#include "simStepping.hh"

simSteppingAction::simSteppingAction(simEventAction *eventAction)
{
    fEventAction = eventAction; // get access to the object we created
}

simSteppingAction::~simSteppingAction() {}

void simSteppingAction::UserSteppingAction(const G4Step *step)
{
    G4LogicalVolume *volume = step->GetPreStepPoint()
                                  ->GetTouchableHandle()
                                  ->GetVolume()
                                  ->GetLogicalVolume();

    // check scoring volume
    const simDetectorConstruction *detectorConstruction =
        static_cast<const simDetectorConstruction *>(
            G4RunManager::GetRunManager()->GetUserDetectorConstruction());

    G4LogicalVolume *fScoringVolume = detectorConstruction->GetScoringVolume();

    // G4cout << volume->GetName() << G4endl;

    if (volume != fScoringVolume)
        return;

    G4double edep = step->GetTotalEnergyDeposit();
    fEventAction->AddEdep(edep); // then accumulate
}
