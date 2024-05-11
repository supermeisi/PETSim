#include "stepping.hh"

MySteppingAction::MySteppingAction(MyEventAction *eventAction)
{
    fEventAction = eventAction;
}

MySteppingAction::~MySteppingAction()
{}

void MySteppingAction::UserSteppingAction(const G4Step *step)
{   
    G4AnalysisManager *man = G4AnalysisManager::Instance();
       
    G4LogicalVolume *volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
    
    const MyDetectorConstruction *detectorConstruction = static_cast<const MyDetectorConstruction*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    
    G4LogicalVolume *fScoringVolume = detectorConstruction->GetScoringVolume();
    
    if(volume != fScoringVolume)
        return;
    
    G4double edep = step->GetTotalEnergyDeposit(); 
    G4StepPoint *preStepPoint = step->GetPreStepPoint();
    G4ThreeVector pos = preStepPoint->GetPosition();

    man->FillH2(0, pos[0], pos[1], edep);

    std::fstream outfile;
    outfile.open("energy.txt", std::ios::app);
    outfile << pos[0] << " " << pos[1] << " " << pos[2] << " " << edep << std::endl;
    outfile.close();

    fEventAction->AddEdep(edep);
    fEventAction->AddPosition(pos);
}
