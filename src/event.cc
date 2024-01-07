#include "event.hh"

MyEventAction::MyEventAction(MyRunAction *)
{
    fEdep = 0.;
}

MyEventAction::~MyEventAction()
{
}

void MyEventAction::BeginOfEventAction(const G4Event *)
{
    fEdep = 0.;
}

void MyEventAction::EndOfEventAction(const G4Event *)
{
    G4cout << "Energy at " << vPos << " deposition: " << fEdep << G4endl;
    
    G4AnalysisManager *man = G4AnalysisManager::Instance();
}
