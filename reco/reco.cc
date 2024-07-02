#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TChain.h"
#include "TH2F.h"
#include "TStyle.h"
#include "TVector3.h"

struct Hits
{
    double time;
    TVector3 pos;
    double energy;
};

// Function for reconstructing the Line of Response (LOR)
TVector3 Reconstruction(Hits hits1, Hits hits2)
{
    if (hits1.energy < 0.4 || hits1.energy > 0.6 ||
        hits2.energy < 0.4 || hits2.energy > 0.6)
    {
        return TVector3(sqrt(-1), sqrt(-1), sqrt(-1));
    }

    double dt = hits1.time - hits2.time;

    if (fabs(dt) > 0.2)
        return TVector3(sqrt(-1), sqrt(-1), sqrt(-1));

    double dist = 0.5 * (299.792458 * dt);

    TVector3 midPoint = 0.5 * (hits2.pos + hits1.pos);

    TVector3 LOR = (hits2.pos - hits1.pos).Unit();

    midPoint.Print();

    std::cout << "Time difference: " << dt << std::endl;
    std::cout << "Distance: " << dist << std::endl;

    TVector3 annihilationPoint = midPoint + LOR * dist;

    annihilationPoint.Print();

    return annihilationPoint;
}

int main(int argv, char **argc)
{
    gStyle->SetPalette(kRainBow);

    const int nRuns = 70;
    const int nFiles = 4;

    TH2F *hPos =
        new TH2F("hPos", "Recostructed Position in XY Plane;X [mm]; Y[mm]", 100,
                 -30, 30, 100, -30, 30);

    hPos->SetMinimum(-1);

    TCanvas *cRun = new TCanvas();

    for (int iRun = 0; iRun < nRuns; iRun++)
    {
        TChain *fHits = new TChain("Hits");

        for (Int_t iFile = 0; iFile < nFiles; iFile++)
        {
            fHits->AddFile(Form("../sim/build/output%d_t%d.root", iRun, iFile));
        }

        Int_t entries = fHits->GetEntries();

        std::cout << "Number of hits: " << entries << std::endl;

        int event;
        double time;
        double x, y, z;
        double energy;

        fHits->SetBranchAddress("fEvent", &event);
        fHits->SetBranchAddress("fTime", &time);
        fHits->SetBranchAddress("fEnergy", &energy);
        fHits->SetBranchAddress("fX", &x);
        fHits->SetBranchAddress("fY", &y);
        fHits->SetBranchAddress("fZ", &z);

        double oldGlobalTime = 0.;
        TVector3 oldPos(0., 0., 0.);

        std::vector<Hits> hits;

        int newEvent, oldEvent;

        for (int iEntry = 0; iEntry < entries; iEntry++)
        {
            fHits->GetEntry(iEntry);

            newEvent = event;

            // Starting reconstruction for each new event
            if (oldEvent != newEvent)
            {
                // Reconstructing LOR for all hit combinations

                for (int iHit = 0; iHit < hits.size(); iHit++)
                {
                    for (int jHit = iHit + 1; jHit < hits.size(); jHit++)
                    {
                        TVector3 annihilationPoint =
                            Reconstruction(hits[iHit], hits[jHit]);

                        // Only fill useful hits into histogram
                        if (annihilationPoint.X() == annihilationPoint.X() &&
                            annihilationPoint.Y() == annihilationPoint.Y() &&
                            annihilationPoint.Z() == annihilationPoint.Z())
                        {
                            annihilationPoint.Print();

                            hPos->Fill(annihilationPoint.X(),
                                       annihilationPoint.Y());
                        }
                    }
                }

                oldEvent = newEvent;
                hits.clear();
            }

            TVector3 pos(x, y, z);

            Hits hit;
            hit.time = time;
            hit.pos = pos;
            hit.energy = energy;

            hits.push_back(hit);
        }

        cRun->cd();
        hPos->Draw("colz");
        cRun->Print("pos.pdf");
        cRun->Print("pos.png");
    }
}
