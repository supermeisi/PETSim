#include <algorithm>
#include <filesystem>
#include <iostream>
#include <queue>
#include <vector>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"

#include "sort.hh"

void check_root_file(int run)
{
    TString input_filename =
            "sorted_output_run" + std::to_string(run) + ".root";
    TString output_filename =
            "sorted_events_run" + std::to_string(run);

    TFile *infile = new TFile(input_filename);
    TTree *sortedTree = (TTree*)infile->Get("sorted_tree");
    
    Int_t event;

    sortedTree->SetBranchAddress("fEvent", &event);

    Int_t nEntries = sortedTree->GetEntries();

    TH1F *hist = new TH1F("hist", "Event Processing", nEntries, 0, nEntries);

    for(int i = 0; i < nEntries; i++)
    {
        sortedTree->GetEntry(i);
        hist->Fill(i, event);
    }

    TCanvas *c1 = new TCanvas();
    hist->Draw();
    c1->Print(output_filename + ".pdf");
    c1->Print(output_filename + ".png");

    infile->Close();
}

int main()
{
    std::string filePattern =
        "output%i_t%j.root"; // Adjusted path based on your file
                                          // naming
    int numRuns = 70;                     // Adjust based on the number of runs
    int numThreads = 8; // Adjust based on the number of threads

    for (int run = 0; run < numRuns; ++run)
    {
        Sort sort;

        // Merge the ROOT files for this run into a single TChain
        TChain *chain = sort.merge_root_files(filePattern, run, numThreads);

        const int chunkSize = 1000; // Adjust this based on available memory
        std::vector<std::string> chunk_files;

        double fTime;
        int fEvent;
        double fX, fY, fZ, fEnergy;

        chain->SetBranchAddress("fTime", &fTime);
        chain->SetBranchAddress("fEvent", &fEvent);
        chain->SetBranchAddress("fX", &fX);
        chain->SetBranchAddress("fY", &fY);
        chain->SetBranchAddress("fZ", &fZ);
        chain->SetBranchAddress("fEnergy", &fEnergy);

        Long64_t nEntries = chain->GetEntries();

        for (Long64_t start = 0; start < nEntries; start += chunkSize)
        {
            std::vector<EventData> entries;
            for (Long64_t i = start; i < std::min(start + chunkSize, nEntries);
                 ++i)
            {
                chain->GetEntry(i);
                entries.push_back({fTime, fEvent, fX, fY, fZ, fEnergy});
            }

            // Sort entries first by event number, then by time within each
            // event number
            std::sort(entries.begin(), entries.end());

            std::string chunk_filename =
                "sorted_chunk_" + std::to_string(run) + "_" +
                std::to_string(start / chunkSize) + ".root";
            sort.write_sorted_chunk(entries, chunk_filename);
            chunk_files.push_back(chunk_filename);
        }

        // Merge sorted chunks into the final output file
        std::string output_filename =
            "sorted_output_run" + std::to_string(run) + ".root";
        sort.merge_sorted_chunks(chunk_files, output_filename);

        // Clean up temporary chunk files
        for (const auto &file : chunk_files)
        {
            fs::remove(file);
        }

        check_root_file(run);
    }

    return 0;
}
