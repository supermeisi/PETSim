#include "sort.hh"

TChain *Sort::merge_root_files(const std::string &filePattern, int run,
                         int numThreads)
{
    TChain *chain = new TChain("Hits"); // The tree name is "Hits"
    for (int j = 0; j < numThreads; ++j)
    {
        std::string filename = filePattern;
        filename.replace(filename.find("%i"), 2, std::to_string(run));
        filename.replace(filename.find("%j"), 2, std::to_string(j));
        if (fs::exists(filename))
        {
            chain->Add(filename.c_str());
        }
    }
    return chain;
}

std::vector<EventData> Sort::collect_entries(TChain *chain)
{
    std::vector<EventData> entries;

    // Set branch addresses
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
    for (Long64_t i = 0; i < nEntries; ++i)
    {
        chain->GetEntry(i);
        entries.push_back({fTime, fEvent, fX, fY, fZ, fEnergy});
    }

    return entries;
}

void Sort::write_sorted_chunk(const std::vector<EventData> &entries,
                        const std::string &filename)
{
    TFile *file = new TFile(filename.c_str(), "RECREATE");
    TTree *tree = new TTree("sorted_tree", "Tree with sorted entries");

    double fTime;
    int fEvent;
    double fX, fY, fZ, fEnergy;

    tree->Branch("fTime", &fTime);
    tree->Branch("fEvent", &fEvent);
    tree->Branch("fX", &fX);
    tree->Branch("fY", &fY);
    tree->Branch("fZ", &fZ);
    tree->Branch("fEnergy", &fEnergy);

    for (const auto &entry : entries)
    {
        fTime = entry.fTime;
        fEvent = entry.fEvent;
        fX = entry.fX;
        fY = entry.fY;
        fZ = entry.fZ;
        fEnergy = entry.fEnergy;
        tree->Fill();
    }

    tree->Write();
    file->Close();
}

void Sort::merge_sorted_chunks(const std::vector<std::string> &chunk_files,
                         const std::string &output_filename)
{
    TChain chain("sorted_tree");
    for (const auto &file : chunk_files)
    {
        chain.Add(file.c_str());
    }

    std::vector<EventData> buffer;
    const int bufferSize = 10000000; // Adjust this based on available memory
    std::priority_queue<EventData, std::vector<EventData>, CompareEventData> pq;

    double fTime;
    int fEvent;
    double fX, fY, fZ, fEnergy;

    chain.SetBranchAddress("fTime", &fTime);
    chain.SetBranchAddress("fEvent", &fEvent);
    chain.SetBranchAddress("fX", &fX);
    chain.SetBranchAddress("fY", &fY);
    chain.SetBranchAddress("fZ", &fZ);
    chain.SetBranchAddress("fEnergy", &fEnergy);

    Long64_t nEntries = chain.GetEntries();

    for (Long64_t i = 0; i < nEntries; ++i)
    {
        chain.GetEntry(i);
        pq.push({fTime, fEvent, fX, fY, fZ, fEnergy});
        if (pq.size() >= bufferSize)
        {
            buffer.push_back(pq.top());
            pq.pop();
        }
    }

    while (!pq.empty())
    {
        buffer.push_back(pq.top());
        pq.pop();
    }

    write_sorted_chunk(buffer, output_filename);
}
