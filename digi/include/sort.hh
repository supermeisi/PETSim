#ifndef SORT_HH
#define SORT_HH

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

namespace fs = std::filesystem;

struct EventData
{
    double fTime;
    int fEvent;
    double fX;
    double fY;
    double fZ;
    double fEnergy;

    bool operator<(const EventData &other) const
    {
        if (fEvent != other.fEvent)
            return fEvent < other.fEvent;
        return fTime < other.fTime;
    }
};

struct CompareEventData
{
    bool operator()(const EventData &a, const EventData &b)
    {
        if (a.fEvent != b.fEvent)
            return a.fEvent >
                   b.fEvent;      // Higher event number has lower priority
        return a.fTime > b.fTime; // Higher time has lower priority
    }
};

class Sort
{
public:
    Sort() {};
    ~Sort() {};

    TChain *merge_root_files(const std::string &filePattern, int run,
                             int numThreads);
    std::vector<EventData> collect_entries(TChain *chain);
    void write_sorted_chunk(const std::vector<EventData> &entries,
                            const std::string &filename);
    void merge_sorted_chunks(const std::vector<std::string> &chunk_files,
                             const std::string &output_filename);
};

#endif
