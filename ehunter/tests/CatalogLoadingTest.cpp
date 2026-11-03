//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
//

#include "io/CatalogLoading.hh"

#include <cstdio>
#include <fstream>
#include <map>
#include <string>

#include "gtest/gtest.h"

#include "io/LocusSpecDecoding.hh"

using namespace ehunter;
using std::string;

namespace
{

// Reference held in memory; chrMasked consists of N only, like a region masked in the reference
class InMemoryReference : public Reference
{
public:
    InMemoryReference()
        : contigInfo_({ { "chr1", 3000 }, { "chrMasked", 3000 } })
    {
        string sequence;
        unsigned state = 12345;
        for (int i = 0; i != 3000; ++i)
        {
            state = state * 1103515245 + 12345;
            sequence += "ACGT"[(state >> 16) % 4];
        }
        sequences_["chr1"] = sequence;
        sequences_["chrMasked"] = string(3000, 'N');
    }

    string getSequence(const string& contigName, int64_t start, int64_t end) const override
    {
        return sequences_.at(contigName).substr(start, end - start);
    }

    string getSequence(const GenomicRegion& region) const override
    {
        return getSequence(contigInfo_.getContigName(region.contigIndex()), region.start(), region.end());
    }

    const ReferenceContigInfo& contigInfo() const override { return contigInfo_; }

private:
    ReferenceContigInfo contigInfo_;
    std::map<string, string> sequences_;
};

string writeCatalog(const string& json)
{
    const string path = ::testing::TempDir() + "eh_catalog_loading_test.json";
    std::ofstream out(path);
    out << json;
    return path;
}

}

TEST(LoadingCatalog, LocusWithNsInFlanks_Skipped)
{
    InMemoryReference reference;
    HeuristicParameters heuristicParams(1000, 20, true, graphtools::AlignerType::DAG_ALIGNER);
    const string path = writeCatalog(R"([
        {"LocusId": "Masked", "LocusStructure": "(CAG)*", "ReferenceRegion": "chrMasked:1200-1230",
         "VariantType": "Repeat"},
        {"LocusId": "Analyzable", "LocusStructure": "(CAG)*", "ReferenceRegion": "chr1:1200-1230",
         "VariantType": "Repeat"}
    ])");

    RegionCatalog catalog = loadLocusCatalogFromDisk(path, heuristicParams, reference);
    std::remove(path.c_str());

    ASSERT_EQ(1u, catalog.size());
    EXPECT_EQ("Analyzable", catalog.front().locusId());
}

TEST(LoadingCatalog, InvalidLocus_Throws)
{
    InMemoryReference reference;
    HeuristicParameters heuristicParams(1000, 20, true, graphtools::AlignerType::DAG_ALIGNER);
    // Two repeats in the structure but only one reference region: a catalog error, not a masked locus
    const string path = writeCatalog(R"([
        {"LocusId": "Invalid", "LocusStructure": "(CAG)*(CCG)*", "ReferenceRegion": "chr1:1200-1230",
         "VariantType": "Repeat"}
    ])");

    EXPECT_ANY_THROW(loadLocusCatalogFromDisk(path, heuristicParams, reference));
    std::remove(path.c_str());
}
