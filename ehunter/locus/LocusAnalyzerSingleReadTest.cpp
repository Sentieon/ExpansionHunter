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

#include "locus/LocusAnalyzer.hh"

#include "gtest/gtest.h"

#include "graphalign/GraphAlignmentOperations.hh"
#include "graphio/AlignmentWriter.hh"
#include "io/GraphBlueprint.hh"
#include "io/RegionGraph.hh"

using namespace ehunter;
using namespace locus;
using graphtools::AlignerType;
using std::vector;

namespace
{
LocusSpecification buildSingleReadStrSpec(const std::string& structure)
{
    auto graph = makeRegionGraph(decodeFeaturesFromRegex(structure));
    vector<GenomicRegion> referenceRegions = { GenomicRegion(1, 1, 2) };
    NodeToRegionAssociation dummyAssociation;
    GenotyperParameters params;
    LocusSpecification locusSpec("region", ChromType::kAutosome, referenceRegions, graph, dummyAssociation, params);
    VariantClassification classification(VariantType::kRepeat, VariantSubtype::kCommonRepeat);
    locusSpec.addVariantSpecification("repeat", classification, GenomicRegion(1, 1, 2), { 1 }, 1);
    return locusSpec;
}
}

TEST(AnalyzingSingleReads, UnpairedSpanningReads_Genotyped)
{
    auto locusSpec = buildSingleReadStrSpec("ATTCGA(C)*ATGTCG");
    HeuristicParameters heuristicParams(1000, 20, true, AlignerType::DAG_ALIGNER, 4, 1, 5, 4, 1);
    heuristicParams.setReadMode(ReadMode::kAny);
    auto writer = std::make_shared<graphtools::BlankAlignmentWriter>();
    graphtools::AlignerSelector selector(heuristicParams.alignerType());

    LocusAnalyzer locusAnalyzer(locusSpec, heuristicParams, writer);

    for (int index = 0; index != 9; ++index)
    {
        Read longAlleleRead(
            ReadId("long" + std::to_string(index), MateNumber::kUnpaired), "TTCGACCCATGTCG", false, 60);
        locusAnalyzer.processMates(longAlleleRead, nullptr, RegionType::kTarget, selector);

        Read shortAlleleRead(
            ReadId("short" + std::to_string(index), MateNumber::kUnpaired), "TTCGACATGTCG", false, 60);
        locusAnalyzer.processMates(shortAlleleRead, nullptr, RegionType::kTarget, selector);
    }

    LocusFindings locusFindings = locusAnalyzer.analyze(Sex::kFemale, boost::none);
    const auto& observed
        = *dynamic_cast<RepeatFindings*>(locusFindings.findingsForEachVariant["repeat"].get());

    ASSERT_EQ(CountTable({ { 1, 9 }, { 3, 9 } }), observed.countsOfSpanningReads());
    ASSERT_TRUE(observed.optionalGenotype());
    ASSERT_EQ(1, observed.optionalGenotype()->shortAlleleSizeInUnits());
    ASSERT_EQ(3, observed.optionalGenotype()->longAlleleSizeInUnits());
    ASSERT_EQ(AlleleCount::kTwo, locusFindings.stats.alleleCount());
    // Mean of the 14 bp and 12 bp reads; without read pairs the fragment length defaults to the read length
    ASSERT_EQ(13, locusFindings.stats.meanReadLength());
    ASSERT_EQ(13, locusFindings.stats.medianFragLength());
}

TEST(AnalyzingSingleReads, PairedModeIgnoresSingletonEvidence)
{
    auto locusSpec = buildSingleReadStrSpec("ATTCGA(C)*ATGTCG");
    HeuristicParameters heuristicParams(1000, 20, true, AlignerType::DAG_ALIGNER, 4, 1, 5, 4, 1);
    auto writer = std::make_shared<graphtools::BlankAlignmentWriter>();
    graphtools::AlignerSelector selector(heuristicParams.alignerType());

    LocusAnalyzer locusAnalyzer(locusSpec, heuristicParams, writer);

    for (int index = 0; index != 9; ++index)
    {
        Read read(ReadId("read" + std::to_string(index), MateNumber::kUnpaired), "TTCGACCCATGTCG", false, 60);
        locusAnalyzer.processMates(read, nullptr, RegionType::kTarget, selector);
    }

    LocusFindings locusFindings = locusAnalyzer.analyze(Sex::kFemale, boost::none);
    const auto& observed
        = *dynamic_cast<RepeatFindings*>(locusFindings.findingsForEachVariant["repeat"].get());

    ASSERT_EQ(CountTable(), observed.countsOfSpanningReads());
}
