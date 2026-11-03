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

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "graphalign/GraphAlignment.hh"
#include "graphalign/GraphAlignmentOperations.hh"
#include "graphcore/Graph.hh"
#include "graphcore/GraphBuilders.hh"

#include "core/LocusStats.hh"
#include "genotyping/AlignMatrix.hh"
#include "genotyping/AlignMatrixFiltering.hh"
#include "genotyping/FragLogliks.hh"
#include "genotyping/StrGenotyper.hh"
#include "io/GraphBlueprint.hh"
#include "io/RegionGraph.hh"

using graphtools::decodeGraphAlignment;
using graphtools::Graph;
using graphtools::GraphAlignment;
using namespace ehunter;
using namespace ehunter::strgt;
using testing::DoubleNear;

TEST(AddingSingleReads, MixedFragments_TrackedByFragment)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGA(C)*ATGTCG"));
    AlignMatrix alignMatrix(1);

    GraphAlignment spanningRead = decodeGraphAlignment(3, "0[3M]1[1M]1[1M]2[4M]", &graph);
    GraphAlignment flankRead = decodeGraphAlignment(0, "0[6M]", &graph);
    GraphAlignment zeroMotifSpanningRead = decodeGraphAlignment(3, "0[3M]2[4M]", &graph);

    alignMatrix.addSingle(spanningRead);
    ASSERT_EQ(1, alignMatrix.numReads());
    ASSERT_EQ(1, alignMatrix.numFrags());
    ASSERT_EQ(1, alignMatrix.fragSize(0));
    ASSERT_EQ(0, alignMatrix.fragFirstRow(0));

    // Reads that neither overlap nor span the repeat are ignored
    alignMatrix.addSingle(flankRead);
    ASSERT_EQ(1, alignMatrix.numReads());
    ASSERT_EQ(1, alignMatrix.numFrags());

    // Reads spanning a zero-length repeat are informative
    alignMatrix.addSingle(zeroMotifSpanningRead);
    ASSERT_EQ(2, alignMatrix.numReads());
    ASSERT_EQ(2, alignMatrix.numFrags());
    ASSERT_EQ(StrAlign('S', 0, 35, 0), alignMatrix.getBestAlign(1));

    alignMatrix.add(spanningRead, flankRead);
    ASSERT_EQ(4, alignMatrix.numReads());
    ASSERT_EQ(3, alignMatrix.numFrags());
    ASSERT_EQ(2, alignMatrix.fragSize(2));
    ASSERT_EQ(2, alignMatrix.fragFirstRow(2));

    alignMatrix.removeFrag(1);
    ASSERT_EQ(3, alignMatrix.numReads());
    ASSERT_EQ(2, alignMatrix.numFrags());
    ASSERT_EQ(1, alignMatrix.fragFirstRow(1));
    ASSERT_EQ(2, alignMatrix.fragSize(1));

    alignMatrix.remove(2);
    ASSERT_EQ(2, alignMatrix.numReads());
    ASSERT_EQ(2, alignMatrix.numFrags());
    ASSERT_EQ(1, alignMatrix.fragSize(1));
}

TEST(AddingSingleReads, PairsOnly_LayoutUnchanged)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGA(C)*ATGTCG"));
    AlignMatrix alignMatrix(1);
    GraphAlignment mate = decodeGraphAlignment(0, "0[6M]", &graph);
    alignMatrix.add(decodeGraphAlignment(3, "0[3M]1[1M]1[1M]2[4M]", &graph), mate);
    alignMatrix.add(decodeGraphAlignment(3, "0[3M]1[1M]1[1M]", &graph), mate);

    ASSERT_EQ(4, alignMatrix.numReads());
    ASSERT_EQ(2, alignMatrix.numFrags());
    ASSERT_EQ(0, alignMatrix.fragFirstRow(0));
    ASSERT_EQ(2, alignMatrix.fragFirstRow(1));

    FragLogliks logliks(1, 8, 20, &alignMatrix);
    ASSERT_EQ(2, logliks.numFrags());
    EXPECT_THAT(logliks.getLoglik(0, 2), DoubleNear(-8.07, 0.1));
}

TEST(SingleReadLogliks, SpanningRead_Computed)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGA(C)*ATGTCG"));
    AlignMatrix alignMatrix(1);
    alignMatrix.addSingle(decodeGraphAlignment(3, "0[3M]1[1M]1[1M]2[4M]", &graph));

    const int motifLen = 1;
    const int readLen = 8;
    const int fragLen = 8;
    FragLogliks logliks(motifLen, readLen, fragLen, &alignMatrix);
    ASSERT_EQ(1, logliks.numFrags());

    // Allele of size 2: spanning alignment with score 45 and 2 * 1 + 8 + 1 = 11 possible read origins
    EXPECT_THAT(logliks.getLoglik(0, 2), DoubleNear(-1.68, 0.05));
    // Allele of size 0: flanking alignment with score 20 and 9 possible read origins
    EXPECT_THAT(logliks.getLoglik(0, 0), DoubleNear(-8.04, 0.05));
}

TEST(GenotypingSingleReads, SpanningReads_Genotyped)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGA(CAG)*ATGTCG"));
    AlignMatrix alignMatrix(1);
    for (int index = 0; index != 10; ++index)
    {
        alignMatrix.addSingle(decodeGraphAlignment(0, "0[6M]1[3M]1[3M]2[6M]", &graph));
        alignMatrix.addSingle(decodeGraphAlignment(0, "0[6M]1[3M]1[3M]1[3M]1[3M]1[3M]2[6M]", &graph));
    }

    const int motifLen = 3;
    const int readLen = 27;
    const int fragLen = 27;
    RepeatGenotype gt = genotype(AlleleCount::kTwo, motifLen, readLen, fragLen, alignMatrix);
    EXPECT_EQ(2, gt.shortAlleleSizeInUnits());
    EXPECT_EQ(5, gt.longAlleleSizeInUnits());

    RepeatGenotype haploidGt = genotype(AlleleCount::kOne, motifLen, readLen, fragLen, alignMatrix);
    EXPECT_TRUE(haploidGt.shortAlleleSizeInUnits() == 2 || haploidGt.shortAlleleSizeInUnits() == 5);
}

TEST(AddingInrepeatReads, LongInrepeatRead_ReadsAdded)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGA(C)*ATGTCG"));
    AlignMatrix alignMatrix(1);
    alignMatrix.addSingle(decodeGraphAlignment(0, "1[1M]1[1M]1[1M]1[1M]", &graph));
    ASSERT_EQ(StrAlign('I', 4, 20, 0), alignMatrix.getBestAlign(0));

    addIrrReadsIfPossibleExpansion(4, alignMatrix, 3, 1);
    ASSERT_EQ(4, alignMatrix.numReads());
    ASSERT_EQ(4, alignMatrix.numFrags());
    ASSERT_EQ(3, alignMatrix.fragFirstRow(3));
    ASSERT_EQ(StrAlign('I', 4, 20, 0), alignMatrix.getBestAlign(3));

    addIrrPairsIfPossibleExpansion(4, alignMatrix, 2);
    ASSERT_EQ(8, alignMatrix.numReads());
    ASSERT_EQ(6, alignMatrix.numFrags());
    ASSERT_EQ(2, alignMatrix.fragSize(5));
}

TEST(AddingInrepeatReads, NoInrepeatReadButLongFlankingRead_PerfectReadsAdded)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGA(C)*ATGTCG"));
    AlignMatrix alignMatrix(1);
    alignMatrix.addSingle(decodeGraphAlignment(3, "0[3M]1[1M]1[1M]1[1M]", &graph));
    ASSERT_EQ(StrAlign('F', 3, 30, 0), alignMatrix.getBestAlign(0));

    // Flanking read with 3 of at most 4 motifs per read: the repeat may extend beyond the read length
    addIrrReadsIfPossibleExpansion(4, alignMatrix, 2, 1);
    ASSERT_EQ(3, alignMatrix.numReads());
    ASSERT_EQ(3, alignMatrix.numFrags());
    ASSERT_EQ(StrAlign('I', 4, 20, 0), alignMatrix.getBestAlign(2));
    ASSERT_EQ(StrAlign('I', 2, 10, 0), alignMatrix.getAlign(2, 2));
    ASSERT_EQ(StrAlign('I', 4, 20, 0), alignMatrix.getAlign(2, 10));

    // Pairs keep the original behaviour: nothing is added without an anchored in-repeat read
    AlignMatrix pairMatrix(1);
    pairMatrix.add(decodeGraphAlignment(3, "0[3M]1[1M]1[1M]1[1M]", &graph), decodeGraphAlignment(0, "0[6M]", &graph));
    addIrrPairsIfPossibleExpansion(4, pairMatrix, 2);
    ASSERT_EQ(2, pairMatrix.numReads());

    // Flanking reads that do not reach deep into the repeat are not enough
    AlignMatrix shortFlankMatrix(1);
    shortFlankMatrix.addSingle(decodeGraphAlignment(3, "0[3M]1[1M]", &graph));
    addIrrReadsIfPossibleExpansion(4, shortFlankMatrix, 2, 1);
    ASSERT_EQ(1, shortFlankMatrix.numReads());
}

TEST(FilteringSingleReads, NoIndels_NothingRemoved)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGA(C)*ATGTCG"));
    AlignMatrix alignMatrix(1);
    alignMatrix.addSingle(decodeGraphAlignment(3, "0[3M]1[1M]1[1M]2[4M]", &graph));
    alignMatrix.add(decodeGraphAlignment(3, "0[3M]1[1M]1[1M]", &graph), decodeGraphAlignment(0, "0[6M]", &graph));

    filter(alignMatrix);
    ASSERT_EQ(3, alignMatrix.numReads());
    ASSERT_EQ(2, alignMatrix.numFrags());
}

TEST(LocusStatsCalculatorForAnyReads, FlankCoverage_DepthEstimated)
{
    graphtools::Graph graph = graphtools::makeStrGraph("TAATG", "CCG", "CCTTATTA");

    LocusStatsCalculator statsCalculator(ChromType::kAutosome, graph, ReadMode::kAny);

    GraphAlignment alignmentStartingOnLeftFlank = decodeGraphAlignment(3, "0[2M]1[2M]", &graph);
    GraphAlignment alignmentStartingInsideRepeat = decodeGraphAlignment(0, "1[3M]", &graph);

    // Flanks have 13 bases in total; each read aligns 2 bases to the left flank
    for (int index = 0; index != 13; ++index)
    {
        statsCalculator.inspectRead(alignmentStartingOnLeftFlank);
        statsCalculator.inspectRead(alignmentStartingInsideRepeat);
    }

    // Without read pairs the fragment length defaults to the read length
    ASSERT_EQ(LocusStats(AlleleCount::kTwo, 4, 4, 2.0), statsCalculator.estimate(Sex::kFemale));
}

TEST(LocusStatsCalculatorForPairedReads, FlankCoverage_OriginalEstimateKept)
{
    graphtools::Graph graph = graphtools::makeStrGraph("TAATG", "CCG", "CCTTATTA");

    LocusStatsCalculator statsCalculator(ChromType::kAutosome, graph);

    GraphAlignment alignmentStartingOnLeftFlank = decodeGraphAlignment(3, "0[2M]1[2M]", &graph);
    for (int index = 0; index != 9; ++index)
    {
        statsCalculator.inspectRead(alignmentStartingOnLeftFlank);
    }

    // 9 reads of length 4 over 13 - 4 = 9 start positions
    ASSERT_EQ(LocusStats(AlleleCount::kTwo, 4, 0, 4.0), statsCalculator.estimate(Sex::kFemale));
}
