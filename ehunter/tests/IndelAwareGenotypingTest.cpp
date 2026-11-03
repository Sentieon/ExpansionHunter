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

#include "genotyping/AlignMatrix.hh"

#include <cmath>

#include "gmock/gmock.h"

#include "graphalign/GraphAlignment.hh"
#include "graphalign/GraphAlignmentOperations.hh"

#include "core/LogSum.hh"
#include "genotyping/FragLogliks.hh"
#include "io/GraphBlueprint.hh"
#include "io/RegionGraph.hh"

using graphtools::Graph;
using graphtools::GraphAlignment;
using testing::DoubleNear;

using namespace ehunter;
using namespace strgt;

// A read whose repeat alignment is free of indels determines the number of motifs of its allele, so smoothing must
// leave it alone
TEST(SmoothingAlignsAroundIndels, ReadWithoutIndels_ScoresUnchanged)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGA(TC)*ATGTCG"));
    GraphAlignment read = decodeGraphAlignment(0, "0[6M]1[2M]1[2M]1[2M]2[6M]", &graph);

    AlignMatrix plainMatrix(1);
    plainMatrix.addSingle(read);
    AlignMatrix smoothMatrix(1, true);
    smoothMatrix.addSingle(read);

    for (int alleleSize = 0; alleleSize != 6; ++alleleSize)
    {
        EXPECT_EQ(plainMatrix.getAlign(0, alleleSize), smoothMatrix.getAlign(0, alleleSize));
    }
}

// A single-base sequencing indel inside the repeat makes the graph aligner traverse one motif more (or fewer) than the
// allele contains; such a read must support both sizes equally instead of only the size that the aligner picked
TEST(SmoothingAlignsAroundIndels, ReadWithIndelInRepeat_ScoresOfNeighboringSizesFlattened)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGA(TC)*ATGTCG"));
    GraphAlignment read = decodeGraphAlignment(0, "0[6M]1[2M]1[2M]1[1M1D]2[6M]", &graph);

    AlignMatrix plainMatrix(1);
    plainMatrix.addSingle(read);
    ASSERT_EQ(1, plainMatrix.getAlign(0, 3).numIndels());
    ASSERT_EQ(StrAlign::Type::kSpanning, plainMatrix.getAlign(0, 3).type());
    // Without smoothing the two sizes that the read cannot distinguish get very different scores
    ASSERT_LT(plainMatrix.getAlign(0, 2).score(), plainMatrix.getAlign(0, 3).score() - 20);

    AlignMatrix smoothMatrix(1, true);
    smoothMatrix.addSingle(read);
    const int topScore = plainMatrix.getAlign(0, 3).score();
    EXPECT_EQ(topScore, smoothMatrix.getAlign(0, 2).score());
    EXPECT_EQ(topScore, smoothMatrix.getAlign(0, 3).score());
    // Sizes further away are not affected
    EXPECT_EQ(plainMatrix.getAlign(0, 1).score(), smoothMatrix.getAlign(0, 1).score());
    // The score used for every allele larger than the matrix row is left alone, so that a read cannot support
    // arbitrarily large expansions
    EXPECT_EQ(plainMatrix.getAlign(0, 20).score(), smoothMatrix.getAlign(0, 20).score());
}

// Without a stutter rate the fragment likelihoods are the ones of the original model
TEST(StutterModel, StutterRateOfZero_LogliksUnchanged)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGA(TC)*ATGTCG"));
    AlignMatrix alignMatrix(1);
    alignMatrix.addSingle(decodeGraphAlignment(0, "0[6M]1[2M]1[2M]1[2M]2[6M]", &graph));

    FragLogliks withoutStutter(2, 18, 18, &alignMatrix);
    FragLogliks withZeroStutter(2, 18, 18, &alignMatrix, 0.0);

    for (int alleleSize = 0; alleleSize != 6; ++alleleSize)
    {
        EXPECT_THAT(
            withZeroStutter.getLoglik(0, alleleSize), DoubleNear(withoutStutter.getLoglik(0, alleleSize), 1e-9));
    }
}

// With a stutter rate the likelihood of a fragment is averaged over the neighboring allele sizes
TEST(StutterModel, NonzeroStutterRate_LogliksAveragedOverNeighboringSizes)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGA(TC)*ATGTCG"));
    AlignMatrix alignMatrix(1);
    alignMatrix.addSingle(decodeGraphAlignment(0, "0[6M]1[2M]1[2M]1[2M]2[6M]", &graph));

    const double stutterRate = 0.05;
    FragLogliks rawLogliks(2, 18, 18, &alignMatrix);
    FragLogliks stutterLogliks(2, 18, 18, &alignMatrix, stutterRate);

    for (int alleleSize = 1; alleleSize != 6; ++alleleSize)
    {
        double expected = std::log(1 - 2 * stutterRate) + rawLogliks.getLoglik(0, alleleSize);
        expected = getLogSum(expected, std::log(stutterRate) + rawLogliks.getLoglik(0, alleleSize - 1));
        expected = getLogSum(expected, std::log(stutterRate) + rawLogliks.getLoglik(0, alleleSize + 1));
        EXPECT_THAT(stutterLogliks.getLoglik(0, alleleSize), DoubleNear(expected, 1e-9));
    }

    // A read matching three motifs still supports three motifs best
    double bestLoglik = stutterLogliks.getLoglik(0, 0);
    int bestSize = 0;
    for (int alleleSize = 1; alleleSize != 8; ++alleleSize)
    {
        if (stutterLogliks.getLoglik(0, alleleSize) > bestLoglik)
        {
            bestLoglik = stutterLogliks.getLoglik(0, alleleSize);
            bestSize = alleleSize;
        }
    }
    EXPECT_EQ(3, bestSize);
}
