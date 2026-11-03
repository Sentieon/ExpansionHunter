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

#include "core/HtsHelpers.hh"

#include <string>
#include <vector>

#include "gtest/gtest.h"

using namespace ehunter;
using namespace ehunter::htshelpers;
using std::string;
using std::vector;

namespace
{
vector<uint32_t> cigar(const vector<std::pair<int, int>>& opsAndLens)
{
    vector<uint32_t> encoding;
    for (const auto& opAndLen : opsAndLens)
    {
        encoding.push_back(bam_cigar_gen(opAndLen.second, opAndLen.first));
    }
    return encoding;
}
}

TEST(ComputingQueryIntervalInRegion, ReadInsideRegion_WholeReadKept)
{
    const auto interval = computeQueryIntervalInRegion(100, cigar({ { BAM_CMATCH, 100 } }), 100, 50, 300);
    EXPECT_EQ(0, interval.start);
    EXPECT_EQ(100, interval.end);
}

TEST(ComputingQueryIntervalInRegion, ReadOverlappingRegionBoundaries_Clipped)
{
    auto interval = computeQueryIntervalInRegion(40, cigar({ { BAM_CMATCH, 100 } }), 100, 50, 300);
    EXPECT_EQ(10, interval.start);
    EXPECT_EQ(100, interval.end);

    interval = computeQueryIntervalInRegion(250, cigar({ { BAM_CMATCH, 100 } }), 100, 50, 300);
    EXPECT_EQ(0, interval.start);
    EXPECT_EQ(50, interval.end);

    interval = computeQueryIntervalInRegion(0, cigar({ { BAM_CMATCH, 1000 } }), 1000, 100, 400);
    EXPECT_EQ(100, interval.start);
    EXPECT_EQ(400, interval.end);
}

TEST(ComputingQueryIntervalInRegion, SoftclipsAtReadEndsInsideRegion_Kept)
{
    // Leading softclip of a read starting inside the region is kept
    auto interval = computeQueryIntervalInRegion(100, cigar({ { BAM_CSOFT_CLIP, 20 }, { BAM_CMATCH, 80 } }), 100, 50, 300);
    EXPECT_EQ(0, interval.start);
    EXPECT_EQ(100, interval.end);

    // Trailing softclip of a read ending inside the region is kept
    interval = computeQueryIntervalInRegion(200, cigar({ { BAM_CMATCH, 80 }, { BAM_CSOFT_CLIP, 20 } }), 100, 50, 300);
    EXPECT_EQ(0, interval.start);
    EXPECT_EQ(100, interval.end);

    // Trailing softclip of a read whose alignment leaves the region is dropped
    interval = computeQueryIntervalInRegion(250, cigar({ { BAM_CMATCH, 80 }, { BAM_CSOFT_CLIP, 20 } }), 100, 50, 300);
    EXPECT_EQ(0, interval.start);
    EXPECT_EQ(50, interval.end);

    // Leading softclip of a read starting before the region is dropped
    interval = computeQueryIntervalInRegion(40, cigar({ { BAM_CSOFT_CLIP, 20 }, { BAM_CMATCH, 80 } }), 100, 50, 300);
    EXPECT_EQ(30, interval.start);
    EXPECT_EQ(100, interval.end);
}

TEST(ComputingQueryIntervalInRegion, IndelsNearBoundaries_Handled)
{
    // Insertion inside the region is kept
    auto interval = computeQueryIntervalInRegion(
        40, cigar({ { BAM_CMATCH, 50 }, { BAM_CINS, 10 }, { BAM_CMATCH, 50 } }), 110, 50, 300);
    EXPECT_EQ(10, interval.start);
    EXPECT_EQ(110, interval.end);

    // Deletion spanning the region start: the read resumes at the first base after the deletion
    interval = computeQueryIntervalInRegion(
        30, cigar({ { BAM_CMATCH, 40 }, { BAM_CDEL, 20 }, { BAM_CMATCH, 40 } }), 80, 80, 300);
    EXPECT_EQ(40, interval.start);
    EXPECT_EQ(80, interval.end);

    // Deletion spanning the region end
    interval = computeQueryIntervalInRegion(
        30, cigar({ { BAM_CMATCH, 40 }, { BAM_CDEL, 20 }, { BAM_CMATCH, 40 } }), 80, 0, 80);
    EXPECT_EQ(0, interval.start);
    EXPECT_EQ(40, interval.end);
}

TEST(ComputingQueryIntervalInRegion, ReadOutsideRegion_EmptyInterval)
{
    auto interval = computeQueryIntervalInRegion(0, cigar({ { BAM_CMATCH, 100 } }), 100, 200, 300);
    EXPECT_EQ(0, interval.length());

    interval = computeQueryIntervalInRegion(400, cigar({ { BAM_CMATCH, 100 } }), 100, 200, 300);
    EXPECT_EQ(0, interval.length());
}

TEST(ClippingReadsToRegion, TypicalReads_Clipped)
{
    const string bases = "ACGTACGTAC";
    Read read(ReadId("frag1", MateNumber::kUnpaired), bases, false, 60);
    GenomicRegion region(1, 100, 105);

    // Alignment overlapping the region is clipped
    auto clipped = clipReadToRegion(read, 1, 98, cigar({ { BAM_CMATCH, 10 } }), region);
    ASSERT_TRUE(clipped);
    EXPECT_EQ("GTACG", clipped->sequence());
    EXPECT_EQ(read.readId(), clipped->readId());
    EXPECT_EQ(60, clipped->mapq());
    EXPECT_TRUE(clipped->isUnpaired());

    // Alignment inside the region is unchanged
    clipped = clipReadToRegion(read, 1, 100, cigar({ { BAM_CMATCH, 10 } }), GenomicRegion(1, 90, 200));
    ASSERT_TRUE(clipped);
    EXPECT_EQ(bases, clipped->sequence());

    // Alignments on another contig, outside of the region, or without a cigar are unchanged
    clipped = clipReadToRegion(read, 2, 98, cigar({ { BAM_CMATCH, 10 } }), region);
    ASSERT_TRUE(clipped);
    EXPECT_EQ(bases, clipped->sequence());

    clipped = clipReadToRegion(read, 1, 200, cigar({ { BAM_CMATCH, 10 } }), region);
    ASSERT_TRUE(clipped);
    EXPECT_EQ(bases, clipped->sequence());

    clipped = clipReadToRegion(read, 1, 98, {}, region);
    ASSERT_TRUE(clipped);
    EXPECT_EQ(bases, clipped->sequence());

    // A read whose overlapping piece is empty is dropped
    clipped = clipReadToRegion(
        read, 1, 95, cigar({ { BAM_CMATCH, 5 }, { BAM_CDEL, 5 }, { BAM_CMATCH, 5 } }), GenomicRegion(1, 100, 105));
    EXPECT_FALSE(clipped);
}

TEST(ClippingReadsToRegion, LongPieces_LengthCapped)
{
    // 30 bp read aligned 100-130 with a 20 bp leading softclip: 50 bp in total
    const string bases(50, 'A');
    Read read(ReadId("frag1", MateNumber::kUnpaired), bases, false, 60);
    const auto readCigar = cigar({ { BAM_CSOFT_CLIP, 20 }, { BAM_CMATCH, 30 } });
    GenomicRegion region(1, 90, 200);

    // Without a cap the softclip of a read starting inside the region is kept
    auto clipped = clipReadToRegion(read, 1, 100, readCigar, region);
    ASSERT_TRUE(clipped);
    EXPECT_EQ(50u, clipped->sequence().length());

    // The softclip is shortened from its outer end first
    Read markedRead(ReadId("frag2", MateNumber::kUnpaired), string(20, 'C') + string(30, 'G'), false, 60);
    clipped = clipReadToRegion(markedRead, 1, 100, readCigar, region, 40);
    ASSERT_TRUE(clipped);
    EXPECT_EQ(string(10, 'C') + string(30, 'G'), clipped->sequence());

    // Then the aligned part is shortened from its end
    clipped = clipReadToRegion(markedRead, 1, 100, readCigar, region, 25);
    ASSERT_TRUE(clipped);
    EXPECT_EQ(string(25, 'G'), clipped->sequence());

    // Reads outside of the region and unaligned reads are capped too
    clipped = clipReadToRegion(markedRead, 1, 1000, readCigar, region, 25);
    ASSERT_TRUE(clipped);
    EXPECT_EQ(string(25, 'G'), clipped->sequence());

    clipped = clipReadToRegion(markedRead, -1, -1, {}, region, 25);
    ASSERT_TRUE(clipped);
    EXPECT_EQ(25u, clipped->sequence().length());

    // A cap of zero means no cap
    clipped = clipReadToRegion(markedRead, 1, 100, readCigar, region, 0);
    ASSERT_TRUE(clipped);
    EXPECT_EQ(50u, clipped->sequence().length());
}
