//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
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

#pragma once

extern "C"
{
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/thread_pool.h"
}

#include <vector>

#include <boost/optional.hpp>

#include "core/GenomicRegion.hh"
#include "core/Read.hh"
#include "core/ReferenceContigInfo.hh"

namespace ehunter
{

namespace htshelpers
{

LinearAlignmentStats decodeAlignmentStats(bam1_t* htsAlignPtr);
bool isPrimaryAlignment(bam1_t* htsAlignPtr);
Read decodeRead(bam1_t* htsAlignPtr);

// Half-open interval of query (read) coordinates
struct QueryInterval
{
    int32_t start = 0;
    int32_t end = 0;
    int32_t length() const { return end - start; }
};

// Computes the piece of a linearly-aligned read that falls inside the reference region [regionStart, regionEnd).
// Soft-clipped read prefixes (suffixes) are included when the alignment starts (ends) inside the region so that
// unaligned repeat sequence at the read ends is retained. The interval is empty if the read does not reach the region.
QueryInterval computeQueryIntervalInRegion(
    int64_t alignStart, const std::vector<uint32_t>& cigar, int32_t queryLen, int64_t regionStart,
    int64_t regionEnd);

std::vector<uint32_t> decodeCigar(bam1_t* htsAlignPtr);

// Restricts a read to the part of its linear alignment that overlaps the region; reads that do not overlap the region
// (or are unaligned) are kept whole and reads whose overlapping piece is empty are dropped (none).
//
// If maxLength > 0 the result is at most maxLength bases long: retained soft clips are shortened from their outer
// ends first, then bases are removed from the end of the read. This bounds the cost of graph alignment, which grows
// with the square of the read length.
boost::optional<Read> clipReadToRegion(
    const Read& read, int32_t contigId, int64_t alignStart, const std::vector<uint32_t>& cigar,
    const GenomicRegion& region, int32_t maxLength = 0);

boost::optional<Read> decodeReadClippedToRegion(bam1_t* htsAlignPtr, const GenomicRegion& region, int32_t maxLength);
ReferenceContigInfo decodeContigInfo(bam_hdr_t* htsHeaderPtr);

} // namespace htshelpers

}
