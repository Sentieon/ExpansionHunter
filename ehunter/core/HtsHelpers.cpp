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

#include "core/HtsHelpers.hh"

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <memory>
#include <string>
#include <utility>

#include "spdlog/spdlog.h"

using std::pair;
using std::string;
using std::vector;

namespace spd = spdlog;

namespace ehunter
{
namespace htshelpers
{

string decodeQuals(bam1_t* htsAlignPtr)
{
    string quals;
    uint8_t* htsQualsPtr = bam_get_qual(htsAlignPtr);
    const int readLength = htsAlignPtr->core.l_qseq;
    quals.resize(readLength);

    for (int index = 0; index < readLength; ++index)
    {
        quals[index] = static_cast<char>(33 + htsQualsPtr[index]);
    }

    return quals;
}

string decodeBases(bam1_t* htsAlignPtr)
{
    string bases;
    uint8_t* htsSeqPtr = bam_get_seq(htsAlignPtr);
    const int32_t readLength = htsAlignPtr->core.l_qseq;
    bases.resize(readLength);

    for (int32_t index = 0; index < readLength; ++index)
    {
        bases[index] = seq_nt16_str[bam_seqi(htsSeqPtr, index)];
    }

    return bases;
}

LinearAlignmentStats decodeAlignmentStats(bam1_t* htsAlignPtr)
{
    LinearAlignmentStats alignmentStats;
    alignmentStats.chromId = htsAlignPtr->core.tid;
    alignmentStats.pos = htsAlignPtr->core.pos;
    alignmentStats.mapq = htsAlignPtr->core.qual;
    alignmentStats.mateChromId = htsAlignPtr->core.mtid;
    alignmentStats.matePos = htsAlignPtr->core.mpos;
    alignmentStats.refEnd = bam_endpos(htsAlignPtr);

    uint32_t samFlag = htsAlignPtr->core.flag;
    alignmentStats.isPaired = samFlag & BAM_FPAIRED;
    alignmentStats.isMapped = !(samFlag & BAM_FUNMAP);
    alignmentStats.isMateMapped = !(samFlag & BAM_FMUNMAP);

    return alignmentStats;
}

bool isPrimaryAlignment(bam1_t* htsAlignPtr)
{
    return !((htsAlignPtr->core.flag & BAM_FSECONDARY) || (htsAlignPtr->core.flag & BAM_FSUPPLEMENTARY));
}

/// Lower-case version of htslib seq_nt16_str table;
class Cache_seq_nt16_str_lc
{
public:
    Cache_seq_nt16_str_lc()
    {
        for (unsigned i(0); i < 16; ++i)
        {
            data[i] = std::tolower(seq_nt16_str[i]);
        }
    }

    char data[16];
};

static Cache_seq_nt16_str_lc seq_nt16_str_lc;

Read decodeRead(bam1_t* htsAlignPtr)
{
    const uint32_t samFlag = htsAlignPtr->core.flag;
    const bool isPaired = samFlag & BAM_FPAIRED;
    const bool isFirstMate = samFlag & BAM_FREAD1;
    const bool isReversed = samFlag & BAM_FREVERSE;

    const char* qname(bam_get_qname(htsAlignPtr));
    MateNumber mateNumber = !isPaired ? MateNumber::kUnpaired
                                      : (isFirstMate ? MateNumber::kFirstMate : MateNumber::kSecondMate);
    ReadId readId(qname, mateNumber);

    string bases;

    {
        // Decode bases and convert low-quality bases to lowercase:
        static const uint8_t lowBaseQualityCutoff(20);

        const uint8_t* htsSeqPtr = bam_get_seq(htsAlignPtr);
        const uint8_t* htsQualPtr = bam_get_qual(htsAlignPtr);
        const int32_t readLength = htsAlignPtr->core.l_qseq;
        bases.resize(readLength);

        for (int32_t index = 0; index < readLength; ++index)
        {
            const char* converter((htsQualPtr[index] <= lowBaseQualityCutoff) ? seq_nt16_str_lc.data : seq_nt16_str);
            bases[index] = converter[bam_seqi(htsSeqPtr, index)];
        }
    }

    return { std::move(readId), std::move(bases), isReversed, static_cast<int>(htsAlignPtr->core.qual) };
}

QueryInterval computeQueryIntervalInRegion(
    int64_t alignStart, const std::vector<uint32_t>& cigar, int32_t queryLen, int64_t regionStart, int64_t regionEnd)
{
    int64_t refPos = alignStart;
    int32_t queryPos = 0;
    int32_t queryStart = -1;
    int32_t queryEnd = -1;

    for (uint32_t cigarOp : cigar)
    {
        const int opType = bam_cigar_op(cigarOp);
        const int64_t opLen = bam_cigar_oplen(cigarOp);
        const bool consumesQuery = bam_cigar_type(opType) & 1;
        const bool consumesReference = bam_cigar_type(opType) & 2;

        if (queryStart < 0 && refPos >= regionStart)
        {
            queryStart = queryPos;
        }
        if (refPos >= regionEnd)
        {
            queryEnd = queryPos;
            break;
        }

        if (consumesReference)
        {
            if (queryStart < 0 && refPos < regionStart && refPos + opLen > regionStart)
            {
                queryStart = consumesQuery ? queryPos + static_cast<int32_t>(regionStart - refPos) : queryPos;
            }
            if (refPos < regionEnd && refPos + opLen > regionEnd)
            {
                queryEnd = consumesQuery ? queryPos + static_cast<int32_t>(regionEnd - refPos) : queryPos;
                break;
            }
        }

        if (consumesQuery)
        {
            queryPos += static_cast<int32_t>(opLen);
        }
        if (consumesReference)
        {
            refPos += opLen;
        }
    }

    if (queryStart < 0)
    {
        // Either the alignment ends before the region or a trailing operation reached the region boundary exactly
        if (queryEnd < 0 && refPos >= regionStart && refPos < regionEnd)
        {
            queryStart = queryPos;
        }
        else
        {
            return QueryInterval();
        }
    }
    if (queryEnd < 0)
    {
        // Alignment ends inside the region: keep the rest of the read, including any trailing soft-clip
        queryEnd = queryLen;
    }

    QueryInterval interval;
    interval.start = std::min(queryStart, queryLen);
    interval.end = std::min(std::max(queryEnd, interval.start), queryLen);
    return interval;
}

std::vector<uint32_t> decodeCigar(bam1_t* htsAlignPtr)
{
    const uint32_t* cigarPtr = bam_get_cigar(htsAlignPtr);
    return std::vector<uint32_t>(cigarPtr, cigarPtr + htsAlignPtr->core.n_cigar);
}

namespace
{
// Number of soft-clipped bases at the start (end) of the read, ignoring hard clips
int32_t leadingSoftclipLength(const std::vector<uint32_t>& cigar)
{
    for (uint32_t cigarOp : cigar)
    {
        if (bam_cigar_op(cigarOp) == BAM_CHARD_CLIP)
        {
            continue;
        }
        return bam_cigar_op(cigarOp) == BAM_CSOFT_CLIP ? static_cast<int32_t>(bam_cigar_oplen(cigarOp)) : 0;
    }
    return 0;
}

int32_t trailingSoftclipLength(const std::vector<uint32_t>& cigar)
{
    for (auto it = cigar.rbegin(); it != cigar.rend(); ++it)
    {
        if (bam_cigar_op(*it) == BAM_CHARD_CLIP)
        {
            continue;
        }
        return bam_cigar_op(*it) == BAM_CSOFT_CLIP ? static_cast<int32_t>(bam_cigar_oplen(*it)) : 0;
    }
    return 0;
}
}

boost::optional<Read> clipReadToRegion(
    const Read& read, int32_t contigId, int64_t alignStart, const std::vector<uint32_t>& cigar,
    const GenomicRegion& region, int32_t maxLength)
{
    const int32_t queryLen = read.sequence().length();
    QueryInterval interval;
    interval.start = 0;
    interval.end = queryLen;

    // Soft-clipped read ends that are part of the interval
    int32_t leadingSoftclip = 0;
    int32_t trailingSoftclip = 0;

    const bool isAlignmentUsable = !cigar.empty() && contigId == region.contigIndex() && alignStart >= 0;
    if (isAlignmentUsable)
    {
        const int64_t alignEnd = alignStart + bam_cigar2rlen(cigar.size(), cigar.data());
        if (alignEnd > region.start() && alignStart < region.end())
        {
            interval = computeQueryIntervalInRegion(alignStart, cigar, queryLen, region.start(), region.end());
            if (interval.length() <= 0)
            {
                return boost::none;
            }
        }

        if (interval.start == 0)
        {
            leadingSoftclip = std::min(leadingSoftclipLength(cigar), interval.length());
        }
        if (interval.end == queryLen)
        {
            trailingSoftclip = std::min(trailingSoftclipLength(cigar), interval.length() - leadingSoftclip);
        }
    }

    // Enforce the length limit: shorten retained soft clips from their outer ends first, then remove bases from the
    // end of the read so that the start of the read (anchored on the left of the locus) is kept
    if (maxLength > 0 && interval.length() > maxLength)
    {
        int32_t excess = interval.length() - maxLength;

        const int32_t leadingTrim = std::min(excess, leadingSoftclip);
        interval.start += leadingTrim;
        excess -= leadingTrim;

        const int32_t trailingTrim = std::min(excess, trailingSoftclip);
        interval.end -= trailingTrim;
        excess -= trailingTrim;

        interval.end -= excess;
    }

    if (interval.start == 0 && interval.end == queryLen)
    {
        return read;
    }

    return Read(read.readId(), read.sequence().substr(interval.start, interval.length()), read.isReversed(), read.mapq());
}

boost::optional<Read> decodeReadClippedToRegion(bam1_t* htsAlignPtr, const GenomicRegion& region, int32_t maxLength)
{
    Read read = decodeRead(htsAlignPtr);
    if (htsAlignPtr->core.flag & BAM_FUNMAP)
    {
        return clipReadToRegion(read, -1, -1, {}, region, maxLength);
    }
    return clipReadToRegion(
        read, htsAlignPtr->core.tid, htsAlignPtr->core.pos, decodeCigar(htsAlignPtr), region, maxLength);
}

ReferenceContigInfo decodeContigInfo(bam_hdr_t* htsHeaderPtr)
{
    vector<pair<string, int64_t>> contigNamesAndSizes;
    const int32_t numContigs = htsHeaderPtr->n_targets;
    contigNamesAndSizes.reserve(numContigs);

    for (int32_t contigIndex = 0; contigIndex != numContigs; ++contigIndex)
    {
        const string contig = htsHeaderPtr->target_name[contigIndex];
        int64_t size = htsHeaderPtr->target_len[contigIndex];
        contigNamesAndSizes.push_back(std::make_pair(contig, size));
    }

    return ReferenceContigInfo(contigNamesAndSizes);
}

} // namespace htshelpers
}
