//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
// Concept: Michael Eberle <meberle@illumina.com>
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

#include "io/SampleStats.hh"

#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

#include <boost/filesystem.hpp>

#include "htslib/hts.h"
#include "htslib/sam.h"

#include "core/HtsHelpers.hh"

using std::pair;
using std::string;
using std::vector;

namespace ehunter
{

int extractReadLength(const string& htsFilePath)
{
    samFile* htsFilePtr = sam_open(htsFilePath.c_str(), "r");
    if (!htsFilePtr)
    {
        throw std::runtime_error("Failed to read " + htsFilePath);
    }
    bam_hdr_t* htsHeaderPtr = sam_hdr_read(htsFilePtr);
    if (!htsHeaderPtr)
    {
        throw std::runtime_error("Failed to read the header of " + htsFilePath);
    }

    enum
    {
        kSupplementaryAlign = 0x800,
        kSecondaryAlign = 0x100
    };

    int readLength = 99;
    bam1_t* htsAlignmentPtr = bam_init1();
    int ret;
    while ((ret = sam_read1(htsFilePtr, htsHeaderPtr, htsAlignmentPtr)) >= 0)
    {
        const bool isSupplementary = htsAlignmentPtr->core.flag & kSupplementaryAlign;
        const bool isSecondary = htsAlignmentPtr->core.flag & kSecondaryAlign;
        const bool isPrimaryAlign = (!isSupplementary) && (!isSecondary);
        if (isPrimaryAlign)
        {
            readLength = htsAlignmentPtr->core.l_qseq;
            break;
        }
    }

    if (ret < 0)
    {
        throw std::runtime_error("Failed to extract a read from " + htsFilePath);
    }

    bam_destroy1(htsAlignmentPtr);
    bam_hdr_destroy(htsHeaderPtr);
    sam_close(htsFilePtr);

    return readLength;
}

ReadMode detectReadMode(const string& htsFilePath, const string& htsReferencePath)
{
    // Reads longer than this are considered long reads even when flagged as paired
    const int kMaxShortReadLength = 1000;
    const int kNumReadsToInspect = 1000;

    std::unique_ptr<samFile, decltype(&hts_close)> htsFilePtr(sam_open(htsFilePath.c_str(), "r"), hts_close);
    if (!htsFilePtr)
    {
        throw std::runtime_error("Failed to read " + htsFilePath);
    }

    // Required step for parsing of some CRAMs
    if (hts_set_fai_filename(htsFilePtr.get(), htsReferencePath.c_str()) != 0)
    {
        throw std::runtime_error("Failed to set index of: " + htsReferencePath);
    }

    std::unique_ptr<bam_hdr_t, decltype(&bam_hdr_destroy)> htsHeaderPtr(sam_hdr_read(htsFilePtr.get()), bam_hdr_destroy);
    if (!htsHeaderPtr)
    {
        throw std::runtime_error("Failed to read the header of " + htsFilePath);
    }

    std::unique_ptr<bam1_t, decltype(&bam_destroy1)> htsAlignmentPtr(bam_init1(), bam_destroy1);

    ReadMode readMode = ReadMode::kPaired;
    int numReadsInspected = 0;
    int ret = 0;
    while (numReadsInspected < kNumReadsToInspect
           && (ret = sam_read1(htsFilePtr.get(), htsHeaderPtr.get(), htsAlignmentPtr.get())) >= 0)
    {
        const uint32_t samFlag = htsAlignmentPtr->core.flag;
        if ((samFlag & BAM_FSECONDARY) || (samFlag & BAM_FSUPPLEMENTARY))
        {
            continue;
        }
        ++numReadsInspected;

        const bool isPaired = samFlag & BAM_FPAIRED;
        const bool isLong = htsAlignmentPtr->core.l_qseq > kMaxShortReadLength;
        if (!isPaired || isLong)
        {
            readMode = ReadMode::kAny;
            break;
        }
    }

    if (ret < -1)
    {
        throw std::runtime_error("Failed to extract a record from " + htsFilePath);
    }

    return readMode;
}

ReferenceContigInfo extractReferenceContigInfo(const std::string& htsFilePath)
{
    std::unique_ptr<samFile, decltype(&hts_close)> htsFilePtr(hts_open(htsFilePath.c_str(), "r"), hts_close);
    if (!htsFilePtr)
    {
        throw std::runtime_error("Failed to read " + htsFilePath);
    }

    std::unique_ptr<bam_hdr_t, decltype(&bam_hdr_destroy)> htsHeaderPtr(
        sam_hdr_read(htsFilePtr.get()), bam_hdr_destroy);
    if (!htsHeaderPtr)
    {
        throw std::runtime_error("Failed to read the header of " + htsFilePath);
    }

    return htshelpers::decodeContigInfo(htsHeaderPtr.get());
}

}
