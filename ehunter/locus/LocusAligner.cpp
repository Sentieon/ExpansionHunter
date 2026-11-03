//
// ExpansionHunter
// Copyright 2016-2021 Illumina, Inc.
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

#include "locus/LocusAligner.hh"

#include "spdlog/spdlog.h"

#include "alignment/AlignmentFilters.hh"
#include "alignment/OperationsOnAlignments.hh"

namespace ehunter
{
namespace locus
{

LocusAligner::LocusAligner(
    std::string locusId, GraphPtr graph, const HeuristicParameters& params, AlignWriterPtr writer)
    : locusId_(std::move(locusId))
    , aligner_(graph, params.kmerLenForAlignment(), params.paddingLength(), params.seedAffixTrimLength())
    , orientationPredictor_(graph, params.orientationPredictorKmerLen(), params.orientationPredictorMinKmerCount())
    , writer_(std::move(writer))
    , readMode_(params.readMode())
    , minSingletonMapq_(params.minSingletonMapq())
{
}

LocusAligner::AlignedPair LocusAligner::align(Read& read, Read* mate, graphtools::AlignerSelector& alignerSelector)
{
    auto readAlign = align(read, alignerSelector);
    auto mateAlign = mate ? align(*mate, alignerSelector) : boost::none;

    int numMatchingBases = static_cast<int>(static_cast<double>(read.sequence().length()) / 7.5);
    numMatchingBases = std::max(numMatchingBases, 10);
    LinearAlignmentParameters parameters;
    const int kMinNonRepeatAlignmentScore = numMatchingBases * parameters.matchScore;

    if (!checkIfLocallyPlacedReadPair(readAlign, mateAlign, kMinNonRepeatAlignmentScore))
    {
        // A read without a mate cannot be anchored by its mate; if it is confidently placed at the locus by the linear
        // aligner it is kept as (in-repeat) evidence provided its graph alignment is of good quality
        const bool keepUnanchoredSingleton = readMode_ == ReadMode::kAny && mate == nullptr && readAlign
            && read.mapq() >= minSingletonMapq_ && readAlign->referenceLength() > 0
            && checkIfPassesAlignmentFilters(*readAlign);

        if (!keepUnanchoredSingleton)
        {
            return { boost::none, boost::none };
        }
    }

    if (readAlign && mateAlign)
    {
        writer_->write(
            locusId_, read.fragmentId(), read.sequence(), read.isFirstMate(), read.isReversed(), read.isReversed(),
            *readAlign);
        writer_->write(
            locusId_, mate->fragmentId(), mate->sequence(), mate->isFirstMate(), mate->isReversed(), mate->isReversed(),
            *mateAlign);
    }
    else if (readMode_ == ReadMode::kAny && mate == nullptr && readAlign)
    {
        if (read.isUnpaired())
        {
            writer_->writeUnpaired(locusId_, read.fragmentId(), read.sequence(), read.isReversed(), *readAlign);
        }
        else
        {
            writer_->write(
                locusId_, read.fragmentId(), read.sequence(), read.isFirstMate(), read.isReversed(), read.isReversed(),
                *readAlign);
        }
    }

    return { readAlign, mateAlign };
}

LocusAligner::OptionalAlign LocusAligner::align(Read& read, graphtools::AlignerSelector& alignerSelector) const
{
    OrientationPrediction predictedOrientation = orientationPredictor_.predict(read.sequence());

    if (predictedOrientation == OrientationPrediction::kAlignsInReverseComplementOrientation)
    {
        read.reverseComplement();
    }
    else if (predictedOrientation == OrientationPrediction::kDoesNotAlign)
    {
        return {};
    }

    if (readMode_ == ReadMode::kPaired)
    {
        auto readAligns = aligner_.align(read.sequence(), alignerSelector);
        if (readAligns.empty())
        {
            return {};
        }

        return computeCanonicalAlignment(readAligns);
    }

    // In read mode "any" (long and unusual reads) a read that cannot be aligned is skipped instead of aborting the run
    try
    {
        auto readAligns = aligner_.align(read.sequence(), alignerSelector);
        if (readAligns.empty())
        {
            return {};
        }

        return computeCanonicalAlignment(readAligns);
    }
    catch (const std::exception& e)
    {
        spdlog::warn(
            "Skipping read {} at locus {} ({} bp): {}", read.fragmentId(), locusId_, read.sequence().length(),
            std::string(e.what()).substr(0, 200));
        return {};
    }
}

}
}
