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

#include <memory>
#include <string>
#include <vector>

#include <boost/optional.hpp>

#include "graphalign/GraphAlignment.hh"
#include "graphcore/Graph.hh"

#include "alignment/AlignmentClassifier.hh"
#include "alignment/GraphVariantAlignmentStats.hh"
#include "core/Read.hh"
#include "genotyping/AlignMatrix.hh"
#include "genotyping/RepeatGenotype.hh"
#include "locus/VariantAnalyzer.hh"

namespace ehunter
{

class RepeatAnalyzer : public VariantAnalyzer
{
public:
    /// \param[in] readMode Reads of indel-prone technologies (single-end short reads and long reads) are only
    /// admitted in ReadMode::kAny, where the motif count of an alignment containing indels is treated as uncertain
    ///
    /// \param[in] stutterRate Probability that a read reports one motif more (or fewer) than the allele it came
    /// from; only used in ReadMode::kAny
    RepeatAnalyzer(
        std::string variantId, const graphtools::Graph& graph, graphtools::NodeId repeatNodeId,
        GenotyperParameters genotyperParams, ReadMode readMode = ReadMode::kPaired, double stutterRate = 0.0)
        : VariantAnalyzer(std::move(variantId), graph, { repeatNodeId }, genotyperParams)
        , repeatUnit_(graph.nodeSeq(repeatNodeId))
        , alignmentStatsCalculator_({ repeatNodeId })
        , countOfInrepeatReadPairs_(0)
        , stutterRate_(readMode == ReadMode::kAny ? stutterRate : 0.0)
        , alignMatrix_(repeatNodeId, readMode == ReadMode::kAny)
    {
    }

    ~RepeatAnalyzer() override = default;

    const std::string& repeatUnit() const { return repeatUnit_; }
    void addInrepeatReadPair() { countOfInrepeatReadPairs_++; }
    void addInrepeatRead() { countOfInrepeatReads_++; }

    void processMates(
        const Read& read, const graphtools::GraphAlignment& readAlignment, const Read& mate,
        const graphtools::GraphAlignment& mateAlignment) override;
    void processRead(const Read& read, const graphtools::GraphAlignment& readAlignment) override;

    std::unique_ptr<VariantFindings> analyze(const LocusStats& stats) override;

private:
    graphtools::NodeId repeatNodeId() const { return nodeIds_.front(); }

    const std::string repeatUnit_;
    GraphVariantAlignmentStatsCalculator alignmentStatsCalculator_;

    int countOfInrepeatReadPairs_;
    int countOfInrepeatReads_ = 0;
    double stutterRate_;
    strgt::AlignMatrix alignMatrix_;
};

}
