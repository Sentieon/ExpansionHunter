//
// ExpansionHunter
// Copyright 2016-2020 Illumina, Inc.
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

#include <limits>
#include <vector>

#include "graphalign/GraphAlignment.hh"

#include "genotyping/StrAlign.hh"

namespace ehunter
{
namespace strgt
{

// Alignments of reads to a repeat, grouped into fragments. A fragment is either a read pair (two consecutive rows)
// or a single read without a mate (one row). Row order is the insertion order.
class AlignMatrix
{
public:
    /// \param[in] smoothIndelAligns Flatten the scores of a read whose alignment to the repeat contains indels over
    /// the allele sizes that this alignment cannot distinguish (see smoothAlignsAroundIndels)
    explicit AlignMatrix(int strNode, bool smoothIndelAligns = false);

    int numReads() const { return alignScoreMatrix_.size(); }
    int numFrags() const { return fragSizes_.size(); }
    int fragSize(int fragIndex) const { return fragSizes_.at(fragIndex); }
    int fragFirstRow(int fragIndex) const { return fragFirstRows_.at(fragIndex); }

    // Adds a read pair (both reads are added if either overlaps the repeat)
    void add(const graphtools::GraphAlignment& read, const graphtools::GraphAlignment& mate);
    // Adds a read without a mate (added if it overlaps or spans the repeat)
    void addSingle(const graphtools::GraphAlignment& read);
    void remove(int readIndex);
    void removeFrag(int fragIndex);
    StrAlign getAlign(int readIndex, int alleleSize) const;
    StrAlign getBestAlign(int readIndex) const;
    const std::vector<std::vector<StrAlign>>& matrix() const { return alignScoreMatrix_; }
    int getMaxMotifCount() const;

    friend void addIrrPairsIfPossibleExpansion(int maxMotifsInRead, AlignMatrix& alignMatrix, int numIrrPairs);
    friend void
    addIrrReadsIfPossibleExpansion(int maxMotifsInRead, AlignMatrix& alignMatrix, int numIrrReads, int motifLen);

private:
    void add(const graphtools::GraphAlignment& graphAlign);
    void smoothAlignsAroundIndels(
        int numMotifsInAlign, const graphtools::GraphAlignment& graphAlign, std::vector<StrAlign>& strAligns) const;
    void rebuildFragFirstRows();
    int findTopInrepeatRow(int maxMotifsInRead) const;

    int strNode_;
    bool smoothIndelAligns_;
    ConsistentAlignmentCalculator alignmentCalculator_;
    std::vector<StrAlign> bestAlignsByRead_;
    std::vector<std::vector<StrAlign>> alignScoreMatrix_;
    std::vector<int> fragSizes_;
    std::vector<int> fragFirstRows_;
};

std::ostream& operator<<(std::ostream& out, const AlignMatrix& matrix);

// Adds copies of the best-scoring long in-repeat read to the matrix as read pairs (paired-end in-repeat read pairs) or
// as single reads (unpaired in-repeat reads) so that they can contribute to the size estimate of long expansions
void addIrrPairsIfPossibleExpansion(int maxMotifsInRead, AlignMatrix& alignMatrix, int numIrrPairs);
void addIrrReadsIfPossibleExpansion(int maxMotifsInRead, AlignMatrix& alignMatrix, int numIrrReads, int motifLen);

}
}
