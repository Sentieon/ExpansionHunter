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

#include "genotyping/AlignMatrix.hh"

#include <algorithm>
#include <cassert>

#include "graphalign/LinearAlignmentParameters.hh"

#include <numeric>

using graphtools::GraphAlignment;
using std::vector;

namespace ehunter
{
namespace strgt
{

AlignMatrix::AlignMatrix(int strNode, bool smoothIndelAligns)
    : strNode_(strNode)
    , smoothIndelAligns_(smoothIndelAligns)
    , alignmentCalculator_(strNode)
{
}

StrAlign AlignMatrix::getAlign(int readIndex, int alleleSize) const
{
    if (readIndex >= static_cast<int>(alignScoreMatrix_.size()))
    {
        throw std::runtime_error("Encountered invalid alignment matrix index " + std::to_string(readIndex));
    }

    if (alleleSize < static_cast<int>(alignScoreMatrix_[readIndex].size()))
    {
        return alignScoreMatrix_[readIndex][alleleSize];
    }
    else
    {
        assert(!alignScoreMatrix_[readIndex].empty());
        return alignScoreMatrix_[readIndex].back();
    }
}

void AlignMatrix::add(const GraphAlignment& read, const GraphAlignment& mate)
{
    const int numMotifsInRead = std::count(read.path().begin(), read.path().end(), strNode_);
    const int numMorifsInMate = std::count(mate.path().begin(), mate.path().end(), strNode_);

    // Motif counts are stored in 16 bits; alignments with more motifs (repeats longer than ~65 kb) are skipped
    if (numMotifsInRead >= StrAlign::kMaxNumMotifs || numMorifsInMate >= StrAlign::kMaxNumMotifs)
    {
        return;
    }

    if (numMotifsInRead != 0 || numMorifsInMate != 0)
    {
        add(read);
        add(mate);
        fragSizes_.push_back(2);
        rebuildFragFirstRows();
    }
}

void AlignMatrix::addSingle(const GraphAlignment& read)
{
    const int numMotifsInRead = std::count(read.path().begin(), read.path().end(), strNode_);
    if (numMotifsInRead >= StrAlign::kMaxNumMotifs)
    {
        return;
    }

    bool overlapsLeftFlank = false;
    bool overlapsRightFlank = false;
    for (const auto nodeId : read.path().nodeIds())
    {
        overlapsLeftFlank = overlapsLeftFlank || (static_cast<int>(nodeId) < strNode_);
        overlapsRightFlank = overlapsRightFlank || (static_cast<int>(nodeId) > strNode_);
    }
    const bool spansRepeat = overlapsLeftFlank && overlapsRightFlank;

    if (numMotifsInRead != 0 || spansRepeat)
    {
        add(read);
        fragSizes_.push_back(1);
        rebuildFragFirstRows();
    }
}

void AlignMatrix::rebuildFragFirstRows()
{
    fragFirstRows_.resize(fragSizes_.size());
    int firstRow = 0;
    for (int fragIndex = 0; fragIndex != static_cast<int>(fragSizes_.size()); ++fragIndex)
    {
        fragFirstRows_[fragIndex] = firstRow;
        firstRow += fragSizes_[fragIndex];
    }
    assert(firstRow == numReads());
}

void AlignMatrix::removeFrag(int fragIndex)
{
    assert(fragIndex < numFrags());
    const int firstRow = fragFirstRows_[fragIndex];
    const int size = fragSizes_[fragIndex];
    bestAlignsByRead_.erase(bestAlignsByRead_.begin() + firstRow, bestAlignsByRead_.begin() + firstRow + size);
    alignScoreMatrix_.erase(alignScoreMatrix_.begin() + firstRow, alignScoreMatrix_.begin() + firstRow + size);
    fragSizes_.erase(fragSizes_.begin() + fragIndex);
    rebuildFragFirstRows();
}

void AlignMatrix::add(const graphtools::GraphAlignment& graphAlign)
{
    vector<StrAlign> strAligns;
    const int numMotifsInAlign = std::count(graphAlign.path().begin(), graphAlign.path().end(), strNode_);
    StrAlign alignToMostConsistentAllele = alignmentCalculator_.findConsistentAlignment(numMotifsInAlign, graphAlign);

    bestAlignsByRead_.push_back(alignToMostConsistentAllele);

    for (int numMotifs = numMotifsInAlign - 1; numMotifs != -1; --numMotifs)
    {
        StrAlign align = alignmentCalculator_.findConsistentAlignment(numMotifs, graphAlign);
        strAligns.emplace_back(align);
    }
    std::reverse(strAligns.begin(), strAligns.end());
    strAligns.emplace_back(alignToMostConsistentAllele);
    StrAlign previousAlign = strAligns.back();
    for (int numMotifs = numMotifsInAlign + 1; numMotifs <= StrAlign::kMaxNumMotifs; ++numMotifs)
    {
        StrAlign align = alignmentCalculator_.findConsistentAlignment(numMotifs, graphAlign);
        if (align.type() == previousAlign.type() && align.score() == previousAlign.score())
        {
            break;
        }
        strAligns.emplace_back(align);
        previousAlign = align;
    }

    if (smoothIndelAligns_)
    {
        smoothAlignsAroundIndels(numMotifsInAlign, graphAlign, strAligns);
    }

    alignScoreMatrix_.emplace_back(strAligns);
}

// The number of motifs of an alignment that contains indels inside the repeat does not determine the size of the
// allele that the read originated from: a single-base sequencing indel inside the repeat typically shifts the motif
// count of the best alignment by one, because aligning the remaining bases to one more (or one fewer) motif is
// cheaper than leaving them unaligned. Such a read supports all allele sizes that its alignment cannot distinguish
// equally well, so its scores are flattened over them; otherwise the alignment of the read to the motif count that
// the aligner happened to pick would be scored as a perfect match and the other sizes would be charged the full
// stutter penalty, turning sequencing indels into support for a spurious second allele.
void AlignMatrix::smoothAlignsAroundIndels(
    int numMotifsInAlign, const graphtools::GraphAlignment& graphAlign, vector<StrAlign>& strAligns) const
{
    if (numMotifsInAlign < 0 || numMotifsInAlign >= static_cast<int>(strAligns.size()))
    {
        return;
    }

    const StrAlign bestAlign = strAligns[numMotifsInAlign];
    if (bestAlign.numIndels() == 0)
    {
        return;
    }

    const int motifLen = graphAlign.path().graphRawPtr()->nodeSeq(strNode_).length();
    if (motifLen <= 0)
    {
        return;
    }

    // The aligned repeat length differs from the length of the allele by at most numIndels bases, so allele sizes
    // within this many motifs of the alignment are indistinguishable
    const int numAmbiguousMotifs = (bestAlign.numIndels() + motifLen - 1) / motifLen;
    const int firstIndex = std::max(0, numMotifsInAlign - numAmbiguousMotifs);

    // The last row entry is also used for every allele larger than the row (see getAlign), so it must keep the score
    // of an allele that the read cannot have come from
    const int lastRowIndex = static_cast<int>(strAligns.size()) - 1;
    const int lastIndex
        = std::min(numMotifsInAlign + numAmbiguousMotifs, std::max(numMotifsInAlign, lastRowIndex - 1));

    for (int index = firstIndex; index <= lastIndex; ++index)
    {
        if (strAligns[index].score() < bestAlign.score())
        {
            strAligns[index] = StrAlign(bestAlign.type(), index, bestAlign.score(), bestAlign.numIndels());
        }
    }
}

int AlignMatrix::getMaxMotifCount() const
{
    int maxMotifCount = 0;
    for (const auto& readAligns : alignScoreMatrix_)
    {
        const int maxMotifsInAlign = static_cast<int>(readAligns.size()) - 1;
        if (maxMotifsInAlign > maxMotifCount)
        {
            maxMotifCount = maxMotifsInAlign;
        }
    }

    return maxMotifCount;
}

StrAlign AlignMatrix::getBestAlign(int readIndex) const { return bestAlignsByRead_[readIndex]; }

std::ostream& operator<<(std::ostream& out, const AlignMatrix& matrix)
{
    const auto& alignsByRead = matrix.matrix();

    const unsigned alignCount(alignsByRead.size());
    std::vector<unsigned> alignIndexes(alignCount);
    std::iota(alignIndexes.begin(), alignIndexes.end(), 0);

    // Sort alignIndexes to canonicalize the alignmatrix output
    std::sort(
        alignIndexes.begin(), alignIndexes.end(),
        [&alignsByRead](const unsigned ai, const unsigned bi) -> bool
        {
            const auto& a(alignsByRead[ai]);
            const auto& b(alignsByRead[bi]);
            if (a.size() < b.size())
            {
                return true;
            }
            if (b.size() < a.size())
            {
                return false;
            }
            for (unsigned i(0); i < a.size(); ++i)
            {
                if (a[i] < b[i])
                {
                    return true;
                }
                if (b[i] < a[i])
                {
                    return false;
                }
            }
            return false;
        });

    auto dumpStrAlign = [&out](const StrAlign& align)
    {
        out << "(";
        switch (align.type())
        {
        case StrAlign::Type::kOutside:
            out << "O";
            break;
        case StrAlign::Type::kInRepeat:
            out << "I";
            break;
        case StrAlign::Type::kSpanning:
            out << "S";
            break;
        case StrAlign::Type::kFlanking:
            out << "F";
            break;
        }
        out << "," << align.numMotifs() << "," << align.score() << "), ";
    };

    for (unsigned origAlignIndex(0); origAlignIndex < alignCount; ++origAlignIndex)
    {
        const unsigned alignIndex(alignIndexes[origAlignIndex]);
        const auto& aligns(alignsByRead[alignIndex]);
        for (const StrAlign& align : aligns)
        {
            dumpStrAlign(align);
        }
        out << "\n";
    }

    return out;
}

int AlignMatrix::findTopInrepeatRow(int maxMotifsInRead) const
{
    // Find highest-scoring in-repeat read
    assert(bestAlignsByRead_.size() == alignScoreMatrix_.size());

    const int longIrrLowerBound = static_cast<int>(0.90 * maxMotifsInRead);

    int topIrrIndex = -1;
    int topIrrScore = -1;
    for (int readIndex = 0; readIndex != numReads(); ++readIndex)
    {
        const StrAlign& align = bestAlignsByRead_[readIndex];
        const bool isLongIrr = align.type() == StrAlign::Type::kInRepeat && align.numMotifs() >= longIrrLowerBound;

        if (isLongIrr && align.score() > topIrrScore)
        {
            topIrrScore = align.score();
            topIrrIndex = readIndex;
        }
    }

    return topIrrIndex;
}

void addIrrPairsIfPossibleExpansion(int maxMotifsInRead, AlignMatrix& alignMatrix, int numIrrPairs)
{
    const int topIrrIndex = alignMatrix.findTopInrepeatRow(maxMotifsInRead);
    if (topIrrIndex == -1)
    {
        return;
    }

    StrAlign irrTopAlign = alignMatrix.bestAlignsByRead_[topIrrIndex];
    vector<StrAlign> irrAligns = alignMatrix.alignScoreMatrix_[topIrrIndex];

    for (int newIrrIndex = 0; newIrrIndex != 2 * numIrrPairs; ++newIrrIndex)
    {
        alignMatrix.bestAlignsByRead_.push_back(irrTopAlign);
        alignMatrix.alignScoreMatrix_.push_back(irrAligns);
    }
    for (int newFragIndex = 0; newFragIndex != numIrrPairs; ++newFragIndex)
    {
        alignMatrix.fragSizes_.push_back(2);
    }
    alignMatrix.rebuildFragFirstRows();
}

void addIrrReadsIfPossibleExpansion(int maxMotifsInRead, AlignMatrix& alignMatrix, int numIrrReads, int motifLen)
{
    StrAlign irrTopAlign(StrAlign::Type::kInRepeat, 0, 0, 0);
    vector<StrAlign> irrAligns;

    const int topIrrIndex = alignMatrix.findTopInrepeatRow(maxMotifsInRead);
    if (topIrrIndex != -1)
    {
        irrTopAlign = alignMatrix.bestAlignsByRead_[topIrrIndex];
        irrAligns = alignMatrix.alignScoreMatrix_[topIrrIndex];
    }
    else
    {
        // Reads without mates cannot be anchored inside the repeat, so in-repeat reads usually reach the matrix only
        // when they are found in off-target regions. Use a perfect in-repeat read as their representative, provided
        // that flanking reads show the repeat extending at least half a read length beyond a flank.
        int longestFlankingSize = 0;
        for (int readIndex = 0; readIndex != alignMatrix.numReads(); ++readIndex)
        {
            const StrAlign& align = alignMatrix.bestAlignsByRead_[readIndex];
            if (align.type() == StrAlign::Type::kFlanking)
            {
                longestFlankingSize = std::max(longestFlankingSize, align.numMotifs());
            }
        }

        if (maxMotifsInRead <= 0 || 2 * longestFlankingSize < maxMotifsInRead)
        {
            return;
        }

        const int matchScore = LinearAlignmentParameters().matchScore;
        for (int numMotifs = 0; numMotifs <= maxMotifsInRead; ++numMotifs)
        {
            irrAligns.emplace_back(StrAlign::Type::kInRepeat, numMotifs, numMotifs * motifLen * matchScore, 0);
        }
        irrTopAlign = irrAligns.back();
    }

    for (int newIrrIndex = 0; newIrrIndex != numIrrReads; ++newIrrIndex)
    {
        alignMatrix.bestAlignsByRead_.push_back(irrTopAlign);
        alignMatrix.alignScoreMatrix_.push_back(irrAligns);
        alignMatrix.fragSizes_.push_back(1);
    }
    alignMatrix.rebuildFragFirstRows();
}

void AlignMatrix::remove(int readIndex)
{
    assert(readIndex < numReads());
    bestAlignsByRead_.erase(bestAlignsByRead_.begin() + readIndex);
    alignScoreMatrix_.erase(alignScoreMatrix_.begin() + readIndex);

    // Shrink the fragment that contained the row
    for (int fragIndex = 0; fragIndex != numFrags(); ++fragIndex)
    {
        if (readIndex < fragFirstRows_[fragIndex] + fragSizes_[fragIndex])
        {
            if (--fragSizes_[fragIndex] == 0)
            {
                fragSizes_.erase(fragSizes_.begin() + fragIndex);
            }
            break;
        }
    }
    rebuildFragFirstRows();
}

}
}
