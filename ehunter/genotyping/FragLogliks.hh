//
// ExpansionHunter
// Copyright 2016-2020 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
//         Konrad Scheffler <kscheffler@illumina.com>,
//         Felix Schlesinger <fschlesinger@illumina.com>
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

#include <unordered_map>
#include <utility>
#include <vector>

#include "genotyping/AlignMatrix.hh"

namespace ehunter
{
namespace strgt
{

class FragLogliks
{
public:
    /// \param[in] stutterRate Probability that a read reports one motif more (or one motif fewer) than the allele
    /// that it originated from; zero disables the stutter model
    FragLogliks(int motifLen, int readLen, int fragLen, const AlignMatrix* alignMatrixPtr, double stutterRate = 0.0)
        : motifLen_(motifLen)
        , readLen_(readLen)
        , fragLen_(fragLen)
        , stutterRate_(stutterRate)
        , alignMatrix_(*alignMatrixPtr)
    {
    }

    // Fragments are read pairs or single reads (see AlignMatrix)
    int numFrags() const { return alignMatrix_.numFrags(); }
    double getLoglik(int fragIndex, int alleleMotifCount);

private:
    // Likelihood of the fragment without the stutter model
    double getRawLoglik(int fragIndex, int alleleMotifCount);

    using FragIndex = int;
    using FragIndexAndNumMotifs = std::pair<FragIndex, int>;

    double computeLoglik(const StrAlign& readAlign, const StrAlign& mateAlign, int alleleMotifCount) const;
    double computeLoglik(const StrAlign& readAlign, int alleleMotifCount) const;

    int motifLen_;
    int readLen_;
    int fragLen_;
    double stutterRate_;
    const AlignMatrix& alignMatrix_;
    std::unordered_map<FragIndexAndNumMotifs, double> fragLogliksBySize_;
    std::unordered_map<FragIndexAndNumMotifs, double> rawFragLogliksBySize_;
};

}
}
