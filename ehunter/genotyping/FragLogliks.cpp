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

#include "genotyping/FragLogliks.hh"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <ostream>
#include <sstream>
#include <stack>
#include <utility>

#include "spdlog/spdlog.h"

#include "core/LogSum.hh"

using std::pair;
using std::vector;

namespace ehunter
{
namespace strgt
{

// Reads do not always report the number of motifs of the allele that they originated from: polymerase slippage and
// sequencing errors inside a repeat change it, typically by one motif. The likelihood of a fragment is therefore
// averaged over the neighboring allele sizes, weighted by the stutter rate of the sequencing technology. With a
// stutter rate of zero this is the original (Illumina paired-end) model.
double FragLogliks::getLoglik(int fragIndex, int alleleMotifCount)
{
    if (stutterRate_ <= 0)
    {
        return getRawLoglik(fragIndex, alleleMotifCount);
    }

    const FragIndexAndNumMotifs index = { fragIndex, alleleMotifCount };
    auto queryResult = fragLogliksBySize_.find(index);
    if (queryResult != fragLogliksBySize_.end())
    {
        return queryResult->second;
    }

    const double neighborWeight = std::log(stutterRate_);
    double loglik = std::log(1.0 - 2 * stutterRate_) + getRawLoglik(fragIndex, alleleMotifCount);
    loglik = getLogSum(loglik, neighborWeight + getRawLoglik(fragIndex, alleleMotifCount + 1));
    if (alleleMotifCount > 0)
    {
        loglik = getLogSum(loglik, neighborWeight + getRawLoglik(fragIndex, alleleMotifCount - 1));
    }

    fragLogliksBySize_[index] = loglik;
    return loglik;
}

double FragLogliks::getRawLoglik(int fragIndex, int alleleMotifCount)
{
    const FragIndexAndNumMotifs index = { fragIndex, alleleMotifCount };
    auto queryResult = rawFragLogliksBySize_.find(index);
    if (queryResult != rawFragLogliksBySize_.end())
    {
        return queryResult->second;
    }

    const int readIndex = alignMatrix_.fragFirstRow(fragIndex);
    const auto& readAlign = alignMatrix_.getAlign(readIndex, alleleMotifCount);

    double loglik = 0;
    if (alignMatrix_.fragSize(fragIndex) == 2)
    {
        const int mateIndex = readIndex + 1;
        assert(mateIndex < alignMatrix_.numReads());
        const auto& mateAlign = alignMatrix_.getAlign(mateIndex, alleleMotifCount);
        loglik = computeLoglik(readAlign, mateAlign, alleleMotifCount);
    }
    else
    {
        loglik = computeLoglik(readAlign, alleleMotifCount);
    }

    rawFragLogliksBySize_[index] = loglik;
    return loglik;
}

// Likelihood of a read without a mate; the fragment is the read itself
double FragLogliks::computeLoglik(const StrAlign& readAlign, int alleleMotifCount) const
{
    const int numPossibleOrigins = alleleMotifCount * motifLen_ + readLen_ + 1;
    int numOriginsForThisFrag = 1;

    if (readAlign.type() == StrAlign::Type::kInRepeat && readAlign.numMotifs() < alleleMotifCount)
    {
        numOriginsForThisFrag += alleleMotifCount - readAlign.numMotifs();
    }

    const double readAlignLoglik = readAlign.score() * std::log(1.3) - 2 * readLen_ * std::log(2.0);
    return std::log(numOriginsForThisFrag) - std::log(numPossibleOrigins) + readAlignLoglik;
}

double FragLogliks::computeLoglik(const StrAlign& readAlign, const StrAlign& mateAlign, int alleleMotifCount) const
{
    const int numPossibleStarts = alleleMotifCount * motifLen_ + fragLen_ + 1;
    const int numPossibleOrigins = numPossibleStarts * numPossibleStarts / 2;
    int numOriginsForThisFrag = 1;

    if (readAlign.type() != StrAlign::Type::kInRepeat && mateAlign.type() == StrAlign::Type::kInRepeat)
    {
        if (mateAlign.numMotifs() < alleleMotifCount)
        {
            numOriginsForThisFrag += alleleMotifCount - mateAlign.numMotifs();
        }
    }
    else if (readAlign.type() == StrAlign::Type::kInRepeat && mateAlign.type() != StrAlign::Type::kInRepeat)
    {
        if (readAlign.numMotifs() < alleleMotifCount)
        {
            numOriginsForThisFrag += alleleMotifCount - readAlign.numMotifs();
        }
    }
    else if (readAlign.type() == StrAlign::Type::kInRepeat && mateAlign.type() == StrAlign::Type::kInRepeat)
    {
        if (readAlign.numMotifs() < alleleMotifCount)
        {
            const int numReadOrigins = alleleMotifCount - readAlign.numMotifs();
            numOriginsForThisFrag += numReadOrigins * numReadOrigins / 2;
        }
    }

    const double readAlignLoglik = readAlign.score() * std::log(1.3) - 2 * readLen_ * std::log(2.0);
    const double mateAlignLoglik = mateAlign.score() * std::log(1.3) - 2 * readLen_ * std::log(2.0);
    return std::log(numOriginsForThisFrag) - std::log(numPossibleOrigins) + readAlignLoglik + mateAlignLoglik;
}

}
}
