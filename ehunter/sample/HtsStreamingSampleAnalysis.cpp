//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
//         Chris Saunders <csaunders@illumina.com>
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

#include "sample/HtsStreamingSampleAnalysis.hh"

#include <memory>

#include "absl/container/flat_hash_map.h"
#include "spdlog/spdlog.h"
#include <boost/optional.hpp>

#include "core/HtsHelpers.hh"
#include "core/ThreadPool.hh"
#include "locus/LocusAnalyzer.hh"
#include "locus/LocusAnalyzerUtil.hh"
#include "sample/GenomeQueryCollection.hh"
#include "sample/HtsFileStreamer.hh"
#include "sample/HtsStreamingReadPairQueue.hh"

using boost::optional;
using ehunter::locus::initializeLocusAnalyzers;
using ehunter::locus::LocusAnalyzer;
using graphtools::AlignmentWriter;
using std::string;
using std::vector;

namespace ehunter
{

namespace
{

/// \brief Mutable data shared by all LocusAnalyzer-processing threads
///
class LocusAnalyzerThreadSharedData
{
public:
    LocusAnalyzerThreadSharedData(const unsigned maxActiveLocusAnalyzerQueues, const unsigned locusAnalyzerCount)
        : isWorkerThreadException(false)
        , readPairQueue(maxActiveLocusAnalyzerQueues, locusAnalyzerCount)
    {
    }

    std::atomic<bool> isWorkerThreadException;
    HtsStreamingReadPairQueue readPairQueue;
    vector<std::unique_ptr<LocusAnalyzer>> locusAnalyzers;
};

/// \brief Data isolated to each LocusAnalyzer-processing thread
///
struct LocusAnalyzerThreadLocalData
{
    std::exception_ptr threadExceptionPtr = nullptr;
    std::shared_ptr<graphtools::AlignerSelector> alignerSelectorPtr;
};

/// Process queue for a single LocusAnalyzer on one thread
void processLocusAnalyzerQueue(
    const int threadIndex, LocusAnalyzerThreadSharedData& locusAnalyzerThreadSharedData,
    std::vector<LocusAnalyzerThreadLocalData>& locusAnalyzerThreadLocalDataPool, const unsigned locusIndex)
{
    if (locusAnalyzerThreadSharedData.isWorkerThreadException.load())
    {
        return;
    }

    LocusAnalyzerThreadLocalData& locusAnalyzerThreadData(locusAnalyzerThreadLocalDataPool[threadIndex]);
    auto& locusAnalyzer(*locusAnalyzerThreadSharedData.locusAnalyzers[locusIndex]);

    boost::optional<HtsStreamingReadPairQueue::ReadPair> readPair;

    try
    {
        while (true)
        {
            locusAnalyzerThreadSharedData.readPairQueue.getNextReadPair(locusIndex, readPair);
            if (not readPair)
            {
                break;
            }
            Read* matePtr = readPair->mate ? &*readPair->mate : nullptr;
            processAnalyzerBundleReadPair(
                locusAnalyzer, readPair->regionType, readPair->inputType, readPair->read, matePtr,
                *locusAnalyzerThreadData.alignerSelectorPtr);
        }
    }
    catch (const std::exception& e)
    {
        locusAnalyzerThreadSharedData.isWorkerThreadException.store(true);
        locusAnalyzerThreadData.threadExceptionPtr = std::current_exception();

        std::ostringstream oss;
        oss << "Exception caught in thread " << threadIndex << " while processing read pair queue for locus: `"
            << locusAnalyzer.locusId() << "`";
        if (readPair)
        {
            oss << " current readPair: `" << readPair->read.fragmentId() << "`";
        }
        oss << ": " << e.what();
        spdlog::error(oss.str());
        throw;
    }
    catch (...)
    {
        locusAnalyzerThreadSharedData.isWorkerThreadException.store(true);
        locusAnalyzerThreadData.threadExceptionPtr = std::current_exception();

        std::ostringstream oss;
        oss << "Exception caught in thread " << threadIndex << " while processing read pair queue for locus: `"
            << locusAnalyzer.locusId() << "`";
        if (readPair)
        {
            oss << " current readPair: `" << readPair->read.fragmentId() << "`";
        }
        spdlog::error(oss.str());
        throw;
    }
}

/// \brief Mutable data shared by all SampleFindings-processing threads
///
class SampleFindingsThreadSharedData
{
public:
    SampleFindingsThreadSharedData()
        : isWorkerThreadException(false)
        , locusIndex(0)
    {
    }

    std::atomic<bool> isWorkerThreadException;
    std::atomic<unsigned> locusIndex;
};

/// \brief Data isolated to each SampleFindings-processing thread
///
struct SampleFindingsThreadLocalData
{
    std::exception_ptr threadExceptionPtr = nullptr;
};

/// \brief Analyze a series of loci on one thread
///
void analyzeLocus(
    const int threadIndex, const Sex sampleSex, vector<std::unique_ptr<LocusAnalyzer>>& locusAnalyzers,
    SampleFindings& sampleFindings, SampleFindingsThreadSharedData& sampleFindingsThreadSharedData,
    std::vector<SampleFindingsThreadLocalData>& sampleFindingsThreadLocalData)
{
    SampleFindingsThreadLocalData& sampleFindingsThreadData(sampleFindingsThreadLocalData[threadIndex]);
    std::string locusId = "Unknown";

    try
    {

        const unsigned size(locusAnalyzers.size());
        while (true)
        {
            if (sampleFindingsThreadSharedData.isWorkerThreadException.load())
            {
                return;
            }
            const auto locusIndex(sampleFindingsThreadSharedData.locusIndex.fetch_add(1));
            if (locusIndex >= size)
            {
                return;
            }

            auto& locusAnalyzer(*locusAnalyzers[locusIndex]);
            locusId = locusAnalyzer.locusId();
            sampleFindings[locusIndex] = locusAnalyzer.analyze(sampleSex, boost::none);
        }
    }
    catch (const std::exception& e)
    {
        sampleFindingsThreadSharedData.isWorkerThreadException = true;
        sampleFindingsThreadData.threadExceptionPtr = std::current_exception();

        spdlog::error("Exception caught in thread {} while analyzing locus: {} : {}", threadIndex, locusId, e.what());
        throw;
    }
    catch (...)
    {
        sampleFindingsThreadSharedData.isWorkerThreadException = true;
        sampleFindingsThreadData.threadExceptionPtr = std::current_exception();

        spdlog::error("Unknown exception caught in thread {} while analyzing locus: {}", threadIndex, locusId);
        throw;
    }
}

/// A read waiting for its mate to be streamed. The linear alignment is retained so that the read can be clipped to
/// each locus it is sent to (only needed in ReadMode::kAny; the cigar is left empty otherwise)
struct PendingRead
{
    Read read;
    int32_t contigId;
    int64_t pos;
    int64_t endPos;
    vector<uint32_t> cigar;
};

}

SampleFindings htsStreamingSampleAnalysis(
    const InputPaths& inputPaths, Sex sampleSex, const HeuristicParameters& heuristicParams, const int threadCount,
    const RegionCatalog& regionCatalog, locus::AlignWriterPtr bamletWriter)
{
    const ReadMode readMode = heuristicParams.readMode();
    const int32_t maxReadLength = heuristicParams.maxLocusReadLength();

    // Setup thread-specific data structures and thread pool
    const unsigned maxActiveLocusAnalyzerQueues(threadCount + 5);
    const unsigned locusAnalyzerCount(regionCatalog.size());
    LocusAnalyzerThreadSharedData locusAnalyzerThreadSharedData(maxActiveLocusAnalyzerQueues, locusAnalyzerCount);
    std::vector<LocusAnalyzerThreadLocalData> locusAnalyzerThreadLocalDataPool(threadCount);
    for (int threadIndex(0); threadIndex < threadCount; ++threadIndex)
    {
        auto& locusAnalyzerThreadData(locusAnalyzerThreadLocalDataPool[threadIndex]);
        locusAnalyzerThreadData.alignerSelectorPtr.reset(
            new graphtools::AlignerSelector(heuristicParams.alignerType()));
    }
    ctpl::thread_pool pool(threadCount);

    spdlog::info("Initializing all loci");
    graphtools::AlignerSelector alignerSelector(heuristicParams.alignerType());
    locusAnalyzerThreadSharedData.locusAnalyzers
        = initializeLocusAnalyzers(regionCatalog, heuristicParams, bamletWriter, threadCount);
    GenomeQueryCollection genomeQuery(locusAnalyzerThreadSharedData.locusAnalyzers, readMode == ReadMode::kPaired);

    // Reference spans that reads are clipped to before being sent to each locus (ReadMode::kAny only)
    vector<optional<GenomicRegion>> clipRegions;
    if (readMode == ReadMode::kAny)
    {
        for (const auto& locusSpec : regionCatalog)
        {
            clipRegions.push_back(locusSpec.readClippingRegion());
        }
    }

    spdlog::info("Streaming reads");

    using PendingReadCatalog = absl::flat_hash_map<string, PendingRead>;
    PendingReadCatalog unpairedReads(1000);

    auto sendReadPair = [&](const AnalyzerBundle& bundle, HtsStreamingReadPairQueue::ReadPair readPair)
    {
        if (locusAnalyzerThreadSharedData.readPairQueue.insertReadPair(bundle.locusIndex, std::move(readPair)))
        {
            pool.push(
                processLocusAnalyzerQueue, std::ref(locusAnalyzerThreadSharedData),
                std::ref(locusAnalyzerThreadLocalDataPool), bundle.locusIndex);
        }
    };

    const unsigned htsDecompressionThreads(std::min(threadCount, 12));
    htshelpers::HtsFileStreamer readStreamer(inputPaths.htsFile(), inputPaths.reference(), htsDecompressionThreads);
    while (readStreamer.trySeekingToNextPrimaryAlignment() && readStreamer.isStreamingAlignedReads())
    {
        // Stop processing reads if an exception is thrown in the worker pool:
        if (locusAnalyzerThreadSharedData.isWorkerThreadException.load())
        {
            break;
        }

        const int32_t readContigId = readStreamer.currentReadContigId();
        const int64_t readPos = readStreamer.currentReadPosition();
        const int64_t readEndPos = std::max<int64_t>(readStreamer.currentReadEndPosition(), readPos + 1);

        const bool isReadNearTargetRegion = readMode == ReadMode::kPaired
            ? genomeQuery.targetRegionMask.query(readContigId, readPos)
            : genomeQuery.targetRegionMask.queryRange(readContigId, readPos, readEndPos);
        const bool isMateNearTargetRegion = genomeQuery.targetRegionMask.query(
            readStreamer.currentMateContigId(), readStreamer.currentMatePosition());
        if (!isReadNearTargetRegion && !isMateNearTargetRegion)
        {
            continue;
        }

        if (not readStreamer.currentIsPaired())
        {
            if (readMode == ReadMode::kPaired)
            {
                continue;
            }

            // Read sequenced without a mate: send it to every locus it overlaps
            vector<AnalyzerBundle> analyzerBundles
                = genomeQuery.analyzerFinder.query(readContigId, readPos, readEndPos);
            for (const auto& bundle : analyzerBundles)
            {
                const auto& clipRegion = clipRegions[bundle.locusIndex];
                optional<Read> read = clipRegion ? readStreamer.decodeReadClippedTo(*clipRegion, maxReadLength)
                                                 : optional<Read>(readStreamer.decodeRead());
                if (read)
                {
                    sendReadPair(
                        bundle, { bundle.regionType, AnalyzerInputType::kReadOnly, std::move(*read), boost::none });
                }
            }
            continue;
        }

        Read read = readStreamer.decodeRead();
        const auto mateIterator = unpairedReads.find(read.fragmentId());
        if (mateIterator == unpairedReads.end())
        {
            const string fragmentId = read.fragmentId();
            vector<uint32_t> cigar = readMode == ReadMode::kAny ? readStreamer.currentCigar() : vector<uint32_t>();
            unpairedReads.emplace(
                fragmentId, PendingRead { std::move(read), readContigId, readPos, readEndPos, std::move(cigar) });
            continue;
        }
        PendingRead mate = std::move(mateIterator->second);
        unpairedReads.erase(mateIterator);

        int64_t readEnd = readPos + read.sequence().length();
        int64_t mateEnd = readStreamer.currentMatePosition() + mate.read.sequence().length();
        if (readMode == ReadMode::kAny)
        {
            readEnd = readEndPos;
            mateEnd = mate.endPos;
        }

        vector<AnalyzerBundle> analyzerBundles = genomeQuery.analyzerFinder.query(
            readContigId, readPos, readEnd, readStreamer.currentMateContigId(), readStreamer.currentMatePosition(),
            mateEnd);

        const unsigned bundleCount(analyzerBundles.size());
        for (unsigned bundleIndex(0); bundleIndex < bundleCount; ++bundleIndex)
        {
            auto& bundle(analyzerBundles[bundleIndex]);

            if (readMode == ReadMode::kPaired)
            {
                if ((bundleIndex + 1) < bundleCount)
                {
                    sendReadPair(bundle, { bundle.regionType, bundle.inputType, read, mate.read });
                }
                else
                {
                    sendReadPair(
                        bundle, { bundle.regionType, bundle.inputType, std::move(read), std::move(mate.read) });
                }
                continue;
            }

            // Clip both mates to the locus before sending them; a mate whose overlapping piece is empty is dropped
            const auto& clipRegion = clipRegions[bundle.locusIndex];
            optional<Read> clippedRead = read;
            optional<Read> clippedMate = mate.read;
            if (clipRegion)
            {
                clippedRead = htshelpers::clipReadToRegion(
                    read, readContigId, readPos, readStreamer.currentCigar(), *clipRegion, maxReadLength);
                clippedMate = htshelpers::clipReadToRegion(
                    mate.read, mate.contigId, mate.pos, mate.cigar, *clipRegion, maxReadLength);
            }

            if (clippedRead && clippedMate)
            {
                sendReadPair(
                    bundle, { bundle.regionType, bundle.inputType, std::move(*clippedRead), std::move(clippedMate) });
            }
            else if (clippedRead && bundle.inputType != AnalyzerInputType::kMateOnly)
            {
                sendReadPair(
                    bundle,
                    { bundle.regionType, AnalyzerInputType::kReadOnly, std::move(*clippedRead), boost::none });
            }
            else if (clippedMate && bundle.inputType != AnalyzerInputType::kReadOnly)
            {
                sendReadPair(
                    bundle,
                    { bundle.regionType, AnalyzerInputType::kReadOnly, std::move(*clippedMate), boost::none });
            }
        }
    }

    pool.stop(true);

    // Rethrow exceptions from the pool in thread order:
    if (locusAnalyzerThreadSharedData.isWorkerThreadException.load())
    {
        for (int threadIndex(0); threadIndex < threadCount; ++threadIndex)
        {
            const LocusAnalyzerThreadLocalData& locusAnalyzerThreadData(locusAnalyzerThreadLocalDataPool[threadIndex]);
            if (locusAnalyzerThreadData.threadExceptionPtr)
            {
                std::rethrow_exception(locusAnalyzerThreadData.threadExceptionPtr);
            }
        }
    }

    spdlog::info("Analyzing read evidence");

    SampleFindingsThreadSharedData sampleFindingsThreadSharedData;
    std::vector<SampleFindingsThreadLocalData> sampleFindingsThreadLocalDataPool(threadCount);

    const unsigned locusCount(locusAnalyzerThreadSharedData.locusAnalyzers.size());
    SampleFindings sampleFindings(locusCount);

    // Start all sampleFindings worker threads
    std::vector<std::thread> sampleFindingsThreads;
    for (int threadIndex(0); threadIndex < threadCount; ++threadIndex)
    {
        sampleFindingsThreads.emplace_back(
            analyzeLocus, threadIndex, sampleSex, std::ref(locusAnalyzerThreadSharedData.locusAnalyzers),
            std::ref(sampleFindings), std::ref(sampleFindingsThreadSharedData),
            std::ref(sampleFindingsThreadLocalDataPool));
    }

    // Rethrow exceptions from worker pool in thread order:
    if (sampleFindingsThreadSharedData.isWorkerThreadException.load())
    {
        for (int threadIndex(0); threadIndex < threadCount; ++threadIndex)
        {
            const auto& sampleFindingsThreadData(sampleFindingsThreadLocalDataPool[threadIndex]);
            if (sampleFindingsThreadData.threadExceptionPtr)
            {
                std::rethrow_exception(sampleFindingsThreadData.threadExceptionPtr);
            }
        }
    }

    for (int threadIndex(0); threadIndex < threadCount; ++threadIndex)
    {
        sampleFindingsThreads[threadIndex].join();
    }

    return sampleFindings;
}

}
