//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
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

// Checks that the linear-time projectAlignmentOntoGraph produces exactly the same graph alignments as the original
// implementation based on repeated Alignment::splitAtReferencePosition calls.

#include <list>
#include <random>
#include <vector>

#include "gtest/gtest.h"

#include "graphalign/GraphAlignment.hh"
#include "graphalign/GraphAlignmentOperations.hh"
#include "graphalign/LinearAlignment.hh"
#include "graphcore/Graph.hh"
#include "graphcore/Path.hh"

#include "io/GraphBlueprint.hh"
#include "io/RegionGraph.hh"

using graphtools::Alignment;
using graphtools::Graph;
using graphtools::GraphAlignment;
using graphtools::NodeId;
using graphtools::Operation;
using graphtools::OperationType;
using graphtools::Path;
using std::list;
using std::vector;

namespace
{
// Original v5.0.0 implementation
GraphAlignment projectAlignmentOntoGraphReference(Alignment linear_alignment, Path path)
{
    vector<Alignment> alignments;

    path.shrinkBy(
        linear_alignment.referenceStart(),
        path.length() - linear_alignment.referenceStart() - linear_alignment.referenceLength());
    linear_alignment.setReferenceStart(0);

    for (size_t node_index = 0; node_index != path.numNodes(); ++node_index)
    {
        const size_t last_position_of_path_on_this_node = path.getNodeOverlapLengthByIndex(node_index);

        const size_t linear_alignment_last_position
            = linear_alignment.referenceStart() + linear_alignment.referenceLength();

        if (linear_alignment_last_position <= last_position_of_path_on_this_node)
        {
            alignments.push_back(linear_alignment);
            break;
        }
        else
        {
            Alignment suffix = linear_alignment.splitAtReferencePosition(last_position_of_path_on_this_node);
            alignments.push_back(linear_alignment);
            linear_alignment = suffix;
            linear_alignment.setReferenceStart(0);
        }
    }

    Alignment& first_alignment = alignments.front();
    first_alignment.setReferenceStart(path.startPosition());

    return GraphAlignment(path, alignments);
}

// Random linear alignment spanning exactly referenceLength bases of the reference
list<Operation> randomOperations(uint32_t referenceLength, std::mt19937& rng)
{
    list<Operation> operations;
    std::uniform_int_distribution<int> lengthDist(1, 7);
    std::uniform_int_distribution<int> typeDist(0, 9);

    if (typeDist(rng) < 3)
    {
        operations.emplace_back(OperationType::kSoftclip, lengthDist(rng));
    }

    uint32_t used = 0;
    OperationType previousType = OperationType::kSoftclip;
    while (used < referenceLength)
    {
        const int typeCode = typeDist(rng);
        OperationType type = OperationType::kMatch;
        if (typeCode == 6 || typeCode == 7)
        {
            type = OperationType::kMismatch;
        }
        else if (typeCode == 8)
        {
            type = OperationType::kInsertionToRef;
        }
        else if (typeCode == 9)
        {
            type = OperationType::kDeletionFromRef;
        }

        // Avoid adjacent operations of the same type and indels next to each other
        const bool isIndel = type == OperationType::kInsertionToRef || type == OperationType::kDeletionFromRef;
        const bool previousIsIndel
            = previousType == OperationType::kInsertionToRef || previousType == OperationType::kDeletionFromRef;
        if (type == previousType || (isIndel && previousIsIndel) || (isIndel && operations.empty()))
        {
            continue;
        }

        uint32_t length = lengthDist(rng);
        if (type != OperationType::kInsertionToRef)
        {
            length = std::min(length, referenceLength - used);
            used += length;
        }
        operations.emplace_back(type, length);
        previousType = type;
    }

    if (typeDist(rng) < 3 && previousType != OperationType::kInsertionToRef)
    {
        operations.emplace_back(OperationType::kSoftclip, lengthDist(rng));
    }

    return operations;
}
}

TEST(ProjectingAlignmentsOntoGraph, RandomAlignments_MatchOriginalImplementation)
{
    Graph graph = ehunter::makeRegionGraph(ehunter::decodeFeaturesFromRegex("ATTCGATCGT(CAG)*CAACAG(CCG)*ATGTCGTACC"));
    // Nodes: 0 left flank, 1 CAG, 2 interruption, 3 CCG, 4 right flank
    std::mt19937 rng(20261019);

    int numCompared = 0;
    for (int trial = 0; trial != 4000; ++trial)
    {
        std::uniform_int_distribution<int> repeatDist(0, 40);
        vector<NodeId> nodes = { 0 };
        for (int index = repeatDist(rng); index != 0; --index)
        {
            nodes.push_back(1);
        }
        nodes.push_back(2);
        for (int index = repeatDist(rng); index != 0; --index)
        {
            nodes.push_back(3);
        }
        nodes.push_back(4);

        std::uniform_int_distribution<int> startDist(0, 9);
        std::uniform_int_distribution<int> endDist(1, 10);
        Path path(&graph, startDist(rng), nodes, endDist(rng));

        std::uniform_int_distribution<uint32_t> offsetDist(0, std::min<uint32_t>(path.length() - 1, 5));
        const uint32_t referenceStart = offsetDist(rng);
        const uint32_t referenceLength = path.length() - referenceStart - offsetDist(rng) / 2;
        if (referenceLength == 0)
        {
            continue;
        }

        Alignment linearAlignment(referenceStart, randomOperations(referenceLength, rng));
        ASSERT_EQ(referenceLength, linearAlignment.referenceLength());

        const GraphAlignment expected = projectAlignmentOntoGraphReference(linearAlignment, path);
        const GraphAlignment observed = graphtools::projectAlignmentOntoGraph(linearAlignment, path);
        ASSERT_EQ(expected, observed) << "trial " << trial << ": " << linearAlignment;
        ++numCompared;
    }

    EXPECT_GT(numCompared, 3000);
}
