// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <iostream>

#include <gtest/gtest.h>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/range/concept.hpp>
#include <seqan3/range/view/to_rank.hpp>
#include <seqan3/std/ranges>

using namespace seqan3;

TEST(view_to_rank, basic)
{
    dna5_vector vec{"ACTTTGATA"_dna5};
    std::vector<uint8_t> cmp{0,1,4,4,4,2,0,4,0};

    // pipe notation
    std::vector<uint8_t> v = vec | view::to_rank | std::ranges::to<std::vector>;
    EXPECT_EQ(cmp, v);

    // function notation
    std::vector<uint8_t> v2(view::to_rank(vec) | std::ranges::to<std::vector>);
    EXPECT_EQ(cmp, v2);

    // combinability
    std::vector<uint8_t> cmp2{0, 4, 0, 2, 4, 4, 4, 1, 0};
    std::vector<uint8_t> v3 = vec | view::to_rank | std::view::reverse | std::ranges::to<std::vector>;
    EXPECT_EQ(cmp2, v3);
}

TEST(view_to_rank, concepts)
{
    dna5_vector vec{"ACTTTGATA"_dna5};
    EXPECT_TRUE(std::ranges::input_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::bidirectional_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::random_access_range<decltype(vec)>);
    EXPECT_FALSE(std::ranges::view<decltype(vec)>);
    EXPECT_TRUE(std::ranges::sized_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::common_range<decltype(vec)>);
    EXPECT_TRUE(const_iterable_range<decltype(vec)>);
    EXPECT_TRUE((std::ranges::output_range<decltype(vec), dna5>));

    auto v1 = vec | view::to_rank;
    EXPECT_TRUE(std::ranges::input_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::bidirectional_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::random_access_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::view<decltype(v1)>);
    EXPECT_TRUE(std::ranges::sized_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::common_range<decltype(v1)>);
    EXPECT_TRUE(const_iterable_range<decltype(v1)>);
    EXPECT_FALSE((std::ranges::output_range<decltype(v1), dna5>));
    EXPECT_FALSE((std::ranges::output_range<decltype(v1), uint8_t>));
}
