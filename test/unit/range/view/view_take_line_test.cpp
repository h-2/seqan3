// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <iostream>

#include <gtest/gtest.h>

#include <range/v3/algorithm/copy.hpp>
#include <range/v3/view/unique.hpp>

#include <seqan3/range/shortcuts.hpp>
#include <seqan3/range/view/take_line.hpp>
#include <seqan3/range/view/single_pass_input.hpp>
#include <seqan3/std/ranges>
#include <seqan3/range/view/to_char.hpp>

using namespace seqan3;

// ============================================================================
//  test templates
// ============================================================================

template <typename adaptor_t>
void do_test(adaptor_t const & adaptor, std::string const & vec)
{
    // pipe notation
    auto v = vec | adaptor;
    EXPECT_EQ("foo", v | std::ranges::to<std::string>);

    // function notation
    std::string v2(adaptor(vec) | std::ranges::to<std::string>);
    EXPECT_EQ("foo", v2);

    // combinability
    auto v3 = vec | adaptor | ranges::view::unique;
    EXPECT_EQ("fo", v3 | std::ranges::to<std::string>);
    std::string v3b = vec | std::view::reverse | adaptor | ranges::view::unique | std::ranges::to<std::string>;
    EXPECT_EQ("rab", v3b);

    // consuming behaviour
    auto v4 = vec | view::single_pass_input;
    auto v5 = std::move(v4) | adaptor;
    EXPECT_EQ("foo", v5 | std::ranges::to<std::string>);
    EXPECT_EQ('b', *begin(v5)); // not newline
}

template <typename adaptor_t>
void do_concepts(adaptor_t const & adaptor)
{
    std::string vec{"foo\nbar"};
    EXPECT_TRUE(std::ranges::input_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::bidirectional_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::random_access_range<decltype(vec)>);
    EXPECT_FALSE(std::ranges::view<decltype(vec)>);
    EXPECT_TRUE(std::ranges::sized_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::common_range<decltype(vec)>);
    EXPECT_TRUE(const_iterable_range<decltype(vec)>);
    EXPECT_TRUE((std::ranges::output_range<decltype(vec), char>));

    auto v1 = vec | adaptor;

    EXPECT_TRUE(std::ranges::input_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::bidirectional_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::random_access_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::view<decltype(v1)>);
    EXPECT_FALSE(std::ranges::sized_range<decltype(v1)>);
    EXPECT_FALSE(std::ranges::common_range<decltype(v1)>);
    EXPECT_TRUE(const_iterable_range<decltype(v1)>);
    EXPECT_TRUE((std::ranges::output_range<decltype(v1), char>));

    auto v2 = vec | view::single_pass_input | adaptor;

    EXPECT_TRUE(std::ranges::input_range<decltype(v2)>);
    EXPECT_FALSE(std::ranges::forward_range<decltype(v2)>);
    EXPECT_FALSE(std::ranges::bidirectional_range<decltype(v2)>);
    EXPECT_FALSE(std::ranges::random_access_range<decltype(v2)>);
    EXPECT_TRUE(std::ranges::view<decltype(v2)>);
    EXPECT_FALSE(std::ranges::sized_range<decltype(v2)>);
    EXPECT_FALSE(std::ranges::common_range<decltype(v2)>);
    EXPECT_FALSE(const_iterable_range<decltype(v2)>);
    EXPECT_TRUE((std::ranges::output_range<decltype(v2), char>));
}

// ============================================================================
//  view_take_line
// ============================================================================

TEST(view_take_line, unix_eol)
{
    do_test(view::take_line, "foo\nbar");
}

TEST(view_take_line, windows_eol)
{
    do_test(view::take_line, "foo\r\nbar");
}

TEST(view_take_line, no_eol)
{
    std::string vec{"foo"};
    std::string v;
    EXPECT_NO_THROW(( v = vec | view::take_line ));
    EXPECT_EQ("foo", v);
}

TEST(view_take_line, eol_at_first_position)
{
    using sbt = std::istreambuf_iterator<char>;
    std::istringstream vec{"\n\nfoo"};
    auto stream_view = std::ranges::subrange<decltype(sbt{vec}), decltype(sbt{})> {sbt{vec}, sbt{}};
    EXPECT_EQ("", stream_view | view::take_line | std::ranges::to<std::string>);
    EXPECT_EQ("foo", stream_view | view::take_line | std::ranges::to<std::string>);
}

TEST(view_take_line, concepts)
{
    do_concepts(view::take_line);
}

// ============================================================================
//  view_take_line_or_throw
// ============================================================================

TEST(view_take_line_or_throw, unix_eol)
{
    do_test(view::take_line_or_throw, "foo\nbar");
}

TEST(view_take_line_or_throw, windows_eol)
{
    do_test(view::take_line_or_throw, "foo\r\nbar");
}

TEST(view_take_line_or_throw, no_eol)
{
    std::string vec{"foo"};
    EXPECT_THROW(std::string v = vec | view::take_line_or_throw,
                 unexpected_end_of_input);
}

TEST(view_take_line_or_throw, concepts)
{
    do_concepts(view::take_line_or_throw);
}

// ============================================================================
//  bug
// ============================================================================

TEST(view_take_line, reverse_bug)
{
    std::string vec{"foo\nbar"};
    auto v1 = vec | view::take_line;
    EXPECT_EQ("foo", std::string(v1));
    EXPECT_TRUE(std::ranges::input_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::bidirectional_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::random_access_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::view<decltype(v1)>);
    EXPECT_FALSE(std::ranges::sized_range<decltype(v1)>);
    EXPECT_FALSE(std::ranges::common_range<decltype(v1)>);
    EXPECT_TRUE(const_iterable_range<decltype(v1)>);
    EXPECT_TRUE((std::ranges::output_range<decltype(v1), char>));

    // No build failure, but wrong results:
//     auto v2 = v1 | std::view::reverse;
//     EXPECT_EQ("oof", std::string(v2));
}
