// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/io/variant_io/header.hpp>

#include <seqan3/test/expect_range_eq.hpp>

#include "../format/vcf_data.hpp"

// using namespace seqan3::literals;
using namespace std::literals;


TEST(var_io_header, spec_from_text)
{
    using svpair = std::pair<std::string_view const, std::string_view>;
    seqan3::var_io::header hdr{std::string{example_from_spec_header}};

    auto const & parsed = hdr.parsed_header();

    EXPECT_EQ(parsed.file_format,       "VCFv4.3");

    // contigs
    ASSERT_EQ(parsed.contigs.size(),                    1);
    EXPECT_EQ(parsed.contigs[0].id,                     "20"sv);
    EXPECT_EQ(parsed.contigs[0].length,                 62435964);
    ASSERT_EQ(parsed.contigs[0].other_fields.size(),    4);
    auto it = parsed.contigs[0].other_fields.begin();
    EXPECT_TRUE(*it   == (svpair{"assembly"sv, "B36"sv}));
    EXPECT_TRUE(*++it == (svpair{"md5"sv, "f126cdf8a6e0c7f379d618ff66beb2da"sv}));
    EXPECT_TRUE(*++it == (svpair{"species"sv, "\"Homo sapiens\""sv}));
    EXPECT_TRUE(*++it == (svpair{"taxonomy"sv, "x"sv}));

    // infos
    ASSERT_EQ(parsed.infos.size(),                      6);

    // info 0
    EXPECT_EQ(parsed.infos[0].id,                       "NS"sv);
    EXPECT_EQ(parsed.infos[0].number,                   1);
    EXPECT_EQ(parsed.infos[0].type,                     seqan3::io_type_id::int32);
    EXPECT_EQ(parsed.infos[0].description,              "\"Number of Samples With Data\""sv);
    EXPECT_EQ(parsed.infos[0].other_fields.size(),      0);

    // info 1
    EXPECT_EQ(parsed.infos[1].id,                       "DP"sv);
    EXPECT_EQ(parsed.infos[1].number,                   1);
    EXPECT_EQ(parsed.infos[1].type,                     seqan3::io_type_id::int32);
    EXPECT_EQ(parsed.infos[1].description,              "\"Total Depth\""sv);
    EXPECT_EQ(parsed.infos[1].other_fields.size(),      0);

    // info 2
    EXPECT_EQ(parsed.infos[2].id,                       "AF"sv);
    EXPECT_EQ(parsed.infos[2].number,                   seqan3::var_io::header_number::A);
    EXPECT_EQ(parsed.infos[2].type,                     seqan3::io_type_id::vector_of_float32);
    EXPECT_EQ(parsed.infos[2].description,              "\"Allele Frequency\""sv);
    EXPECT_EQ(parsed.infos[2].other_fields.size(),      0);

    // info 3
    EXPECT_EQ(parsed.infos[3].id,                       "AA"sv);
    EXPECT_EQ(parsed.infos[3].number,                   1);
    EXPECT_EQ(parsed.infos[3].type,                     seqan3::io_type_id::string);
    EXPECT_EQ(parsed.infos[3].description,              "\"Ancestral Allele\""sv);
    EXPECT_EQ(parsed.infos[3].other_fields.size(),      0);

    // info 4
    EXPECT_EQ(parsed.infos[4].id,                       "DB"sv);
    EXPECT_EQ(parsed.infos[4].number,                   0);
    EXPECT_EQ(parsed.infos[4].type,                     seqan3::io_type_id::flag);
    EXPECT_EQ(parsed.infos[4].description,              "\"dbSNP membership, build 129\""sv);
    EXPECT_EQ(parsed.infos[4].other_fields.size(),      0);

    // info 5
    EXPECT_EQ(parsed.infos[5].id,                       "H2"sv);
    EXPECT_EQ(parsed.infos[5].number,                   0);
    EXPECT_EQ(parsed.infos[5].type,                     seqan3::io_type_id::flag);
    EXPECT_EQ(parsed.infos[5].description,              "\"HapMap2 membership\""sv);
    EXPECT_EQ(parsed.infos[5].other_fields.size(),      0);

    //TODO rest

}
