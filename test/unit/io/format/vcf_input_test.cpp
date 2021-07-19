// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/io/format/format_vcf_input_handler.hpp>
#include <seqan3/io/variant_io/reader.hpp>

#include "vcf_input_data.hpp"

using namespace seqan3::literals;

template <bool shallow>
void main_test()
{
    std::istringstream istr{std::string{example_from_spec}};

    using record_t = std::conditional_t<shallow,
                                        seqan3::record<std::remove_cvref_t<decltype(seqan3::var_io::field_types_shallow)>,
                                                       std::remove_cvref_t<decltype(seqan3::var_io::default_field_ids)>>,
                                        seqan3::record<std::remove_cvref_t<decltype(seqan3::var_io::field_types_deep)>,
                                                       std::remove_cvref_t<decltype(seqan3::var_io::default_field_ids)>>>;
    using string_t = std::conditional_t<shallow, std::string_view, std::string>;

    seqan3::input_format_handler<seqan3::format_vcf> handler{istr, seqan3::var_io::reader_options{}};

    record_t rec;

    /* FIRST RECORD */
    handler.parse_next_record_into(rec);

    EXPECT_EQ(get<seqan3::field::chrom>(rec),               0);
    EXPECT_EQ(get<seqan3::field::pos>(rec),                 14370);
    EXPECT_EQ(get<seqan3::field::id>(rec),                  "rs6054257");
    EXPECT_RANGE_EQ(get<seqan3::field::ref>(rec),           "G"_dna5);

    EXPECT_EQ(get<seqan3::field::alt>(rec).size(),          1);
    ASSERT_EQ(get<seqan3::field::alt>(rec)[0],              "A");

    EXPECT_EQ(get<seqan3::field::qual>(rec),                29.0);
    ASSERT_EQ(get<seqan3::field::filter>(rec).size(),       1);
    EXPECT_EQ(get<seqan3::field::filter>(rec)[0],           0);

    ASSERT_EQ(get<seqan3::field::info>(rec).size(),         5ULL);
    EXPECT_EQ(get<seqan3::field::info>(rec)[0].first,       0);
    ASSERT_TRUE(std::holds_alternative<int32_t>(get<seqan3::field::info>(rec)[0].second));
    EXPECT_EQ(get<int32_t>(get<seqan3::field::info>(rec)[0].second),      3);

    EXPECT_EQ(get<seqan3::field::info>(rec)[1].first,       1);
    ASSERT_TRUE(std::holds_alternative<int32_t>(get<seqan3::field::info>(rec)[1].second));
    EXPECT_EQ(get<int32_t>(get<seqan3::field::info>(rec)[1].second),      14);

    EXPECT_EQ(get<seqan3::field::info>(rec)[2].first,       2);
    ASSERT_TRUE(std::holds_alternative<std::vector<float>>(get<seqan3::field::info>(rec)[2].second));
    EXPECT_EQ(get<std::vector<float>>(get<seqan3::field::info>(rec)[2].second).size(), 1);
    EXPECT_FLOAT_EQ(get<std::vector<float>>(get<seqan3::field::info>(rec)[2].second)[0], 0.5);

    EXPECT_EQ(get<seqan3::field::info>(rec)[3].first,       4);
    ASSERT_TRUE(std::holds_alternative<bool>(get<seqan3::field::info>(rec)[3].second));
    EXPECT_EQ(get<bool>(get<seqan3::field::info>(rec)[3].second), 1);

    EXPECT_EQ(get<seqan3::field::info>(rec)[4].first,       5);
    ASSERT_TRUE(std::holds_alternative<bool>(get<seqan3::field::info>(rec)[4].second));
    EXPECT_EQ(get<bool>(get<seqan3::field::info>(rec)[4].second), 1);


    ASSERT_EQ(get<seqan3::field::genotypes>(rec).size(),    4ULL);
    EXPECT_EQ(get<seqan3::field::genotypes>(rec)[0].first,  0);
    ASSERT_TRUE(std::holds_alternative<std::vector<string_t>>(get<seqan3::field::genotypes>(rec)[0].second));
    ASSERT_EQ(get<std::vector<string_t>>(get<seqan3::field::genotypes>(rec)[0].second).size(), 3);
    EXPECT_EQ(get<std::vector<string_t>>(get<seqan3::field::genotypes>(rec)[0].second)[0], "0|0");
    EXPECT_EQ(get<std::vector<string_t>>(get<seqan3::field::genotypes>(rec)[0].second)[1], "1|0");
    EXPECT_EQ(get<std::vector<string_t>>(get<seqan3::field::genotypes>(rec)[0].second)[2], "1/1");

    EXPECT_EQ(get<seqan3::field::genotypes>(rec)[1].first,  1);
    ASSERT_TRUE(std::holds_alternative<std::vector<int32_t>>(get<seqan3::field::genotypes>(rec)[1].second));
    ASSERT_EQ(get<std::vector<int32_t>>(get<seqan3::field::genotypes>(rec)[1].second).size(), 3);
    EXPECT_EQ(get<std::vector<int32_t>>(get<seqan3::field::genotypes>(rec)[1].second)[0], 48);
    EXPECT_EQ(get<std::vector<int32_t>>(get<seqan3::field::genotypes>(rec)[1].second)[1], 48);
    EXPECT_EQ(get<std::vector<int32_t>>(get<seqan3::field::genotypes>(rec)[1].second)[2], 43);

    EXPECT_EQ(get<seqan3::field::genotypes>(rec)[2].first,  2);
    ASSERT_TRUE(std::holds_alternative<std::vector<int32_t>>(get<seqan3::field::genotypes>(rec)[2].second));
    ASSERT_EQ(get<std::vector<int32_t>>(get<seqan3::field::genotypes>(rec)[2].second).size(), 3);
    EXPECT_EQ(get<std::vector<int32_t>>(get<seqan3::field::genotypes>(rec)[2].second)[0], 1);
    EXPECT_EQ(get<std::vector<int32_t>>(get<seqan3::field::genotypes>(rec)[2].second)[1], 8);
    EXPECT_EQ(get<std::vector<int32_t>>(get<seqan3::field::genotypes>(rec)[2].second)[2], 5);

    EXPECT_EQ(get<seqan3::field::genotypes>(rec)[3].first,  3);
    ASSERT_TRUE(std::holds_alternative<std::vector<std::vector<int32_t>>>(get<seqan3::field::genotypes>(rec)[3].second));
    ASSERT_EQ(get<std::vector<std::vector<int32_t>>>(get<seqan3::field::genotypes>(rec)[3].second).size(), 3);
    ASSERT_EQ(get<std::vector<std::vector<int32_t>>>(get<seqan3::field::genotypes>(rec)[3].second)[0].size(), 2);
    EXPECT_EQ(get<std::vector<std::vector<int32_t>>>(get<seqan3::field::genotypes>(rec)[3].second)[0][0], 51);
    EXPECT_EQ(get<std::vector<std::vector<int32_t>>>(get<seqan3::field::genotypes>(rec)[3].second)[0][1], 51);
    ASSERT_EQ(get<std::vector<std::vector<int32_t>>>(get<seqan3::field::genotypes>(rec)[3].second)[1].size(), 2);
    EXPECT_EQ(get<std::vector<std::vector<int32_t>>>(get<seqan3::field::genotypes>(rec)[3].second)[1][0], 51);
    EXPECT_EQ(get<std::vector<std::vector<int32_t>>>(get<seqan3::field::genotypes>(rec)[3].second)[1][1], 51);
    ASSERT_EQ(get<std::vector<std::vector<int32_t>>>(get<seqan3::field::genotypes>(rec)[3].second)[2].size(), 2);
    EXPECT_EQ(get<std::vector<std::vector<int32_t>>>(get<seqan3::field::genotypes>(rec)[3].second)[2][0],
              seqan3::var_io::missing_value<int32_t>);
    EXPECT_EQ(get<std::vector<std::vector<int32_t>>>(get<seqan3::field::genotypes>(rec)[3].second)[2][1],
              seqan3::var_io::missing_value<int32_t>);

    EXPECT_EQ(get<seqan3::field::_private>(rec).header_ptr->raw_header(), example_from_spec_header);

    /* SECOND RECORD */
    handler.parse_next_record_into(rec);

    EXPECT_EQ(get<seqan3::field::chrom>(rec),               0);
    EXPECT_EQ(get<seqan3::field::pos>(rec),                 17330);
    EXPECT_EQ(get<seqan3::field::id>(rec),                  ".");
    EXPECT_RANGE_EQ(get<seqan3::field::ref>(rec),           "T"_dna5);

    EXPECT_EQ(get<seqan3::field::alt>(rec).size(),          1);
    ASSERT_EQ(get<seqan3::field::alt>(rec)[0],              "A");

    EXPECT_EQ(get<seqan3::field::qual>(rec),                3.0);
    ASSERT_EQ(get<seqan3::field::filter>(rec).size(),       1);
    EXPECT_EQ(get<seqan3::field::filter>(rec)[0],           1);

    ASSERT_EQ(get<seqan3::field::info>(rec).size(),         3ULL);
    EXPECT_EQ(get<seqan3::field::info>(rec)[0].first,       0);
    ASSERT_TRUE(std::holds_alternative<int32_t>(get<seqan3::field::info>(rec)[0].second));
    EXPECT_EQ(get<int32_t>(get<seqan3::field::info>(rec)[0].second),      3);

    EXPECT_EQ(get<seqan3::field::info>(rec)[1].first,       1);
    ASSERT_TRUE(std::holds_alternative<int32_t>(get<seqan3::field::info>(rec)[1].second));
    EXPECT_EQ(get<int32_t>(get<seqan3::field::info>(rec)[1].second),      11);

    EXPECT_EQ(get<seqan3::field::info>(rec)[2].first,       2);
    ASSERT_TRUE(std::holds_alternative<std::vector<float>>(get<seqan3::field::info>(rec)[2].second));
    EXPECT_EQ(get<std::vector<float>>(get<seqan3::field::info>(rec)[2].second).size(), 1);
    EXPECT_FLOAT_EQ(get<std::vector<float>>(get<seqan3::field::info>(rec)[2].second)[0], 0.017);

    ASSERT_EQ(get<seqan3::field::genotypes>(rec).size(),    4ULL);
    EXPECT_EQ(get<seqan3::field::genotypes>(rec)[0].first,  0);
    ASSERT_TRUE(std::holds_alternative<std::vector<string_t>>(get<seqan3::field::genotypes>(rec)[0].second));
    ASSERT_EQ(get<std::vector<string_t>>(get<seqan3::field::genotypes>(rec)[0].second).size(), 3);
    EXPECT_EQ(get<std::vector<string_t>>(get<seqan3::field::genotypes>(rec)[0].second)[0], "0|0");
    EXPECT_EQ(get<std::vector<string_t>>(get<seqan3::field::genotypes>(rec)[0].second)[1], "0|1");
    EXPECT_EQ(get<std::vector<string_t>>(get<seqan3::field::genotypes>(rec)[0].second)[2], "0/0");

    EXPECT_EQ(get<seqan3::field::genotypes>(rec)[1].first,  1);
    ASSERT_TRUE(std::holds_alternative<std::vector<int32_t>>(get<seqan3::field::genotypes>(rec)[1].second));
    ASSERT_EQ(get<std::vector<int32_t>>(get<seqan3::field::genotypes>(rec)[1].second).size(), 3);
    EXPECT_EQ(get<std::vector<int32_t>>(get<seqan3::field::genotypes>(rec)[1].second)[0], 49);
    EXPECT_EQ(get<std::vector<int32_t>>(get<seqan3::field::genotypes>(rec)[1].second)[1], 3);
    EXPECT_EQ(get<std::vector<int32_t>>(get<seqan3::field::genotypes>(rec)[1].second)[2], 41);

    EXPECT_EQ(get<seqan3::field::genotypes>(rec)[2].first,  2);
    ASSERT_TRUE(std::holds_alternative<std::vector<int32_t>>(get<seqan3::field::genotypes>(rec)[2].second));
    ASSERT_EQ(get<std::vector<int32_t>>(get<seqan3::field::genotypes>(rec)[2].second).size(), 3);
    EXPECT_EQ(get<std::vector<int32_t>>(get<seqan3::field::genotypes>(rec)[2].second)[0], 3);
    EXPECT_EQ(get<std::vector<int32_t>>(get<seqan3::field::genotypes>(rec)[2].second)[1], 5);
    EXPECT_EQ(get<std::vector<int32_t>>(get<seqan3::field::genotypes>(rec)[2].second)[2], 3);

    EXPECT_EQ(get<seqan3::field::genotypes>(rec)[3].first,  3);
    ASSERT_TRUE(std::holds_alternative<std::vector<std::vector<int32_t>>>(get<seqan3::field::genotypes>(rec)[3].second));
    ASSERT_EQ(get<std::vector<std::vector<int32_t>>>(get<seqan3::field::genotypes>(rec)[3].second).size(), 3);
    ASSERT_EQ(get<std::vector<std::vector<int32_t>>>(get<seqan3::field::genotypes>(rec)[3].second)[0].size(), 2);
    EXPECT_EQ(get<std::vector<std::vector<int32_t>>>(get<seqan3::field::genotypes>(rec)[3].second)[0][0], 58);
    EXPECT_EQ(get<std::vector<std::vector<int32_t>>>(get<seqan3::field::genotypes>(rec)[3].second)[0][1], 50);
    ASSERT_EQ(get<std::vector<std::vector<int32_t>>>(get<seqan3::field::genotypes>(rec)[3].second)[1].size(), 2);
    EXPECT_EQ(get<std::vector<std::vector<int32_t>>>(get<seqan3::field::genotypes>(rec)[3].second)[1][0], 65);
    EXPECT_EQ(get<std::vector<std::vector<int32_t>>>(get<seqan3::field::genotypes>(rec)[3].second)[1][1], 3);
    ASSERT_EQ(get<std::vector<std::vector<int32_t>>>(get<seqan3::field::genotypes>(rec)[3].second)[2].size(), 0);

    EXPECT_EQ(get<seqan3::field::_private>(rec).header_ptr->raw_header(), example_from_spec_header);

    /* THIRD RECORD */
    handler.parse_next_record_into(rec);

    EXPECT_EQ(get<seqan3::field::chrom>(rec),               0);
    EXPECT_EQ(get<seqan3::field::pos>(rec),                 1110696);
    EXPECT_EQ(get<seqan3::field::id>(rec),                  "rs6040355");
    EXPECT_RANGE_EQ(get<seqan3::field::ref>(rec),           "A"_dna5);

    EXPECT_EQ(get<seqan3::field::alt>(rec).size(),          2);
    ASSERT_EQ(get<seqan3::field::alt>(rec)[0],              "G");
    ASSERT_EQ(get<seqan3::field::alt>(rec)[1],              "T");

    EXPECT_FLOAT_EQ(get<seqan3::field::qual>(rec),          67.0);

    ASSERT_EQ(get<seqan3::field::filter>(rec).size(),       1);
    EXPECT_EQ(get<seqan3::field::filter>(rec)[0],           0);

    ASSERT_EQ(get<seqan3::field::info>(rec).size(),         5ULL);
    EXPECT_EQ(get<seqan3::field::info>(rec)[0].first,       0);
    ASSERT_TRUE(std::holds_alternative<int32_t>(get<seqan3::field::info>(rec)[0].second));
    EXPECT_EQ(get<int32_t>(get<seqan3::field::info>(rec)[0].second),      2);

    EXPECT_EQ(get<seqan3::field::info>(rec)[1].first,       1);
    ASSERT_TRUE(std::holds_alternative<int32_t>(get<seqan3::field::info>(rec)[1].second));
    EXPECT_EQ(get<int32_t>(get<seqan3::field::info>(rec)[1].second),      10);

    EXPECT_EQ(get<seqan3::field::info>(rec)[2].first,       2);
    ASSERT_TRUE(std::holds_alternative<std::vector<float>>(get<seqan3::field::info>(rec)[2].second));
    EXPECT_EQ(get<std::vector<float>>(get<seqan3::field::info>(rec)[2].second).size(), 2);
    EXPECT_FLOAT_EQ(get<std::vector<float>>(get<seqan3::field::info>(rec)[2].second)[0], 0.333);
    EXPECT_FLOAT_EQ(get<std::vector<float>>(get<seqan3::field::info>(rec)[2].second)[1], 0.667);

    EXPECT_EQ(get<seqan3::field::info>(rec)[3].first,       3);
    ASSERT_TRUE(std::holds_alternative<string_t>(get<seqan3::field::info>(rec)[3].second));
    EXPECT_EQ(get<string_t>(get<seqan3::field::info>(rec)[3].second), "T");

    EXPECT_EQ(get<seqan3::field::info>(rec)[4].first,       4);
    ASSERT_TRUE(std::holds_alternative<bool>(get<seqan3::field::info>(rec)[4].second));
    EXPECT_EQ(get<bool>(get<seqan3::field::info>(rec)[4].second), 1);


    ASSERT_EQ(get<seqan3::field::genotypes>(rec).size(),    4ULL);
    EXPECT_EQ(get<seqan3::field::genotypes>(rec)[0].first,  0);
    ASSERT_TRUE(std::holds_alternative<std::vector<string_t>>(get<seqan3::field::genotypes>(rec)[0].second));
    ASSERT_EQ(get<std::vector<string_t>>(get<seqan3::field::genotypes>(rec)[0].second).size(), 3);
    EXPECT_EQ(get<std::vector<string_t>>(get<seqan3::field::genotypes>(rec)[0].second)[0], "1|2");
    EXPECT_EQ(get<std::vector<string_t>>(get<seqan3::field::genotypes>(rec)[0].second)[1], "2|1");
    EXPECT_EQ(get<std::vector<string_t>>(get<seqan3::field::genotypes>(rec)[0].second)[2], "2/2");

    EXPECT_EQ(get<seqan3::field::genotypes>(rec)[1].first,  1);
    ASSERT_TRUE(std::holds_alternative<std::vector<int32_t>>(get<seqan3::field::genotypes>(rec)[1].second));
    ASSERT_EQ(get<std::vector<int32_t>>(get<seqan3::field::genotypes>(rec)[1].second).size(), 3);
    EXPECT_EQ(get<std::vector<int32_t>>(get<seqan3::field::genotypes>(rec)[1].second)[0], 21);
    EXPECT_EQ(get<std::vector<int32_t>>(get<seqan3::field::genotypes>(rec)[1].second)[1], 2);
    EXPECT_EQ(get<std::vector<int32_t>>(get<seqan3::field::genotypes>(rec)[1].second)[2], 35);

    EXPECT_EQ(get<seqan3::field::genotypes>(rec)[2].first,  2);
    ASSERT_TRUE(std::holds_alternative<std::vector<int32_t>>(get<seqan3::field::genotypes>(rec)[2].second));
    ASSERT_EQ(get<std::vector<int32_t>>(get<seqan3::field::genotypes>(rec)[2].second).size(), 3);
    EXPECT_EQ(get<std::vector<int32_t>>(get<seqan3::field::genotypes>(rec)[2].second)[0], 6);
    EXPECT_EQ(get<std::vector<int32_t>>(get<seqan3::field::genotypes>(rec)[2].second)[1], 0);
    EXPECT_EQ(get<std::vector<int32_t>>(get<seqan3::field::genotypes>(rec)[2].second)[2], 4);

    EXPECT_EQ(get<seqan3::field::genotypes>(rec)[3].first,  3);
    ASSERT_TRUE(std::holds_alternative<std::vector<std::vector<int32_t>>>(get<seqan3::field::genotypes>(rec)[3].second));
    ASSERT_EQ(get<std::vector<std::vector<int32_t>>>(get<seqan3::field::genotypes>(rec)[3].second).size(), 3);
    ASSERT_EQ(get<std::vector<std::vector<int32_t>>>(get<seqan3::field::genotypes>(rec)[3].second)[0].size(), 2);
    EXPECT_EQ(get<std::vector<std::vector<int32_t>>>(get<seqan3::field::genotypes>(rec)[3].second)[0][0], 23);
    EXPECT_EQ(get<std::vector<std::vector<int32_t>>>(get<seqan3::field::genotypes>(rec)[3].second)[0][1], 27);
    ASSERT_EQ(get<std::vector<std::vector<int32_t>>>(get<seqan3::field::genotypes>(rec)[3].second)[1].size(), 2);
    EXPECT_EQ(get<std::vector<std::vector<int32_t>>>(get<seqan3::field::genotypes>(rec)[3].second)[1][0], 18);
    EXPECT_EQ(get<std::vector<std::vector<int32_t>>>(get<seqan3::field::genotypes>(rec)[3].second)[1][1], 2);
    ASSERT_EQ(get<std::vector<std::vector<int32_t>>>(get<seqan3::field::genotypes>(rec)[3].second)[2].size(), 0);

    EXPECT_EQ(get<seqan3::field::_private>(rec).header_ptr->raw_header(), example_from_spec_header);

    /* FOURTH RECORD */
    handler.parse_next_record_into(rec);

    EXPECT_EQ(get<seqan3::field::chrom>(rec),               0);
    EXPECT_EQ(get<seqan3::field::pos>(rec),                 1230237);
    EXPECT_EQ(get<seqan3::field::id>(rec),                  ".");
    EXPECT_RANGE_EQ(get<seqan3::field::ref>(rec),           "T"_dna5);

    EXPECT_EQ(get<seqan3::field::alt>(rec).size(),          0);

    EXPECT_FLOAT_EQ(get<seqan3::field::qual>(rec),          47.0);

    ASSERT_EQ(get<seqan3::field::filter>(rec).size(),       1);
    EXPECT_EQ(get<seqan3::field::filter>(rec)[0],           0);

    ASSERT_EQ(get<seqan3::field::info>(rec).size(),         3ULL);
    EXPECT_EQ(get<seqan3::field::info>(rec)[0].first,       0);
    ASSERT_TRUE(std::holds_alternative<int32_t>(get<seqan3::field::info>(rec)[0].second));
    EXPECT_EQ(get<int32_t>(get<seqan3::field::info>(rec)[0].second),      3);

    EXPECT_EQ(get<seqan3::field::info>(rec)[1].first,       1);
    ASSERT_TRUE(std::holds_alternative<int32_t>(get<seqan3::field::info>(rec)[1].second));
    EXPECT_EQ(get<int32_t>(get<seqan3::field::info>(rec)[1].second),      13);

    EXPECT_EQ(get<seqan3::field::info>(rec)[2].first,       3);
    ASSERT_TRUE(std::holds_alternative<string_t>(get<seqan3::field::info>(rec)[2].second));
    EXPECT_EQ(get<string_t>(get<seqan3::field::info>(rec)[2].second), "T");

    ASSERT_EQ(get<seqan3::field::genotypes>(rec).size(),    4ULL);
    EXPECT_EQ(get<seqan3::field::genotypes>(rec)[0].first,  0);
    ASSERT_TRUE(std::holds_alternative<std::vector<string_t>>(get<seqan3::field::genotypes>(rec)[0].second));
    ASSERT_EQ(get<std::vector<string_t>>(get<seqan3::field::genotypes>(rec)[0].second).size(), 3);
    EXPECT_EQ(get<std::vector<string_t>>(get<seqan3::field::genotypes>(rec)[0].second)[0], "0|0");
    EXPECT_EQ(get<std::vector<string_t>>(get<seqan3::field::genotypes>(rec)[0].second)[1], "0|0");
    EXPECT_EQ(get<std::vector<string_t>>(get<seqan3::field::genotypes>(rec)[0].second)[2], "0/0");

    EXPECT_EQ(get<seqan3::field::genotypes>(rec)[1].first,  1);
    ASSERT_TRUE(std::holds_alternative<std::vector<int32_t>>(get<seqan3::field::genotypes>(rec)[1].second));
    ASSERT_EQ(get<std::vector<int32_t>>(get<seqan3::field::genotypes>(rec)[1].second).size(), 3);
    EXPECT_EQ(get<std::vector<int32_t>>(get<seqan3::field::genotypes>(rec)[1].second)[0], 54);
    EXPECT_EQ(get<std::vector<int32_t>>(get<seqan3::field::genotypes>(rec)[1].second)[1], 48);
    EXPECT_EQ(get<std::vector<int32_t>>(get<seqan3::field::genotypes>(rec)[1].second)[2], 61);

    EXPECT_EQ(get<seqan3::field::genotypes>(rec)[2].first,  2);
    ASSERT_TRUE(std::holds_alternative<std::vector<int32_t>>(get<seqan3::field::genotypes>(rec)[2].second));
    ASSERT_EQ(get<std::vector<int32_t>>(get<seqan3::field::genotypes>(rec)[2].second).size(), 3);
    EXPECT_EQ(get<std::vector<int32_t>>(get<seqan3::field::genotypes>(rec)[2].second)[0], 7);
    EXPECT_EQ(get<std::vector<int32_t>>(get<seqan3::field::genotypes>(rec)[2].second)[1], 4);
    EXPECT_EQ(get<std::vector<int32_t>>(get<seqan3::field::genotypes>(rec)[2].second)[2], 2);

    EXPECT_EQ(get<seqan3::field::genotypes>(rec)[3].first,  3);
    ASSERT_TRUE(std::holds_alternative<std::vector<std::vector<int32_t>>>(get<seqan3::field::genotypes>(rec)[3].second));
    ASSERT_EQ(get<std::vector<std::vector<int32_t>>>(get<seqan3::field::genotypes>(rec)[3].second).size(), 3);
    ASSERT_EQ(get<std::vector<std::vector<int32_t>>>(get<seqan3::field::genotypes>(rec)[3].second)[0].size(), 2);
    EXPECT_EQ(get<std::vector<std::vector<int32_t>>>(get<seqan3::field::genotypes>(rec)[3].second)[0][0], 56);
    EXPECT_EQ(get<std::vector<std::vector<int32_t>>>(get<seqan3::field::genotypes>(rec)[3].second)[0][1], 60);
    ASSERT_EQ(get<std::vector<std::vector<int32_t>>>(get<seqan3::field::genotypes>(rec)[3].second)[1].size(), 2);
    EXPECT_EQ(get<std::vector<std::vector<int32_t>>>(get<seqan3::field::genotypes>(rec)[3].second)[1][0], 51);
    EXPECT_EQ(get<std::vector<std::vector<int32_t>>>(get<seqan3::field::genotypes>(rec)[3].second)[1][1], 51);
    ASSERT_EQ(get<std::vector<std::vector<int32_t>>>(get<seqan3::field::genotypes>(rec)[3].second)[2].size(), 0);

    EXPECT_EQ(get<seqan3::field::_private>(rec).header_ptr->raw_header(), example_from_spec_header);

    /* FIFTH RECORD */
    handler.parse_next_record_into(rec);

    EXPECT_EQ(get<seqan3::field::chrom>(rec),               0);
    EXPECT_EQ(get<seqan3::field::pos>(rec),                 1234567);
    EXPECT_EQ(get<seqan3::field::id>(rec),                  "microsat1");
    EXPECT_RANGE_EQ(get<seqan3::field::ref>(rec),           "GTC"_dna5);

    ASSERT_EQ(get<seqan3::field::alt>(rec).size(),          2);
    EXPECT_EQ(get<seqan3::field::alt>(rec)[0],              "G");
    EXPECT_EQ(get<seqan3::field::alt>(rec)[1],              "GTCT");

    EXPECT_FLOAT_EQ(get<seqan3::field::qual>(rec),          50.0);

    ASSERT_EQ(get<seqan3::field::filter>(rec).size(),       1);
    EXPECT_EQ(get<seqan3::field::filter>(rec)[0],           0);

    ASSERT_EQ(get<seqan3::field::info>(rec).size(),         3ULL);
    EXPECT_EQ(get<seqan3::field::info>(rec)[0].first,       0);
    ASSERT_TRUE(std::holds_alternative<int32_t>(get<seqan3::field::info>(rec)[0].second));
    EXPECT_EQ(get<int32_t>(get<seqan3::field::info>(rec)[0].second),      3);

    EXPECT_EQ(get<seqan3::field::info>(rec)[1].first,       1);
    ASSERT_TRUE(std::holds_alternative<int32_t>(get<seqan3::field::info>(rec)[1].second));
    EXPECT_EQ(get<int32_t>(get<seqan3::field::info>(rec)[1].second),      9);

    EXPECT_EQ(get<seqan3::field::info>(rec)[2].first,       3);
    ASSERT_TRUE(std::holds_alternative<string_t>(get<seqan3::field::info>(rec)[2].second));
    EXPECT_EQ(get<string_t>(get<seqan3::field::info>(rec)[2].second), "G");

    ASSERT_EQ(get<seqan3::field::genotypes>(rec).size(),    3ULL);
    EXPECT_EQ(get<seqan3::field::genotypes>(rec)[0].first,  0);
    ASSERT_TRUE(std::holds_alternative<std::vector<string_t>>(get<seqan3::field::genotypes>(rec)[0].second));
    ASSERT_EQ(get<std::vector<string_t>>(get<seqan3::field::genotypes>(rec)[0].second).size(), 3);
    EXPECT_EQ(get<std::vector<string_t>>(get<seqan3::field::genotypes>(rec)[0].second)[0], "0/1");
    EXPECT_EQ(get<std::vector<string_t>>(get<seqan3::field::genotypes>(rec)[0].second)[1], "0/2");
    EXPECT_EQ(get<std::vector<string_t>>(get<seqan3::field::genotypes>(rec)[0].second)[2], "1/1");

    EXPECT_EQ(get<seqan3::field::genotypes>(rec)[1].first,  1);
    ASSERT_TRUE(std::holds_alternative<std::vector<int32_t>>(get<seqan3::field::genotypes>(rec)[1].second));
    ASSERT_EQ(get<std::vector<int32_t>>(get<seqan3::field::genotypes>(rec)[1].second).size(), 3);
    EXPECT_EQ(get<std::vector<int32_t>>(get<seqan3::field::genotypes>(rec)[1].second)[0], 35);
    EXPECT_EQ(get<std::vector<int32_t>>(get<seqan3::field::genotypes>(rec)[1].second)[1], 17);
    EXPECT_EQ(get<std::vector<int32_t>>(get<seqan3::field::genotypes>(rec)[1].second)[2], 40);

    EXPECT_EQ(get<seqan3::field::genotypes>(rec)[2].first,  2);
    ASSERT_TRUE(std::holds_alternative<std::vector<int32_t>>(get<seqan3::field::genotypes>(rec)[2].second));
    ASSERT_EQ(get<std::vector<int32_t>>(get<seqan3::field::genotypes>(rec)[2].second).size(), 3);
    EXPECT_EQ(get<std::vector<int32_t>>(get<seqan3::field::genotypes>(rec)[2].second)[0], 4);
    EXPECT_EQ(get<std::vector<int32_t>>(get<seqan3::field::genotypes>(rec)[2].second)[1], 2);
    EXPECT_EQ(get<std::vector<int32_t>>(get<seqan3::field::genotypes>(rec)[2].second)[2], 3);

    EXPECT_EQ(get<seqan3::field::_private>(rec).header_ptr->raw_header(), example_from_spec_header);

    /* END */

    EXPECT_THROW(handler.parse_next_record_into(rec), std::runtime_error);
}

TEST(vcf, main_shallow)
{
    main_test<true>();
}

TEST(vcf, main_deep)
{
    main_test<false>();
}
