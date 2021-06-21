// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * brief Provides the seqan3::input_format_handler<seqan3::format_vcf> .
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <iostream>
#include <regex>
#include <string>
#include <string_view>
#include <vector>

#include <seqan3/alphabet/views/char_strictly_to.hpp>
#include <seqan3/alphabet/views/char_to.hpp>
#include <seqan3/core/range/type_traits.hpp>
#include <seqan3/io/detail/ignore_output_iterator.hpp>
#include <seqan3/io/detail/misc.hpp>
#include <seqan3/io/format/input_format_handler_base.hpp>
#include <seqan3/io/format/format_vcf.hpp>
#include <seqan3/io/plaintext_io/reader.hpp>
#include <seqan3/io/stream/iterator.hpp>
#include <seqan3/io/variant_io/misc.hpp>
#include <seqan3/range/detail/misc.hpp>
#include <seqan3/std/algorithm>
#include <seqan3/std/ranges>
#include <seqan3/utility/type_list/traits.hpp>
#include <seqan3/utility/views/to.hpp>

namespace seqan3
{

template <>
class input_format_handler<format_vcf> : public input_format_handler_base<input_format_handler<format_vcf>>
{
private:
    /* CRTP */
    using base_t = input_format_handler_base<input_format_handler<format_vcf>>;
    using base_t::stream;

    friend base_t;

    /* RAW RECORD HANDLING*/
    using format_fields     = decltype(var_io::default_field_ids);
    using raw_record_type   = seqan3::record<list_traits::repeat<format_fields::size, std::string_view>,
                                             format_fields>;

    raw_record_type             raw_record;
    plain_io::detail::plaintext_input_iterator<char, std::char_traits<char>, plain_io::record_kind::line_and_fields> file_it;

    void read_raw_record()
    {
        ++file_it;

        get<field::chrom> (raw_record) = (*file_it).fields[0];
        get<field::pos>   (raw_record) = (*file_it).fields[1];
        get<field::id>    (raw_record) = (*file_it).fields[2];
        get<field::ref>   (raw_record) = (*file_it).fields[3];
        get<field::alt>   (raw_record) = (*file_it).fields[4];
        get<field::qual>  (raw_record) = (*file_it).fields[5];
        get<field::filter>(raw_record) = (*file_it).fields[6];
        get<field::info>  (raw_record) = (*file_it).fields[7];

        // fields[7].end() that is guaranteed to be char*
        char const * end_qual   = (*file_it).fields[7].data() + (*file_it).fields[7].size();
        // line.end() that is guaranteed to be char*
        char const * end_line   = (*file_it).line.data()      + (*file_it).line.size();
        // genotypes fo from end of qual til end of line (possibly empty)
        get<field::genotypes>  (raw_record) = std::string_view{end_qual, static_cast<size_t>(end_line - end_qual)};

        // header does not need to be reset
    }

    /* PARSED RECORD HANDLING */
    var_io::header header;

    // override the general case to handle MISSING value
//     template <field field_id, typename parsed_field_t>
//     void parse_field(tag_t<field_id> const & /**/, parsed_field_t & parsed_field) const
//     {
//         if (get<field_id>(raw_record) != var_io::special_value::missing)
//             static_cast<base_t const &>(*this).parse_field(tag<field_id>, parsed_field);
//     }

    using base_t::parse_field;

    void parse_field(tag_t<field::ref> const & /**/, var_io::allele & parsed_field) const
    {
        std::string_view raw_field = get<field::alt>(raw_record);
        parsed_field = parse_alt_impl(raw_field);
    }

    // TODO: make this back_insertable instead vector
    void parse_field(tag_t<field::alt> const & /**/, std::vector<var_io::allele> & parsed_field) const
    {
        std::string_view raw_field = get<field::alt>(raw_record);

        size_t last_begin = 0;
        for (size_t i = 0; i <= raw_field.size(); ++i)
        {
            if (i == raw_field.size() || raw_field[i] == ',')
            {
                parsed_field.push_back(parse_alt_impl(raw_field.substr(last_begin, i)));
                last_begin = i + 1;
            }
        }
    }

    template <typename parsed_field_t>
    void parse_field(tag_t<field::genotypes> const & /**/, parsed_field_t & parsed_field) const
    {
        if ((*file_it).fields.size() > 10) // there are optional fields
        {
            // TODO store into dictionary
        }
    }

    void parse_field(tag_t<field::header> const & /**/, var_io::header const * & parsed_field) const
    {
        parsed_field = &header;
    }


    /* TODO move to mixin: */

    static var_io::allele parse_alt_impl(std::string_view in)
    {
        std::regex is_seq{"^[ACGTNacgtn]+$"};

        if (in == "*")
            return var_io::special_value::unknown;
        else if (in == ".")
            return var_io::special_value::missing;
        else if (std::regex_match(in.begin(), in.end(), is_seq))
            return in | views::char_to<dna5> | views::to<std::vector<dna5>>;
        else
            return std::string{in};
    }

public:

    input_format_handler()                                              = default;
    input_format_handler(input_format_handler const &)                  = delete;
    input_format_handler(input_format_handler && )                      = default;
    input_format_handler & operator=(input_format_handler const &)      = delete;
    input_format_handler & operator=(input_format_handler &&)           = default;

    template <typename config_t>
    input_format_handler(std::istream & str, config_t const & /*cfg*/) :
        base_t{str}, file_it{str, false /*no_init!*/}

    {
        /* potentially extract useful runtime options from config and store as members */

        bool file_format_read = false;
        while (file_it != std::default_sentinel && file_it.peak() == '#')
        {
            ++file_it;
            std::string_view l = file_it->line;
            header.add_raw_line(l);
        }

        get<field::header>(raw_record)  = header.raw_header();
    }
};

} // namespace seqan
