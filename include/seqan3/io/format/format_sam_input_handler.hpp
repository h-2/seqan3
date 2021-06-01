// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * brief Provides the seqan3::format_fasta.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <iostream>
#include <string>
#include <string_view>
#include <vector>

#include <seqan3/core/range/type_traits.hpp>
#include <seqan3/io/detail/ignore_output_iterator.hpp>
#include <seqan3/io/detail/misc.hpp>
#include <seqan3/io/format/input_format_handler_base.hpp>
#include <seqan3/io/stream/iterator.hpp>
#include <seqan3/range/detail/misc.hpp>
#include <seqan3/range/views/char_strictly_to.hpp>
#include <seqan3/range/views/istreambuf.hpp>
#include <seqan3/range/views/to_char.hpp>
#include <seqan3/range/views/take_line.hpp>
#include <seqan3/range/views/take_until.hpp>
#include <seqan3/std/algorithm>
#include <seqan3/std/ranges>
#include <seqan3/utility/char_operations/predicate.hpp>

namespace seqan3
{

template <>
class input_format_handler<format_sam> : public input_format_handler_base<input_format_handler<format_sam>>
{
private:
    /* CRTP */
    using base_t = input_format_handler_base<input_format_handler<format_sam>>;
    using base_t::stream;

    friend base_t;

    /* RAW RECORD HANDLING*/
    using format_fields = fields<field::qname,
                                 field::flag,
                                 field::rname,
                                 field::pos,
                                 field::mapq,
                                 field::cigar,
                                 field::rnext,
                                 field::pnext,
                                 field::tlen,
                                 field::seq,
                                 field::qual,
                                 field::optionals,
                                 field::header>;
    using raw_record_type = seqan3::record<list_traits::repeat<format_fields::as_array.size(), std::string_view>,
                                           format_fields>;

    raw_record_type raw_record;
    plain_file_iterator<char, std::char_traits<char>, true> file_it;
    std::string raw_header;

    void read_raw_record()
    {
        ++file_it;

        get<field::qname>(raw_record) = (*file_it).fields[0];
        get<field::flag> (raw_record) = (*file_it).fields[1];
        get<field::rname>(raw_record) = (*file_it).fields[2];
        get<field::pos>  (raw_record) = (*file_it).fields[3];
        get<field::mapq> (raw_record) = (*file_it).fields[4];
        get<field::cigar>(raw_record) = (*file_it).fields[5];
        get<field::rnext>(raw_record) = (*file_it).fields[6];
        get<field::pnext>(raw_record) = (*file_it).fields[7];
        get<field::tlen> (raw_record) = (*file_it).fields[8];
        get<field::seq>  (raw_record) = (*file_it).fields[9];
        get<field::qual> (raw_record) = (*file_it).fields[10];

        // == .end() that is guaranteed to be char*
        char * end_qual   = (*file_it).fields[10].data() + (*file_it).fields[10].size();
        char * end_line   = (*file_it).line.data() +       (*file_it).line.size();
        // == string_view{} if there are no optional fields
        get<field::optionals>  (raw_record) = std::string_view{end_line, end_line - end_qual};

        // header_ptr does not need to be reset
    }

    /* PARSED RECORD HANDLING */
    // override the general case to handle "*"
    template <field field_id, typename parsed_field_t>
    void parse_field(tag_t<field_id> const & /**/, parsed_field_t & parsed_field)
    {
        std::string_view raw_field = get<field_id>(raw_record);
        if (raw_field != "*")
            static_cast<base_t &>(*this).parse_field(tag<field_id>, parsed_field);
    }

    void parse_field(tag_t<flag> const & /**/, sam_flag & parsed_field)
    {
        std::string_view raw_field = get<flag>(raw_record);
        uint16_t tmp = 0;

        if (auto r = std::from_chars(raw_field.begin(), raw_field.end(), parsed_field); r.ec != std::errc{})
            throw format_error{std::string{"Failed to convert \""} + std::string{raw_field} + "\" into a number."};

        parsed_field = sam_flag{tmp};
    }

    template <typename parsed_field_t>
        requires std::ranges::output_range<parsed_field_t, cigar>
    void parse_field(tag_t<field::cigar> const & /**/, parsed_field_t & parsed_field)
    {
        //TODO parse cigar
    }

    template <typename parsed_field_t>
    void parse_field(tag_t<field::tags> const & /**/, parsed_field_t & parsed_field)
    {
        if ((*file_it).fields.size() > 10) // there are optional fields
        {
            // TODO store into dictionary
        }
    }

    void parse_field(tag_t<field::header_ptr> const & /**/, std::string * & parsed_field)
    {
        parsed_field = &header;
    }



public:

    input_format_handler()                                              = default;
    input_format_handler(input_format_handler const &)                  = delete;
    input_format_handler(input_format_handler && )                      = default;
    input_format_handler & operator=(input_format_handler const &)      = delete;
    input_format_handler & operator=(input_format_handler &&)           = default;

    template <typename config_t>
    input_format_handler(std::istream & str, config_t const & /*cfg*/) :
        base_t{str}, file_it{*str.rdbuf(), false /*no_init!*/}

    {
        /* potentially extract useful runtime options from config and store as members */

        while (file_it != std::default_sentinel && file_it.peak() == '#')
        {
            ++file_it;
            header += (*file_it).line;
        }

        get<field::header>(raw_record) = std::string_view{raw_header};
    }
};

} // namespace seqan
