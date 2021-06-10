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
#include <seqan3/io/alignment_map_io/misc.hpp>
#include <seqan3/io/plaintext_io/reader.hpp>
#include <seqan3/io/format/format_sam.hpp>
#include <seqan3/io/format/input_format_handler_base.hpp>
#include <seqan3/io/stream/iterator.hpp>
#include <seqan3/utility/char_operations/predicate.hpp>
#include <seqan3/utility/type_list/traits.hpp>

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
    using format_fields = tag_t<field::qname,
                                field::flag,
                                field::ref_id,
                                field::pos,
                                field::mapq,
                                field::cigar,
                                field::next_ref_id,
                                field::next_pos,
                                field::tlen,
                                field::seq,
                                field::qual,
                                field::optionals,
                                field::header>;
    using raw_record_type = seqan3::record<list_traits::repeat<format_fields::size, std::string_view>,
                                           format_fields>;

    raw_record_type raw_record;
    plain_io::detail::plaintext_input_iterator<char, std::char_traits<char>, plain_io::record_kind::line_and_fields> file_it;
    std::string raw_header;
    am_io::header parsed_header;

    void read_raw_record()
    {
        ++file_it;

        if (file_it->fields.size() < 11)
            throw format_error{"Encountered line with less than the 11 required columns."};

        get<field::qname>      (raw_record) = (*file_it).fields[0];
        get<field::flag>       (raw_record) = (*file_it).fields[1];
        get<field::ref_id>     (raw_record) = (*file_it).fields[2];
        get<field::pos>        (raw_record) = (*file_it).fields[3];
        get<field::mapq>       (raw_record) = (*file_it).fields[4];
        get<field::cigar>      (raw_record) = (*file_it).fields[5];
        get<field::next_ref_id>(raw_record) = (*file_it).fields[6];
        get<field::next_pos>   (raw_record) = (*file_it).fields[7];
        get<field::tlen>       (raw_record) = (*file_it).fields[8];
        get<field::seq>        (raw_record) = (*file_it).fields[9];
        get<field::qual>       (raw_record) = (*file_it).fields[10];
        get<field::header>     (raw_record) = raw_header;

        // == .end() that is guaranteed to be char*
        char const * end_qual   = (*file_it).fields[10].data() + (*file_it).fields[10].size();
        char const * end_line   = (*file_it).line.data() +       (*file_it).line.size();
        // == string_view{} if there are no optional fields
        get<field::optionals>  (raw_record) = std::string_view{end_qual, static_cast<size_t>(end_line - end_qual)};
    }

    /* PARSED RECORD HANDLING */
    //TODO fix this in base class
    static constexpr record field_parsers = record<type_list<>, tag_t<>>{};

    // override the general case to handle "*"
    template <field field_id, typename parsed_field_t>
    void parse_field(tag_t<field_id> const & /**/, parsed_field_t & parsed_field) const
    {
        std::string_view raw_field = get<field_id>(raw_record);
        if (raw_field != "*")
            static_cast<base_t const *>(this)->parse_field(tag<field_id>, parsed_field);
    }

    void parse_field(tag_t<field::flag> const & /**/, am_io::flag & parsed_field) const
    {
        std::string_view raw_field = get<field::flag>(raw_record);
        uint16_t tmp = 0;

        if (auto r = std::from_chars(raw_field.data(), raw_field.data() + raw_field.size(), tmp);
            r.ec != std::errc{})
        {
            throw format_error{std::string{"Failed to convert \""} + std::string{raw_field} + "\" into a number."};
        }

        parsed_field = am_io::flag{tmp};
    }

    template <typename parsed_field_t>
        requires std::ranges::output_range<parsed_field_t, cigar>
    void parse_field(tag_t<field::cigar> const & /**/, parsed_field_t & parsed_field) const
    {
        std::string_view raw_field = get<field::cigar>(raw_record);
        char const * b = std::ranges::data(raw_field);
        char const * e = b + std::ranges::size(raw_field);

        while (b != e)
        {
            uint64_t i = 0;
            std::from_chars_result res = std::from_chars(b, e, i);
            if (res.ec != std::errc{} || res.ptr == e)
                throw format_error{"Corrupted cigar string encountered"};

            b = res.ptr;
            parsed_field.emplace_back(i, assign_char_to(*b, exposition_only::cigar_operation{}));
            ++b;
        }
    }

    template <typename parsed_field_t>
    void parse_field(tag_t<field::optionals> const & /**/, parsed_field_t & parsed_field) const
    {
        if ((*file_it).fields.size() > 10) // there are optional fields
        {
            // TODO store into dictionary
        }
    }

    void parse_field(tag_t<field::header> const & /**/, am_io::header const * & parsed_field) const
    {
        parsed_field = &parsed_header;
    }

    static void parse_header(std::string_view raw_header, am_io::header & parsed_header)
    {
        //TODO implement
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

        while (file_it != std::default_sentinel && file_it.peak() == '@')
        {
            ++file_it;
            raw_header += (*file_it).line;
        }

        parse_header(raw_header, parsed_header);
    }
};

} // namespace seqan
