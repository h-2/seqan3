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
class input_format_handler<format_fasta> : public input_format_handler_base<input_format_handler<format_fasta>>
{
private:
    /* CRTP */
    using base_t = input_format_handler_base<input_format_handler<format_fasta>>;
    using base_t::stream;

    friend base_t;

    /* RAW RECORD HANDLING*/
    using format_fields = fields<field::id, field::seq>;
    using raw_record_type = seqan3::record<list_traits::repeat<format_fields::as_array.size(), std::string_view>,
                                           format_fields>;

    raw_record_type raw_record;

    // will be externalised later
    std::string buffer_buffer;
    void read_raw_record()
    {
        buffer_buffer.clear();
        raw_record.clear();

        auto stream_view = views::istreambuf(*stream);
        constexpr auto is_id = is_char<'>'> || is_char<';'>;

        std::ranges::copy(stream_view | views::take_line_or_throw                    // read line
                                      | std::views::drop_while(is_id || is_blank),    // skip leading >
                          std::cpp20::back_inserter(buffer_buffer));

        size_t sep_pos = buffer_buffer.size();


        std::ranges::copy(stream_view | views::take_until(is_id),                   // until next header (or end)
                          std::cpp20::back_inserter(buffer_buffer));

        get<field::id>(raw_record) = std::string_view{buffer_buffer.data(), sep_pos};
        get<field::seq>(raw_record) = std::string_view{buffer_buffer.data() + sep_pos, buffer_buffer.size() - sep_pos};
    }

    /* PARSED RECORD HANDLING */
    static constexpr record field_parsers = detail::make_record<fields<field::seq>>(
        std::views::filter(!(is_space || is_digit))     // digits and spaces filtered from seq
        );

public:

    input_format_handler()                                              = default;
    input_format_handler(input_format_handler const &)                  = delete;
    input_format_handler(input_format_handler && )                      = default;
    input_format_handler & operator=(input_format_handler const &)      = delete;
    input_format_handler & operator=(input_format_handler &&)           = default;

    template <typename config_t>
    input_format_handler(std::istream & str, config_t const & /*cfg*/) :
        base_t{str}
    {
        /* potentially extract useful runtime options from config and store as members */
    }
};

} // namespace seqan
