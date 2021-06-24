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

#include <string_view>
#include <vector>

#include <seqan3/alphabet/views/char_strictly_to.hpp>
#include <seqan3/io/detail/misc.hpp>
#include <seqan3/io/detail/record.hpp>
#include <seqan3/io/utility.hpp>
#include <seqan3/io/record.hpp>
#include <seqan3/std/charconv>
#include <seqan3/std/span>

namespace seqan3::detail
{

template <field field_id, typename input_format_handler_t>
    requires (requires { input_format_handler_t::field_parser_views; })
auto get_or_view_all(input_format_handler_t && handler)
{
    return detail::get_or<field_id>(handler, std::views::all /*NOOP*/);
}

template <field field_id, typename input_format_handler_t>
auto get_or_view_all(input_format_handler_t && handler)
{
    return std::views::all;
}

} // namespace seqan3::detail

namespace seqan3
{

template <typename format_t>
class input_format_handler;


template <typename derived_t>
class input_format_handler_base
{
private:

    /* CRTP STUFF */
    friend derived_t;

    derived_t * to_derived()
    {
        return static_cast<derived_t *>(this);
    }

    derived_t const * to_derived() const
    {
        return static_cast<derived_t const *>(this);
    }

    /* members */
    std::istream * stream = nullptr;


    /* stuff for turning raw record into parsed record */
    template <field ... field_ids, typename parsed_record_t>
    void parse_record(tag_t<field_ids...> const & /**/, parsed_record_t & parsed_record)
    {
        (to_derived()->parse_field_wrapper(tag<field_ids>, parsed_record), ...);
    }

    template <field field_id, typename parsed_record_t>
    void parse_field_wrapper(tag_t<field_id> const & /**/, parsed_record_t & parsed_record)
    {
        if constexpr (parsed_record_t::field_ids::contains(field_id))
        {
            std::string_view raw_field = get<field_id>(to_derived()->raw_record);
            auto & parsed_field        = get<field_id>(parsed_record);

            using parsed_field_t = std::remove_reference_t<decltype(parsed_field)>;

            if constexpr (std::same_as<parsed_field_t, std::span<std::byte>>)
            {
                parsed_field = std::span<std::byte>{reinterpret_cast<std::byte*>(raw_field.data()),
                                                    raw_field.size()};
            }
            else
            {
                to_derived()->parse_field(tag<field_id>, parsed_field);
            }
        }
        // fields that are not in format or not in target record are simply ignored
    }

    // sane default for handling most sequences and arithmetic fields; override this for special behaviour
    template <field field_id, typename parsed_field_t>
    void parse_field(tag_t<field_id> const & /**/, parsed_field_t & parsed_field)
    {
        std::string_view raw_field = get<field_id>(to_derived()->raw_record);

        if constexpr (std::ranges::range<parsed_field_t>) //TODO output_range
        {
            using target_alph_type = std::ranges::range_value_t<parsed_field_t>;

            auto adaptor = detail::get_or_view_all<field_id>(to_derived());

            if constexpr (std::constructible_from<target_alph_type, char>) // no alphabet conversion
            {
                detail::sized_range_copy(raw_field | adaptor,
                                         parsed_field);
            }
            else if constexpr (alphabet<target_alph_type>)
            {
                detail::sized_range_copy(raw_field | adaptor | views::char_strictly_to<target_alph_type>,
                                         parsed_field);
            }
            else
            {
                static_assert(arithmetic<parsed_field_t>/*always false*/,
                          "Format X does not know how to convert field Y into type Z. Provide different traits or a "
                          "custom format handler.");
                //TODO replace X Y and Z with actual strings generated from types.
            }
        }
        else if constexpr (arithmetic<parsed_field_t>)
        {
            if (auto r = std::from_chars(raw_field.data(), raw_field.data() + raw_field.size(), parsed_field);
                r.ec != std::errc{})
            {
                throw format_error{std::string{"Failed to convert \""} + std::string{raw_field} + "\" into a number."};
            }
        }
        else
        {
            static_assert(arithmetic<parsed_field_t>/*always false*/,
                          "Format X does not know how to convert field Y into type Z. Provide different traits or a "
                          "custom format handler.");
            //TODO replace X Y and Z with actual strings generated from types.
        }
    }


    // private to prevent wrong derivation
    input_format_handler_base()                                                   = default;
    input_format_handler_base(input_format_handler_base const &)                  = delete;
    input_format_handler_base(input_format_handler_base && )                      = default;
    input_format_handler_base & operator=(input_format_handler_base const &)      = delete;
    input_format_handler_base & operator=(input_format_handler_base &&)           = default;

    input_format_handler_base(std::istream & str)
        : stream{&str}
    {}

public:

    template <typename parsed_record_t>
    void parse_next_record_into(parsed_record_t & parsed_record)
    {
        if (std::istreambuf_iterator<char>{*stream} == std::istreambuf_iterator<char>{})
            throw 42;

        // create new raw record
        to_derived()->read_raw_record();

        // create new parsed record
        parsed_record.clear();
        to_derived()->parse_record(typename derived_t::format_fields{}, parsed_record);
    }

};

} // namespace seqan
