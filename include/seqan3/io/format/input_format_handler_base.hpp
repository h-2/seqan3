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

    /*!\name Parsing individual fields - defaults (step 3)
     * \{
     */
    //!\brief Not parsing at all / *raw IO*.
    static void parse_field_impl(std::string_view const in, std::span<std::byte> & parsed_field)
    {
        parsed_field = std::span<std::byte>{reinterpret_cast<std::byte*>(in.data()), in.size()};
    }

    //!\brief Parsing into string views. NOTE: binary formats may want to = delete override this.
    static void parse_field_impl(std::string_view const in, std::string_view & parsed_field)
    {
        parsed_field = in;
    }

    //!\brief Parsing into transformed string views. NOTE: binary formats may want to = delete override this.
    template <typename fun_t>
    static void parse_field_impl(std::string_view const in,
                                 std::ranges::transform_view<std::string_view, fun_t> & parsed_field)
    {
        parsed_field = std::ranges::transform_view<std::string_view, fun_t>{in, fun_t{}};
    }

    //!\brief Parse into string-like types.
    template <typename parsed_field_t>
        requires std::ranges::range<parsed_field_t> &&
                 detail::back_insertable_with<parsed_field_t, char>
    static void parse_field_impl(std::string_view const in, parsed_field_t & parsed_field)
    {
//         auto adaptor = detail::get_or_view_all<field_id>(to_derived());
//         detail::sized_range_copy(raw_field | adaptor, parsed_field);
        detail::sized_range_copy(in, parsed_field);
    }

    //!\brief Parse into containers of alphabets.
    template <typename parsed_field_t>
        requires (std::ranges::range<parsed_field_t> &&
                  !detail::back_insertable_with<parsed_field_t, char> &&
                  alphabet<std::ranges::range_reference_t<parsed_field_t>> &&
                  detail::back_insertable_with<parsed_field_t, std::ranges::range_reference_t<parsed_field_t>>)
    static void parse_field_impl(std::string_view const in, parsed_field_t & parsed_field)
    {
        using target_alph_type = std::ranges::range_value_t<parsed_field_t>;
        detail::sized_range_copy(in | views::char_strictly_to<target_alph_type>,
                                 parsed_field);
//         auto adaptor = detail::get_or_view_all<field_id>(to_derived());
//         detail::sized_range_copy(raw_field | adaptor | views::char_strictly_to<target_alph_type,
//                                  parsed_field);
    }


    //!\brief Parse into a a numerical type
    template <typename parsed_field_t>
        requires arithmetic<parsed_field_t>
    static void parse_field_impl(std::string_view const in, parsed_field_t & parsed_field)
    {
        detail::string_to_number(in, parsed_field);
    }
    //!\}

    /*!\name Parsing individual fields (step 2)
     * \{
     */
    //!\brief Default is no handler.
    template <typename tag_type, typename parsed_field_t>
    void parse_field(tag_type const & /**/, parsed_field_t & parsed_field)
    {
            static_assert(arithmetic<parsed_field_t>/*always false*/,
                          "Format X does not know how to convert field Y into type Z. Provide different traits or a "
                          "custom format handler.");
            //TODO replace X Y and Z with actual strings generated from types.
    }

    //!\brief Various target types have sane default implementations.
    template <field field_id, typename parsed_field_t>
    void parse_field(tag_t<field_id> const & /**/, parsed_field_t & parsed_field)
        requires (requires { derived_t::parse_field_impl(std::string_view{}, parsed_field); })
    {
        to_derived()->parse_field_impl(get<field_id>(to_derived()->raw_record), parsed_field);
    }
    //!\}

    /*!\name Parsing record (step 1)
     * \{
     */
    template <field field_id, typename parsed_record_t>
    void parse_record_impl(tag_t<field_id> const & /**/, parsed_record_t & parsed_record)
    {
        if constexpr (parsed_record_t::field_ids::contains(field_id))
        {
            auto & parsed_field = get<field_id>(parsed_record);
            to_derived()->parse_field(tag<field_id>, parsed_field);

        }
        // fields that are not in format or not in target record are simply ignored
    }

    template <field ... field_ids, typename parsed_record_t>
    void parse_record(tag_t<field_ids...> const & /**/, parsed_record_t & parsed_record)
    {
        (to_derived()->parse_record_impl(tag<field_ids>, parsed_record), ...);
    }
    //!\}


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
