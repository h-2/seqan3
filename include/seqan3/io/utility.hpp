// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the seqan3::record template and the seqan3::field enum.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <tuple>

#include <seqan3/utility/type_list/type_list.hpp>
#include <seqan3/utility/views/eager_split.hpp>

namespace seqan3
{


template <typename t, typename ... ts>
concept is_one_of = (std::same_as<t, ts> || ...);

template <auto ... more_vs>
struct tag_t
{
    static constexpr size_t size = 0;

    static constexpr auto as_tuple = std::tuple{};

    static constexpr bool contains(auto &&)
    {
        return false;
    }

    static constexpr size_t index_of(auto &&)
    {
        return static_cast<size_t>(-1ULL);
    }
};

template <auto v, auto ... more_vs>
struct tag_t<v, more_vs...>
{
    static constexpr auto first_value = v;

    static constexpr size_t size = sizeof...(more_vs) + 1;

    static constexpr auto as_tuple = std::tuple{v, more_vs...};

    static constexpr auto as_array = [] ()
    {
        if constexpr ((std::same_as<decltype(v), decltype(more_vs)> && ...))
        {
            return std::array<decltype(v), size>{v, more_vs...};
        }
        else
        {
            return;
        }
    } ();

    static constexpr bool unique_values = ((v != more_vs) && ...);

    static constexpr bool contains(auto && s)
    {
        return s == v || ((s == more_vs) || ...);
    }

    static constexpr size_t index_of(auto && s)
    {
        for (size_t i = 0; i < size; ++i)
            if (as_array[i] == s)
                return i;
        return static_cast<size_t>(-1ULL);
    }

};

template <auto ... more_vs>
inline constexpr tag_t<more_vs...> tag{};


template <typename type, typename ... more_types>
inline constexpr type_list<type, more_types...> type_tag{};

//!\brief Enumerator to ease "dynamic typing" in alignment map IO and variant IO.
enum class io_type_id
{
    flag,
    char8,
    int32,
    float32,
    string,
    vector_of_char8,
    vector_of_int8,
    vector_of_uint8,
    vector_of_int16,
    vector_of_uint16,
    vector_of_int32,
    vector_of_uint32,
    vector_of_float32,
    vector_of_string
};

//!\brief Variant to handle "dynamic typing" in alignment map IO and variant IO.
template <bool shallow = true>
using io_type_variant = std::variant<int8_t,                    // var_io only
                                     char,                      // am_io only
                                     int32_t,
                                     float,
                                     std::conditional_t<shallow, std::string_view, std::string>,
                                     std::vector<char>,         // var_io only
                                     std::vector<int8_t>,       // am_io only
                                     std::vector<uint8_t>,      // am_io only
                                     std::vector<int16_t>,      // am_io only
                                     std::vector<uint16_t>,     // am_io only
                                     std::vector<int32_t>,
                                     std::vector<uint32_t>,     // am_io only
                                     std::vector<float>,
                                     std::vector<std::conditional_t<shallow, std::string_view, std::string>>>; // var_io only

template <bool shallow = true>
using io_type_vector_variant = std::variant<std::vector<int8_t>,                    // var_io only
                                            std::vector<char>,                      // am_io only
                                            std::vector<int32_t>,
                                            std::vector<float>,
                                            std::vector<std::conditional_t<shallow, std::string_view, std::string>>,
                                            std::vector<std::vector<char>>,         // var_io only
                                            std::vector<std::vector<int8_t>>,       // am_io only
                                            std::vector<std::vector<uint8_t>>,      // am_io only
                                            std::vector<std::vector<int16_t>>,      // am_io only
                                            std::vector<std::vector<uint16_t>>,     // am_io only
                                            std::vector<std::vector<int32_t>>,
                                            std::vector<std::vector<uint32_t>>,     // am_io only
                                            std::vector<std::vector<float>>,
                                            std::vector<std::vector<std::conditional_t<shallow, std::string_view, std::string>>>>; // var_io only


} // namespace seqan3

namespace seqan3::detail
{

void string_to_number(std::string_view input, arithmetic auto & number)
{
    std::from_chars_result res = std::from_chars(input.data(), input.data() + input.size(), number);
    if (res.ec != std::errc{} || res.ptr != input.data() + input.size())
        throw std::runtime_error{std::string{"Could not convert \""} + std::string{input} + "\" into a number."};
}

template <typename t>
concept is_io_type_variant = is_one_of<t, io_type_variant<true>, io_type_variant<false>>;

template <typename t>
concept is_io_type_vector_variant = is_one_of<t, io_type_vector_variant<true>, io_type_vector_variant<false>>;

template <typename t>
    requires is_io_type_variant<t> || is_io_type_vector_variant<t>
inline void init_io_type_variant(io_type_id const id_, t & output)
{
    switch (id_)
    {
        case io_type_id::flag:
        {
            constexpr size_t id = static_cast<size_t>(io_type_id::flag);
            output.template emplace<id>();
            return;
        }
        case io_type_id::char8:
        {
            constexpr size_t id = static_cast<size_t>(io_type_id::char8);
            output.template emplace<id>();
            return;
        }
        case io_type_id::int32:
        {
            constexpr size_t id = static_cast<size_t>(io_type_id::int32);
            output.template emplace<id>();
            return;
        }
        case io_type_id::float32:
        {
            constexpr size_t id = static_cast<size_t>(io_type_id::float32);
            output.template emplace<id>();
            return;
        }
        case io_type_id::string:
        {
            constexpr size_t id = static_cast<size_t>(io_type_id::string);
            output.template emplace<id>();
            return;
        }
        case io_type_id::vector_of_char8:
        {
            constexpr size_t id = static_cast<size_t>(io_type_id::vector_of_char8);
            output.template emplace<id>();
            return;
        }
        case io_type_id::vector_of_int8:
        {
            constexpr size_t id = static_cast<size_t>(io_type_id::vector_of_int8);
            output.template emplace<id>();
            return;
        }
        case io_type_id::vector_of_uint8:
        {
            constexpr size_t id = static_cast<size_t>(io_type_id::vector_of_uint8);
            output.template emplace<id>();
            return;
        }
        case io_type_id::vector_of_int16:
        {
            constexpr size_t id = static_cast<size_t>(io_type_id::vector_of_int16);
            output.template emplace<id>();
            return;
        }
        case io_type_id::vector_of_uint16:
        {
            constexpr size_t id = static_cast<size_t>(io_type_id::vector_of_uint16);
            output.template emplace<id>();
            return;
        }
        case io_type_id::vector_of_int32:
        {
            constexpr size_t id = static_cast<size_t>(io_type_id::vector_of_int32);
            output.template emplace<id>();
            return;
        }
        case io_type_id::vector_of_uint32:
        {
            constexpr size_t id = static_cast<size_t>(io_type_id::vector_of_uint32);
            output.template emplace<id>();
            return;
        }
        case io_type_id::vector_of_float32:
        {
            constexpr size_t id = static_cast<size_t>(io_type_id::vector_of_float32);
            output.template emplace<id>();
            return;
        }
        case io_type_id::vector_of_string:
        {
            constexpr size_t id = static_cast<size_t>(io_type_id::vector_of_string);
            output.template emplace<id>();
            return;
        }
    }
}

} // namespace seqan3::detail
