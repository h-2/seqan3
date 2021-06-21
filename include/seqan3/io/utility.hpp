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

namespace seqan3
{

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
    vector_of_int8,
    vector_of_uint8,
    vector_of_int16,
    vector_of_uint16,
    vector_of_int32,
    vector_of_uint32,
    vector_of_float32
};

//!\brief Variant to handle "dynamic typing" in alignment map IO and variant IO.
using io_type_variant = std::variant<bool,
                                     char,
                                     int32_t,
                                     float,
                                     std::string,
                                     std::vector<int8_t>,
                                     std::vector<uint8_t>,
                                     std::vector<int16_t>,
                                     std::vector<uint16_t>,
                                     std::vector<int32_t>,
                                     std::vector<uint32_t>,
                                     std::vector<float>>;


} // namespace seqan3

