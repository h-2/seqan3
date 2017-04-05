// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2017, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2017, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ============================================================================
// Author: Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
// ============================================================================

#pragma once

#include <type_traits>

namespace seqan3::detail
{

/* A tuple class that is an aggregate POD type.
 * It serves as the basis of alphabet_composition.
 */

template <typename ...types>
struct pod_tuple
{
    static_assert(sizeof...(types) > 0, "Need at least one element in pod_tuple.");
    static_assert(sizeof...(types) < 6, "Only up to 5 elements currently supported by pod_tuple.");
};

#define SEQAN_NOT_POD "If you are not going to insert a POD type, use std::tuple instead."

template <typename type0>
struct pod_tuple<type0>
{
    static_assert(std::is_pod_v<type0>, SEQAN_NOT_POD);
    type0 _val0;
};

template <typename type0, typename type1>
struct pod_tuple<type0, type1>
{
    static_assert(std::is_pod_v<type0>, SEQAN_NOT_POD);
    static_assert(std::is_pod_v<type1>, SEQAN_NOT_POD);
    type0 _val0;
    type1 _val1;
};

template <typename type0, typename type1, typename type2>
struct pod_tuple<type0, type1, type2>
{
    static_assert(std::is_pod_v<type0>, SEQAN_NOT_POD);
    static_assert(std::is_pod_v<type1>, SEQAN_NOT_POD);
    static_assert(std::is_pod_v<type2>, SEQAN_NOT_POD);
    type0 _val0;
    type1 _val1;
    type2 _val2;
};

template <typename type0, typename type1, typename type2, typename type3>
struct pod_tuple<type0, type1, type2, type3>
{
    static_assert(std::is_pod_v<type0>, SEQAN_NOT_POD);
    static_assert(std::is_pod_v<type1>, SEQAN_NOT_POD);
    static_assert(std::is_pod_v<type2>, SEQAN_NOT_POD);
    static_assert(std::is_pod_v<type3>, SEQAN_NOT_POD);
    type0 _val0;
    type1 _val1;
    type2 _val2;
    type3 _val3;
};

template <typename type0, typename type1, typename type2, typename type3, typename type4>
struct pod_tuple<type0, type1, type2, type3, type4>
{
    static_assert(std::is_pod_v<type0>, SEQAN_NOT_POD);
    static_assert(std::is_pod_v<type1>, SEQAN_NOT_POD);
    static_assert(std::is_pod_v<type2>, SEQAN_NOT_POD);
    static_assert(std::is_pod_v<type3>, SEQAN_NOT_POD);
    static_assert(std::is_pod_v<type4>, SEQAN_NOT_POD);
    type0 _val0;
    type1 _val1;
    type2 _val2;
    type3 _val3;
    type4 _val4;
};

#undef SEQAN_NOT_POD

// TODO move this somewhere in core
template <std::size_t i, typename head_t, typename ...tail_types >
struct get_ith_type :
    get_ith_type<i - 1, tail_types...>
{
//TODO fix this
//     static_assert(i > sizeof...(tail_types), "Trying to access a type behind end of pack.");
};

template <typename head_t, typename ...tail_types >
struct get_ith_type<0, head_t, tail_types...>
{
   using type = head_t;
};

template <std::size_t i, typename ...types>
using get_ith_type_t = typename get_ith_type<i, types...>::type;

} // namespace seqan3::detail

namespace std
{

template <std::size_t i, typename ...types>
    requires i < sizeof...(types)
constexpr auto & get(seqan3::detail::pod_tuple<types...> & t)
{
    if constexpr (i == 0)
        return t._val0;
    else if constexpr (i == 1)
        return t._val1;
    else if constexpr (i == 2)
        return t._val2;
    else if constexpr (i == 3)
        return t._val3;
    else
        return t._val4;
}

template <std::size_t i, typename ...types>
    requires i < sizeof...(types)
constexpr auto const & get(seqan3::detail::pod_tuple<types...> const & t) noexcept
{
    if constexpr (i == 0)
        return t._val0;
    else if constexpr (i == 1)
        return t._val1;
    else if constexpr (i == 2)
        return t._val2;
    else if constexpr (i == 3)
        return t._val3;
    else
        return t._val4;
}

// extra overloads for temporaries required, because members of temporaries may only be returned as temporaries

template <std::size_t i, typename ...types>
    requires i < sizeof...(types)
constexpr auto && get(seqan3::detail::pod_tuple<types...> && t)
{
    if constexpr (i == 0)
        return std::move(t._val0);
    else if constexpr (i == 1)
        return std::move(t._val1);
    else if constexpr (i == 2)
        return std::move(t._val2);
    else if constexpr (i == 3)
        return std::move(t._val3);
    else
        return std::move(t._val4);
}

template <std::size_t i, typename ...types>
    requires i < sizeof...(types)
constexpr auto const && get(seqan3::detail::pod_tuple<types...> const && t) noexcept
{
    if constexpr (i == 0)
        return std::move(t._val0);
    else if constexpr (i == 1)
        return std::move(t._val1);
    else if constexpr (i == 2)
        return std::move(t._val2);
    else if constexpr (i == 3)
        return std::move(t._val3);
    else
        return std::move(t._val4);
}

//TODO type based std::get

// element type access
template <std::size_t i, typename ...types >
struct tuple_element<i, seqan3::detail::pod_tuple<types...>>
{
    using type = seqan3::detail::get_ith_type_t<i, types...>;
};

} // namespace std
