// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\cond DEV
 * \file
 * \brief Provides auxiliary data structures and functions for seqan3::record and seqan3::fields.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \endcond
 */

#pragma once

#include <seqan3/std/ranges>

#include <seqan3/core/concept/tuple.hpp>
#include <seqan3/io/record.hpp>
#include <seqan3/utility/type_list/traits.hpp>

namespace seqan3::detail
{

// ----------------------------------------------------------------------------
// get_or_ignore
// ----------------------------------------------------------------------------

/*!\addtogroup io
 *!\{
 */
//!\brief Access an element in a std::tuple or seqan3::record; return reference to std::ignore if not contained.
template <field f, typename field_types, typename field_ids>
auto & get_or_ignore(record<field_types, field_ids> & r)
{
    if constexpr (field_ids::contains(f))
        return std::get<field_ids::index_of(f)>(r);
    else
        return std::ignore;
}

//!\copydoc seqan3::detail::get_or_ignore
template <field f, typename field_types, typename field_ids>
auto const & get_or_ignore(record<field_types, field_ids> const & r)
{
    if constexpr (field_ids::contains(f))
        return std::get<field_ids::index_of(f)>(r);
    else
        return std::ignore;
}

//!\copydoc seqan3::detail::get_or_ignore
template <size_t i, template <tuple_like ...types_> typename tuple_like_t, typename ...types>
auto & get_or_ignore(tuple_like_t<types...> & t)
{
    if constexpr (i < sizeof...(types))
        return std::get<i>(t);
    else
        return std::ignore;
}

//!\copydoc seqan3::detail::get_or_ignore
template <size_t i, template <tuple_like ...types_> typename tuple_like_t, typename ...types>
auto const & get_or_ignore(tuple_like_t<types...> const & t)
{
    if constexpr (i < sizeof...(types))
        return std::get<i>(t);
    else
        return std::ignore;
}
//!\}


// ----------------------------------------------------------------------------
// get_or
// ----------------------------------------------------------------------------

/*!\addtogroup io
 *!\{
 */
//!\brief Access an element in a std::tuple or seqan3::record; return or_value if not contained.
template <field f, typename field_types, typename field_ids, typename or_type>
decltype(auto) get_or(record<field_types, field_ids> & r, or_type && or_value)
{
    if constexpr (field_ids::contains(f))
        return std::get<field_ids::index_of(f)>(r);
    else
        return std::forward<or_type>(or_value);
}

//!\copydoc seqan3::detail::get_or
template <field f, typename field_types, typename field_ids, typename or_type>
decltype(auto) get_or(record<field_types, field_ids> const & r, or_type && or_value)
{
    if constexpr (field_ids::contains(f))
        return std::get<field_ids::index_of(f)>(r);
    else
        return std::forward<or_type>(or_value);
}

//!\copydoc seqan3::detail::get_or
template <size_t i, typename or_type, typename ...types>
decltype(auto) get_or(std::tuple<types...> & t, or_type && or_value)
{
    if constexpr (i < sizeof...(types))
        return std::get<i>(t);
    else
        return std::forward<or_type>(or_value);
}

//!\copydoc seqan3::detail::get_or
template <size_t i, typename or_type, typename ...types>
decltype(auto) get_or(std::tuple<types...> const & t, or_type && or_value)
{
    if constexpr (i < sizeof...(types))
        return std::get<i>(t);
    else
        return std::forward<or_type>(or_value);
}
//!\}

} // namespace seqan3::detail
