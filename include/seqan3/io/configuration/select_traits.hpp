// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the seqan3::io_cfg::select_traits configuration object.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <seqan3/io/configuration/detail.hpp>

namespace seqan3::io_cfg
{

template <typename type_traits_t>
struct select_traits_t : public pipeable_config_element<select_traits_t<type_traits_t>>
{

    constexpr select_traits_t() = default;

    using type = type_traits_t;

    //!\privatesection
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr detail::io_config_id id{detail::io_config_id::traits};
};

template <typename type_traits_t>
select_traits_t<type_traits_t> select_traits{};

} // namespace seqan3::io_cfg
