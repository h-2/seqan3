// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the seqan3::io_cfg::select_fields configuration object.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <seqan3/io/configuration/detail.hpp>
#include <seqan3/io/detail/record.hpp>

namespace seqan3::io_cfg
{

template <typename fields_type>
struct select_fields_t : public pipeable_config_element<select_fields_t<fields_type>>
{
    constexpr select_fields_t() = default;

    using type = fields_type;

    //!\privatesection
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr detail::io_config_id id{detail::io_config_id::fields};
};

template <field ...field_args>
select_fields_t<fields<field_args...>> select_fields{};

} // namespace seqan3::io_cfg
