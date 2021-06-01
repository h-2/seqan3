// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief  Provides some utility functions for the io configurations.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/algorithm/all.hpp>

namespace seqan3::detail
{

/*!\brief An internal enum to check for a consistent configuration object.
 * \ingroup alignment_configuration
 */
enum struct io_config_id : uint8_t
{
    fields,
    formats,
    traits,
    ref_info, // only alignment
    SIZE          //!< Represents the number of configuration elements.
};

// ----------------------------------------------------------------------------
// compatibility_table
// ----------------------------------------------------------------------------

/*!\brief Declaration of algorithm specific compatibility table.
 * \ingroup algorithm
 * \tparam algorithm_id_type The type of the algorithm specific id. Algorithm configurations must maintain
 *                           this table to allow validation checks.
 */
template <>
inline constexpr std::array<std::array<bool, static_cast<uint8_t>(io_config_id::SIZE)>,
                            static_cast<uint8_t>(io_config_id::SIZE)> compatibility_table<io_config_id>
{
    {
        { 0, 1 , 1, 1}, //  0: fields
        { 1, 0 , 1, 1}, //  1: formats
        { 1, 1 , 0, 1},  //  2: traits
        { 1, 1 , 1, 0}  //  2: ref_info
    }
};

} // namespace seqan3::detail
