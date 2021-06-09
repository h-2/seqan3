// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::seq_io::reader and corresponding traits classes.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <seqan3/io/utility.hpp>


namespace seqan3::seq_io
{

//!\brief Default fields for seqan3::am_io::reader_options.
inline constexpr auto default_field_ids = tag<field::id,
                                              field::seq,
                                              field::qual>;
} // namespace seqan3::seq_io
