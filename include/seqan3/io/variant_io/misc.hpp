// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the seqan3::var_io::tag_dictionary class and auxiliaries.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <string>
#include <tuple>
#include <variant>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/io/record.hpp>
#include <seqan3/io/utility.hpp>

namespace seqan3::var_io
{

//!\brief An enumerator denoting field special states in a variant file.
//!\ingroup variant_io
enum class special_value
{
    missing,    //!< "."
    unknown,    //!< "*"
};

/*!\brief A type representing variant file alleles.
 * \ingroup variant_io
 *
 * \details
 *
 * Alleles in variant files are encoded as one of the following
 *
 *  1. A seqan3::var_io::special_value if missing/absent.
 *  2. A std::vector<seqan3::dna5> if a single character or sequence of DNA.
 *  3. A std::string if they are anything else (imprecise structural variant, breakpoint-string etc).
 */
using allele = std::variant<special_value, std::vector<dna5>, std::string>;

/*!\brief A type representing an variant file INFO field.
 * \ingroup variant_io
 */
using info_entry = std::pair<int32_t, io_type_variant>;

//!\brief Default fields for seqan3::var_io::reader_options.
//!\ingroup variant_io
inline constexpr auto default_field_ids = tag<field::chrom,
                                              field::pos,
                                              field::id,
                                              field::ref,
                                              field::alt,
                                              field::qual,
                                              field::filter,
                                              field::info,
                                              field::genotypes,
                                              field::header>;

} // namespace seqan3::var_io
