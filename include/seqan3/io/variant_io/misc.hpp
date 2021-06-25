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
#include <seqan3/core/detail/debug_stream_type.hpp>
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


constexpr bool is_missing(char const c)
{
    return c == '.';
}

constexpr bool is_missing(std::string_view const s)
{
    return s == ".";
}

template <typename int_t>
    requires (std::same_as<int_t, int8_t> || std::same_as<int_t, int16_t> || std::same_as<int_t, int32_t>)
constexpr bool is_missing(int8_t const i)
{
    return i == std::numeric_limits<int_t>::lowest();
}


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

/*!\brief A type representing an variant file INFO field [index of the INFO in header, value].
 * \ingroup variant_io
 */
using info_element = std::pair<int32_t, io_type_variant>;

/*!\brief A type representing an variant file QUAL field [missing or a float value].
 * \ingroup variant_io
 */
using qual = std::variant<special_value, float>;

/*!\brief A type representing an genotype.
 * \ingroup variant_io
 *
 * \details
 *
 * Genotypes / samples are represented as decribed in the BCF specification, i.e. information is grouped by FORMAT
 * identifier, not by sample.
 *
 * This entry consists of the FORMAT index in the file's header and a vector of values. The size of the vector is:
 *
 *   * equal to the number of samples; or
 *   * 0 -- if the field is missing from all samples.
 *
 * The variant vector is guaranteed to be over the type defined in the header. Note that this is a vector over such
 * types (one element per sample!), so seqan3::io_type_id::vector_of_int32 corresponds to std::vector<std::vector<int32_t>>.
 * See seqan3::io_type_vector_variant for more details.
 *
 * If fields are missing from some samples but not others, the vector will have full size but the respective values
 * will be set to the missing value (see seqan3::var_io::is_missing()) or be the empty vector (in case the element type
 * is a vector).
 */
using genotype_element = std::pair<int32_t, io_type_vector_variant>;

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

namespace seqan3
{

template <typename char_t, typename special_value_t>
    requires std::same_as<std::remove_cvref_t<special_value_t>, var_io::special_value>
inline debug_stream_type<char_t> & operator<<(debug_stream_type<char_t> & s, special_value_t && v)
{
    s << (v == var_io::special_value::missing ? '.' : '*');
    return s;
}

} // namespace seqan3
