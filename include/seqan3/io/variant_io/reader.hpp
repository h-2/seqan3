// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::reader and corresponding traits classes.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <string>
#include <vector>

#include <seqan3/alphabet/cigar/all.hpp>
#include <seqan3/alphabet/nucleotide/sam_dna16.hpp>
#include <seqan3/alphabet/quality/phred63.hpp>
#include <seqan3/io/variant_io/header.hpp>
#include <seqan3/io/variant_io/misc.hpp>
#include <seqan3/io/detail/reader_base.hpp>
#include <seqan3/io/format/format_vcf.hpp>
#include <seqan3/io/format/format_vcf_input_handler.hpp>
#include <seqan3/io/utility.hpp>


namespace seqan3::var_io
{

/*!\brief Shallow field types for variant io.
 *!\ingroup variant_io
 *
 * \details
 *
 * These traits define a record type with minimal memory allocations. It is the recommended record type when
 * streaming over the input file.
 *
 * Where possible numeric IDs are used (BCF-style).
 *
 * \warning Shallow types
 *
 * These records are not self-contained, i.e. they depend on caches and will become invalid when the file
 * is moved to the next record.
 * Since some elements in the record are views, it may not be trivial to change all values.
 */
inline constexpr auto shallow_field_types = type_tag<int32_t,                             // field::chrom,
                                                     int32_t,                             // field::pos,
                                                     std::string_view,                    // field::id,
                                                     decltype(std::string_view{} | views::char_to<seqan3::dna5>),// field::ref,
                                                     std::vector<std::string_view>,       // field::alt,
                                                     float,                               // field::qual,
                                                     std::vector<int32_t>,                // field::filter,
                                                     std::vector<info_element<true>>,     // field::info,
                                                     std::vector<genotype_element<true>>, // field::genotypes,
//                                                      std::string_view, // field::genotypes,
                                                     header const *>;                     // field::header>;


/*!\brief Deep field types for variant io.
 *!\ingroup variant_io
 *
 * \details
 *
 * TODO
 */
inline constexpr auto    deep_field_types = type_tag<int32_t,                             // field::chrom,
                                                     int32_t,                             // field::pos,
                                                     std::string,                         // field::id,
                                                     std::vector<seqan3::dna5>,           // field::ref,
                                                     std::vector<std::string>,            // field::alt,
                                                     float,                               // field::qual,
                                                     std::vector<int32_t>,                // field::filter,
                                                     std::vector<info_element<false>>,    // field::info,
                                                     std::vector<genotype_element<false>>,// field::genotypes,
                                                     header const *>;

//!\brief Every field is configured as a std::span of std::byte (this enables "raw" io).
//!\ingroup variant_io
inline constexpr auto raw_field_types = list_traits::repeat<default_field_ids.size, std::span<std::byte>>{};

//!\brief Every field is configured as a std::string.
//!\ingroup variant_io
inline constexpr auto string_field_types = list_traits::repeat<default_field_ids.size, std::string>{};



/*!\brief Options that can be used to configure the behaviour of seqan3::var_io::reader.
 * \tparam field_ids_t   Type of the field_ids member (usually deduced).
 * \tparam field_types_t Type of the field_types member (usually deduced).
 * \tparam formats_t     Type of the formats member (usually deduced).
 * \ingroup variant_io
 *
 * \details
 *
 * TODO describe how to easily initialise this
 */
template <typename field_ids_t = decltype(default_field_ids),
          typename field_types_t = decltype(shallow_field_types),
          typename formats_t = type_list<format_vcf>>
struct reader_options
{
    //!\brief The fields that shall be contained in each record; a seqan3::tag over seqan3::field.
    field_ids_t field_ids{};

    /*!\brief The types corresponding to each field; a seqan3::type_tag over the types.
     *
     * \details
     *
     * See seqan3::var_io::reader for an overview of the supported field/type combinations.
     */
    field_types_t field_types{};

    /*!\brief The formats that input files can take; a seqan3::type_tag over the types.
     *
     * \details
     *
     * See seqan3::var_io::reader for an overview of the the supported formats.
     */
    formats_t formats{};

    //!\brief Options that are passed on to the internal stream oject.
    transparent_istream_options stream_options{};

    //TODO static_assert

};

// ----------------------------------------------------------------------------
// reader
// ----------------------------------------------------------------------------

/*!\brief A class for reading alignment map files, e.g. SAM, BAM, CRAM.
 * \tparam options_t A specialisation of seqan3::var_io::reader_options.
 * \ingroup variant_io
 *
 * \details
 *
 * TODO
 */
template <typename options_t = reader_options<>>
class reader : public reader_base<options_t>
{
    using base_t = reader_base<options_t>;
    using typename base_t::format_type;
public:

    // need these for CTAD
    explicit reader(std::filesystem::path const & filename, options_t const & opt = options_t{}) :
        base_t{filename, opt}
    {}

    reader(std::filesystem::path const & filename, format_type const & fmt, options_t const & opt = options_t{}) :
        base_t{filename, fmt, opt}
    {}

    reader(std::istream & stream, format_type const & fmt, options_t const & opt) :
        base_t{stream, fmt, opt}
    {}

    reader(std::istream && stream, format_type const & fmt, options_t const & opt) :
        base_t{std::move(stream), fmt, opt}
    {}
};

} // namespace seqan3
