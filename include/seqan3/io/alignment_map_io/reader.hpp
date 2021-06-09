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
#include <seqan3/io/alignment_map_io/misc.hpp>
#include <seqan3/io/detail/reader_base.hpp>
#include <seqan3/io/format/format_sam.hpp>
#include <seqan3/io/format/format_sam_input_handler.hpp>
#include <seqan3/io/utility.hpp>


namespace seqan3::am_io
{

//!\brief The default types corresponding to seqan3::am_io::default_field_ids.
//!\ingroup alignment_map_io
inline constexpr auto default_field_types = type_tag<std::string,
                                                     flag,
                                                     int32_t,
                                                     int32_t,
                                                     int8_t,
                                                     std::vector<seqan3::cigar>,
                                                     int32_t,
                                                     int32_t,
                                                     int32_t,
                                                     std::vector<seqan3::sam_dna16>,
                                                     std::vector<seqan3::phred63>,
                                                     tag_dictionary,
                                                     header const *>;

//!\brief Every field is configured as a std::span of std::byte (this enables "raw" io).
//!\ingroup alignment_map_io
inline constexpr auto raw_field_types = list_traits::repeat<default_field_ids.size, std::span<std::byte>>{};

//!\brief Every field is configured as a std::string (this enables "raw" io).
//!\ingroup alignment_map_io
inline constexpr auto string_field_types = list_traits::repeat<default_field_ids.size, std::string>{};



/*!\brief Options that can be used to configure the behaviour of seqan3::am_io::reader.
 * \tparam field_ids_t   Type of the field_ids member (usually deduced).
 * \tparam field_types_t Type of the field_types member (usually deduced).
 * \tparam formats_t     Type of the formats member (usually deduced).
 * \ingroup alignment_map_io
 *
 * \details
 *
 * TODO describe how to easily initialise this
 */
template <typename field_ids_t = decltype(default_field_ids),
          typename field_types_t = decltype(default_field_types),
          typename formats_t = type_list<format_sam>>
struct reader_options
{
    //!\brief The fields that shall be contained in each record; a seqan3::tag over seqan3::field.
    field_ids_t field_ids{};

    /*!\brief The types corresponding to each field; a seqan3::type_tag over the types.
     *
     * \details
     *
     * See seqan3::am_io::reader for an overview of the supported field/type combinations.
     */
    field_types_t field_types{};

    /*!\brief The formats that input files can take; a seqan3::type_tag over the types.
     *
     * \details
     *
     * See seqan3::am_io::reader for an overview of the the supported formats.
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
 * \tparam options_t A specialisation of seqan3::am_io::reader_options.
 * \ingroup alignment_map_io
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
