// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the seqan3::var_io::tag_dictionary class and auxiliaries.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <deque>
#include <map>
#include <variant>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/core/concept/core_language.hpp>
#include <seqan3/core/detail/template_inspection.hpp>
#include <seqan3/io/record.hpp>
#include <seqan3/io/utility.hpp>
#include <seqan3/std/concepts>
#include <seqan3/utility/char_operations/predicate.hpp>

namespace seqan3::var_io
{

//!\brief An enumerator denoting variant file special states.
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
 *  2. A std::vector<seqan3::dna5> if a simple character or sequence of DNA.
 *  3. A std::string if they are anything else (imprecise structural variant, breakpoint-string etc).
 */
using allele = std::variant<special_value, std::vector<dna5>, std::string>;

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



//!\brief Stores the header information of alignment files.
//!\ingroup variant_io
struct header
{
    std::string file_format;         //!< The file format version. Note: this is overwritten by our formats on output.

    struct contig_t
    {
        std::string                        id;
        int64_t                            length = -1;
        std::map<std::string, std::string> other_fields;
    };

    std::deque<contig_t> contigs;
    std::unordered_map<std::string_view, size_t> contig_name_to_index;

    struct info_t
    {
        std::string id;
        uint32_t    number; // TODO enum with special values
        uint32_t    type;   // TODO enum with special values
        std::string description;
        std::map<std::string, std::string> other_fields;
    };

    std::deque<info_t> infos;
    std::unordered_map<std::string_view, size_t> info_name_to_index;

    struct filter_t
    {
        std::string id;
        std::string description;
        std::map<std::string, std::string> other_fields;
    };

    std::deque<info_t> filters;
    std::unordered_map<std::string_view, size_t> filter_name_to_index;

    struct format_t
    {
        std::string id;
        uint32_t    number; // TODO enum with special values
        uint32_t    type;   // TODO enum with special values
        std::string description;
        std::map<std::string, std::string> other_fields;
    };

    std::deque<info_t> formats;
    std::unordered_map<std::string_view, size_t> format_name_to_index;

    std::vector<std::string> other_lines;         //!< The list of comments.
    //TODO read group
};

} // namespace seqan3::var_io
