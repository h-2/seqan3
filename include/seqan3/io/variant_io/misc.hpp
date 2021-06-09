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

#include <seqan3/core/concept/core_language.hpp>
#include <seqan3/core/detail/template_inspection.hpp>
#include <seqan3/io/record.hpp>
#include <seqan3/io/utility.hpp>
#include <seqan3/range/container/small_string.hpp>
#include <seqan3/std/concepts>
#include <seqan3/utility/char_operations/predicate.hpp>

namespace seqan3::var_io
{

enum class special_value
{
    missing,    //!< "."
    unknown,    //!< "*"
};

using allele = std::variant<special_value,
                            std::vector<dna5>,       // regular allele
                            std::string>;            // anything else, e.g. breakpoint-string




//!\brief Default fields for seqan3::var_io::reader_options.
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
    //!\brief Stores information of the program/tool that was used to create the file.
    struct program_info_t
    {
        std::string id;                //!< A unique (file scope) id.
        std::string name;              //!< The official name.
        std::string command_line_call; //!< The command line call that produces the file.
        std::string previous;          //!< The id of the previous program if program calls were chained.
        std::string description;       //!< A description of the program and/or program call.
        std::string version;           //!< The program/tool version.
    };

    std::string format_version; //!< The file format version. Note: this is overwritten by our formats on output.
    std::string sorting;        //!< The sorting of the file. SAM: [unknown, unsorted, queryname, coordinate].
    std::string subsorting;     //!< The sub-sorting of the file. SAM: [unknown, unsorted, queryname, coordinate](:[A-Za-z0-9_-]+)+.
    std::string grouping;       //!< The grouping of the file. SAM: [none, query, reference].

    std::vector<program_info_t> program_infos; //!< The list of program information.

    std::vector<std::string> comments;         //!< The list of comments.

    std::deque<std::string> ref_names;
    std::vector<uint32_t> ref_lengths;

    std::unordered_map<std::string_view, int32_t> ref_name2ref_id;

    //TODO read group
};

} // namespace seqan3::var_io
