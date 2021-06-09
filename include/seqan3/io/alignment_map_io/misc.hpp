// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the seqan3::am_io::tag_dictionary class and auxiliaries.
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

namespace seqan3::am_io
{

//!\brief Default fields for seqan3::am_io::reader_options.
inline constexpr auto default_field_ids = tag<field::qname,
                                              field::flag,
                                              field::ref_id,
                                              field::pos,
                                              field::mapq,
                                              field::cigar,
                                              field::next_ref_id,
                                              field::next_pos,
                                              field::tlen,
                                              field::seq,
                                              field::qual,
                                              field::optionals,
                                              field::header>;


//!\brief Stores the header information of alignment files.
//!\ingroup alignment_map_io
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


/*!\brief An enum flag that describes the properties of an aligned read (given as a SAM record).
 * \ingroup alignment_map_io
 *
 * The SAM flag is a bitwise flags, which means that each value corresponds to a specific bit that is set and that they
 * can be combined and tested using binary operations.
 * See this [tutorial](https://www.codeproject.com/Articles/13740/The-Beginner-s-Guide-to-Using-Enum-Flags) for an
 * introduction on bitwise operations on enum flags.
 *
 * Example:
 *
 * \include test/snippet/io/alignment_file/flags.cpp
 *
 * Adapted from the [SAM specifications](https://samtools.github.io/hts-specs/SAMv1.pdf) are the following additional
 * information to some flag values:
 * * For each read/contig in a SAM file, it is required that one and only one line associated with the read
 *   has neither the seqan3::am_io::flag::secondary_alignment nor the seqan3::am_io::flag::supplementary_alignment flag value
 *   set (satisfies `FLAG & 0x900 == 0 `). This line is called the **primary alignment** of the read.
 * * seqan3::am_io::flag::secondary_alignment (bit `0x100`) marks the alignment not to be used in certain analyses when
 *   the tools in use are aware of this bit. It is typically used to flag alternative mappings when multiple mappings
 *   are presented in a SAM.
 * * seqan3::am_io::flag::supplementary_alignment (bit `0x800`) indicates that the corresponding alignment line is part
 *   of a chimeric alignment. If the SAM/BAM file corresponds to long reads (nanopore/pacbio) this happens when
 *   reads are split before being aligned and the best matching part is marked as primary, while all other aligned
 *   parts are marked supplementary.
 * * seqan3::am_io::flag::unmapped (bit `0x4`) is the only reliable place to tell whether the read is unmapped.
 *   If seqan3::am_io::flag::unmapped is set, no assumptions can be made about RNAME, POS, CIGAR, MAPQ, and
 *   seqan3::am_io::flag::proper_pair, seqan3::am_io::flag::secondary_alignment, and seqan3::am_io::flag::supplementary_alignment
 *   (bits `0x2`, `0x100`, and `0x800`).
 * * seqan3::am_io::flag::on_reverse_strand (bit `0x10`) indicates  whether the read sequence has been reverse complemented
 *   and the quality string is reversed.  When bit seqan3::am_io::flag::unmapped (`0x4`) is unset, this
 *   corresponds to the strand to which the segment has been mapped: seqan3::am_io::flag::on_reverse_strand (bit `0x10`)
 *   unset indicates the forward strand, while set indicates the reverse strand. When seqan3::am_io::flag::unmapped (`0x4`)
 *   is set, this indicates whether the unmapped read is stored in its original orientation as it came off the
 *   sequencing machine.
 * * seqan3::am_io::flag::first_in_pair and seqan3::am_io::flag::second_in_pair (bits `0x40` and `0x80`) reflect the read
 *   ordering within each template inherent in the sequencing technology used. If seqan3::am_io::flag::first_in_pair and
 *   seqan3::am_io::flag::second_in_pair (`0x40` and `0x80`) are both set, the read is part of a linear template, but it
 *   is neither the first nor the last read. If both are unset, the index of the read in the template is unknown.
 *   This may happen for a non-linear template or when this information is lost during data processing.
 * * If seqan3::am_io::flag::paired (bit `0x1`) is unset, no assumptions can be made about seqan3::am_io::flag::proper_pair,
 *   seqan3::am_io::flag::mate_unmapped, seqan3::am_io::flag::mate_on_reverse_strand, seqan3::am_io::flag::first_in_pair and
 *   seqan3::am_io::flag::second_in_pair (bits `0x2`, `0x8`, `0x20`, `0x40` and `0x80`).
 *
 * \sa https://broadinstitute.github.io/picard/explain-flags.html
 */
enum class flag : uint16_t
{
   none                    = 0,     //!< None of the flags below are set.
   paired                  = 0x1,   //!< The aligned read is paired (paired-end sequencing).
   proper_pair             = 0x2,   //!< The two aligned reads in a pair have a proper distance between each other.
   unmapped                = 0x4,   //!< The read is not mapped to a reference (unaligned).
   mate_unmapped           = 0x8,   //!< The mate of this read is not mapped to a reference (unaligned).
   on_reverse_strand       = 0x10,  //!< The read sequence has been reverse complemented before being mapped (aligned).
   mate_on_reverse_strand  = 0x20,  //!< The mate sequence has been reverse complemented before being mapped (aligned).
   first_in_pair           = 0x40,  //!< Indicates the ordering (see details in the seqan3::flag description).
   second_in_pair          = 0x80,  //!< Indicates the ordering (see details in the seqan3::flag description).
   secondary_alignment     = 0x100, //!< This read alignment is an alternative (possibly suboptimal) to the primary.
   failed_filter           = 0x200, //!< The read alignment failed a filter, e.g. quality controls.
   duplicate               = 0x400, //!< The read is marked as a PCR duplicate or optical duplicate.
   supplementary_alignment = 0x800  //!< This sequence is part of a split alignment and is not the primary alignment.
};


/*!\brief Overload for the seqan3::flags.
 * \tparam char_t Type char type of the debug_stream.
 * \param stream The seqan3::debug_stream.
 * \param flag The flag to print.
 * \relates seqan3::debug_stream_type
 */
template <typename char_t>
inline debug_stream_type<char_t> & operator<<(debug_stream_type<char_t> & stream, flag const flag)
{
    return stream << static_cast<int16_t>(flag);
}

} // namespace seqan3::am_io

namespace seqan3
{

//!\brief Enables bitwise operations for seqan3::flags.
template <>
constexpr bool add_enum_bitwise_operators<seqan3::am_io::flag> = true;

} // namespace seqan3

namespace seqan3::literals
{

/*!\brief The SAM tag literal, such that tags can be used in constant expressions.
 * \ingroup alignment_map_io
 * \tparam char_t The char type. Usually `char`. Parameter pack `...s` must be of
 *                length 2, since SAM tags consist of two letters (char0 and char1).
 * \returns The unique identifier of the SAM tag computed by char0 * 128 + char1.
 *
 * \details
 *
 * A SAM tag consists of two letters, initialized via the string literal ""_sam_tag,
 * which delegate to its unique id.
 * e.g.
 *
 * \snippet test/snippet/io/alignment_file/tag_dictionary/tag_dictionary.cpp tag
 *
 * The purpose of those tags is to fill or query the seqan3::am_io::tag_dictionary
 * for a specific key (tag_id) and retrieve the corresponding value.
 *
 * \sa seqan3::am_io::tag_dictionary
 */

#ifdef __cpp_nontype_template_parameter_class
template <small_string<2> str> // TODO: better handling if too large string is provided?
constexpr uint16_t operator""_sam_tag()
{
#else // GCC/Clang extension
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
template <typename char_t, char_t ...s>
constexpr uint16_t operator""_sam_tag()
{
    static_assert(std::same_as<char_t, char>, "Illegal SAM tag: Type must be char.");
    constexpr small_string<sizeof...(s)> str{std::array<char, sizeof...(s)>{s...}};
#pragma GCC diagnostic pop
#endif

    static_assert(str.size() == 2, "Illegal SAM tag: Exactly two characters must be given.");

    char constexpr char0 = str[0];
    char constexpr char1 = str[1];

    static_assert((is_alpha(char0) && is_alnum(char1)),
                  "Illegal SAM tag: a SAM tag must match /[A-Za-z][A-Za-z0-9]/.");

    return static_cast<uint16_t>(char0) * 256 + static_cast<uint16_t>(char1);
}

} // namespae seqan3::literals

namespace seqan3::am_io
{

//!\brief std::variant of allowed types for optional tag fields of the SAM format.
//!\ingroup alignment_map_io
using sam_tag_variant = std::variant<char, int32_t, float, std::string,
                                        std::vector<int8_t>, std::vector<uint8_t>,
                                        std::vector<int16_t>, std::vector<uint16_t>,
                                        std::vector<int32_t>, std::vector<uint32_t>,
                                        std::vector<float>>;

//!\brief Each SAM tag type char identifier. Index corresponds to the seqan3::am_io::sam_tag_variant types.
// char constexpr sam_tag_trait_char[11]       = {'A',  'i',  'f',  'Z',  'B', 'B', 'B', 'B', 'B', 'B', 'B'};
//!\brief Each types SAM tag type extra char id. Index corresponds to the seqan3::am_io::sam_tag_variant types.
// char constexpr sam_tag_trait_char_extra[11] = {'\0', '\0', '\0', '\0', 'c', 'C', 's', 'S', 'i', 'I', 'f'};

/*!\brief The generic base class.
 * \ingroup alignment_map_io
 *
 * \attention This is a pure base class that needs to be specialized in order to
 *            be used.
 *
 * ### How to specialize the type for your custom tag
 *
 * All known tags of the SAM specifications already have a pre-defined type.
 * If you want to specify the type of your custom tag (the SAM specifications
 * recommend to use X?, Y? or Z?) you need to overload the seqan3::sam_tag_trait
 * struct in the following way: (take tag "XX" as an example)
 *
 * \snippet test/snippet/io/alignment_file/tag_dictionary/tag_dictionary.cpp type_overload
 *
 * Everything else, like the get and set functions and correct SAM output
 * (XX:i:? in this case) is handled by the seqan3::am_io::tag_dictionary.
 *
 * The seqan3::sam_tag_trait is overloaded the following SAM tags:
 *
 * | Tag Name | SeqAn Type Implementation |
 * | -------- | --------------------- |
 * | "AM"_sam_tag | int32_t               |
 * | "AS"_sam_tag | int32_t               |
 * | "BC"_sam_tag | std::string           |
 * | "BQ"_sam_tag | std::string           |
 * | "BZ"_sam_tag | std::string           |
 * | "CB"_sam_tag | std::string           |
 * | "CC"_sam_tag | std::string           |
 * | "CG"_sam_tag | std::vector<int32_t>  |
 * | "CM"_sam_tag | int32_t               |
 * | "CO"_sam_tag | std::string           |
 * | "CP"_sam_tag | int32_t               |
 * | "CQ"_sam_tag | std::string           |
 * | "CR"_sam_tag | std::string           |
 * | "CS"_sam_tag | std::string           |
 * | "CT"_sam_tag | std::string           |
 * | "CY"_sam_tag | std::string           |
 * | "E2"_sam_tag | std::string           |
 * | "FI"_sam_tag | int32_t               |
 * | "FS"_sam_tag | std::string           |
 * | "FZ"_sam_tag | std::vector<uint16_t> |
 * | "H0"_sam_tag | int32_t               |
 * | "H1"_sam_tag | int32_t               |
 * | "H2"_sam_tag | int32_t               |
 * | "HI"_sam_tag | int32_t               |
 * | "IH"_sam_tag | int32_t               |
 * | "LB"_sam_tag | std::string           |
 * | "MC"_sam_tag | std::string           |
 * | "MD"_sam_tag | std::string           |
 * | "MI"_sam_tag | std::string           |
 * | "MQ"_sam_tag | int32_t               |
 * | "NH"_sam_tag | int32_t               |
 * | "NM"_sam_tag | int32_t               |
 * | "OC"_sam_tag | std::string           |
 * | "OP"_sam_tag | int32_t               |
 * | "OQ"_sam_tag | std::string           |
 * | "OX"_sam_tag | std::string           |
 * | "PG"_sam_tag | std::string           |
 * | "PQ"_sam_tag | int32_t               |
 * | "PT"_sam_tag | std::string           |
 * | "PU"_sam_tag | std::string           |
 * | "Q2"_sam_tag | std::string           |
 * | "QT"_sam_tag | std::string           |
 * | "QX"_sam_tag | std::string           |
 * | "R2"_sam_tag | std::string           |
 * | "RG"_sam_tag | std::string           |
 * | "RT"_sam_tag | std::string           |
 * | "RX"_sam_tag | std::string           |
 * | "SA"_sam_tag | std::string           |
 * | "SM"_sam_tag | int32_t               |
 * | "TC"_sam_tag | int32_t               |
 * | "U2"_sam_tag | std::string           |
 * | "UQ"_sam_tag | int32_t               |
 *
 */
template <uint16_t tag_value>
struct sam_tag_trait
{
    //!\brief The type for all unknown tags with no extra overload defaults to an std::variant.
    using type = sam_tag_variant;
};

//!\brief Short cut helper for seqan3::sam_tag_trait::type.
//!\ingroup alignment_map_io
template <uint16_t tag_value>
using sam_tag_trait_t = typename sam_tag_trait<tag_value>::type;

//!\cond
template <> struct sam_tag_trait<"AM"_sam_tag> { using type = int32_t; };
template <> struct sam_tag_trait<"AS"_sam_tag> { using type = int32_t; };
template <> struct sam_tag_trait<"BC"_sam_tag> { using type = std::string; };
template <> struct sam_tag_trait<"BQ"_sam_tag> { using type = std::string; };
template <> struct sam_tag_trait<"BZ"_sam_tag> { using type = std::string; };
template <> struct sam_tag_trait<"CB"_sam_tag> { using type = std::string; };
template <> struct sam_tag_trait<"CC"_sam_tag> { using type = std::string; };
template <> struct sam_tag_trait<"CG"_sam_tag> { using type = std::vector<int32_t>; };
template <> struct sam_tag_trait<"CM"_sam_tag> { using type = int32_t; };
template <> struct sam_tag_trait<"CO"_sam_tag> { using type = std::string; };
template <> struct sam_tag_trait<"CP"_sam_tag> { using type = int32_t; };
template <> struct sam_tag_trait<"CQ"_sam_tag> { using type = std::string; };
template <> struct sam_tag_trait<"CR"_sam_tag> { using type = std::string; };
template <> struct sam_tag_trait<"CS"_sam_tag> { using type = std::string; };
template <> struct sam_tag_trait<"CT"_sam_tag> { using type = std::string; };
template <> struct sam_tag_trait<"CY"_sam_tag> { using type = std::string; };
template <> struct sam_tag_trait<"E2"_sam_tag> { using type = std::string; };
template <> struct sam_tag_trait<"FI"_sam_tag> { using type = int32_t; };
template <> struct sam_tag_trait<"FS"_sam_tag> { using type = std::string; };
template <> struct sam_tag_trait<"FZ"_sam_tag> { using type = std::vector<uint16_t>; };

// template <> struct sam_tag_trait<"GC"_sam_tag> {};
// template <> struct sam_tag_trait<"GQ"_sam_tag> {};
// template <> struct sam_tag_trait<"GS"_sam_tag> {};

template <> struct sam_tag_trait<"H0"_sam_tag> { using type = int32_t; };
template <> struct sam_tag_trait<"H1"_sam_tag> { using type = int32_t; };
template <> struct sam_tag_trait<"H2"_sam_tag> { using type = int32_t; };
template <> struct sam_tag_trait<"HI"_sam_tag> { using type = int32_t; };
template <> struct sam_tag_trait<"IH"_sam_tag> { using type = int32_t; };
template <> struct sam_tag_trait<"LB"_sam_tag> { using type = std::string; };
template <> struct sam_tag_trait<"MC"_sam_tag> { using type = std::string; };
template <> struct sam_tag_trait<"MD"_sam_tag> { using type = std::string; };

// template <> struct sam_tag_trait<"MF"_sam_tag> {};

template <> struct sam_tag_trait<"MI"_sam_tag> { using type = std::string; };
template <> struct sam_tag_trait<"MQ"_sam_tag> { using type = int32_t; };
template <> struct sam_tag_trait<"NH"_sam_tag> { using type = int32_t; };
template <> struct sam_tag_trait<"NM"_sam_tag> { using type = int32_t; };
template <> struct sam_tag_trait<"OC"_sam_tag> { using type = std::string; };
template <> struct sam_tag_trait<"OP"_sam_tag> { using type = int32_t; };
template <> struct sam_tag_trait<"OQ"_sam_tag> { using type = std::string; };
template <> struct sam_tag_trait<"OX"_sam_tag> { using type = std::string; };
template <> struct sam_tag_trait<"PG"_sam_tag> { using type = std::string; };
template <> struct sam_tag_trait<"PQ"_sam_tag> { using type = int32_t; };
template <> struct sam_tag_trait<"PT"_sam_tag> { using type = std::string; };
template <> struct sam_tag_trait<"PU"_sam_tag> { using type = std::string; };
template <> struct sam_tag_trait<"Q2"_sam_tag> { using type = std::string; };
template <> struct sam_tag_trait<"QT"_sam_tag> { using type = std::string; };
template <> struct sam_tag_trait<"QX"_sam_tag> { using type = std::string; };
template <> struct sam_tag_trait<"R2"_sam_tag> { using type = std::string; };
template <> struct sam_tag_trait<"RG"_sam_tag> { using type = std::string; };
template <> struct sam_tag_trait<"RT"_sam_tag> { using type = std::string; };
template <> struct sam_tag_trait<"RX"_sam_tag> { using type = std::string; };

// template <> struct sam_tag_trait<"S2"_sam_tag> {};

template <> struct sam_tag_trait<"SA"_sam_tag> { using type = std::string; };
template <> struct sam_tag_trait<"SM"_sam_tag> { using type = int32_t; };

// template <> struct sam_tag_trait<"SQ"_sam_tag> {};

template <> struct sam_tag_trait<"TC"_sam_tag> { using type = int32_t; };
template <> struct sam_tag_trait<"U2"_sam_tag> { using type = std::string; };
template <> struct sam_tag_trait<"UQ"_sam_tag> { using type = int32_t; };
//!\endcond

/*!\brief The SAM tag dictionary class that stores all optional SAM fields.
 * \ingroup alignment_map_io
 *
 * \details
 *
 * ### SAM tags
 *
 * A SAM tag consists of two letters, initialized via the string literal ""_sam_tag,
 * which delegates to its unique id (type uint16_t).
 * Example:
 *
 * \snippet test/snippet/io/alignment_file/tag_dictionary/tag_dictionary.cpp tag
 *
 * The purpose of those tags is to fill or query the seqan3::am_io::tag_dictionary
 * for a specific key (tag_id) and retrieve the corresponding value.
 *
 * ### SAM tag types
 *
 * Note that a SAM tag is always associated with a specific type.
 * In the SAM format, the type is indicated in the second argument of the
 * TAG:TYPE:VALUE field. For example "NM:i:3" specifies the NM tag of an integer
 * type with value 3.
 * In seqan3, the types for
 * [known](https://samtools.github.io/hts-specs/SAMtags.pdf) SAM tags
 * are pre-defined by a type trait called seqan3::sam_tag_trait. You can access
 * the type via:
 *
 * \snippet test/snippet/io/alignment_file/tag_dictionary/tag_dictionary.cpp tag_type_t
 *
 * which is the short cut for:
 *
 * \snippet test/snippet/io/alignment_file/tag_dictionary/tag_dictionary.cpp tag_type
 *
 * The following types are allowed by the
 * [SAM specifications](https://samtools.github.io/hts-specs/SAMtags.pdf):
 *
 * |Type | Regexp matching VALUE                  | Description                             | SeqAn Type           |
 * |-----|----------------------------------------|-----------------------------------------|----------------------|
 * | A   | [!-~]                                  |  Printable character                    | char                 |
 * | i   | [-+]?[0-9]+                            |  Signed integer                         | int32_t              |
 * | f   | [-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)? |  Single-precision floating number       | float                |
 * | Z   | [ !-~]*                                |  Printable string, including space      | std::string          |
 * | H   | ([0-9A-F][0-9A-F])*                    |  Byte array in the Hex format           | std::vector<uint8_t> |
 * | B   | [cCsSiIf]\(,[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)+ |  Integer or numeric array | std::vector<T>       |
 *
 * For  an  integer  or  numeric  array  (type  ‘B’),  the second  letter can be
 * one of ‘cCsSiIf’, corresponding to type **T** = int8_t, uint8_t, int16_t,
 * uint16_t, int32_t, uint32_t and float, respectively.
 *
 * ### Using the tag_dictionary
 *
 * The dictionary can be accessed via the functions seqan3::am_io::tag_dictionary::get() and
 * seqan3::am_io::tag_dictionary::set(). Every time the SAM tag you wish to query
 * for must be given as a template argument to the functions.
 *
 * Example:
 *
 * \include test/snippet/io/alignment_file/tag_dictionary/general_usage.cpp
 *
 * \attention You can get any SAM_tag out of the dictionary, even if the tag is
 *            user defined, but note that for unknown tags the return type is an
 *            [std::variant](https://en.cppreference.com/w/cpp/utility/variant).
 *            If you want specify the return type of your custom tag, you need
 *            to overload the seqan3::sam_tag_trait type trait.
 *
 * Unknown Tag Example:
 *
 * \include test/snippet/io/alignment_file/tag_dictionary/unknown_tag.cpp
 *
 * As mentioned before you can either overload the type trait seqan3::sam_tag_trait
 * for the tag "XZ" or learn more about an std::variant at
 * https://en.cppreference.com/w/cpp/utility/variant.
 *
 * \sa seqan3::sam_tag_trait
 * \sa https://en.cppreference.com/w/cpp/utility/variant
 * \sa https://samtools.github.io/hts-specs/SAMv1.pdf
 * \sa https://samtools.github.io/hts-specs/SAMtags.pdf
 */
class tag_dictionary : public std::map<uint16_t, sam_tag_variant>
{
private:
    //!\brief The base type.
    using base_type = std::map<uint16_t, sam_tag_variant>;

public:
    //!\brief The variant type defining all valid SAM tag field types.
    using variant_type = sam_tag_variant;

    /*!\name Getter function for the seqan3::am_io::tag_dictionary.
     *\brief Gets the value of known SAM tags by its correct type instead of the std::variant.
     * \tparam tag The unique tag id of a SAM tag.
     * \returns The value corresponding to the key `tag` of type seqan3::sam_tag_trait<tag>::type.
     *
     * \details
     *
     * See the seqan3::am_io::tag_dictionary detailed documentation below for an example.
     *
     * \attention This function is only available for tags that have an
     *            seqan3::sam_tag_trait<tag>::type overload. See the type trait
     *            documentation for further details.
     * \{
     */

    //!\brief Uses std::map::operator[] for access and default initializes new keys.
    template <uint16_t tag>
    //!\cond
        requires (!std::same_as<sam_tag_trait_t<tag>, variant_type>)
    //!\endcond
    auto & get() &
    {
        if ((*this).count(tag) == 0)
            (*this)[tag] = sam_tag_trait_t<tag>{}; // set correct type if tag is not set yet on

        return std::get<sam_tag_trait_t<tag>>((*this)[tag]);
    }

    //!\brief Uses std::map::operator[] for access and default initializes new keys.
    template <uint16_t tag>
    //!\cond
        requires (!std::same_as<sam_tag_trait_t<tag>, variant_type>)
    //!\endcond
    auto && get() &&
    {
        if ((*this).count(tag) == 0)
            (*this)[tag] = sam_tag_trait_t<tag>{}; // set correct type if tag is not set yet on

        return std::get<sam_tag_trait_t<tag>>(std::move((*this)[tag]));
    }

    //!\brief Uses std::map::at() for access and throws when the key is unknown.
    //!\throws std::out_of_range if map has no key `tag`.
    template <uint16_t tag>
    //!\cond
        requires (!std::same_as<sam_tag_trait_t<tag>, variant_type>)
    //!\endcond
    auto const & get() const &
    {
        return std::get<sam_tag_trait_t<tag>>((*this).at(tag));
    }

    //!\brief Uses std::map::at() for access and throws when the key is unknown.
    //!\throws std::out_of_range if map has no key `tag`.
    template <uint16_t tag>
    //!\cond
        requires (!std::same_as<sam_tag_trait_t<tag>, variant_type>)
    //!\endcond
    auto const && get() const &&
    {
        return std::get<sam_tag_trait_t<tag>>(std::move((*this).at(tag)));
    }
    //!\}
};

} // namespace seqan3::am_io
