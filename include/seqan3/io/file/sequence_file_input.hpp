// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::sequence_file_input and corresponding traits classes.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <string>
#include <vector>

#include <seqan3/io/file/input_file_base.hpp>
#include <seqan3/io/format/format_fasta.hpp>


namespace seqan3
{

struct sequence_file_input_default_traits_dna
{
    /*!\name Member types
     * \brief Definitions to satisfy seqan3::sequence_file_input_traits.
     * \{
     */

    using types = type_list<std::string, std::vector<phred42>, std::vector<dna5>, std::vector<dna5q>>;
    using ids   = fields<     field::id,          field::qual,        field::seq,    field::seq_qual>;

    //!\}
};

struct sequence_file_input_default_traits_aa
{
    /*!\name Member types
     * \brief Definitions to satisfy seqan3::sequence_file_input_traits.
     * \{
     */

    using types = type_list<std::string, std::vector<phred42>, std::vector<aa27>>;
    using ids   = fields<     field::id,          field::qual,        field::seq>;

    //!\}
};

} // namespace seqan3

namespace seqan3::io_cfg
{

inline constexpr configuration sequence_file_default_configuration =
    select_fields<field::id, field::qual, field::seq>       |
    select_traits<sequence_file_input_default_traits_dna>  |
    select_formats<format_fasta>;

} // namespace seqan3::io_cfg


namespace seqan3
{
// ----------------------------------------------------------------------------
// sequence_file_input
// ----------------------------------------------------------------------------

/*!\brief A class for reading sequence files, e.g. FASTA, FASTQ ...
 * \ingroup sequence
 * \tparam traits_type          An auxiliary type that defines certain member types and constants, must satisfy
 *                              seqan3::sequence_file_input_traits.
 * \tparam selected_field_ids   A seqan3::fields type with the list and order of desired record entries; all fields
 *                              must be in seqan3::sequence_file_input::field_ids.
 * \tparam valid_formats        A seqan3::type_list of the selectable formats (each must meet
 *                              seqan3::sequence_file_input_format).
 *
 * \details
 *
 * ### Introduction
 *
 * Sequence files are the most generic and common biological files. Well-known formats include
 * FastA and FastQ, but some may also be interested in treating SAM or BAM files as sequence
 * files, discarding the alignment.
 *
 * The Sequence file abstraction supports reading four different fields:
 *
 *   1. seqan3::field::seq
 *   2. seqan3::field::id
 *   3. seqan3::field::qual
 *   4. seqan3::field::seq_qual (sequence and qualities in one range)
 *
 * The first three fields are retrieved by default (and in that order). The last field may be selected to have
 * sequence and qualities directly stored in a more memory-efficient combined container. If you select the last
 * field you may not select seqan3::field::seq or seqan3::field::qual.
 *
 * ### Construction and specialisation
 *
 * This class comes with two constructors, one for construction from a file name and one for construction from
 * an existing stream and a known format. The first one automatically picks the format based on the extension
 * of the file name. The second can be used if you have a non-file stream, like std::cin or std::istringstream,
 * that you want to read from and/or if you cannot use file-extension based detection, but know that your input
 * file has a certain format.
 *
 * In most cases the template parameters are deduced completely automatically:
 *
 * \include test/snippet/io/sequence_file/sequence_file_input_template_deduction.cpp
 * Reading from an std::istringstream:
 * \include test/snippet/io/sequence_file/sequence_file_input_istringstream.cpp
 *
 * Note that this is not the same as writing `sequence_file_input<>` (with angle brackets). In the latter case they are
 * explicitly set to their default values, in the former case
 * [automatic deduction](https://en.cppreference.com/w/cpp/language/class_template_argument_deduction) happens which
 * chooses different parameters depending on the constructor arguments. For opening from file, `sequence_file_input<>`
 * would have also worked, but for opening from stream it would not have.
 *
 * In some cases, you do need to specify the arguments, e.g. if you want to read amino acids:
 *
 * \include test/snippet/io/sequence_file/sequence_file_input_aminoacid.cpp
 *
 * You can define your own traits type to further customise the types used by and returned by this class, see
 * seqan3::sequence_file_default_traits_dna for more details. As mentioned above, specifying at least one
 * template parameter yourself means that you loose automatic deduction so if you want to read amino acids **and**
 * want to read from a string stream you need to give all types yourself:
 *
 * \include test/snippet/io/sequence_file/sequence_file_input_template_specification.cpp
 *
 * ### Reading record-wise
 *
 * You can iterate over this file record-wise:
 *
 * \include test/snippet/io/sequence_file/sequence_file_input_record_iter.cpp
 *
 * In the above example, rec has the type \ref record_type which is a specialisation of seqan3::record and behaves
 * like an std::tuple (that's why we can access it via get). Instead of using the seqan3::field based interface on
 * the record, you could also use `std::get<0>` or even `std::get<dna4_vector>` to retrieve the sequence, but it is
 * not recommended, because it is more error-prone.
 *
 * *Note:* It is important to write `auto &` and not just `auto`, otherwise you will copy the record on every iteration.
 * Since the buffer gets "refilled" on every iteration, you can also move the data out of the record if you want
 * to store it somewhere without copying:
 *
 * \include test/snippet/io/sequence_file/sequence_file_input_auto_ref.cpp
 *
 * ### Reading record-wise (decomposed records)
 *
 * Instead of using `get` on the record, you can also use
 * [structured bindings](https://en.cppreference.com/w/cpp/language/structured_binding)
 * to decompose the record into its elements:
 *
 * \include test/snippet/io/sequence_file/sequence_file_input_decomposed.cpp
 *
 * In this case you immediately get the two elements of the tuple: `seq` of \ref sequence_type and `id` of
 * \ref id_type. **But beware: with structured bindings you do need to get the order of elements correctly!**
 *
 * ### Reading record-wise (custom fields)
 *
 * If you want to skip specific fields from the record you can pass a non-empty fields trait object to the
 * sequence_file_input constructor to select the fields that should be read from the input. For example to choose a
 * combined field for SEQ and QUAL (see above). Or to never actually read the QUAL, if you don't need it.
 * The following snippets demonstrate the usage of such a fields trait object.
 *
 * \include test/snippet/io/sequence_file/sequence_file_input_custom_fields.cpp
 *
 * When reading a file, all fields not present in the file (but requested implicitly or via the `selected_field_ids`
 * parameter) are ignored.
 *
 * ### Views on files
 *
 * Since SeqAn files are ranges, you can also create views over files. A useful example is to filter the records
 * based on certain criteria, e.g. minimum length of the sequence field:
 *
 * \include test/snippet/io/sequence_file/sequence_file_input_file_view.cpp
 *
 * ### End of file
 *
 * You can check whether a file is at end by comparing begin() and end() (if they are the same, the file is at end).
 *
 * ### Formats
 *
 * We currently support reading the following formats:
 *   * seqan3::format_fasta
 *   * seqan3::format_fastq
 *   * seqan3::format_embl
 *   * seqan3::format_genbank
 *   * seqan3::format_sam
 */
template <typename config_t = std::remove_cvref_t<decltype(io_cfg::sequence_file2_default_configuration)>>
class sequence_file_input : public input_file_base<config_t>
{
    using base_t = input_file_base<config_t>;
public:

    // need these for CTAD
    sequence_file_input(std::filesystem::path filename, config_t const & cfg = config_t{}) :
        base_t{filename, cfg}
    {}

    sequence_file_input(std::istream & stream, config_t const & cfg) :
        base_t{stream, cfg}
    {}

    sequence_file_input(std::istream && stream, config_t const & cfg) :
        base_t{std::move(stream), cfg}
    {}
};

} // namespace seqan3
