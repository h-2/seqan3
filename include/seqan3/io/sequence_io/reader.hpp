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

#include <string>
#include <vector>

#include <seqan3/io/detail/reader_base.hpp>
#include <seqan3/io/format/format_fasta.hpp>

namespace seqan3::seq_io
{

//!\brief The default types for reading DNA data.
//!\ingroup sequence_io
inline constexpr auto dna_field_types = type_tag<std::string, std::vector<dna15>, std::vector<phred63>>;

//!\brief The default types corresponding to seqan3::am_io::default_field_ids.
//!\ingroup sequence_io
inline constexpr auto protein_field_types = type_tag<std::string, std::vector<aa27>, std::vector<phred63>>;

//!\brief Every field is configured as a std::span of std::byte (this enables "raw" io).
//!\ingroup sequence_io
inline constexpr auto raw_field_types = list_traits::repeat<default_field_ids.size, std::span<std::byte>>{};

//!\brief Every field is configured as a std::string (this enables "raw" io).
//!\ingroup sequence_io
inline constexpr auto string_field_types = list_traits::repeat<default_field_ids.size, std::string>{};

/*!\brief Options that can be used to configure the behaviour of seqan3::am_io::reader.
 * \tparam field_ids_t   Type of the field_ids member (usually deduced).
 * \tparam field_types_t Type of the field_types member (usually deduced).
 * \tparam formats_t     Type of the formats member (usually deduced).
 * \ingroup sequence_io
 *
 * \details
 *
 * By default, the reader options assume DNA data. You select seqan3::seq_io::protein_field_types to
 * read protein data or seqan3::seq_io::string_field_types to store in an agnostic type.
 *
 * TODO describe how to easily initialise this
 */
template <typename field_ids_t = decltype(default_field_ids),
          typename field_types_t = decltype(dna_field_types),
          typename formats_t = type_list<format_fasta>>
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

/*!\brief A class for reading sequence files, e.g. FASTA, FASTQ ...
 * \ingroup sequence
 * \tparam traits_type          An auxiliary type that defines certain member types and constants, must satisfy
 *                              seqan3::seq_io::reader_traits.
 * \tparam selected_field_ids   A seqan3::fields type with the list and order of desired record entries; all fields
 *                              must be in seqan3::seq_io::reader::field_ids.
 * \tparam valid_formats        A seqan3::type_list of the selectable formats (each must meet
 *                              seqan3::seq_io::reader_format).
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
 * \include test/snippet/io/sequence_file/reader_template_deduction.cpp
 * Reading from an std::istringstream:
 * \include test/snippet/io/sequence_file/reader_istringstream.cpp
 *
 * Note that this is not the same as writing `reader<>` (with angle brackets). In the latter case they are
 * explicitly set to their default values, in the former case
 * [automatic deduction](https://en.cppreference.com/w/cpp/language/class_template_argument_deduction) happens which
 * chooses different parameters depending on the constructor arguments. For opening from file, `reader<>`
 * would have also worked, but for opening from stream it would not have.
 *
 * In some cases, you do need to specify the arguments, e.g. if you want to read amino acids:
 *
 * \include test/snippet/io/sequence_file/reader_aminoacid.cpp
 *
 * You can define your own traits type to further customise the types used by and returned by this class, see
 * seqan3::sequence_file_default_traits_dna for more details. As mentioned above, specifying at least one
 * template parameter yourself means that you loose automatic deduction so if you want to read amino acids **and**
 * want to read from a string stream you need to give all types yourself:
 *
 * \include test/snippet/io/sequence_file/reader_template_specification.cpp
 *
 * ### Reading record-wise
 *
 * You can iterate over this file record-wise:
 *
 * \include test/snippet/io/sequence_file/reader_record_iter.cpp
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
 * \include test/snippet/io/sequence_file/reader_auto_ref.cpp
 *
 * ### Reading record-wise (decomposed records)
 *
 * Instead of using `get` on the record, you can also use
 * [structured bindings](https://en.cppreference.com/w/cpp/language/structured_binding)
 * to decompose the record into its elements:
 *
 * \include test/snippet/io/sequence_file/reader_decomposed.cpp
 *
 * In this case you immediately get the two elements of the tuple: `seq` of \ref sequence_type and `id` of
 * \ref id_type. **But beware: with structured bindings you do need to get the order of elements correctly!**
 *
 * ### Reading record-wise (custom fields)
 *
 * If you want to skip specific fields from the record you can pass a non-empty fields trait object to the
 * reader constructor to select the fields that should be read from the input. For example to choose a
 * combined field for SEQ and QUAL (see above). Or to never actually read the QUAL, if you don't need it.
 * The following snippets demonstrate the usage of such a fields trait object.
 *
 * \include test/snippet/io/sequence_file/reader_custom_fields.cpp
 *
 * When reading a file, all fields not present in the file (but requested implicitly or via the `selected_field_ids`
 * parameter) are ignored.
 *
 * ### Views on files
 *
 * Since SeqAn files are ranges, you can also create views over files. A useful example is to filter the records
 * based on certain criteria, e.g. minimum length of the sequence field:
 *
 * \include test/snippet/io/sequence_file/reader_file_view.cpp
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
template <typename options_t = reader_options<>>
class reader : public reader_base<options_t>
{
    using base_t = input_file_base<options_t>;
public:

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
