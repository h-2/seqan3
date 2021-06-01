// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::input_file_base and corresponding traits classes.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <cassert>
#include <seqan3/std/filesystem>
#include <fstream>
#include <string>
#include <variant>
#include <vector>

#include <seqan3/alphabet/adaptation/char.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/nucleotide/dna15.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/alphabet/quality/qualified.hpp>
#include <seqan3/core/detail/pack_algorithm.hpp>
#include <seqan3/io/exception.hpp>
#include <seqan3/io/record.hpp>
#include <seqan3/io/detail/in_file_iterator.hpp>
#include <seqan3/io/detail/misc_input.hpp>
#include <seqan3/io/detail/record.hpp>
#include <seqan3/io/format/format_fasta.hpp>
#include <seqan3/io/format/format_fasta_input_handler.hpp>
#include <seqan3/io/sequence_file/default_configuration.hpp>
// #include <seqan3/io/sequence_file/input_format_concept.hpp>
// #include <seqan3/io/sequence_file/format_embl.hpp>
// #include <seqan3/io/sequence_file/format_fastq.hpp>
// #include <seqan3/io/sequence_file/format_genbank.hpp>
// #include <seqan3/io/alignment_file/format_sam.hpp>
#include <seqan3/utility/type_list/traits.hpp>


namespace seqan3
{



template <typename config_t>
constexpr bool valid_config()
{
    //TODO implement

//     static_assert([] () constexpr
//                   {
//                       for (field f : selected_field_ids::as_array)
//                           if (!field_ids::contains(f))
//                               return false;
//                       return true;
//                   }(),
//                   "You selected a field that is not valid for sequence files, please refer to the documentation "
//                   "of input_file_base::field_ids for the accepted values.");
//
//     static_assert([] () constexpr
//                   {
//                       return !(selected_field_ids::contains(field::seq_qual) &&
//                                (selected_field_ids::contains(field::seq) ||
//                                (selected_field_ids::contains(field::qual))));
//                   }(),
//                   "You may not select field::seq_qual and either of field::seq and field::qual at the same time.");


//     template <typename t>
// SEQAN3_CONCEPT input_file_base_traits = requires (t v)
// {
//     requires writable_alphabet<typename t::sequence_alphabet>;
//     requires writable_alphabet<typename t::sequence_legal_alphabet>;
//     requires explicitly_convertible_to<typename t::sequence_legal_alphabet, typename t::sequence_alphabet>;
//     requires sequence_container<typename t::template sequence_container<typename t::sequence_alphabet>>;
//
//     requires writable_alphabet<typename t::id_alphabet>;
//     requires sequence_container<typename t::template id_container<typename t::id_alphabet>>;
//
//     requires writable_quality_alphabet<typename t::quality_alphabet>;
//     requires sequence_container<typename t::template quality_container<typename t::quality_alphabet>>;
// };

    return true;
}

// ----------------------------------------------------------------------------
// input_file_base
// ----------------------------------------------------------------------------

template <typename config_t>
class input_file_base
{
    static_assert(valid_config<config_t>(), "Invalid configuration.");

    /*!\name Template arguments
     * \brief Exposed as member types for public access.
     * \{
     */
    //!\brief A traits type that defines aliases and template for storage of the fields.
    using traits_type = typename get_type_t<config_t, io_cfg::select_traits_t>::type;
    //!\brief A seqan3::fields list with the fields selected for the record.
    using selected_field_ids = typename get_type_t<config_t, io_cfg::select_fields_t>::type;
    //!\brief A seqan3::type_list with the possible formats.
    using valid_formats = typename get_type_t<config_t, io_cfg::select_formats_t>::type;
    //!\brief Character type of the stream(s).
    using stream_char_type      = char;
    //!\}
public:


    /*!\name Field types and record type
     * \brief These types are relevant for record/row-based reading; they may be manipulated via the \ref traits_type
     * to achieve different storage behaviour.
     * \{
     */

    //!\brief The type of the record, a specialisation of seqan3::record; acts as a tuple of the selected field types.
    using record_type           = record<detail::select_types_with_ids_t<typename traits_type::types,
                                                                         typename traits_type::ids,
                                                                         selected_field_ids>,
                                         selected_field_ids>;
    //!\}

    /*!\name Range associated types
     * \brief The types necessary to facilitate the behaviour of an input range (used in record-wise reading).
     * \{
     */
    //!\brief The value_type is the \ref record_type.
    using value_type        = record_type;
    //!\brief The reference type.
    using reference         = record_type &;
    //!\brief The const_reference type is void, because files are not const-iterable.
    using const_reference   = void;
    //!\brief An unsigned integer type, usually std::size_t.
    using size_type         = size_t;
    //!\brief A signed integer type, usually std::ptrdiff_t.
    using difference_type   = std::make_signed_t<size_t>;
    //!\brief The iterator type of this view (an input iterator).
    using iterator          = detail::in_file_iterator<input_file_base>;
    //!\brief The const iterator type is void, because files are not const-iterable.
    using const_iterator    = void;
    //!\brief The type returned by end().
    using sentinel          = std::default_sentinel_t;
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    //!\brief Default constructor is explicitly deleted, you need to give a stream or file name.
    input_file_base() = delete;
    //!\brief Copy construction is explicitly deleted, because you can't have multiple access to the same file.
    input_file_base(input_file_base const &) = delete;
    //!\brief Copy assignment is explicitly deleted, because you can't have multiple access to the same file.
    input_file_base & operator=(input_file_base const &) = delete;
    //!\brief Move construction is defaulted.
    input_file_base(input_file_base &&) = default;
    //!\brief Move assignment is defaulted.
    input_file_base & operator=(input_file_base &&) = default;
    //!\brief Destructor is defaulted.
    ~input_file_base() = default;

    /*!\brief Construct from filename.
     * \param[in] filename      Path to the file you wish to open.
     * \param[in] fields_tag    A seqan3::fields tag. [optional]
     * \throws seqan3::file_open_error If the file could not be opened, e.g. non-existant, non-readable, unknown format.
     *
     * \details
     *
     * In addition to the file name, you may specify a custom seqan3::fields type which may be easier than
     * defining all the template parameters.
     *
     * ### Decompression
     *
     * This constructor transparently applies a decompression stream on top of the file stream in case
     * the file is detected as being compressed.
     * See the section on \link io_compression compression and decompression \endlink for more information.
     */
    input_file_base(std::filesystem::path filename,
                    config_t const & cfg = config_t{}) :
        primary_stream{new std::ifstream{}, stream_deleter_default},
        config{cfg}
    {
        primary_stream->rdbuf()->pubsetbuf(stream_buffer.data(), stream_buffer.size());
        static_cast<std::basic_ifstream<char> *>(primary_stream.get())->open(filename,
                                                                             std::ios_base::in | std::ios::binary);

        if (!primary_stream->good())
            throw file_open_error{"Could not open file " + filename.string() + " for reading."};

        // possibly add intermediate compression stream
        secondary_stream = detail::make_secondary_istream(*primary_stream, filename);

        // initialise format handler or throw if format is not found
        detail::set_format(format, filename);


    }

    /*!\brief Construct from an existing stream and with specified format.
     * \tparam file_format   The format of the file in the stream, must satisfy seqan3::input_file_base_format.
     * \param[in] stream     The stream to operate on; must be derived of std::basic_istream.
     * \param[in] format_tag The file format tag.
     * \param[in] fields_tag A seqan3::fields tag. [optional]
     *
     * \details
     *
     * ### Decompression
     *
     * This constructor transparently applies a decompression stream on top of the stream in case
     * it is detected as being compressed.
     * See the section on \link io_compression compression and decompression \endlink for more information.
     */
    input_file_base(std::istream & stream, config_t const & cfg) :
        primary_stream{&stream, stream_deleter_noop},
        config{cfg}
    {
        static_assert(config_t::template exists<io_cfg::select_formats_t>(), "Select format of stream!");
         // possibly add intermediate compression stream
        secondary_stream = detail::make_secondary_istream(*primary_stream);

        // assert length == 1
        using file_format = list_traits::at<0, valid_formats>;

        format = file_format{};
    }

    //!\overload
    input_file_base(std::istream && stream, config_t const & cfg) :
        primary_stream{new std::istream{std::move(stream)}, stream_deleter_default},
        config{cfg}
    {
        static_assert(config_t::template exists<io_cfg::select_formats_t>(), "Select format of stream!");
         // possibly add intermediate compression stream
        secondary_stream = detail::make_secondary_istream(*primary_stream);

        // assert length == 1
        using file_format = list_traits::at<0, valid_formats>;

        format = file_format{};
    }
    //!\}

    /*!\name Range interface
     * \brief Provides functions for record based reading of the file.
     * \{
     */
    /*!\brief Returns an iterator to current position in the file.
     * \returns An iterator pointing to the current position in the file.
     * \throws seqan3::format_error
     *
     * Equals end() if the file is at end.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * Throws seqan3::format_error if the first record could not be read into the buffer.
     */
    iterator begin()
    {
        // buffer first record
        if (!first_record_was_read)
        {
            // set format-handler
            std::visit([&] (auto f)
                {
                    format_handler = input_format_handler<decltype(f)>{*secondary_stream, config};
                }, format);

            // read first record
            read_next_record();
            first_record_was_read = true;
        }

        return {*this};
    }

    /*!\brief Returns a sentinel for comparison with iterator.
     * \returns Iterator to the first element.
     *
     * This element acts as a placeholder; attempting to dereference it results in undefined behaviour.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    sentinel end() noexcept
    {
        return {};
    }

    /*!\brief Return the record we are currently at in the file.
     * \returns A reference to the currently buffered record.
     *
     * This function returns a reference to the currently buffered record, it is identical to dereferencing begin(),
     * but begin also always points to the current record on single pass input ranges:
     *
     * \include test/snippet/io/sequence_file/input_file_base_return_record.cpp
     *
     * It most situations using the iterator interface or a range-based for-loop are preferable to using front(),
     * because you can only move to the next record via the iterator.
     *
     * In any case, don't forget the reference! If you want to save the data from the record elsewhere, use move:
     *
     * \include test/snippet/io/sequence_file/input_file_base_record_move.cpp
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    reference front() noexcept
    {
        return *begin();
    }
    //!\}

protected:
    //!\privatesection

    config_t config;

    /*!\name Data buffers
     * \{
     */
    //!\brief Buffer for a single record.
    record_type record_buffer;
    //!\brief A larger (compared to stl default) stream buffer to use when reading from a file.
    std::vector<char> stream_buffer{std::vector<char>(1'000'000)};
    //!\}

    /*!\name Stream / file access
     * \{
     */
    //!\brief The type of the internal stream pointers. Allows dynamically setting ownership management.
    using stream_ptr_t = std::unique_ptr<std::basic_istream<stream_char_type>,
                                         std::function<void(std::basic_istream<stream_char_type>*)>>;
    //!\brief Stream deleter that does nothing (no ownership assumed).
    static void stream_deleter_noop(std::basic_istream<stream_char_type> *) {}
    //!\brief Stream deleter with default behaviour (ownership assumed).
    static void stream_deleter_default(std::basic_istream<stream_char_type> * ptr) { delete ptr; }

    //!\brief The primary stream is the user provided stream or the file stream if constructed from filename.
    stream_ptr_t primary_stream{nullptr, stream_deleter_noop};
    //!\brief The secondary stream is a compression layer on the primary or just points to the primary (no compression).
    stream_ptr_t secondary_stream{nullptr, stream_deleter_noop};

    //!\brief Tracks whether the very first record is buffered when calling begin().
    bool first_record_was_read{false};
    //!\brief File is at position 1 behind the last record.
    bool at_end{false};

    //!\brief Type of the format, an std::variant over the `valid_formats`.
    using format_type           = detail::transfer_template_args_onto_t<valid_formats, std::variant>;
//     template <typename format_t>
//     using format_handler_alias  = input_format_handler<format_t, config_t>;
    using format_handler_type   = detail::transfer_template_args_onto_t<list_traits::transform<input_format_handler,
                                                                                               valid_formats>,
                                                                        std::variant>;

    //!\brief The actual std::variant holding a pointer to the detected/selected format.
    format_type format;
    format_handler_type format_handler;
    //!\}

    //!\brief Tell the format to move to the next record and update the buffer.
    void read_next_record()
    {
        // at end if we could not read further
        if ((std::istreambuf_iterator<stream_char_type>{*secondary_stream} ==
             std::istreambuf_iterator<stream_char_type>{}))
        {
            at_end = true;
            return;
        }

        if (at_end)
            return;

        assert(!format_handler.valueless_by_exception());
        std::visit([&] (auto & f) { f.parse_next_record_into(record_buffer); },
                   format_handler);
    }

    //!\brief Befriend iterator so it can access the buffers.
    friend iterator;
};

/*!\name Type deduction guides
 * \relates seqan3::input_file_base
 * \{
 */

//!\brief Deduces the sequence input file type from the stream and the format.
// template <typename config_type>
// input_file_base(std::filesystem::path, config_type)
//     -> input_file_base<config_type>;
//
// //!\brief Deduces the sequence input file type from the stream, the format and the field ids.
// template <typename config_type>
// input_file_base(std::istream && stream, config_type)
//     -> input_file_base<config_type>;
//
// //!\overload
// template <typename config_type>
// input_file_base(std::istream & stream, config_type)
//     -> input_file_base<config_type>;
//!\}

} // namespace seqan3
