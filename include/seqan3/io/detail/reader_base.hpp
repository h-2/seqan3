// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::reader_base and corresponding traits classes.
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
// #include <seqan3/io/sequence_file/default_optionsuration.hpp>
// #include <seqan3/io/sequence_file/input_format_concept.hpp>
// #include <seqan3/io/sequence_file/format_embl.hpp>
// #include <seqan3/io/sequence_file/format_fastq.hpp>
// #include <seqan3/io/sequence_file/format_genbank.hpp>
// #include <seqan3/io/alignment_file/format_sam.hpp>
#include <seqan3/io/stream/transparent_istream.hpp>
#include <seqan3/utility/type_list/traits.hpp>


namespace seqan3
{

// ----------------------------------------------------------------------------
// reader_base
// ----------------------------------------------------------------------------


//TODO document format handling
template <typename options_t>
class reader_base
{
public:
    /*!\name Field types and record type
     * \brief These types size_Tare relevant for record/row-based reading; they may be manipulated via the \ref traits_type
     * to achieve different storage behaviour.
     * \{
     */
    //!\brief The type of the record, a specialisation of seqan3::record; acts as a tuple of the selected field types.
    using record_type = record<std::remove_cvref_t<decltype(options_t::field_types)>,
                               std::remove_cvref_t<decltype(options_t::field_ids)>>;
    //!\brief The iterator type of this view (an input iterator).
    using iterator          = detail::in_file_iterator<reader_base>;
    //!\brief The type returned by end().
    using sentinel          = std::default_sentinel_t;
    //!\}

protected:
    //!\privatesection
    /*!\name Format handling
     * \{
     */
    //!\brief A seqan3::type_list with the possible formats.
    using valid_formats = std::remove_cvref_t<decltype(options_t::formats)>;

    //!\brief Type of the format, an std::variant over the `valid_formats`.
    using format_type           = detail::transfer_template_args_onto_t<valid_formats, std::variant>;
    using format_handler_type   = detail::transfer_template_args_onto_t<list_traits::transform<input_format_handler,
                                                                                               valid_formats>,
                                                                        std::variant>;
    //!\}

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    //!\brief Default constructor is explicitly deleted, you need to give a stream or file name.
    reader_base() = delete;
    //!\brief Copy construction is explicitly deleted, because you can't have multiple access to the same file.
    reader_base(reader_base const &) = delete;
    //!\brief Copy assignment is explicitly deleted, because you can't have multiple access to the same file.
    reader_base & operator=(reader_base const &) = delete;
    //!\brief Move construction is defaulted.
    reader_base(reader_base &&) = default;
    //!\brief Move assignment is defaulted.
    reader_base & operator=(reader_base &&) = default;
    //!\brief Destructor is defaulted.
    ~reader_base() = default;

    /*!\brief Construct from filename.
     * \param[in] filename  Path to the file you wish to open.
     * \param[in] frmt      The file format given as e.g. `format_fasta{}` [optional]
     * \param[in] opt       Reader options (exact type depends on specialisation). [optional]
     * \throws seqan3::file_open_error If the file could not be opened, e.g. non-existant, non-readable, unknown format.
     *
     * \details
     *
     * In addition to the file name, you may specify a custom seqan3::fields type which may be easier than
     * defining all the template parameters.
     *
     *
     * ### Decompression
     *
     * This constructor transparently applies a decompression stream on top of the file stream in case
     * the file is detected as being compressed.
     * See the section on \link io_compression compression and decompression \endlink for more information.
     */
    reader_base(std::filesystem::path const & filename, format_type const & fmt, options_t const & opt = options_t{}) :
        options{opt},
        stream{filename, opt.stream_options},
        format{fmt}
    {}

    //!\overload
    explicit reader_base(std::filesystem::path const & filename, options_t const & opt = options_t{}) :
        options{opt},
        stream{filename, opt.stream_options}
    {
        // initialise format handler or throw if format is not found
        detail::set_format(format, stream.truncated_filename());
    }

    /*!\brief Construct from an existing stream and with specified format.
     * \param[in] str  The stream to operate on.
     * \param[in] frmt The file format given as e.g. `format_fasta{}`.
     * \param[in] opt  Reader options (exact type depends on specialisation). [optional]
     *
     * \details
     *
     * ### Decompression
     *
     * This constructor transparently applies a decompression stream on top of the stream in case
     * it is detected as being compressed.
     * See the section on \link io_compression compression and decompression \endlink for more information.
     */
    reader_base(std::istream & str, format_type const & frmt, options_t const & opt = options_t{}) :
        options{opt},
        stream{str, opt.stream_options},
        format{frmt}
    {}

    //!\overload
    reader_base(std::istream && str, format_type const & frmt, options_t const & opt = options_t{}) :
        options{opt},
        stream{std::move(str), opt.stream_options},
        format{frmt}
    {}
    //!\}

    /*!\name Range interface
     * \brief Provides functions for record based reading of the file.
     * \{
     */
    /*!\brief Returns an iterator to current position in the file.
     * \returns An iterator pointing to the current position in the file.
     * \throws seqan3::format_error
     *
     * It is safe to call this function repeatedly, but it will always return an iterator pointing to the current
     * record in the file (and not seek back to the beginning).
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
                format_handler = input_format_handler<decltype(f)>{stream, options};
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
     * \include test/snippet/io/sequence_file/reader_base_return_record.cpp
     *
     * It most situations using the iterator interface or a range-based for-loop are preferable to using front(),
     * because you can only move to the next record via the iterator.
     *
     * In any case, don't forget the reference! If you want to save the data from the record elsewhere, use move:
     *
     * \include test/snippet/io/sequence_file/reader_base_record_move.cpp
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    record_type & front() noexcept
    {
        return *begin();
    }
    //!\}

protected:
    //!\privatesection


    /*!\name Data buffers
     * \{
     */
    options_t                 options;
    transparent_istream<char> stream;

    //!\brief Buffer for a single record.
    record_type record_buffer;
    //!\}

    //!\brief Tracks whether the very first record is buffered when calling begin().
    bool first_record_was_read{false};
    //!\brief File is at position 1 behind the last record.
    bool at_end{false};



    //!\brief The actual std::variant holding a pointer to the detected/selected format.
    format_type format;
    format_handler_type format_handler;
    //!\}

    //!\brief Tell the format to move to the next record and update the buffer.
    void read_next_record()
    {
        // at end if we could not read further
        if ((std::istreambuf_iterator<char>{stream} ==
             std::istreambuf_iterator<char>{}))
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
 * \relates seqan3::reader_base
 * \{
 */

//!\brief Deduces the sequence input file type from the stream and the format.
// template <typename options_type>
// reader_base(std::filesystem::path, options_type)
//     -> reader_base<options_type>;
//
// //!\brief Deduces the sequence input file type from the stream, the format and the field ids.
// template <typename options_type>
// reader_base(std::istream && stream, options_type)
//     -> reader_base<options_type>;
//
// //!\overload
// template <typename options_type>
// reader_base(std::istream & stream, options_type)
//     -> reader_base<options_type>;
//!\}

} // namespace seqan3
