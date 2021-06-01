// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::plaintext_file_input.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <seqan3/io/stream/iterator.hpp>
#include <seqan3/std/iterator>

#include <iostream> //DEBUG
#include <variant>

namespace seqan3
{

//!\brief The value type of seqan3::plaintext_file_input if every line is split into fields.
struct plaintext_record
{
    //!\brief The entire line (exluding EOL characters but including delimiters).
    std::string_view line;
    //!\brief A range of the individual fields (without delimiters or EOL characters).
    std::vector<std::string_view> fields;
};


enum class plaintext_record_kind
{
    line,
    line_and_fields
};

} // namespace seqan3

namespace seqan3::detail
{


template <typename char_t,
          typename traits_t = std::char_traits<char_t>,
          plaintext_record_kind record_kind = plaintext_record_kind::line>
class plaintext_iterator
{
private:
    //!\brief Down-cast pointer to the stream-buffer.
    stream_buffer_exposer<char_t, traits_t> * stream_buf = nullptr;

    std::string             overflow_buffer;
    std::vector<size_t> field_end_positions;

    plaintext_record record;

    bool at_end = false;
    char_t field_sep = '\t';
    char_t rec_sep = '\n';

public:
    /*!\name Associated types
     * \{
     */
    using difference_type   = ptrdiff_t;                //!< Defaults to ptrdiff_t.
    //!\brief The char type of the stream.
    using value_type        = std::conditional_t<record_kind == plaintext_record_kind::line, std::string_view, plaintext_record>;
    using reference         = value_type &;             //!< The char type of the stream.
    using pointer           = void;                     //!< Has no pointer type.
    using iterator_category = std::input_iterator_tag;  //!< Pure input iterator.
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    plaintext_iterator()                                        noexcept = default; //!< Defaulted.
    plaintext_iterator(plaintext_iterator const &)              noexcept = default; //!< Defaulted.
    plaintext_iterator(plaintext_iterator &&)                   noexcept = default; //!< Defaulted.
    plaintext_iterator & operator=(plaintext_iterator const &)  noexcept = default; //!< Defaulted.
    plaintext_iterator & operator=(plaintext_iterator &&)       noexcept = default; //!< Defaulted.
    ~plaintext_iterator()                                       noexcept = default; //!< Defaulted.

    //!\brief Construct from a stream buffer.
    explicit plaintext_iterator(std::basic_streambuf<char_t, traits_t> & ibuf, bool const init = true) :
        stream_buf{reinterpret_cast<stream_buffer_exposer<char_t, traits_t> *>(&ibuf)}
    {
        if (init) // read first record
        {
            operator++();
        }
        else // only ensure streambuffer is readable
        {
            assert(stream_buf != nullptr);
            stream_buf->underflow(); // ensure the stream buffer has content on construction
            if (stream_buf->gptr() == stream_buf->egptr())
                at_end = true;
        }
    }

    //!\brief Construct from a stream.
    explicit plaintext_iterator(std::basic_istream<char_t, traits_t> & istr, bool const init = true) :
        plaintext_iterator{*istr.rdbuf(), init}
    {}


    plaintext_iterator(std::basic_streambuf<char_t, traits_t> & ibuf, char const sep, bool const init = true)
        requires (record_kind == plaintext_record_kind::line_and_fields)
        : plaintext_iterator{ibuf, init}
    {
        field_sep = sep;
    }

    plaintext_iterator(std::basic_istream<char_t, traits_t> & istr, char const sep, bool const init = true)
        requires (record_kind == plaintext_record_kind::line_and_fields)
        : plaintext_iterator{*istr.rdbuf(), sep, init}
    {}

    //!\}

    /*!\name Arithmetic operators
     * \{
     */
    //!\brief Advance by one and rebuffer if necessary (vtable lookup iff rebuffering).
    plaintext_iterator & operator++()
    {
        assert(stream_buf != nullptr);

        if (at_end)
            return *this;;

        if (stream_buf->gptr() == stream_buf->egptr()) // possible to be on empty buffer
        {
            stream_buf->underflow();
            if (stream_buf->gptr() == stream_buf->egptr())
            {
                at_end = true;
                return *this;
            }
        }

        overflow_buffer.clear();
        if constexpr (record_kind == plaintext_record_kind::line_and_fields)
            field_end_positions.clear();

        bool rec_end_found  = false;
        bool has_overflowed = false;
        size_t count        = 0;
        size_t old_count    = 0;
        char * data_begin   = stream_buf->gptr(); // point into stream buffer by default

        //TODO: does this handle premature EOF?
        while (!rec_end_found)
        {
            for (count = 0; count < (stream_buf->egptr() - stream_buf->gptr()); ++count)
            {
                if (stream_buf->gptr()[count] == rec_sep)
                {
                    rec_end_found = true;
                    break;
                }
                else
                {
                    if constexpr (record_kind == plaintext_record_kind::line_and_fields)
                        if (stream_buf->gptr()[count] == field_sep)
                            field_end_positions.push_back(old_count + count);
                }
            }

            if (!rec_end_found)
            {
                has_overflowed = true;
                overflow_buffer.resize(old_count + count);
                std::ranges::copy(stream_buf->gptr(), stream_buf->egptr(), overflow_buffer.data() + old_count);

                old_count += count;
                stream_buf->gbump(count);
                stream_buf->underflow();
            }
        }

        if (has_overflowed)
        {
            // need to copy last data
            overflow_buffer.resize(old_count + count);
            std::ranges::copy(stream_buf->gptr(),
                              stream_buf->gptr() + count,
                              overflow_buffer.data() + old_count);

            // make data pointer point into overflow
            data_begin = overflow_buffer.data();
        }

        size_t end_of_record = old_count + count;
        // dirty hack for CR: skip it in the buffer but don't add to output
        if (count > 0 && (stream_buf->gptr()[count - 1] == '\r'))
            --end_of_record;

        // move get pointer to point BEHIND current record / beginning of next. This may make gptr == egptr!
        stream_buf->gbump(count + 1);

        /* create the record */
        record.line = std::string_view{data_begin, end_of_record};
        if constexpr (record_kind == plaintext_record_kind::line_and_fields)
        {
            // add last end position
            field_end_positions.push_back(end_of_record);

            record.fields.clear();
            for (size_t i = 0; i < field_end_positions.size(); ++i)
            {
                if (i == 0)
                {
                    record.fields.emplace_back(data_begin,               // ptr
                                               field_end_positions[0]);  // size
                }
                else
                {
                    record.fields.emplace_back(data_begin + field_end_positions[i - 1] + 1,                 // ptr
                                               field_end_positions[i] -  field_end_positions[i - 1] - 1);   // size
                }
            }
        }

        return *this;
    }

    //!\overload
    void operator++(int)
    {
        ++(*this);
    }
    //!\}

    //!\brief Read current value from buffer (no vtable lookup, safe even at end).
    reference operator*()
    {
        if constexpr (record_kind == plaintext_record_kind::line_and_fields)
            return record;
        else
            return record.line;
    }

    //!\brief Show the character behind the current record.
    // ATTENTION: calling this function may invalidate the current record!
    // UB if stream is at_end
    char_t peak()
    {
        assert(stream_buf != nullptr);
        if (stream_buf->gptr() == stream_buf->egptr())
            stream_buf->underflow();
        assert(stream_buf->gptr() != stream_buf->egptr());
        return *stream_buf->gptr();
    }


    /*!\name Comparison operators
     * \brief We define comparison only against the sentinel.
     * \{
     */
    //!\brief True if the read buffer is not empty; involves no vtable lookup.
    friend bool operator==(plaintext_iterator const & lhs, std::default_sentinel_t const &) noexcept
    {
        return lhs.at_end;
    }

    //!\brief True if the read buffer is empty; involves no vtable lookup.
    friend bool operator!=(plaintext_iterator const & lhs, std::default_sentinel_t const &) noexcept
    {
        return !(lhs == std::default_sentinel);
    }

    //!\brief True if the read buffer is not empty; involves no vtable lookup.
    friend bool operator==(std::default_sentinel_t const &, plaintext_iterator const & rhs) noexcept
    {
        return rhs == std::default_sentinel;
    }

    //!\brief True if the read buffer is empty; involves no vtable lookup.
    friend bool operator!=(std::default_sentinel_t const &, plaintext_iterator const & rhs) noexcept
    {
        return !(rhs == std::default_sentinel);
    }
    //!\}
};

} // namespace seqan3::detail

#if 0

namespace seqan3
{

struct plaintext_file_header
{
    struct {} none;
    struct {} first_line;

    struct starts_with
    {
        char c;
    };

    using _t = std::variant<decltype(none), decltype(first_line), starts_with>;
};

enum class plaintext_record_kind
{
    line,
    line_and_fields
};

template <plaintext_record_kind record>
class plaintext_file_input
{
public:
    /*!\name Range associated types
     * \brief The types necessary to facilitate the behaviour of an input range (used in record-wise reading).
     * \{
     */
    //!\brief The iterator type of this view (an input iterator).
    using iterator          = detail::plaintext_iterator<char, std::char_traits<char_t>, record>;
    //!\brief The const iterator type is void, because files are not const-iterable.
    using const_iterator    = void;
    //!\brief The type returned by end().
    using sentinel          = std::default_sentinel_t;
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    //!\brief Default constructor is explicitly deleted, you need to give a stream or file name.
    plaintext_file_input() = delete;
    //!\brief Copy construction is explicitly deleted, because you can't have multiple access to the same file.
    plaintext_file_input(plaintext_file_input const &) = delete;
    //!\brief Copy assignment is explicitly deleted, because you can't have multiple access to the same file.
    plaintext_file_input & operator=(plaintext_file_input const &) = delete;
    //!\brief Move construction is defaulted.
    plaintext_file_input(plaintext_file_input &&) = default;
    //!\brief Move assignment is defaulted.
    plaintext_file_input & operator=(plaintext_file_input &&) = default;
    //!\brief Destructor is defaulted.
    ~plaintext_file_input() = default;

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
    plaintext_file_input(std::filesystem::path filename, plaintext_file_header::_t = plaintext_file_header::none) :
        primary_stream{new std::ifstream{}, stream_deleter_default}
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
     * \tparam file_format   The format of the file in the stream, must satisfy seqan3::plaintext_file_input_format.
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
    plaintext_file_input(std::istream & stream, config_t const & cfg) :
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
    plaintext_file_input(std::istream && stream, config_t const & cfg) :
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
     * \include test/snippet/io/sequence_file/plaintext_file_input_return_record.cpp
     *
     * It most situations using the iterator interface or a range-based for-loop are preferable to using front(),
     * because you can only move to the next record via the iterator.
     *
     * In any case, don't forget the reference! If you want to save the data from the record elsewhere, use move:
     *
     * \include test/snippet/io/sequence_file/plaintext_file_input_record_move.cpp
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

    iterator it;

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
    //!\}

};

} // namespace seqan3

#endif
