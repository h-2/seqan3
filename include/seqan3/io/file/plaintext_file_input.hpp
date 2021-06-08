// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::reader.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <string_view>
#include <variant>

#include <seqan3/io/file/plaintext_file_util.hpp>
#include <seqan3/io/stream/iterator.hpp>
#include <seqan3/io/stream/transparent_istream.hpp>

namespace seqan3::plain_io::detail
{

template <typename char_t,
          typename traits_t = std::char_traits<char_t>,
          record_kind record_kind_ = record_kind::line>
class plaintext_input_iterator
{
private:
    //!\brief Down-cast pointer to the stream-buffer.
    seqan3::detail::stream_buffer_exposer<char_t, traits_t> * stream_buf = nullptr;

    //!\brief Place to store lines that overlap buffer boundaries.
    std::basic_string<char_t, traits_t> overflow_buffer;
    //!\brief Temporary storage for field delimiter positions.
    std::vector<size_t>                 field_end_positions;

    //!\brief The record.
    record record_;

    //!\brief Whether iterator is at end.
    bool   at_end    = false;
    //!\brief Delimiter between fields.
    char_t field_sep = '\t';
    //!\brief Delimiter between records [not exposed to modification ATM].
    char_t rec_sep   = '\n';

public:
    /*!\name Associated types
     * \{
     */
    using difference_type   = ptrdiff_t;                //!< Defaults to ptrdiff_t.
    //!\brief The record type.
    using value_type        = std::conditional_t<record_kind_ == record_kind::line,
                                                 std::string_view,
                                                 record>;
    using reference         = std::conditional_t<record_kind_ == record_kind::line,
                                                 std::string_view,
                                                 record const &>;
    using pointer           = value_type const *;       //!< Has no pointer type.
    using iterator_category = std::input_iterator_tag;  //!< Pure input iterator.
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    plaintext_input_iterator()                                        noexcept = default; //!< Defaulted.
    plaintext_input_iterator(plaintext_input_iterator const &)              noexcept = default; //!< Defaulted.
    plaintext_input_iterator(plaintext_input_iterator &&)                   noexcept = default; //!< Defaulted.
    plaintext_input_iterator & operator=(plaintext_input_iterator const &)  noexcept = default; //!< Defaulted.
    plaintext_input_iterator & operator=(plaintext_input_iterator &&)       noexcept = default; //!< Defaulted.
    ~plaintext_input_iterator()                                       noexcept = default; //!< Defaulted.

    //!\brief Construct from a stream buffer.
    explicit plaintext_input_iterator(std::basic_streambuf<char_t, traits_t> & ibuf, bool const init = true) :
        stream_buf{reinterpret_cast<seqan3::detail::stream_buffer_exposer<char_t, traits_t> *>(&ibuf)}
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

    plaintext_input_iterator(std::basic_streambuf<char_t, traits_t> & ibuf, char_t const sep, bool const init = true)
        requires (record_kind_ == record_kind::line_and_fields)
        : plaintext_input_iterator{ibuf, init}
    {
        field_sep = sep;
    }

    //!\brief Construct from a stream.
    explicit plaintext_input_iterator(std::basic_istream<char_t, traits_t> & istr, bool const init = true) :
        plaintext_input_iterator{*istr.rdbuf(), init}
    {}


    plaintext_input_iterator(std::basic_istream<char_t, traits_t> & istr, char_t const sep, bool const init = true)
        requires (record_kind_ == record_kind::line_and_fields)
        : plaintext_input_iterator{*istr.rdbuf(), sep, init}
    {}

    //!\}

    /*!\name Arithmetic operators
     * \{
     */
    //!\brief Advance by one and rebuffer if necessary (vtable lookup iff rebuffering).
    plaintext_input_iterator & operator++()
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
        if constexpr (record_kind_ == record_kind::line_and_fields)
            field_end_positions.clear();

        bool rec_end_found  = false;
        bool has_overflowed = false;
        size_t count        = 0;
        size_t old_count    = 0;
        char_t * data_begin   = stream_buf->gptr(); // point into stream buffer by default

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
                    if constexpr (record_kind_ == record_kind::line_and_fields)
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
        record_.line = std::string_view{data_begin, end_of_record};
        if constexpr (record_kind_ == record_kind::line_and_fields)
        {
            // add last end position
            field_end_positions.push_back(end_of_record);

            record_.fields.clear();
            for (size_t i = 0; i < field_end_positions.size(); ++i)
            {
                if (i == 0)
                {
                    record_.fields.emplace_back(data_begin,               // ptr
                                               field_end_positions[0]);  // size
                }
                else
                {
                    record_.fields.emplace_back(data_begin + field_end_positions[i - 1] + 1,                 // ptr
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
    reference operator*() const
    {
        if constexpr (record_kind_ == record_kind::line_and_fields)
            return record_;
        else
            return record_.line;
    }

    pointer operator->() const
    {
        return &**this;
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
    friend bool operator==(plaintext_input_iterator const & lhs, std::default_sentinel_t const &) noexcept
    {
        return lhs.at_end;
    }

    //!\brief True if the read buffer is empty; involves no vtable lookup.
    friend bool operator!=(plaintext_input_iterator const & lhs, std::default_sentinel_t const &) noexcept
    {
        return !(lhs == std::default_sentinel);
    }

    //!\brief True if the read buffer is not empty; involves no vtable lookup.
    friend bool operator==(std::default_sentinel_t const &, plaintext_input_iterator const & rhs) noexcept
    {
        return rhs == std::default_sentinel;
    }

    //!\brief True if the read buffer is empty; involves no vtable lookup.
    friend bool operator!=(std::default_sentinel_t const &, plaintext_input_iterator const & rhs) noexcept
    {
        return !(rhs == std::default_sentinel);
    }
    //!\}
};

} // namespace seqan3::detail



namespace seqan3::plain_io
{


/*!\brief Line-wise reader of plaintext files; supports transparent decompression.
 * \tparam record_kind_ Whether to split lines on delimiter (e.g. TSV files) or not.
 *
 * \details
 *
 * Main features:
 *   * Reads a file line-by-line. Optionally splits lines on the provided delimiter.
 *   * Very fast: the lines/fields are provided as string-views into the read buffer.
 *   * No dynamic memory allocations happen while iterating over the file.
 *   * Different header formats are supported ("none", "first line", "lines starting with #"...).
 *   * Supports opening files from filenames and wrapping existing streams (e.g. std::cin).
 *   * Automatically detects compressed files/streams and transparently decompresses.
 *
 * Lines read never include the end-of-line character, except in multi-line headers where lines are separated by
 * `'\n'`. Legacy Windows line-endings (including carriage-return) are supported when reading but whether or not any
 * were present in the file is not exposed to the user (they are not preserved even in the header!).
 *
 * This class is particularly well-suited for fast lowlevel parsing of plaintext files, TSV/CSV files, SAM files,
 * VCF files et cetera.
 *
 * ### Attention
 *
 * This file performs line-wise buffering internally. If the file that you are attempting to read contains unreasonably
 * long lines (or is in-fact binary), performance will degrade severely! Line-lengths up 10.000 or even 100.000 might
 * work, but files with e.g. full chromosomes or genomes in a single line are not supported.
 */
template <record_kind record_kind_>
class reader
{
public:
    /*!\name Range associated types
     * \brief The types necessary to facilitate the behaviour of an input range (used in record-wise reading).
     * \{
     */
    //!\brief The iterator type of this view (an input iterator).
    using iterator          = detail::plaintext_input_iterator<char, std::char_traits<char>, record_kind_>;
    //!\brief The const iterator type is void, because files are not const-iterable.
    using const_iterator    = void;
    //!\brief The type returned by end().
    using sentinel          = std::default_sentinel_t;
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    //!\brief Default constructor is explicitly deleted, you need to give a stream or file name.
    reader()                                              = delete;
    //!\brief Copy construction is explicitly deleted, because you can't have multiple access to the same file.
    reader(reader const &)                  = delete;
    //!\brief Move construction is defaulted.
    reader(reader &&)                       = default;
    //!\brief Copy assignment is explicitly deleted, because you can't have multiple access to the same file.
    reader & operator=(reader const &)      = delete;
    //!\brief Move assignment is defaulted.
    reader & operator=(reader &&)           = default;
    //!\brief Destructor is defaulted.
    ~reader()                                             = default;

    /*!\brief Construct from filename.
     * \param[in] filename        Path to the file you wish to open.
     * \param[in] field_separator Delimiter between fields in a line. [optional]
     * \param[in] header          Whether to treat certain lines as header; see seqan3::plain_io::header_kind. [optional]
     * \param[in] istream_options Options passed to the underlying stream; see seqan3::transparent_istream_options. [optional]
     * \throws seqan3::file_open_error If the file could not be opened, e.g. non-existant or non-readable.
     *
     * \details
     *
     * TODO example
     *
     * ### Decompression
     *
     * This constructor transparently applies a decompression stream on top of the file stream in case
     * the file is detected as being compressed.
     * See the section on \link io_compression compression and decompression \endlink for more information.
     */
    explicit reader(std::filesystem::path const & filename,
                    char const field_separator,
                    header_kind header = header_kind::none,
                    transparent_istream_options const & istream_options = transparent_istream_options{})
        requires (record_kind_ == record_kind::line_and_fields)
        : stream{filename, istream_options}, it{stream, field_separator}
    {
        read_header(header);
    }

    //!\overload
    explicit reader(std::filesystem::path const & filename,
                    header_kind header = header_kind::none,
                    transparent_istream_options const & istream_options = transparent_istream_options{})
        requires (record_kind_ == record_kind::line)
        : stream{filename, istream_options}, it{stream}
    {
        read_header(header);
    }

    /*!\brief Construct from an existing stream and with specified format.
     * \param[in] str             The stream to open from.
     * \param[in] field_separator Delimiter between fields in a line. [optional]
     * \param[in] header          Whether to treat certain lines as header; see seqan3::plain_io::header_kind. [optional]
     * \param[in] istream_options Options passed to the underlying stream; see seqan3::transparent_istream_options. [optional]
     * \throws seqan3::file_open_error If the file could not be opened, e.g. non-existant or non-readable.
     *
     * \details
     *
     * TODO example
     *
     * ### Decompression
     *
     * This constructor transparently applies a decompression stream on top of the stream in case
     * it is detected as being compressed.
     * See the section on \link io_compression compression and decompression \endlink for more information.
     */
    explicit reader(std::istream & str,
                    char const field_separator,
                    header_kind header = header_kind::none,
                    transparent_istream_options const & istream_options = transparent_istream_options{})
        requires (record_kind_ == record_kind::line_and_fields)
        : stream{str, istream_options}, it{stream, field_separator}
    {
        read_header(header);
    }

    //!\overload
    explicit reader(std::istream & str,
                    header_kind header = header_kind::none,
                    transparent_istream_options const & istream_options = transparent_istream_options{})
        requires (record_kind_ == record_kind::line)
        : stream{str, istream_options}, it{stream}
    {
        read_header(header);
    }

    //!\overload
    explicit reader(std::istream && str,
                    char const field_separator,
                    header_kind header = header_kind::none,
                    transparent_istream_options const & istream_options = transparent_istream_options{})
        requires (record_kind_ == record_kind::line_and_fields)
        : stream{std::move(str), istream_options}, it{stream, field_separator}
    {
        read_header(header);
    }

    //!\overload
    explicit reader(std::istream && str,
                    header_kind header = header_kind::none,
                    transparent_istream_options const & istream_options = transparent_istream_options{})
        requires (record_kind_ == record_kind::line)
        : stream{std::move(str), istream_options}, it{stream}
    {
        read_header(header);
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

        return it;
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
     * \returns A std::string_view of the current line or a reference to seqan3::record.
     *
     * This function is identical to calling *begin().
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    decltype(auto) front()
    {
        return *begin();
    }
    //!\}

    //!\brief The header from the file.
    std::string_view header() noexcept
    {
        return headr;
    }

protected:
    //!\brief Process the header (if requested).
    void read_header(header_kind _header)
    {
        if (_header.is_none())
        {
            return;
        }
        else if (_header.is_first_line())
        {
            if (it == std::default_sentinel)
                return;

            headr += current_line();
            headr += "\n";
            ++it;
        }
        else
        {
            while (it != std::default_sentinel &&
                   current_line().size() > 0   &&
                   current_line()[0] == _header.get_starts_with())
            {
                headr += current_line();
                headr += "\n";
                ++it;
            }
        }
    }

    //!\brief Return the current line (indepent of record_kind_).
    std::string_view current_line()
    {
        if constexpr (record_kind_ == record_kind::line)
            return *it;
        else
            return it->line;
    }

    //!\brief The underlying stream object.
    transparent_istream<char> stream;
    //!\brief The stream iterator.
    iterator it;

    //!\brief The stored header.
    std::string headr;
};

template <typename t>
reader(t &&,
                     char const,
                     header_kind header = header_kind::none,
                     transparent_istream_options const & istream_options = transparent_istream_options{})
-> reader<record_kind::line_and_fields>;

template <typename t>
reader(t &&,
                     header_kind header = header_kind::none,
                     transparent_istream_options const & istream_options = transparent_istream_options{})
-> reader<record_kind::line>;

// reader(std::istream &,
//                      char const,
//                      header_kind header = header_kind::none,
//                      transparent_istream_options const & istream_options = transparent_istream_options{})
// -> reader<record_kind::line_and_fields>;
//
// reader(std::istream &,
//                      header_kind header = header_kind::none,
//                      transparent_istream_options const & istream_options = transparent_istream_options{})
// -> reader<record_kind::line>;
//
// reader(std::istream &&,
//                      char const,
//                      header_kind header = header_kind::none,
//                      transparent_istream_options const & istream_options = transparent_istream_options{})
// -> reader<record_kind::line_and_fields>;
//
// reader(std::istream &&,
//                      header_kind header = header_kind::none,
//                      transparent_istream_options const & istream_options = transparent_istream_options{})
// -> reader<record_kind::line>;


} // namespace seqan3

