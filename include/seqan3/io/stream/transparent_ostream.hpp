// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::transparent_ostream.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <fstream>
#include <iostream>
#include <thread>

#include <seqan3/io/stream/compression.hpp>
#include <seqan3/io/detail/misc_input.hpp>
#include <seqan3/std/filesystem>
#include <seqan3/std/span>

namespace seqan3
{

//!\brief Options that can be provided to seqan3::transparent_ostream.
struct transparent_ostream_options
{
    //!\brief Size of the buffer used when opening a file from a filename.
    size_t buffer1_size = 1024 * 1024;
    //!\brief Size of the buffer used for the compression stream.
    size_t buffer2_size = 1024 * 1024 * 4;

    /*!\brief Which compressor to use.
     *
     * \details
     *
     * For ostream opened from filenames, the default is to detect the desired compression from the file's extension.
     * But you can also specifiy a desired compression format manually.
     */
    compression_format compression = compression_format::detect;

    /*!\brief The compression level to use by the algorithm.
     *
     * \details
     *
     * The default value is -1 which maps to the default value of the respective algorithm (6 for GZ/BGZF and 9 for
     * BZip2). ZLIB macros and numeric values between -1 and 9 are supported.
     */
    int compression_level = -1;

    /*!\brief Maximum number of threads to use for compression.
     *
     * \details
     *
     * This value is currently only relevant for BGZF compressed streams/files. Note that these threads refer to the
     * total number of used threads, i.e. a value of 4 means that three extra threads are spawned.
     *
     * The default value for this "all available CPUs" as more threads typically improve performance (albeit not
     * linearly).
     *
     * ### Attention
     *
     * A value of 1 is currently not supported by the BGZF implementation!
     */
    size_t threads = std::thread::hardware_concurrency();
};

/*!\brief A std::ostream that automatically detects compressed streams and transparently decompresses them.
 * \tparam char_t The character type of the stream.
 */
template <typename char_t = char>
class transparent_ostream : public std::basic_ostream<char_t>
{
private:
    /*TODO
     * evaluate whether to use custom sized buffer on both streams
     * evaluate whether to use custom sized buffer strings passed in (what happens to stuff in old buffer?)
     * tie in threads
     */

    //!\brief The options.
    transparent_ostream_options options_;
    //!\brief The stream buffer.
    std::vector<char>           stream1_buffer;
    std::vector<char>           stream2_buffer;
    //!\brief Filename (if stream was opened from path).
    std::filesystem::path       filename_;
    //!\brief Filename after possible compression extensions have been removed.
    std::filesystem::path       truncated_filename_;


    //!\brief The type of the internal stream pointers. Allows dynamically setting ownership management.
    using stream_ptr_t = std::unique_ptr<std::basic_ostream<char_t>,
                                         std::function<void(std::basic_ostream<char_t>*)>>;
    //!\brief Stream deleter that does nothing (no ownership assumed).
    static void stream_deleter_noop(std::basic_ostream<char_t> *) {}
    //!\brief Stream deleter with default behaviour (ownership assumed).
    static void stream_deleter_default(std::basic_ostream<char_t> * ptr) { delete ptr; }

    //!\brief The primary stream is the user provided stream or the file stream if constructed from filename.
    stream_ptr_t primary_stream{nullptr, stream_deleter_noop};
    //!\brief The secondary stream is a compression layer on the primary or just points to the primary (no compression).
    stream_ptr_t secondary_stream{nullptr, stream_deleter_noop};


    //!\brief This function reads the magic bytes from the stream and adds a decompression layer if necessary.
    void set_secondary_stream()
    {
        assert(primary_stream->good());

        /* detect compression format */
        if (options_.compression == compression_format::detect)
        {
            if (filename_.empty())
            {
                throw file_open_error{"Cannot auto-detect compression type from arbitrary streams."
                                      "Please select \"none\" or a specific compression format."};
            }
            options_.compression = detail::detect_format_from_filename(filename_);
        }

        // Thread handling
        if (options_.compression == compression_format::bgzf)
        {

            if (options_.threads == 1) // TODO this needs a real resolution
                throw file_open_error{"BGZF compression with only one thread is currently not supported."};
            else
                --options_.threads; // bgzf spawns **additional** threads, but user sets total
        }

        std::span<std::string> file_extensions{};
        std::ostream * sec = nullptr;
        switch (options_.compression)
        {
            case compression_format::bgzf:
                sec = detail::make_ostream<compression_format::bgzf>(*primary_stream,
                                                                     options_.threads,
                                                                     static_cast<size_t>(8ul),
                                                                     options_.compression_level);
                file_extensions = compression_traits<compression_format::bgzf>::file_extensions;
                break;
            case compression_format::gz:
                sec = detail::make_ostream<compression_format::gz>(*primary_stream, options_.compression_level);
                file_extensions = compression_traits<compression_format::gz>::file_extensions;
                break;
            case compression_format::bz2:
                sec = detail::make_ostream<compression_format::bz2>(*primary_stream, options_.compression_level);
                file_extensions = compression_traits<compression_format::bz2>::file_extensions;
                break;
            case compression_format::zstd:
                sec = detail::make_ostream<compression_format::zstd>(*primary_stream, options_.compression_level);
                file_extensions = compression_traits<compression_format::zstd>::file_extensions;
                break;
            default:
                break;
        }

        if (sec == nullptr)
            secondary_stream = stream_ptr_t{&*primary_stream, stream_deleter_noop};
        else
            secondary_stream = stream_ptr_t{sec, stream_deleter_default};

        // truncate the filename in truncated_filename_ to show that compression has taken place
        if (filename_.has_extension())
        {
            std::string extension = filename_.extension().string().substr(1);
            if (std::ranges::find(file_extensions, extension) !=  std::ranges::end(file_extensions))
                truncated_filename_.replace_extension();
        }
    }

    void init()
    {
        truncated_filename_ = filename_;

        // possibly add intermediate compression stream
        set_secondary_stream();
        assert(secondary_stream != nullptr);

        // make this behave like the secondary stream
        this->rdbuf(secondary_stream->rdbuf());
    }

public:

    /*!\name Constructors, destructor and assignment
     * \{
     */
    //!\brief Manually defined default constructor that behaves as expected.
    transparent_ostream() : std::basic_ostream<char_t>{} {}
    transparent_ostream(transparent_ostream const &)                = delete;   //!< Defaulted.
    transparent_ostream & operator=(transparent_ostream const &)    = delete;   //!< Defaulted.
    //TODO double check that this works:
    transparent_ostream & operator=(transparent_ostream &&)         = default;  //!< Defaulted.

    //!\brief Manually defined move constructor that behaves as expected.
    transparent_ostream(transparent_ostream && rhs)
    {
        std::swap(options_,             rhs.options_);
        std::swap(stream1_buffer,       rhs.stream1_buffer);
        std::swap(stream2_buffer,       rhs.stream2_buffer);
        std::swap(filename_,            rhs.filename_);
        std::swap(truncated_filename_,  rhs.truncated_filename_);
        std::swap(primary_stream,       rhs.primary_stream);
        std::swap(secondary_stream,     rhs.secondary_stream);

        this->set_rdbuf(secondary_stream->rdbuf());
    }
    /*!\brief Construct from a filename.
     * \param[in] filename  The filename to open.
     * \param[in] options   See seqan3::transparent_ostream_options.
     *
     * \details
     *
     * The compression format is auto-detected from the filename by default. It can manually be selected via the
     * options.
     *
     * The stream is opened in binary mode and provided with a buffer the size of options.buffer1_size.
     */
    explicit transparent_ostream(std::filesystem::path filename,
                                 transparent_ostream_options options = transparent_ostream_options{}) :
        primary_stream{new std::ofstream{}, stream_deleter_default},
        options_{std::move(options)},
        filename_{std::move(filename)}
    {
        stream1_buffer.resize(options.buffer1_size);

        primary_stream->rdbuf()->pubsetbuf(stream1_buffer.data(), stream1_buffer.size());
        static_cast<std::basic_ofstream<char> *>(primary_stream.get())->open(filename_,
                                                                             std::ios_base::out | std::ios::binary);

        if (!primary_stream->good())
            throw file_open_error{"Could not open file " + filename_.string() + " for reading."};

        init();
    }

    /*!\brief Construct from a stream.
     * \param[in] stream  The stream to wrap.
     * \param[in] options See seqan3::transparent_ostream_options.
     *
     * \details
     *
     * The compression format is "none" by default. It can manually be selected via the options.
     *
     *
     */
    explicit transparent_ostream(std::ostream & stream,
                                 transparent_ostream_options options = transparent_ostream_options{.compression = compression_format::none}) :
        primary_stream{&stream, stream_deleter_noop},
        options_{std::move(options)}
    {
        init();
    }

    //!\overload
    explicit transparent_ostream(std::ostream && stream,
                                 transparent_ostream_options options = transparent_ostream_options{.compression = compression_format::none}) :
        primary_stream{new std::ostream{std::move(stream)}, stream_deleter_default},
        options_{std::move(options)}
    {
        init();
    }

    //!\brief The filename this object was created from; empty if this object was not created from a file.
    std::filesystem::path const & filename()
    {
        return filename_;
    }

    /*!\brief The filename this object was created from without compression-specific suffix.
     *
     * \details
     *
     * If this object was created from e.g. "foo.fasta.gz", #filename() will return the full name, but this function
     * will return only "foo.fasta". This is useful for determining the format "inside" the compression.
     *
     * If this object was not created from a file, an empty path is returned.
     */
    std::filesystem::path const & truncated_filename()
    {
        return truncated_filename_;
    }
};

} // namespace seqan3
