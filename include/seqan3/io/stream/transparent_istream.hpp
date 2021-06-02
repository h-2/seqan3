// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::transparent_istream.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <fstream>
#include <iostream>
#include <thread>

#include <seqan3/core/platform.hpp>
#include <seqan3/io/stream/compression.hpp>
#include <seqan3/io/detail/misc_input.hpp>
#include <seqan3/std/filesystem>
#include <seqan3/std/span>

namespace seqan3
{

//!\brief Options that can be provided to seqan3::transparent_istream.
struct transparent_istream_options
{
    //!\brief Size of the buffer used when opening a file from a filename.
    size_t buffer1_size = 1024 * 1024;
    //!\brief Size of the buffer used for the decompression stream.
    size_t buffer2_size = 1024 * 1024 * 4;

    /*!\brief Which decompressor to use.
     *
     * \details
     *
     * In almost all cases, you will want to have this auto-detected. You can explicitly set this to
     * seqan3::compression_format::gz when opening BGZF-compressed files to use the regular, single-threaded
     * GZ-decompressor, but it is simpler to just set threads to 1.
     */
    compression_format compression = compression_format::detect;

    /*!\brief Maximum number of threads to use for decompression.
     *
     * \details
     *
     * This value is currently only relevant for BGZF compressed streams/files. Note that these threads refer to the
     * total number of used threads, i.e. a value of 4 means that three extra threads are spawned. A value of 1 will
     * result in the regular GZ decompressor being used.
     *
     * The default value for this is 4 or "available CPUs - 1" if that is less than 4. The reason is that for default
     * the compression levels of the GZip blocks, performance degrades with more than 4 threads. If you use files
     * compressed with stronger settings, more threads might be useful.
     *
     * **4 threads can provide up to 3.5x speed-up.**
     */
    size_t threads  = std::max<size_t>(1, std::min<size_t>(4, std::thread::hardware_concurrency() - 1));
};

template <typename char_t>
class transparent_istream : public std::basic_istream<char_t>
{
private:
    /*TODO
     * evaluate whether to use custom sized buffer on both streams
     * evaluate whether to use custom sized buffer strings passed in (what happens to stuff in old buffer?)
     * tie in threads
     */

    //!\brief The type of the internal stream pointers. Allows dynamically setting ownership management.
    using stream_ptr_t = std::unique_ptr<std::basic_istream<char_t>,
                                         std::function<void(std::basic_istream<char_t>*)>>;
    //!\brief Stream deleter that does nothing (no ownership assumed).
    static void stream_deleter_noop(std::basic_istream<char_t> *) {}
    //!\brief Stream deleter with default behaviour (ownership assumed).
    static void stream_deleter_default(std::basic_istream<char_t> * ptr) { delete ptr; }

    //!\brief The primary stream is the user provided stream or the file stream if constructed from filename.
    stream_ptr_t primary_stream{nullptr, stream_deleter_noop};
    //!\brief The secondary stream is a compression layer on the primary or just points to the primary (no compression).
    stream_ptr_t secondary_stream{nullptr, stream_deleter_noop};


    //!\brief The options.
    transparent_istream_options _options;
    //!\brief The stream buffer.
    std::vector<char>           stream1_buffer;
    std::vector<char>           stream2_buffer;
    //!\brief Filename (if stream was opened from path).
    std::filesystem::path       _filename;
    //!\brief Filename after possible compression extensions have been removed.
    std::filesystem::path       _truncated_filename;

    //!\brief This function reads the magic bytes from the stream and adds a decompression layer if necessary.
    void set_secondary_stream()
    {
        assert(primary_stream->good());

        /* detect compression format */
        std::string magic_header = detail::read_magic_header(*primary_stream);
        compression_format selected_compression{};

        if (_options.compression == compression_format::detect)
        {
            selected_compression = detail::detect_format_from_magic_header(magic_header);
        }
        else
        {
            selected_compression = _options.compression;
            if (!detail::header_matches_dyn(selected_compression, magic_header))
                throw file_open_error{"The file has a different compression format than the one selected."};
        }

        // Thread handling
        if (selected_compression == compression_format::bgzf)
        {
            if (_options.threads == 1)
                selected_compression = compression_format::gz;
            else
                --_options.threads; // bgzf spawns **additional** threads, but user sets total
        }

        std::span<std::string> file_extensions{};
        switch (selected_compression)
        {
            case compression_format::bgzf:
                secondary_stream = stream_ptr_t{detail::make_istream<compression_format::bgzf>(*primary_stream,
                                                                                              _options.threads),
                                                stream_deleter_default};
                file_extensions = compression_traits<compression_format::bgzf>::file_extensions;
                break;
            case compression_format::gz:
                secondary_stream = stream_ptr_t{detail::make_istream<compression_format::gz>(*primary_stream),
                                                stream_deleter_default};
                file_extensions = compression_traits<compression_format::gz>::file_extensions;
                break;
            case compression_format::bz2:
                secondary_stream = stream_ptr_t{detail::make_istream<compression_format::bz2>(*primary_stream),
                                                stream_deleter_default};
                file_extensions = compression_traits<compression_format::bz2>::file_extensions;
                break;
            case compression_format::zstd:
                secondary_stream = stream_ptr_t{detail::make_istream<compression_format::zstd>(*primary_stream),
                                                stream_deleter_default};
                file_extensions = compression_traits<compression_format::zstd>::file_extensions;
                break;
            default:
                secondary_stream = stream_ptr_t{&*primary_stream, stream_deleter_noop};
                return;
        }

        // truncate the filename in _truncated_filename to show that decompression has taken place
        if (_filename.has_extension())
        {
            std::string extension = _filename.extension().string().substr(1);
            if (std::ranges::find(file_extensions, extension) !=  std::ranges::end(file_extensions))
                _truncated_filename.replace_extension();
        }
    }

    void init()
    {
        _truncated_filename = _filename;

        // possibly add intermediate compression stream
        set_secondary_stream();

        // make this behave like the secondary stream
        this->rdbuf(secondary_stream->rdbuf());
    }

public:

    /*!\name Constructors, destructor and assignment
     * \{
     */
    transparent_istream()                                           = delete;   //!< Defaulted.
    transparent_istream(transparent_istream const &)                = delete;   //!< Defaulted.
    transparent_istream(transparent_istream &&)                     = default;  //!< Defaulted.
    transparent_istream & operator=(transparent_istream const &)    = delete;   //!< Defaulted.
    transparent_istream & operator=(transparent_istream &&)         = default;  //!< Defaulted.

    /*!\brief Construct from a filename.
     * \param[in] filename  The filename to open.
     * \param[in] options   See seqan3::transparent_istream_options.
     */
    explicit transparent_istream(std::filesystem::path filename,
                                 transparent_istream_options options = transparent_istream_options{}) :
        primary_stream{new std::ifstream{}, stream_deleter_default},
        _options{std::move(options)},
        _filename{std::move(filename)}
    {
        stream1_buffer.resize(options.buffer1_size);

        primary_stream->rdbuf()->pubsetbuf(stream1_buffer.data(), stream1_buffer.size());
        static_cast<std::basic_ifstream<char> *>(primary_stream.get())->open(_filename,
                                                                             std::ios_base::in | std::ios::binary);

        if (!primary_stream->good())
            throw file_open_error{"Could not open file " + _filename.string() + " for reading."};

        init();
    }

    /*!\brief Construct from a stream.
     * \param[in] stream    The stream to wrap.
     * \param[in] options   See seqan3::transparent_istream_options.
     */
    explicit transparent_istream(std::istream & stream,
                                 transparent_istream_options options = transparent_istream_options{}) :
        primary_stream{&stream, stream_deleter_noop},
        _options{std::move(options)}
    {
        init();
    }

    //!\overload
    explicit transparent_istream(std::istream && stream,
                                 transparent_istream_options options = transparent_istream_options{}) :
        primary_stream{new std::istream{std::move(stream)}, stream_deleter_default},
        _options{std::move(options)}
    {
        init();
    }
    //!\}

    //!\brief The filename this object was created from; empty if this object was not created from a file.
    std::filesystem::path const & filename()
    {
        return _filename;
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
        return _truncated_filename;
    }
};

} // namespace seqan3
