// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::compression_format and seqan3::compression_traits.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <string_view>

#include <seqan3/core/platform.hpp>
#ifdef SEQAN3_HAS_BZIP2
    #include <seqan3/contrib/stream/bz2_istream.hpp>
    #include <seqan3/contrib/stream/bz2_ostream.hpp>
#endif
#ifdef SEQAN3_HAS_ZLIB
    #include <seqan3/contrib/stream/bgzf_istream.hpp>
    #include <seqan3/contrib/stream/bgzf_ostream.hpp>
    #include <seqan3/contrib/stream/bgzf_stream_util.hpp>
    #include <seqan3/contrib/stream/gz_istream.hpp>
    #include <seqan3/contrib/stream/gz_ostream.hpp>
#endif
#include <seqan3/std/filesystem>

namespace seqan3
{

//!\brief Possible formats for stream (de-)compression.
//!\ingroup io
enum class compression_format
{
    none,   //!< No compression.
    detect, //!< Detect compression format automatically.
    bgzf,   //!< Blocked GZip format TODO link.
    gz,     //!< GZip format.
    bz2,    //!< BZip2 format.
    zstd    //!< ZStd format.
};

/*!\brief Traits of the compression formats.
 * \tparam format The seqan3::compression_format whose traits are provided.
 * \ingroup io
 *
 * \details
 *
 * Note that this is **not** a customisation point. You may not specialise it.
 */
template <compression_format format>
struct compression_traits
{
    //!\brief The compression format as a human-readable string.
    static constexpr std::string_view as_string = "";

    /*!\brief The valid file extension for the compression format.
     *
     * \details
     *
     * Note that this static member is not `const` and may be modified to allow
     * user-provided extensions.
     */
    static inline std::vector<std::string> file_extensions
    {
        {"gz", "bgz", "bgzf"}
    };

    //!\brief The magic byte sequence to disambiguate the compression format.
    static constexpr std::string_view magic_header = "";

    //!\brief Whether this compression format was available at build-time.
    static constexpr bool available = false;

    //!\brief The default input stream corresponding to the format.
    template <typename char_t>
    using basic_istream = std::void_t<char_t>;

    //!\brief The default output stream corresponding to the format.
    template <typename char_t>
    using basic_ostream = std::void_t<char_t>;
};

/*!\brief Traits for the seqan3::compression_format::bgzf.
 * \see seqan3::compression_traits
 */
template <>
struct compression_traits<compression_format::bgzf> : compression_traits<compression_format::none>
{
    //!\copydoc seqan3::compression_traits<compression_format::none>::as_string
    static constexpr std::string_view as_string = "BGZF";

    //!\copydoc seqan3::compression_traits<compression_format::none>::file_extensions
    static inline std::vector<std::string> file_extensions
    {
        {"gz", "bgz", "bgzf"}
    };

    //!\copydoc seqan3::compression_traits<compression_format::none>::magic_header
    static constexpr std::string_view magic_header{
        // GZip header
        "\x1f\x8b\x08"
        //FLG[MTIME         ] XFL OS [XLEN  ]
        "\x04\x00\x00\x00\x00\x00\xff\x06\x00"
        // B   C [SLEN  ][BSIZE ]
        "\x42\x43\x02\x00\x00\x00",
        18};

#ifdef SEQAN3_HAS_ZLIB
    //!\copydoc seqan3::compression_traits<compression_format::none>::available
    static constexpr bool available = true;

    //!\copydoc seqan3::compression_traits<compression_format::none>::basic_istream
    template <typename char_t>
    using basic_istream = contrib::basic_bgzf_istream<char_t>;

    //!\copydoc seqan3::compression_traits<compression_format::none>::basic_ostream
    template <typename char_t>
    using basic_ostream = contrib::basic_bgzf_ostream<char_t>;
#endif
};


/*!\brief Traits for the seqan3::compression_format::gz.
 * \see seqan3::compression_traits
 */
template <>
struct compression_traits<compression_format::gz> : compression_traits<compression_format::none>
{
    //!\copydoc seqan3::compression_traits<compression_format::none>::as_string
    static constexpr std::string_view as_string = "GZip";

    //!\copydoc seqan3::compression_traits<compression_format::none>::file_extensions
    static inline std::vector<std::string> file_extensions
    {
        {"gz"}
    };

    //!\copydoc seqan3::compression_traits<compression_format::none>::magic_header
    static constexpr std::string_view magic_header{"\x1f\x8b\x08", 3};

#ifdef SEQAN3_HAS_ZLIB
    //!\copydoc seqan3::compression_traits<compression_format::none>::available
    static constexpr bool available = true;

    //!\copydoc seqan3::compression_traits<compression_format::none>::basic_istream
    template <typename char_t>
    using basic_istream = contrib::basic_gz_istream<char_t>;

    //!\copydoc seqan3::compression_traits<compression_format::none>::basic_ostream
    template <typename char_t>
    using basic_ostream = contrib::basic_gz_ostream<char_t>;
#endif
};

/*!\brief Traits for the seqan3::compression_format::bz2.
 * \see seqan3::compression_traits
 */
template <>
struct compression_traits<compression_format::bz2> : compression_traits<compression_format::none>
{
    //!\copydoc seqan3::compression_traits<compression_format::none>::as_string
    static constexpr std::string_view as_string = "BZip2";

    //!\copydoc seqan3::compression_traits<compression_format::none>::file_extensions
    static inline std::vector<std::string> file_extensions
    {
        {"bz2"}
    };

    //!\copydoc seqan3::compression_traits<compression_format::none>::magic_header
    static constexpr std::string_view magic_header{"\x42\x5a\x68", 3};

#ifdef SEQAN3_HAS_BZIP2
    //!\copydoc seqan3::compression_traits<compression_format::none>::available
    static constexpr bool available = true;

    //!\copydoc seqan3::compression_traits<compression_format::none>::basic_istream
    template <typename char_t>
    using basic_istream = contrib::basic_bz2_istream<char_t>;

    //!\copydoc seqan3::compression_traits<compression_format::none>::basic_ostream
    template <typename char_t>
    using basic_ostream = contrib::basic_bz2_ostream<char_t>;
#endif
};

/*!\brief Traits for the seqan3::compression_format::zstd.
 * \see seqan3::compression_traits
 */
template <>
struct compression_traits<compression_format::zstd> : compression_traits<compression_format::none>
{
    //!\copydoc seqan3::compression_traits<compression_format::none>::as_string
    static constexpr std::string_view as_string = "ZStandard";

    //!\copydoc seqan3::compression_traits<compression_format::none>::file_extensions
    static inline std::vector<std::string> file_extensions
    {
        {"zstd"}
    };

    //!\copydoc seqan3::compression_traits<compression_format::none>::magic_header
    static constexpr std::string_view magic_header{"\x28\xb5\x2f\xfd", 4};

    // NOT YET SUPPORTED BY SEQAN
};

} // namespace seqan3

namespace seqan3::detail
{

//-------------------------------------------------------------------------------
// header_matches
//-------------------------------------------------------------------------------

//!\brief By default, the given argument is just compared with the stored magic header.
template <compression_format format>
constexpr bool header_matches(std::string_view to_compare)
{
    return to_compare.starts_with(compression_traits<format>::magic_header);
}

//!\brief For seqan3::compression_format::bgzf not all values are compared.
template <>
constexpr bool header_matches<compression_format::bgzf>(std::string_view to_compare)
{
    std::string_view magic_header = compression_traits<compression_format::bgzf>::magic_header;

    return (to_compare[0] == magic_header[0] &&                                                          // GZ_ID1
            to_compare[1] == magic_header[1] &&                                                          // GZ_ID2
            to_compare[2] == magic_header[2] &&                                                          // GZ_CM
            (to_compare[3] & magic_header[3]) != 0 &&                                                    // FLG_FEXTRA
            to_little_endian(*reinterpret_cast<uint16_t const *>(&to_compare[10])) == magic_header[10] &&// BGZF_ID1
            to_compare[12] == magic_header[12] &&                                                        // BGZF_ID2
            to_compare[13] == magic_header[13] &&                                                        // BGZF_SLEN
            to_little_endian(*reinterpret_cast<uint16_t const *>(&to_compare[14])) == magic_header[14]); // BGZF_XLEN
}

//!\brief Header matches "none" if it doesn't match any known compression.
template <>
constexpr bool header_matches<compression_format::none>(std::string_view to_compare)
{
    return !(header_matches<compression_format::bgzf>(to_compare)  ||
             header_matches<compression_format::gz>(to_compare)    ||
             header_matches<compression_format::bz2>(to_compare)   ||
             header_matches<compression_format::zstd>(to_compare));
}

//!\brief Header matches "none" if it doesn't match any known compression.
template <>
constexpr bool header_matches<compression_format::detect>(std::string_view)
{
    return true;
}

//-------------------------------------------------------------------------------
// header_matches_dyn
//-------------------------------------------------------------------------------

inline bool header_matches_dyn(compression_format format, std::string_view to_compare)
{
    switch (format)
    {
        case compression_format::none:   return header_matches<compression_format::none>(to_compare);
        case compression_format::detect: return header_matches<compression_format::detect>(to_compare);
        case compression_format::bgzf:   return header_matches<compression_format::bgzf>(to_compare);
        case compression_format::gz:     return header_matches<compression_format::gz>(to_compare);
        case compression_format::bz2:    return header_matches<compression_format::bz2>(to_compare);
        case compression_format::zstd:   return header_matches<compression_format::zstd>(to_compare);
    }

    return false;
}

//-------------------------------------------------------------------------------
// read_magic_header
//-------------------------------------------------------------------------------

inline std::string read_magic_header(std::istream & istr)
{
    std::string ret;

    assert(istr.good());
    std::istreambuf_iterator<char> it{istr};

    size_t max_chars = compression_traits<compression_format::bgzf>::magic_header.size();
    while (ret.size() < max_chars && it != std::istreambuf_iterator<char>{})
    {
        ret.push_back(*it);
        ++it;
    }

    // unget all read chars.
    for (size_t i = 0 ; i < ret.size(); ++i)
        istr.unget();

    return ret;
}

//-------------------------------------------------------------------------------
// detect_format_from_magic_header
//-------------------------------------------------------------------------------

/*!\brief Deduce seqan3::compression_format from a magic header string.
 *
 * \details
 *
 * Note that this function checks BGZF before GZ, since the latter's magic header is a prefix of the former's.
 */
inline compression_format detect_format_from_magic_header(std::string_view magic_header)
{
    if (header_matches<compression_format::bgzf>(magic_header))
        return compression_format::bgzf;
    else if (header_matches<compression_format::gz>(magic_header))
        return compression_format::gz;
    else if (header_matches<compression_format::bz2>(magic_header))
        return compression_format::bz2;
    else if (header_matches<compression_format::zstd>(magic_header))
        return compression_format::zstd;
    else
        return compression_format::none;
}

//-------------------------------------------------------------------------------
// detect_format_from_filename
//-------------------------------------------------------------------------------

/*!\brief Deduce seqan3::compression_format from a filename extension.
 *
 * \details
 *
 * Note that this function checks BGZF before GZ which means that it always selects BGZF for the extension ".gz".
 * This is desired, because many biological formats expect this.
 */
inline compression_format detect_format_from_filename(std::filesystem::path const & path)
{
    auto ext = path.extension().string();

    if (!ext.starts_with('.'))
        return compression_format::none;
    else
        ext = ext.substr(1);

    for (std::string_view ext2 : compression_traits<compression_format::bgzf>::file_extensions)
        if (ext == ext2)
            return compression_format::bgzf;
    for (std::string_view ext2 : compression_traits<compression_format::gz>::file_extensions)
        if (ext == ext2)
            return compression_format::gz;
    for (std::string_view ext2 : compression_traits<compression_format::bz2>::file_extensions)
        if (ext == ext2)
            return compression_format::bz2;
    for (std::string_view ext2 : compression_traits<compression_format::zstd>::file_extensions)
        if (ext == ext2)
            return compression_format::zstd;
    for (std::string_view ext2 : compression_traits<compression_format::bgzf>::file_extensions)
        if (ext == ext2)
            return compression_format::bgzf;

    return compression_format::none;
}


//-------------------------------------------------------------------------------
// make_istream
//-------------------------------------------------------------------------------

template <compression_format format, typename ... args_t>
std::istream * make_istream(args_t && ... args)
{
    if constexpr (compression_traits<format>::available)
    {
        return new typename compression_traits<format>::template basic_istream<char>{std::forward<args_t>(args)...};
    }
    else
    {
        throw file_open_error{std::string{"The file is "} + std::string{compression_traits<format>::as_string} +
                              "-compressed, but SeqAn3 is built without support for this format."};
        return nullptr;
    }
}

//-------------------------------------------------------------------------------
// make_ostream
//-------------------------------------------------------------------------------

template <compression_format format, typename ... args_t>
std::ostream * make_ostream(args_t && ... args)
{
    if constexpr (compression_traits<format>::available)
    {
        return new typename compression_traits<format>::template basic_ostream<char>{std::forward<args_t>(args)...};
    }
    else
    {
        throw file_open_error{std::string{compression_traits<format>::as_string} + "-compression was selected, "
                              "but SeqAn3 is built without support for this format."};
        return nullptr;
    }
}

} //namespace seqan3::detail
