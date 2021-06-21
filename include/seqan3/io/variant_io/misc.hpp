// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the seqan3::var_io::tag_dictionary class and auxiliaries.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <deque>
#include <map>
#include <variant>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/core/concept/core_language.hpp>
#include <seqan3/core/detail/template_inspection.hpp>
#include <seqan3/io/record.hpp>
#include <seqan3/io/utility.hpp>
#include <seqan3/std/concepts>
#include <seqan3/utility/char_operations/predicate.hpp>
#include <seqan3/utility/views/eager_split.hpp>

namespace seqan3::var_io
{

//!\brief An enumerator denoting variant file special states.
//!\ingroup variant_io
enum class special_value
{
    missing,    //!< "."
    unknown,    //!< "*"
};

/*!\brief A type representing variant file alleles.
 * \ingroup variant_io
 *
 * \details
 *
 * Alleles in variant files are encoded as one of the following
 *
 *  1. A seqan3::var_io::special_value if missing/absent.
 *  2. A std::vector<seqan3::dna5> if a single character or sequence of DNA.
 *  3. A std::string if they are anything else (imprecise structural variant, breakpoint-string etc).
 */
using allele = std::variant<special_value, std::vector<dna5>, std::string>;

//!\brief Default fields for seqan3::var_io::reader_options.
//!\ingroup variant_io
inline constexpr auto default_field_ids = tag<field::chrom,
                                              field::pos,
                                              field::id,
                                              field::ref,
                                              field::alt,
                                              field::qual,
                                              field::filter,
                                              field::info,
                                              field::genotypes,
                                              field::header>;

//!\brief Scoped (but weakly typed) enum for "Number" special values in seqan3::var_io::header INFO fields.
//!\ingroup variant_io
struct header_number
{
    enum
    {
        A   = -1,   //!< One value per alternate allele.
        R   = -2,   //!< One value for each possible allele (including ref) -> A + 1.
        G   = -3,   //!< One value per Genotype.
        dot = -4    //!< Unknown, unspecified or unbounded.
    };
};

/*!\brief Stores the header information of alignment files.
 * \ingroup variant_io
 */
class header
{
public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    header()
    {
        raw_data.reserve(100*1024); // reserve 100KB
    }
    header(header const &)                = default;
    header(header &&)                     = default;
    header & operator=(header const &)    = default;
    header & operator=(header &&)         = default;
    ~header()                             = default;

    header(std::string plaintext_header)
        : raw_data{std::move(plaintext_header)}
    {
        for (std::string_view line : raw_data | views::eager_split('\n'))
            parse_line(line);
    }

    //!\}

    /*!\name Member types
     * \{
     */
    //!\brief Type of the contig field header line.
    struct contig_t
    {
        std::string_view                             id;                //!< The ID.
        int64_t                                      length = -1;       //!< Length of the contig (-1 if absent).
        std::map<std::string_view, std::string_view> other_fields;      //!< Other entries.
    };

    //!\brief Type of a INFO field header line.
    struct info_t
    {
        std::string_view                             id;                //!< The ID.
        //!\brief Number of values, see also seqan3::var_io::header_number.
        int32_t                                      number;
        io_type_id                                   type;              //!< Type of the field.
        std::string_view                             description;       //!< Description.
        std::map<std::string_view, std::string_view> other_fields;      //!< Other entries.
    };

    //!\brief Type of a FILTER field header line.
    struct filter_t
    {
        std::string_view                             id;                //!< The ID.
        std::string_view                             description;       //!< Description.
        std::map<std::string_view, std::string_view> other_fields;      //!< Other entries.
    };

    using format_t = info_t;                                            //!< Type of a FORMAT field header line.

    //!\brief A datastructure that contains all parsed data fields as views into the raw header.
    struct parsed_data_t
    {
        std::string_view                             file_format;           //!< The file format version.
        std::vector<contig_t>                        contigs;               //!< Header lines describing contigs.
        std::unordered_map<std::string_view, size_t> contig_id_to_index;    //!< ID-string to position in #contigs.
        std::vector<info_t>                          infos;                 //!< Header lines describing INFO fields.
        std::unordered_map<std::string_view, size_t> info_id_to_index;      //!< ID-string to position in #infos.
        std::vector<filter_t>                        filters;               //!< Header lines describing FILTER fields.
        std::unordered_map<std::string_view, size_t> filter_id_to_index;    //!< ID-string to position in #filters.
        std::vector<format_t>                        formats;               //!< Header lines describing FORMAT fields.
        std::unordered_map<std::string_view, size_t> format_id_to_index;    //!< ID-string to position in #formats.
        std::vector<std::string_view>                other_lines;           //!< Any other lines in the header.
    };

    //!\}

    /*!\name Functions for adding data
     * \brief It is recommended to only use these functions when creating a header "from scratch".
     * \{
     */

    //!\brief Add a a line to the header and parse it.
    void add_raw_line(std::string_view l)
    {
        if (!append_to_header(l)) // if hash-maps are regenerated, the line is already parsed
            parse_line(std::string_view{raw_data.data() + raw_data.size() - l.size() - 1, l.size()});
    }

//     //!\brief Add a contig.
//     void set_file_format(std::string_view f)
//     {
//         if (!data.empty())
//             throw std::runtime_error{"The file format must be set before everything else."}
//
//         unparse_file_format(f);
//         file_format_read = true;
//     }
//
//     //!\brief Add INFO.
//     void add_info(info_t c)
//     {
//         infos.push_back(std::move(c));
//         info_id_to_index[infos.back().id] = infos.size() - 1;
//     }
//
//     //!\brief Add a contig.
//     void add_contig(contig_t c)
//     {
//         contigs.push_back(std::move(c));
//         contig_id_to_index[contigs.back().id] = contigs.size() - 1;
//     }



//     void clear()
//     {
//         raw_data.clear();
//         file_format.clear();
//         contigs.clear();
//         contig_id_to_index.clear();
//         infos.clear();
//         info_id_to_index.clear();
//         filters.clear();
//         filter_name_to_index.clear();
//         formats.clear();
//         format_name_to_index.clear();
//         other_lines.clear();
//
//         file_format_read = false;
//     }

    std::string_view raw_header() const
    {
        return raw_data;
    }

    parsed_data_t const & parsed_header() const
    {
        return parsed_data;
    }

protected:
    std::string   raw_data;     //!< The header in plaintext-form.
    parsed_data_t parsed_data;  //!< The parsed data

private:
    bool file_format_read = false;

    //!\brief Append to the raw_data and regenerate hash-maps if necessary.
    //!\returns whether hash-maps where regenerated.
    bool append_to_header(auto && ... args)
    {
        size_t added_size = (0 + ... + args.size());
        size_t old_capacity = raw_data.capacity();

        ((raw_data += args), ...);
        raw_data += '\n';

        if (old_capacity != raw_data.capacity()) // string_views were invalidated
        {
            std::cerr << "REGENERATING HEADER\n"; // TODO DEBUG

            *this = header{std::move(raw_data)}; // regenerate hash-maps
            return true;
        }
        else
        {
            return false;
        }
    }

//     void unparse_file_format(std::string_view new_value)
//     {
//         append_to_header(std::string_view{"##fileformat="}, new_value);
//     }


    void parse_line(std::string_view l)
    {
        if (file_format_read == false)
        {
            if (l.starts_with("##fileformat="))
            {
                parsed_data.file_format = l.substr(13);
                file_format_read = true;
            }
            else
            {
                throw format_error{"File does not begin with \"##fileformat\"."};
            }
        }
        else if (l.starts_with("##fileformat="))
        {
            throw format_error{"File has two lines that begin with \"##fileformat\"."};
        }
        else if (l.starts_with("##INFO="))
        {
            parse_info_or_format_line(strip_angular_brackets(l.substr(7)), true);
        }
        else if (l.starts_with("##FILTER="))
        {
            parse_filter_line(strip_angular_brackets(l.substr(9)));
        }
        else if (l.starts_with("##FORMAT="))
        {
            parse_info_or_format_line(strip_angular_brackets(l.substr(9)), false);
        }
        else if (l.starts_with("##contig="))
        {
            parse_contig_line(strip_angular_brackets(l.substr(9)));
        }
        else
        {
            parsed_data.other_lines.push_back(l);
        }

    }

    void parse_info_or_format_line(std::string_view l, bool is_info)
    {
        info_t new_info;
        new_info.other_fields = to_dictionary(l);

        /* ID */
        auto id = new_info.other_fields.extract("ID");
        if (id.key() != "ID")
            throw format_error{"INFO or FORMAT line does not contain ID field."};
        else
            new_info.id = id.mapped();

        /* Number */
        auto number = new_info.other_fields.extract("Number");
        if (number.key() != "Number")
            throw format_error{"INFO or FORMAT line does not contain Number field."};
        else
            new_info.number = parse_number(number.mapped());

        /* Type */
        auto type = new_info.other_fields.extract("Type");
        if (type.key() != "Type")
            throw format_error{"INFO or FORMAT line does not contain Type field."};
        else
            new_info.type = parse_type(type.mapped());

        /* Description */
        auto description = new_info.other_fields.extract("Description");
        if (description.key() != "Description")
            throw format_error{"INFO or FORMAT line does not contain Description field."};
        else
            new_info.description = description.mapped();

        if (is_info)
            parsed_data.infos.push_back(std::move(new_info));
        else
            parsed_data.formats.push_back(std::move(new_info));
    }

    void parse_filter_line(std::string_view l)
    {
        filter_t new_filter;
        new_filter.other_fields = to_dictionary(l);

        /* ID */
        auto id = new_filter.other_fields.extract("ID");
        if (id.key() != "ID")
            throw format_error{"FILTER line does not contain ID field."};
        else
            new_filter.id = id.mapped();

        /* Description */
        auto description = new_filter.other_fields.extract("Description");
        if (description.key() != "Description")
            throw format_error{"FILTER line does not contain Description field."};
        else
            new_filter.description = description.mapped();

        parsed_data.filters.push_back(std::move(new_filter));
    }

    void parse_contig_line(std::string_view l)
    {
        contig_t new_contig;
        new_contig.other_fields = to_dictionary(l);

        /* ID */
        auto id = new_contig.other_fields.extract("ID");
        if (id.key() != "ID")
            throw format_error{"FILTER line does not contain ID field."};
        else
            new_contig.id = id.mapped();

        /* Length */
        auto description = new_contig.other_fields.extract("length");
        if (description.key() == "length")
        {
            int64_t ret = -1;
            std::string_view s = description.mapped();

            if (auto tmp = std::from_chars(s.data(), s.data() + s.size(), ret);
                tmp.ec != std::errc{} || tmp.ptr != s.data() + s.size())
            {
                throw format_error{std::string{"Cannot convert the following string to a number: "} + std::string{s}};
            }

            new_contig.length = ret;
        }

        parsed_data.contigs.push_back(std::move(new_contig));
    }

    static inline std::string_view strip_angular_brackets(std::string_view in)
    {
        if (in.size() < 2 || in.front() != '<' || in.back() != '>')
            throw format_error{"Structured line does not contain \"<\" and \">\" at right places."};
        return in.substr(1, in.size() - 2);
    }

    static int32_t parse_number(std::string_view in)
    {

        switch (in[0])
        {
            case 'A': return header_number::A;
            case 'R': return header_number::R;
            case 'G': return header_number::G;
            case '.': return header_number::dot;
            default:
            {
                int32_t ret = 0;
                auto tmp = std::from_chars(in.data(), in.data() + in.size(), ret);
                if (tmp.ec != std::errc{} || tmp.ptr != in.data() + in.size())
                {
                    throw format_error{std::string{"Cannot convert the following string to a number: "} +
                                       std::string{in}};
                }
                return ret;
            }
        }
        return header_number::dot;
    }

    static io_type_id parse_type(std::string_view in)
    {
        io_type_id ret{};

        if (in == "Integer")
        {
            ret = io_type_id::int32;
        }
        else if (in == "Float")
        {
            ret = io_type_id::float32;
        }
        else if (in == "Flag")
        {
            ret = io_type_id::flag;
        }
        else if (in == "Character")
        {
            ret = io_type_id::char8;
        }
        else if (in == "String")
        {
            ret = io_type_id::string;
        }
        else
        {
            throw format_error{std::string{"Cannot convert the following string to a type identifier: "} +
                               std::string{in}};
        }

        return ret;
    }

    static std::map<std::string_view, std::string_view> to_dictionary(std::string_view value_pairs)
    {
        std::map<std::string_view, std::string_view> ret;

        for (std::string_view pair : value_pairs | views::eager_split(',', true))
        {
            auto pair_split = pair | views::eager_split('=');
            auto it1 = pair_split.begin();
            auto it2 = std::ranges::next(it1);
            auto it3 = std::ranges::next(it2);

            if (it1 == std::default_sentinel || it2 == std::default_sentinel || it3 != std::default_sentinel)
            {
                throw format_error{std::string{"Could not parse the following string into a dictionary: "} +
                                   std::string{pair}};
            }

            ret[*it1] = *it2;
        }

        return ret;
    }



};

} // namespace seqan3::var_io
