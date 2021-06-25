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

#include <string>
#include <string_view>
#include <map>
#include <vector>

#include <seqan3/io/utility.hpp>
#include <seqan3/utility/views/eager_split.hpp>

namespace seqan3::var_io
{
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
        init();
    }
    header(header const &)                = default;
    header(header &&)                     = default;
    header & operator=(header const &)    = default;
    header & operator=(header &&)         = default;
    ~header()                             = default;

    header(std::string plaintext_header)
        : raw_data{std::move(plaintext_header)}
    {
        init();
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

        std::vector<std::string_view>                samples;               //!< IDs of samples (genotyping).
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
        append_to_header(l);
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
    //!\brief Add a contig.
    void add_contig(contig_t const & contig)
    {
        unparse_contig(contig);
    }

    //!\overload
    void add_contig(std::string_view const contig_name)
    {
        unparse_contig({ .id = contig_name });
    }

    //!\brief Add INFO.
    void add_info(info_t info)
    {
        unparse_info(info);
    }

    //!\overload
    void add_info(std::string_view const info_name)
    {
        unparse_info({ .id = info_name });
    }

    //!\brief Add a FILTER.
    void add_filter(filter_t const & filter)
    {
        unparse_filter(filter);
    }

    //!\overload
    void add_filter(std::string_view const filter_name)
    {
        unparse_filter({ .id = filter_name });
    }

    //!\brief Add FORMAT.
    void add_format(format_t format)
    {
        unparse_format(format);
    }

    //!\overload
    void add_format(std::string_view const format_name)
    {
        unparse_format({ .id = format_name });
    }

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

    void init()
    {
        /* needs to be first entry! */
        filter_t pass_filter{"PASS", "\"All filters passed\""};
        parsed_data.filter_id_to_index["PASS"] = 0;
        parsed_data.filters.push_back(std::move(pass_filter));
    }

//     bool append_to_header(auto && ... args)
//     {
//         size_t added_size = (0 + ... + args.size());
//         size_t old_capacity = raw_data.capacity();
//
//         ((raw_data += args), ...);
//         raw_data += '\n';
//
//         if (old_capacity != raw_data.capacity()) // string_views were invalidated
//         {
//             std::cerr << "REGENERATING HEADER\n"; // TODO DEBUG
//
//             *this = header{std::move(raw_data)}; // regenerate hash-maps
//             return true;
//         }
//         else
//         {
//             return false;
//         }
//     }

    //!\brief Append to the raw_data and regenerate hash-maps if necessary. String may not contain newlines!
    //!\tparam param_t Either std::string_view/std::string or a range therof.
    template <typename param_t>
    void append_to_header(param_t const & string)
    {
        size_t old_size     = raw_data.size();
        size_t old_capacity = raw_data.capacity();

        if constexpr (std::same_as<param_t, std::string_view> || std::same_as<param_t, std::string>)
        {
            raw_data += string;
        }
        else // assuming range
        {
            for (std::string_view s : string)
                raw_data += s;
        }
        raw_data += '\n';

        size_t new_size     = raw_data.size();
        size_t new_capacity = raw_data.capacity();


        if (old_capacity != new_capacity) // string_views were invalidated
        {
            std::cerr << "REGENERATING HEADER\n"; // TODO DEBUG

            *this = header{std::move(raw_data)}; // regenerate hash-maps
        }
        else
        {
            parse_line(std::string_view{raw_data.data() + old_size, new_size - old_size - 1});
        }
    }

    void unparse_contig(contig_t const & contig)
    {
        std::vector<std::string_view> strings;
        std::string length_buffer;

        strings.push_back("##contig=<ID=");
        strings.push_back(contig.id);
        if (contig.length != -1)
        {
            strings.push_back(",length=");
            length_buffer = std::to_string(contig.length); //TODO replace with std::to_chars
            strings.push_back(length_buffer); // this would dangling, but append_to_header is called in-scope and copies anyway
        }

        for (auto [ key, value ] : contig.other_fields)
        {
            strings.push_back(",");
            strings.push_back(key);
            strings.push_back("=");
            strings.push_back(value);
        }
        strings.push_back(">");

        append_to_header(strings);
    }

    void unparse_info(info_t const & info)
    {
        // TODO implement
    }

    void unparse_filter(filter_t const & info)
    {
        // TODO implement
    }

    void unparse_format(format_t const & info)
    {
        // TODO implement
    }


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
        if (id.empty())
            throw format_error{"INFO or FORMAT line does not contain ID field."};
        else
            new_info.id = id.mapped();

        /* Number */
        auto number = new_info.other_fields.extract("Number");
        if (number.empty())
            throw format_error{"INFO or FORMAT line does not contain Number field."};
        else
            new_info.number = parse_number(number.mapped());

        /* Type */
        auto type = new_info.other_fields.extract("Type");
        if (type.empty())
            throw format_error{"INFO or FORMAT line does not contain Type field."};
        else
            new_info.type = parse_type(type.mapped(), new_info.number);

        /* Description */
        auto description = new_info.other_fields.extract("Description");
        if (description.empty())
            throw format_error{"INFO or FORMAT line does not contain Description field."};
        else
            new_info.description = description.mapped();

        if (is_info)
        {
            if (parsed_data.info_id_to_index.contains(new_info.id))
                throw format_error{std::string{"Duplicate INFO ID \""} + std::string{new_info.id} + "\" in HEADER."};

            parsed_data.info_id_to_index[new_info.id] = parsed_data.infos.size();
            parsed_data.infos.push_back(std::move(new_info));
        }
        else
        {
            if (parsed_data.format_id_to_index.contains(new_info.id))
                throw format_error{std::string{"Duplicate FORMAT ID \""} + std::string{new_info.id} + "\" in HEADER."};

            parsed_data.format_id_to_index[new_info.id] = parsed_data.formats.size();
            parsed_data.formats.push_back(std::move(new_info));
        }
    }

    void parse_filter_line(std::string_view l)
    {
        // Should never be 0, since we always add PASS entry in construction
        assert(parsed_data.filters.size() != 0);
        assert(parsed_data.filters[0].id == "PASS");

        filter_t new_filter;
        new_filter.other_fields = to_dictionary(l);

        /* ID */
        auto id = new_filter.other_fields.extract("ID");
        if (id.empty())
            throw format_error{"FILTER line does not contain ID field."};
        else
            new_filter.id = id.mapped();

        /* Description */
        auto description = new_filter.other_fields.extract("Description");
        if (description.empty())
            throw format_error{"FILTER line does not contain Description field."};
        else
            new_filter.description = description.mapped();

        // PASS line was added by us before and is now swapped with user-provided
        if (new_filter.id == "PASS")
        {
            std::swap(parsed_data.filters[0], new_filter);
        }
        else
        {
            if (parsed_data.filter_id_to_index.contains(new_filter.id))
                throw format_error{std::string{"Duplicate FILTER ID \""} + std::string{new_filter.id} + "\" in HEADER."};

            parsed_data.filter_id_to_index[new_filter.id] = parsed_data.filters.size();
            parsed_data.filters.push_back(std::move(new_filter));
        }
    }

    void parse_contig_line(std::string_view l)
    {
        contig_t new_contig;
        new_contig.other_fields = to_dictionary(l);

        /* ID */
        auto id = new_contig.other_fields.extract("ID");
        if (id.empty())
            throw format_error{"FILTER line does not contain ID field."};
        else
            new_contig.id = id.mapped();

        /* Length */
        auto length = new_contig.other_fields.extract("length");
        if (!length.empty())
        {
            int64_t ret = -1;
            std::string_view s = length.mapped();
            detail::string_to_number(s, ret);
            new_contig.length = ret;
        }

        if (parsed_data.contig_id_to_index.contains(new_contig.id))
                throw format_error{std::string{"Duplicate CONTIG ID \""} + std::string{new_contig.id} + "\" in HEADER."};

        parsed_data.contig_id_to_index[new_contig.id] = parsed_data.contigs.size();
        parsed_data.contigs.push_back(std::move(new_contig));
    }

    static inline std::string_view strip_angular_brackets(std::string_view const in)
    {
        if (in.size() < 2 || in.front() != '<' || in.back() != '>')
            throw format_error{"Structured line does not contain \"<\" and \">\" at right places."};
        return in.substr(1, in.size() - 2);
    }

    static int32_t parse_number(std::string_view const in)
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
                detail::string_to_number(in, ret);
                return ret;
            }
        }
        return header_number::dot;
    }

    static io_type_id parse_type(std::string_view const in, int32_t const number)
    {
        io_type_id ret{};

        if (in == "Flag")
        {
            ret = io_type_id::flag;
            if (number != 0)
                throw format_error{std::string{"Flags must always have number 0 in header."}};
            return ret;
        }

        if (number == 0)
            throw format_error{std::string{"Only flags may have number 0 in header."}};

        if (in == "Integer")
        {
            if (number == 1)
                ret = io_type_id::int32;
            else
                ret = io_type_id::vector_of_int32;
        }
        else if (in == "Float")
        {
            if (number == 1)
                ret = io_type_id::float32;
            else
                ret = io_type_id::vector_of_float32;
        }
        else if (in == "Character")
        {
            if (number == 1)
                ret = io_type_id::char8;
            else
                ret = io_type_id::vector_of_char8;
        }
        else if (in == "String")
        {
            if (number == 1)
                ret = io_type_id::string;
            else
                ret = io_type_id::vector_of_string;
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
