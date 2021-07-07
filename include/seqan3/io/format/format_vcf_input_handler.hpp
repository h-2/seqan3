// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * brief Provides the seqan3::input_format_handler<seqan3::format_vcf> .
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <iostream>
#include <regex>
#include <string>
#include <string_view>
#include <vector>

#include <seqan3/alphabet/views/char_strictly_to.hpp>
#include <seqan3/alphabet/views/char_to.hpp>
#include <seqan3/core/range/type_traits.hpp>
#include <seqan3/core/debug_stream/detail/to_string.hpp>
#include <seqan3/core/debug_stream.hpp>                 //TODO evaluate if there is a besser solution
#include <seqan3/io/detail/ignore_output_iterator.hpp>
#include <seqan3/io/detail/misc.hpp>
#include <seqan3/io/format/input_format_handler_base.hpp>
#include <seqan3/io/format/format_vcf.hpp>
#include <seqan3/io/plaintext_io/reader.hpp>
#include <seqan3/io/stream/iterator.hpp>
#include <seqan3/io/variant_io/header.hpp>
#include <seqan3/io/variant_io/misc.hpp>
#include <seqan3/range/detail/misc.hpp>
#include <seqan3/std/algorithm>
#include <seqan3/std/ranges>
#include <seqan3/utility/type_list/traits.hpp>
#include <seqan3/utility/views/to.hpp>

namespace seqan3
{

template <>
class input_format_handler<format_vcf> : public input_format_handler_base<input_format_handler<format_vcf>>
{
private:
    /* CRTP */
    using base_t = input_format_handler_base<input_format_handler<format_vcf>>;
    using base_t::stream;

    friend base_t;

    /* RAW RECORD HANDLING*/
    using format_fields     = decltype(var_io::default_field_ids);
    using raw_record_type   = seqan3::record<list_traits::repeat<format_fields::size, std::string_view>,
                                             format_fields>;

    raw_record_type raw_record;
    plain_io::detail::plaintext_input_iterator<char, std::char_traits<char>, plain_io::record_kind::line_and_fields> file_it;
    std::string_view last_chrom;
    size_t           last_chrom_index = -1;
    size_t           line = 0;



    [[noreturn]] void error(auto const & ... messages)
    {
        std::string message = "[SeqAn3 VCF format error in line " + detail::to_string(line) + "]";
        message += detail::to_string(messages...);

        throw format_error{message};
    }

    void read_raw_record()
    {
        ++line;
        ++file_it;

        get<field::chrom> (raw_record) = (*file_it).fields[0];
        get<field::pos>   (raw_record) = (*file_it).fields[1];
        get<field::id>    (raw_record) = (*file_it).fields[2];
        get<field::ref>   (raw_record) = (*file_it).fields[3];
        get<field::alt>   (raw_record) = (*file_it).fields[4];
        get<field::qual>  (raw_record) = (*file_it).fields[5];
        get<field::filter>(raw_record) = (*file_it).fields[6];
        get<field::info>  (raw_record) = (*file_it).fields[7];

        // fields[7].end() that is guaranteed to be char*
        char const * end_qual   = (*file_it).fields[7].data() + (*file_it).fields[7].size();
        // line.end() that is guaranteed to be char*
        char const * end_line   = (*file_it).line.data()      + (*file_it).line.size();
        // genotypes fo from end of qual til end of line (possibly empty)
        get<field::genotypes>  (raw_record) = std::string_view{end_qual, static_cast<size_t>(end_line - end_qual)};

        // header does not need to be reset
    }



    /* OPTIONS */

    bool warn_on_missing_header_entries = false;

    /* PARSED RECORD HANDLING */
    var_io::header header;

    //!\brief Default field handlers.
    using base_t::parse_field;
    using base_t::parse_field_impl;

    //!\brief Parse the CHROM field. Reading chrom as number means getting the index (not converting string to number).
    template <std::integral parsed_field_t>
    void parse_field(tag_t<field::chrom> const & /**/,  parsed_field_t & parsed_field)
    {
        std::string_view raw_field = get<field::chrom>(raw_record);

        if (raw_field == last_chrom) // reuse cached index
        {
            parsed_field = static_cast<parsed_field_t>(last_chrom_index);
        }
        else
        {
            if (auto it = header.parsed_header().contig_id_to_index.find(raw_field);
                it == header.parsed_header().contig_id_to_index.end()) // contig name was not in header, insert!
            {
                header.add_contig(raw_field);
                last_chrom_index = header.parsed_header().contigs.size() - 1;
                parsed_field = static_cast<parsed_field_t>(last_chrom_index);

                if (warn_on_missing_header_entries)
                {
                    debug_stream << "[seqan3::var_io::warning] The contig name \"" << raw_field << "\" found on line "
                                 << line << "was not present in the header.\n";
                }
            }
            else
            {
                parsed_field = static_cast<parsed_field_t>(it->second);
                last_chrom_index = static_cast<size_t>(it->second);
            }

            last_chrom = raw_field;
        }
    }

    /* POS, ID, REF are handled correctly by default */

    //!\brief Overload for parsing ALT.
    template <typename parsed_field_t>
        requires detail::back_insertable<parsed_field_t> &&
                 std::ranges::range<std::ranges::range_reference_t<parsed_field_t>>
    void parse_field(tag_t<field::alt> const & /**/, parsed_field_t & parsed_field)
    {
        std::string_view raw_field = get<field::alt>(raw_record);

        if (raw_field != ".")
        {
            for (std::string_view subfield : raw_field | views::eager_split(','))
            {
                parsed_field.emplace_back();
                // delegate parsing of element to base
                parse_field_impl(subfield, parsed_field.back());
            }
        }
    }

    //!\brief Overload for parsing QUAL.
    template <arithmetic parsed_field_t>
    void parse_field(tag_t<field::qual> const & /**/, parsed_field_t & parsed_field)
    {
        std::string_view raw_field = get<field::qual>(raw_record);

        if (raw_field == ".")
        {
            parsed_field = var_io::missing_value<parsed_field_t>;
        }
        else
        {
            parsed_field = parsed_field_t{};
            parse_field_impl(raw_field, parsed_field); // arithmetic parsing
        }
    }

    //!\brief Overload for parsing FILTER as strings.
    template <typename parsed_field_t>
        requires (detail::back_insertable<parsed_field_t> &&
                  is_one_of<std::ranges::range_value_t<parsed_field_t>, std::string, std::string_view>)
    void parse_field(tag_t<field::filter> const & /**/, parsed_field_t & parsed_field)
    {
        using elem_t = std::ranges::range_value_t<parsed_field_t>;
        std::string_view raw_field = get<field::filter>(raw_record);

        if (raw_field == ".")
            return;

        for (std::string_view subfield : raw_field | views::eager_split(';'))
            parsed_field.push_back(elem_t{subfield});
    }

    //!\brief Overload for parsing FILTER as indexes.
    template <typename parsed_field_t>
    requires (detail::back_insertable<parsed_field_t> &&
              std::integral<std::ranges::range_value_t<parsed_field_t>>)
    void parse_field(tag_t<field::filter> const & /**/, parsed_field_t & parsed_field)
    {
        using elem_t = std::ranges::range_value_t<parsed_field_t>;
        std::string_view raw_field = get<field::filter>(raw_record);

        if (raw_field == ".")
            return;

        for (std::string_view subfield : raw_field | views::eager_split(';'))
        {
            if (auto it = header.parsed_header().filter_id_to_index.find(raw_field);
                     it == header.parsed_header().filter_id_to_index.end()) // filter name was not in header, insert!
            {
                header.add_filter(subfield);
                parsed_field.push_back(static_cast<elem_t>(header.parsed_header().filters.size()));

                if (warn_on_missing_header_entries)
                {
                    debug_stream << "[seqan3::var_io::warning] The filter name \"" << subfield << "\" found on line "
                                 << line << "was not present in the header.\n";
                }

            }
            else
            {
                parsed_field.push_back(static_cast<elem_t>(it->second));
            }
        }
    }

    //!\brief Overload for parsing INFO.
    //TODO make this back_insertable
    template <typename key_t, typename value_t>
    void parse_field(tag_t<field::info> const & /**/, std::vector<std::pair<key_t, value_t>> & parsed_field)
    {
        static_assert(is_one_of<key_t, int32_t, std::string, std::string_view>,
                      "The key type for reading INFO fields must be int32_t, std::string or std::string_view.");
        static_assert(is_one_of<value_t, io_type_variant<true>, io_type_variant<false>, std::string, std::string_view>,
                      "The value type for reading INFO fields must be int32_t, std::string or std::string_view.");
        static_assert(!is_one_of<value_t, io_type_variant<true>, io_type_variant<false>> || std::same_as<key_t, int32_t>,
                      "The key type for reading INFO fields may only be io_type_variant if the key type is int32_t.");

        std::string_view raw_field = get<field::info>(raw_record);

        if (raw_field == ".")
            return;

        for (std::string_view subfield : raw_field | views::eager_split(';'))
        {
            std::pair<key_t, value_t> ret{};

            auto key_value_view = subfield | views::eager_split('=');
            auto key_it = key_value_view.begin();
            auto val_it = std::ranges::next(key_it);
            auto post_it = std::ranges::next(val_it);

            if (key_it == key_value_view.end() || (val_it != key_value_view.end() && post_it != key_value_view.end()))
                error("Could not read INFO fields from the following string: ", subfield);

            /* PARSE KEY */
            std::string_view key = *key_it;
            if constexpr (std::same_as<key_t, int32_t>)
            {
                if (auto it = header.parsed_header().info_id_to_index.find(key);
                    it == header.parsed_header().filter_id_to_index.end()) // info name was not in header, insert!
                {
                    var_io::header::info_t info;
                    info.id = key;
                    if (val_it == key_value_view.end()) // no "=" → flag
                    {
                        info.type   = io_type_id::flag;
                        info.number = 0;
                    }
                    else if (val_it->find(',') != std::string_view::npos) // found comma → assume vector-of-strings
                    {
                        info.type   = io_type_id::vector_of_string;
                        info.number = var_io::header_number::dot;
                    }
                    else // assume string as type
                    {
                        info.type   = io_type_id::string;
                        info.number = 1;
                    }

                    header.add_info(info);

                    ret.first = static_cast<key_t>(header.parsed_header().infos.size());

                    if (warn_on_missing_header_entries)
                    {
                        debug_stream << "[seqan3::var_io::warning] The INFO name \"" << raw_field << "\" found on line "
                                     << line << "was not present in the header.\n";
                    }
                }
                else
                {
                    ret.first = static_cast<key_t>(it->second);
                }
            }
            else // key_t is std::string or std::string_view
            {
                ret.first = key_t{key};
            }

            /* PARSE VALUE */
            if (val_it == key_value_view.end()) // no "=" → flag
            {
                if constexpr (is_one_of<value_t, io_type_variant<true>, io_type_variant<false>>)
                {
                    if (header.parsed_header().infos.back().type != io_type_id::flag ||
                        header.parsed_header().infos.back().number != 0)
                    {
                        error("INFO field \"", key, "\" is not a flag and should come with a value -- but does not.");
                    }

                    ret.second = true; // set variant to boolean/flag state
                }
                // no-op for std::string
            }
            else // any other type than flag
            {
                if constexpr (is_one_of<value_t, io_type_variant<true>, io_type_variant<false>>)
                {
                    size_t num_val = detail::parse_io_type_variant(header.parsed_header().infos[ret.first].type,
                                                                   *val_it,
                                                                   ret.second);
                    if (size_t exp_val = header.parsed_header().infos[ret.first].number;
                        warn_on_missing_header_entries && num_val != exp_val)
                    {
                        debug_stream << "[seqan3::var_io::warning] Expected to find "
                                     << exp_val
                                     << "values for the INFO field "
                                     << raw_field
                                     << "but found: "
                                     << num_val
                                     << "\n";
                    }
                }
                else // val_t == std::string or std::string_view
                {
                    ret.second = value_t{*val_it};
                }
            }

            parsed_field.push_back(std::move(ret));
        }
    }

//     template <typename key_t, typename value_t>
//     void parse_field(tag_t<field::genotypes> const & /**/, std::vector<std::pair<key_t, value_t>> & parsed_field)
//     {
//         static_assert(is_one_of<key_t, int32_t, std::string, std::string_view>,
//                       "The key type for reading INFO fields must be int32_t, std::string or std::string_view.");
//         static_assert(detail::is_io_type_vector_variant<value_t> ||
//                       is_one_of<value_t, std::vector<std::string>, std::vector<std::string_view>>,
//                       "The value type for reading INFO fields must be io_type_vector or a vector of std::string or "
//                       "vector of std::string_view.");
//         static_assert(!detail::is_io_type_vector_variant<key_t> || std::same_as<key_t, int32_t>,
//                       "The key type for reading INFO fields may only be io_type_variant if the key type is int32_t.");
//
//         size_t column_number          = file_it->fields.size();
//         size_t expected_column_number = header.parsed_header().samples.size() + 8;
//
//         if (column_number > 8) // there are genotypes
//         {
//             if (column_number != expected_column_number)
//                 error("Expected ", expected_column_number, " columns in line but found", column_number, ".");
//
//             std::string_view format_names = file_it->fields[8];
//
//             /* parse keys */
//             for (std::string_view format_name : format_names | views::eager_split(':'))
//             {
//                 if constexpr (std::same_as<key_t, int32_t>)
//                 {
//                     int32_t format_index = -1;
//                     if (auto it = header.parsed_header().format_id_to_index.find(format_name);
//                         it == header.parsed_header().format_id_to_index.end()) // format name was not in header, insert!
//                     {
//                         //TODO
//                     }
//                     else
//                     {
//                         format_index = it->second;
//                     }
//
//                     parsed_field.push_back({format_index});
//
//                     header::format_t format const & header.parsed_header().formats[format_index];
//
//                     if constexpr (detail::is_io_type_vector_variant<value_t>)
//                     {
//                         detail::init_io_type_variant(format.type, parsed_field.back().second);
//                         auto reserve = [s = column_number - 8] (auto & vec) { vec.reserve(s); };
//                         std::visit(reserve, parsed_field.back().second);
//                     }
//                     else
//                     {
//                         parsed_field.back().second.reserve(column_number - 8);
//                     }
//
//                 }
//                 else
//                 {
//                     parsed_field.push_back(key_t{format_name});
//                 }
//             }
//
//             /* parse values/samples */
//             for (size_t i = 9; i < column_number; ++i)
//             {
//                 std::string_view sample = file_it->fields[i];
//
//                 if constexpr (std::same_as<key_t, int32_t>)
//                 {
//
//                 }
//
//                 size_t field_no = 0;
//                 for (std::string_view field : sample | views::eager_split(':'))
//                 {
//
//             }
//
//                 header::format_t format const & header.parsed_header().formats[format_index];
//
//
//                 auto resize = [s = header.parsed_header().samples.count()] (auto & vec) { vec.resize(s); };
//
//
//         }
//
//     }

    template <typename value_t>
        requires detail::is_io_type_vector_variant<value_t>
    void parse_field(tag_t<field::genotypes> const & /**/, std::vector<std::pair<int32_t, value_t>> & parsed_field)
    {
        size_t column_number          = file_it->fields.size();
        size_t expected_column_number = header.parsed_header().samples.size() ?
                                        header.parsed_header().samples.size() + 9 : 8;

        if (column_number > 8) // there are genotypes
        {
            if (column_number != expected_column_number)
                error("Expected ", expected_column_number, " columns in line but found", column_number, ".");

            std::string_view format_names = file_it->fields[8];

            /* parse keys */
            size_t formats = 0;
            for (std::string_view format_name : format_names | views::eager_split(':'))
            {
                int32_t format_index = -1;
                if (auto it = header.parsed_header().format_id_to_index.find(format_name);
                    it == header.parsed_header().format_id_to_index.end()) // format name was not in header, insert!
                {
                    //TODO
                }
                else
                {
                    format_index = it->second;
                }

                parsed_field.push_back({format_index, {}});

                var_io::header::format_t const & format = header.parsed_header().formats[format_index];

                detail::init_io_type_variant(format.type, parsed_field.back().second);
                auto reserve = [s = column_number - 8] (auto & vec) { vec.reserve(s); };
                std::visit(reserve, parsed_field.back().second);

                ++formats;
            }

            /* parse values/samples */
            for (size_t i = 9; i < column_number; ++i)
            {
                std::string_view sample = file_it->fields[i];

                auto fields_view = sample | views::eager_split(':');
                auto fields_it = std::ranges::begin(fields_view);

                for (size_t i = 0; i < formats; ++i)
                {
                    std::string_view field{};

                    if (fields_it != std::default_sentinel)
                    {
                        field = *fields_it;
                        ++fields_it;
                    }
                    else // this handles trailing dropped fields
                    {
                        field = ".";
                    }

                    auto parse_and_append = [field, i] (auto & variant)
                    {
                        if constexpr (std::same_as<std::ranges::range_value_t<decltype(variant)>, bool>)
                            variant.push_back(true);
                        else
                            variant.emplace_back();

                        detail::parse_io_type_data_t{field}(variant.back());
                    };
                    std::visit(parse_and_append, parsed_field[i].second);
                }
            }
        }
    }


    void parse_field(tag_t<field::header> const & /**/, var_io::header const * & parsed_field)
    {
        parsed_field = &header;
    }

public:

    input_format_handler()                                              = default;
    input_format_handler(input_format_handler const &)                  = delete;
    input_format_handler(input_format_handler && )                      = default;
    input_format_handler & operator=(input_format_handler const &)      = delete;
    input_format_handler & operator=(input_format_handler &&)           = default;

    template <typename config_t>
    input_format_handler(std::istream & str, config_t const & /*cfg*/) :
        base_t{str}, file_it{str, false /*no_init!*/}

    {
        /* potentially extract useful runtime options from config and store as members */

        while (file_it != std::default_sentinel && file_it.peak() == '#')
        {
            ++file_it;
            std::string_view l = file_it->line;
            header.add_raw_line(l);
            ++line;
        }

        get<field::header>(raw_record)  = header.raw_header();
    }
};

} // namespace seqan
