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
#include <seqan3/io/detail/ignore_output_iterator.hpp>
#include <seqan3/io/detail/misc.hpp>
#include <seqan3/io/format/input_format_handler_base.hpp>
#include <seqan3/io/format/format_vcf.hpp>
#include <seqan3/io/plaintext_io/reader.hpp>
#include <seqan3/io/stream/iterator.hpp>
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

    /* pos is handled correctly by default */

    //!\brief Read ID as string, but read "." as empty string.
    template <typename parsed_field_t>
    void parse_field(tag_t<field::id> const & /**/, parsed_field_t & parsed_field)
    {
        if (get<field::id>(raw_record) != ".") // dot means empty, otherwise delegate to base
            static_cast<base_t &>(*this).parse_field(tag<field::id>, parsed_field);
    }

    //!\brief Overload for parsing REF.
    void parse_field(tag_t<field::ref> const & /**/, var_io::allele & parsed_field)
    {
        std::string_view raw_field = get<field::ref>(raw_record);
        parsed_field = parse_allele_impl(raw_field);
    }

    //!\brief Overload for parsing ALT.
    // TODO: make this back_insertable instead vector
    void parse_field(tag_t<field::alt> const & /**/, std::vector<var_io::allele> & parsed_field)
    {
        std::string_view raw_field = get<field::alt>(raw_record);

        for (std::string_view subfield : raw_field | views::eager_split(','))
            parsed_field.push_back(parse_allele_impl(subfield));
    }

    //!\brief Overload for parsing QUAL.
    template <arithmetic parsed_field_t>
    void parse_field(tag_t<field::qual> const & /**/,
                     std::variant<var_io::special_value, parsed_field_t> & parsed_field)
    {
        std::string_view raw_field = get<field::qual>(raw_record);

        if (raw_field == ".")
        {
            parsed_field = var_io::special_value::missing;
        }
        else if (raw_field == "*")
        {
            parsed_field = var_io::special_value::unknown;
        }
        else
        {
            parsed_field = parsed_field_t{};
            static_cast<base_t &>(*this).parse_field(tag<field::qual>, std::get<1>(parsed_field)); // arithmetic parsing
        }
    }

    //!\brief Overload for parsing FILTER as strings.
    void parse_field(tag_t<field::filter> const & /**/, std::vector<std::string> & parsed_field)
    {
        std::string_view raw_field = get<field::filter>(raw_record);

        if (raw_field == ".")
            return;

        for (std::string_view subfield : raw_field | views::eager_split(';'))
            parsed_field.push_back(std::string{subfield});
    }

    //!\brief Overload for parsing FILTER as indexes.
    template <std::integral parsed_field_t>
    void parse_field(tag_t<field::filter> const & /**/,  std::vector<parsed_field_t> & parsed_field)
    {
        std::string_view raw_field = get<field::filter>(raw_record);

        if (raw_field == ".")
            return;

        for (std::string_view subfield : raw_field | views::eager_split(';'))
        {
            if (auto it = header.parsed_header().filter_id_to_index.find(raw_field);
                     it == header.parsed_header().filter_id_to_index.end()) // filter name was not in header, insert!
            {
                header.add_filter(raw_field);
                parsed_field.push_back(static_cast<parsed_field_t>(header.parsed_header().filters.size()));

                if (warn_on_missing_header_entries)
                {
                    debug_stream << "[seqan3::var_io::warning] The filter name \"" << raw_field << "\" found on line "
                                 << line << "was not present in the header.\n";
                }

            }
            else
            {
                parsed_field.push_back(static_cast<parsed_field_t>(it->second));
            }
        }
    }

    //!\brief Overload for parsing INFO.
    template <typename key_t, typename value_t>
    void parse_field(tag_t<field::info> const & /**/, std::vector<std::pair<key_t, value_t>> & parsed_field)
    {
        static_assert((std::same_as<key_t, int32_t>     && std::same_as<value_t, io_type_variant>) ||
                      (std::same_as<key_t, std::string> && std::same_as<value_t, std::string>),
                      "Unsupported key/value types for reading INFO field. See documentation of var_io::reader.");

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
            else // std::string
            {
                ret.first = std::string{key};
            }

            if (val_it == key_value_view.end()) // no "=" → flag
            {
                if constexpr (std::same_as<value_t, io_type_variant>)
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
            else
            {
                if constexpr (std::same_as<value_t, io_type_variant>)
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
                else // val_t == std::string
                {
                    ret.second = std::string{*val_it};
                }
            }

            parsed_field.push_back(std::move(ret));
        }
    }

    void parse_field(tag_t<field::genotypes> const & /**/, std::vector<genotype_element> & parsed_field)
    {
        size_t column_number          = file_it->fields.size();
        size_t expected_column_number = header.parsed_header().samples.count() + 8;

        if (column_number > 8) // there are genotypes
        {
            if (column_number != expected_column_number)
                error("Expected ", expected_column_number, " columns in line but found", column_number, ".");

            std::string_view format_names = file_it->fields[8];

            for (std::string_view format_name : format_names | views::eager_split(':'))
            {
                size_t format_index = -1;
                if (auto it = header.parsed_header().format_id_to_index.find(format_name);
                    it == header.parsed_header().format_id_to_index.end()) // format name was not in header, insert!
                {
                    //TODO
                }
                else
                {
                    format_index = it->second;
                }

                header::format_t format const & header.parsed_header().formats[format_index];

                auto resize = [s = header.parsed_header().samples.count()] (auto & vec) { vec.resize(s); };


        }
    }

    void parse_field(tag_t<field::header> const & /**/, var_io::header const * & parsed_field)
    {
        parsed_field = &header;
    }


    /* TODO move to mixin: */

    static var_io::allele parse_allele_impl(std::string_view in)
    {
        std::regex is_seq{"^[ACGTNacgtn]+$"};

        if (in == "*")
            return var_io::special_value::unknown;
        else if (in == ".")
            return var_io::special_value::missing;
        else if (std::regex_match(in.begin(), in.end(), is_seq))
            return in | views::char_to<dna5> | views::to<std::vector<dna5>>;
        else
            return std::string{in};
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

        bool file_format_read = false;
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
