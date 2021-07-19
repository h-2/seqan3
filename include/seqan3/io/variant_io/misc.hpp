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
#include <tuple>
#include <variant>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/core/detail/debug_stream_type.hpp>
#include <seqan3/io/record.hpp>
#include <seqan3/io/utility.hpp>

namespace seqan3::var_io
{

template <typename t>
inline t missing_value = t{};


template <>
inline std::string missing_value<std::string> = "*";

template <>
inline constexpr std::string_view missing_value<std::string_view> = "*";


template <typename int_t>
    requires (std::same_as<int_t, int8_t> || std::same_as<int_t, int16_t> || std::same_as<int_t, int32_t>)
inline constexpr int_t missing_value<int_t> = std::numeric_limits<int_t>::lowest();

template <>
inline float missing_value<float> = [] ()
{
    uint32_t i = 0x7F800001U;
    return *reinterpret_cast<float*>(&i);
} ();


/*!\brief A type representing an variant file INFO field [index of the INFO in header, value].
 * \ingroup variant_io
 */
template <bool shallow = true>
using info_element = std::pair<int32_t, io_type_variant<shallow>>;

/*!\brief A type representing an genotype.
 * \ingroup variant_io
 *
 * \details
 *
 * Genotypes / samples are represented as decribed in the BCF specification, i.e. information is grouped by FORMAT
 * identifier, not by sample.
 *
 * This entry consists of the FORMAT index in the file's header and a vector of values. The size of the vector is:
 *
 *   * equal to the number of samples; or
 *   * 0 -- if the field is missing from all samples.
 *
 * The variant vector is guaranteed to be over the type defined in the header. Note that this is a vector over such
 * types (one element per sample!), so seqan3::io_type_id::vector_of_int32 corresponds to std::vector<std::vector<int32_t>>.
 * See seqan3::io_type_vector_variant for more details.
 *
 * If fields are missing from some samples but not others, the vector will have full size but the respective values
 * will be set to the missing value (see seqan3::var_io::is_missing()) or be the empty vector (in case the element type
 * is a vector).
 */
template <bool shallow = true>
using genotype_element = std::pair<int32_t, io_type_vector_variant<shallow>>;

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
                                              field::_private>;
} // namespace seqan3::var_io

namespace seqan3::detail
{

struct parse_io_type_data_t
{
    std::string_view input;
    static constexpr std::string_view missing = ".";

    constexpr size_t operator()(bool & output) const
    {
        output = true;
        return 0;
    }

    inline size_t operator()(std::vector<bool>::reference output) const
    {
        output = true;
        return 0;
    }

    constexpr size_t operator()(char & output) const
    {
        assert(input.size() == 1);
        output = input[0];
        return 1;
    }

    template <arithmetic arith_t>
    inline size_t operator()(arith_t & output) const
    {
        if (input == missing)
            output = var_io::missing_value<arith_t>;
        else
            string_to_number(input, output);
        return 1;
    }

    inline size_t operator()(std::string & output) const
    {
        if (input != missing)
            output = std::string{input};
        return 1;
    }

    inline size_t operator()(std::string_view & output) const
    {
        output = input;
        return 1;
    }

    template <typename elem_t>
    //TODO back_insertable
    inline size_t operator()(std::vector<elem_t> & vec) const
    {
        if (input != missing)
        {
            for (std::string_view s : input | views::eager_split(','))
            {
                vec.emplace_back();
                parse_io_type_data_t{s}(vec.back());
            }
        }

        return vec.size();
    }
};

/*!\brief Create an seqan3::io_type_variant from a string and a known seqan3::io_type_id.
 * \ingroup io
 * \param[in]  id           ID of the type that shall be read.
 * \param[in]  input_string The string data to read from.
 * \param[out] output       The object to store the result in.
 * \returns The number of elements stored in the output in case ID is one of the "vector_of_"-types; 1 otherwise.
 */
inline size_t parse_io_type_variant(io_type_id const id,
                                    std::string_view const input_string,
                                    is_io_type_variant auto & output)
{
    init_io_type_variant(id, output);
    return std::visit(parse_io_type_data_t{input_string}, output);
}

} // namespace seqan3::detail
