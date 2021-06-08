// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// Copyright (c) 2020-2021, deCODE Genetics
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides utility types for seqan3::plaintext_file_input and seqan3::plaintext_file_output.
 * \author Hannes Hauswedell <hannes.hauswedell AT decode.is>
 */

#pragma once

#include <stdexcept>
#include <string_view>
#include <vector>

#include <seqan3/core/platform.hpp>

namespace seqan3::plain_io
{

//!\brief The value type of seqan3::plaintext_file_input if every line is split into fields.
struct record
{
    //!\brief The entire line (exluding EOL characters but including delimiters).
    std::string_view line;
    //!\brief A range of the individual fields (without delimiters or EOL characters).
    std::vector<std::string_view> fields;
};

//!\brief Option to switch between reading-by-line and splitting a line into fields.
enum class record_kind
{
    line,               //!< Only the line is provided.
    line_and_fields     //!< The line is provided and also individual fields (seqan3::plaintext_record).
};


/*!\brief A helper for specifying the header of a seqan3::plaintext_file_input.
 * \tparam record_kind Whether to split lines on delimiter (e.g. TSV files) or not.
 */
class header_kind
{
private:
    //!\brief The state of the variable, encoded as the value of the char or the special values -200 and -300.
    int state = -200;

public:
    //!\brief The state representing "no header".
    static constexpr struct {} none{};
    //!\brief The state representing "first line is header".
    static constexpr struct {} first_line{};

    //!\brief Type of the state representing "all lines that start with character X".
    struct starts_with
    {
        //!\privatesection
        char c;
    };

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr header_kind()                                           noexcept = default;
    constexpr header_kind(header_kind const & )             noexcept = default;
    constexpr header_kind(header_kind && )                  noexcept = default;
    constexpr header_kind & operator=(header_kind const & ) noexcept = default;
    constexpr header_kind & operator=(header_kind && )      noexcept = default;

    constexpr header_kind(decltype(none))         noexcept : state{-200} {}
    constexpr header_kind(decltype(first_line))   noexcept : state{-300} {}
    constexpr header_kind(starts_with s)          noexcept : state{s.c}  {}
    //!\}

    /*!\name Functions for retrieving the state.
     * \{
     */
    constexpr bool is_none()        const noexcept { return state == -200; }
    constexpr bool is_first_line()  const noexcept { return state == -300; }
    constexpr bool is_starts_with() const noexcept { return state != -200 && state != -300; }

    char get_starts_with() const
    {
        if (!is_starts_with())
        {
            throw std::logic_error{"Tried to read starts_with from header_kind but it was in a "
                                   "different state."};
        }

        return static_cast<char>(state);
    }
    //!\}
};

} // namespace seqan3
