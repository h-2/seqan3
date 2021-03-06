// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
 * \brief Provides seqan3::NucleotideAlphabet.
 */

#pragma once

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/std/concepts>

// ============================================================================
// forwards
// ============================================================================

//!\cond
namespace seqan3::custom
{

void complement();

} // namespace seqan3::custom
//!\endcond

// ============================================================================
// complement()
// ============================================================================

namespace seqan3::detail::adl::only
{

//!\brief Functor definition for seqan3::complement.
struct complement_fn
{
private:
    SEQAN3_CPO_IMPL(2, complement(v)                     ) // ADL
    SEQAN3_CPO_IMPL(1, seqan3::custom::complement(v)     ) // customisation namespace
    SEQAN3_CPO_IMPL(0, v.complement()                    ) // member

public:
    //!\brief Operator definition.
    template <typename nucleotide_t>
    //!\cond
        requires requires (nucleotide_t const nucl) { { impl(priority_tag<2>{}, nucl) }; }
    //!\endcond
    constexpr auto operator()(nucleotide_t const nucl) const noexcept
    {
        static_assert(noexcept(impl(priority_tag<2>{}, nucl)),
            "Only overloads that are marked noexcept are picked up by seqan3::complement().");
        static_assert(std::Same<nucleotide_t, decltype(impl(priority_tag<2>{}, nucl))>,
            "The return type of your complement() implementation must be 'nucleotide_t'.");

        return impl(priority_tag<2>{}, nucl);
    }
};

} // namespace seqan3::detail::adl::only

namespace seqan3
{

/*!\name Function objects (Nucleotide)
 * \{
 */

/*!\brief Return the complement of a nucleotide object.
 * \tparam your_type Type of the argument.
 * \param  nucl      The nucleotide object for which you want to receive the complement.
 * \returns The complement character of `nucl`, e.g. 'C' for 'G'.
 * \ingroup nucleotide
 * \details
 *
 * This is a function object. Invoke it with the parameter(s) specified above.
 *
 * It acts as a wrapper and looks for three possible implementations (in this order):
 *
 *   1. A free function `complement(your_type const a)` in the namespace of your type (or as `friend`).
 *      The function must be marked `noexcept` (`constexpr` is not required, but recommended) and the
 *      return type be `your_type`.
 *   2. A free function `complement(your_type const a)` in `namespace seqan3::custom`.
 *      The same restrictions apply as above.
 *   3. A member function called `complement()`.
 *      It must be marked `noexcept` (`constexpr` is not required, but recommended) and the return type be
 *      `your_type`.
 *
 * Every nucleotide alphabet type must provide one of the above.
 *
 * ### Example
 *
 * \include test/snippet/alphabet/nucleotide/complement_fn.cpp
 *
 * ### Customisation point
 *
 * This is a customisation point (see \ref about_customisation). To specify the behaviour for your own alphabet type,
 * simply provide one of the three functions specified above.
 */
inline constexpr auto complement = detail::adl::only::complement_fn{};
//!\}

// ============================================================================
// NucleotideAlphabet concept
// ============================================================================

/*!\interface seqan3::NucleotideAlphabet <>
 * \extends seqan3::Alphabet
 * \brief A concept that indicates whether an alphabet represents nucleotides.
 * \ingroup nucleotide
 *
 * \details
 *
 * In addition to the requirements for seqan3::Alphabet, the NucleotideAlphabet introduces
 * a requirement for a complement function: seqan3::NucleotideAlphabet::complement.
 *
 * ### Requirements
 *
 *   1. `t` shall model seqan3::Alphabet
 *   2. seqan3::complement needs to be defined for objects of type `t`
 *
 * See the documentation pages for the respective requirements.
 *
 * ### Related types
 *
 * If a given type `t` models this concept, the following types typically do so, as well:
 *
 *   * `t &`
 *   * `t const`
 *   * `t const &`
 */
//!\cond
template <typename t>
SEQAN3_CONCEPT NucleotideAlphabet = Alphabet<t> && requires (t val)
{
    { seqan3::complement(val) };
};
//!\endcond

} // namespace seqan3
