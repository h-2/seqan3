// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Sara Hetzel <sara.hetzel AT fu-berlin.de>
 * \brief Provides seqan3::view::translate and seqan3::view::translate_single.
 */

#pragma once

#include <vector>
#include <stdexcept>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/alphabet/aminoacid/translation.hpp>
#include <seqan3/core/add_enum_bitwise_operators.hpp>
#include <seqan3/core/metafunction/range.hpp>
#include <seqan3/range/container/small_string.hpp>
#include <seqan3/range/detail/random_access_iterator.hpp>
#include <seqan3/range/view/deep.hpp>
#include <seqan3/range/view/detail.hpp>
#include <seqan3/std/concepts>
#include <seqan3/std/ranges>
#include <seqan3/range/container/concept.hpp>

// ============================================================================
//  forwards
// ============================================================================

namespace seqan3::detail
{

template <std::ranges::View urng_t>
//!\cond
    requires std::ranges::SizedRange<urng_t> &&
             std::ranges::RandomAccessRange<urng_t> &&
             NucleotideAlphabet<reference_t<urng_t>>
//!\endcond
class view_translate;

template <std::ranges::View urng_t>
//!\cond
    requires std::ranges::SizedRange<urng_t> &&
             std::ranges::RandomAccessRange<urng_t> &&
             NucleotideAlphabet<reference_t<urng_t>>
//!\endcond
class view_translate_single;

} // namespace seqan3::detail

// ============================================================================
//  translation_frames
// ============================================================================

namespace seqan3
{

//!\brief Specialisation values for single and multiple translation frames.
enum class translation_frames : uint8_t
{
    FWD_FRAME_0 = 1,                                    //!< The first forward frame starting at position 0
    FWD_FRAME_1 = 1 << 1,                               //!< The second forward frame starting at position 1
    FWD_FRAME_2 = 1 << 2,                               //!< The third forward frame starting at position 2
    REV_FRAME_0 = 1 << 3,                               //!< The first reverse frame starting at position 0
    REV_FRAME_1 = 1 << 4,                               //!< The second reverse frame starting at position 1
    REV_FRAME_2 = 1 << 5,                               //!< The third reverse frame starting at position 2
    FWD_REV_0 = FWD_FRAME_0 | REV_FRAME_0,              //!< The first forward and first reverse frame
    FWD_REV_1 = FWD_FRAME_1 | REV_FRAME_1,              //!< The second forward and second reverse frame
    FWD_REV_2 = FWD_FRAME_2 | REV_FRAME_2,              //!< The first third and third reverse frame
    FWD = FWD_FRAME_0 | FWD_FRAME_1 | FWD_FRAME_2,      //!< All forward frames
    REV = REV_FRAME_0 | REV_FRAME_1 | REV_FRAME_2,      //!< All reverse frames
    SIX_FRAME = FWD | REV                               //!< All frames
};

//!\brief Enable bitwise operators for enum translation_frames.
template<>
constexpr bool add_enum_bitwise_operators<translation_frames> = true;

}

namespace seqan3::detail
{

// ============================================================================
//  translate_fn (adaptor definition for both views)
// ============================================================================

/*!\brief Definition of the range adaptor object type for seqan3::view::translate and seqan3::view::translate_single.
 * \tparam single Switch between seqan3::view::translate and seqan3::view::translate_single.
 */
template <bool single>
struct translate_fn
{
    //!\brief The default frames parameter for the translation view adaptors.
    static constexpr translation_frames default_frames = single ?
                                                         translation_frames::FWD_FRAME_0 :
                                                         translation_frames::SIX_FRAME;

    //!\brief Store the argument and return a range adaptor closure object.
    constexpr auto operator()(translation_frames const tf = default_frames) const
    {
        return detail::adaptor_from_functor{*this, tf};
    }

    //!\brief Directly return an instance of the view, initialised with the given parameters.
    template <std::ranges::Range urng_t>
    constexpr auto operator()(urng_t && urange, translation_frames const tf = default_frames) const;
//     {
//         static_assert(std::ranges::ViewableRange<urng_t>,
//             "The range parameter to view::translate[_single] cannot be a temporary of a non-view range.");
//         static_assert(std::ranges::SizedRange<urng_t>,
//             "The range parameter to view::translate[_single] must model std::ranges::SizedRange.");
//         static_assert(std::ranges::RandomAccessRange<urng_t>,
//             "The range parameter to view::translate[_single] must model std::ranges::RandomAccessRange.");
//         static_assert(NucleotideAlphabet<reference_t<urng_t>>,
//             "The range parameter to view::translate[_single] must be over elements of seqan3::NucleotideAlphabet.");
//
//         if constexpr (single)
//         {
//             static_assert(std::ranges::View<decltype(detail::view_translate_single{std::forward<urng_t>(urange), tf})>);
//             return detail::view_translate_single{std::forward<urng_t>(urange), tf};
//         }
//         else
//             return detail::view_translate{std::forward<urng_t>(urange), tf};
//     }

    //!\brief This adaptor is usuable without setting the frames parameter in which case the default is chosen.
    template <std::ranges::Range urng_t>
    constexpr friend auto operator|(urng_t && urange, translate_fn const & me)
    {
        return me(std::forward<urng_t>(urange));
    }
};

// ============================================================================
//  view_translate_single (range definition)
// ============================================================================

/*!\brief The return type of seqan3::view::translate_single.
 * \implements std::ranges::View
 * \implements std::ranges::SizedRange
 * \implements std::ranges::RandomAccessRange
 * \ingroup view
 */
template <std::ranges::View urng_t>
//!\cond
    requires std::ranges::SizedRange<urng_t> &&
             std::ranges::RandomAccessRange<urng_t> &&
             NucleotideAlphabet<reference_t<urng_t>>
//!\endcond
class view_translate_single : public ranges::view_base
{
private:
    //!\brief The input range (of ranges).
    urng_t urange;
    //!\brief The frame that should be used for translation.
    translation_frames tf;
    //!\brief Error thrown if tried to be used with multiple frames.
    static inline small_string multiple_frame_error{"Error: Invalid type of frame. Choose one out of FWD_FRAME_0, "
                                                       "REV_FRAME_0, FWD_FRAME_1, REV_FRAME_1, FWD_FRAME_2 and "
                                                       "REV_FRAME_2."};
public:

    /*!\name Member types
     * \{
     */
    //!\brief The reference_type.
    using reference         = aa27;
    //!\brief The const_reference type.
    using const_reference   = aa27;
    //!\brief The value_type (which equals the reference_type with any references removed).
    using value_type        = aa27;
    //!\brief The size_type.
    using size_type         = size_type_t<urng_t>;
    //!\brief A signed integer type, usually std::ptrdiff_t.
    using difference_type   = difference_type_t<urng_t>;
    //!\brief The iterator type of this view (a random access iterator).
    using iterator          = detail::random_access_iterator<view_translate_single>;
    //!\brief The const_iterator type is equal to the iterator type.
    using const_iterator    = detail::random_access_iterator<view_translate_single const>;
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    view_translate_single()                                                        noexcept = default; //!< Defaulted.
    constexpr view_translate_single(view_translate_single const & rhs)             noexcept = default; //!< Defaulted.
    constexpr view_translate_single(view_translate_single && rhs)                  noexcept = default; //!< Defaulted.
    constexpr view_translate_single & operator=(view_translate_single const & rhs) noexcept = default; //!< Defaulted.
    constexpr view_translate_single & operator=(view_translate_single && rhs)      noexcept = default; //!< Defaulted.
    ~view_translate_single()                                                       noexcept = default; //!< Defaulted.


    /*!\brief Construct from another View.
     * \param[in] _urange The underlying range.
     * \param[in] _tf The frame that should be used for translation.
     *
     * ### Exceptions
     *
     * Throws if multiple frames are given as _tf input argument.
     */
    view_translate_single(urng_t _urange, translation_frames const _tf = translation_frames::FWD_FRAME_0)
        : urange{std::move(_urange)}, tf{_tf}
    {
        if (__builtin_popcount(static_cast<uint8_t>(_tf)) > 1)
        {
            throw std::invalid_argument(multiple_frame_error.c_str());
        }
    }

    /*!\brief Construct from another range.
     * \param[in] _urange The underlying range.
     * \param[in] _tf The frame that should be used for translation.
     *
     * ### Exceptions
     *
     * Throws if multiple frames are given as _tf input argument.
     */
    template<std::ranges::ViewableRange rng_t>
    //!\cond
        requires !std::Same<std::remove_reference_t<rng_t>, view_translate_single> &&
                 !std::Same<std::remove_reference_t<rng_t>, urng_t> &&
                 std::Constructible<urng_t, std::ranges::all_view<rng_t>>
    //!\endcond
    view_translate_single(rng_t && _urange, translation_frames const _tf = translation_frames::FWD_FRAME_0)
        : view_translate_single{std::view::all(std::forward<rng_t>(_urange)), _tf}
    {}
    //!\}

    /*!\name Iterators
     * \{
     */
    /*!\brief Returns an iterator to the first element of the container.
     * \returns Iterator to the first element.
     *
     * If the container is empty, the returned iterator will be equal to end().
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    iterator begin() noexcept
    {
        return {*this, 0};
    }

    //!\overload
    const_iterator begin() const noexcept
    {
        return {*this, 0};
    }

    //!\overload
    const_iterator cbegin() const noexcept
    {
        return begin();
    }

    /*!\brief Returns an iterator to the element following the last element of the container.
     * \returns Iterator to the first element.
     *
     * This element acts as a placeholder; attempting to dereference it results in undefined behaviour.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    iterator end() noexcept
    {
        return {*this, size()};
    }

    //!\overload
    const_iterator end() const noexcept
    {
        return {*this, size()};
    }

    //!\overload
    const_iterator cend() const noexcept
    {
        return end();
    }
    //!\}

    /*!\brief Returns the number of elements in the view.
     * \returns The number of elements in the container.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * Strong exception guarantee (never modifies data).
     */
    size_type size()
    {
        switch (tf)
        {
            case translation_frames::FWD_FRAME_0:
                [[fallthrough]];
            case translation_frames::REV_FRAME_0:
                return seqan3::size(urange) / 3;
                break;
            case translation_frames::FWD_FRAME_1:
                [[fallthrough]];
            case translation_frames::REV_FRAME_1:
                return (seqan3::size(urange) - 1) / 3;
                break;
            case translation_frames::FWD_FRAME_2:
                [[fallthrough]];
            case translation_frames::REV_FRAME_2:
                return (seqan3::size(urange) - 2) / 3;
                break;
            default:
                throw std::invalid_argument(multiple_frame_error.c_str());
                break;
        }
    }

    //!\overload
    size_type size() const
    {
        switch (tf)
        {
            case translation_frames::FWD_FRAME_0:
                [[fallthrough]];
            case translation_frames::REV_FRAME_0:
                return seqan3::size(urange) / 3;
                break;
            case translation_frames::FWD_FRAME_1:
                [[fallthrough]];
            case translation_frames::REV_FRAME_1:
                return (seqan3::size(urange) - 1) / 3;
                break;
            case translation_frames::FWD_FRAME_2:
                [[fallthrough]];
            case translation_frames::REV_FRAME_2:
                return (seqan3::size(urange) - 2) / 3;
                break;
            default:
                throw std::invalid_argument(multiple_frame_error.c_str());
                break;
        }
    }

    /*!\name Element access
     * \{
     */
    /*!\brief Return the n-th element.
     * \param[in] n The element to retrieve.
     *
     * Accessing an element behind the last causes undefined behaviour. In debug mode an assertion checks the size of
     * the container.
     *
     * ### Exceptions
     *
     * Strong exception guarantee (never modifies data).
     *
     * ### Complexity
     *
     * Constant.
     */
    reference operator[](size_type const n)
    {
     assert(n < size());
     switch (tf)
     {
         case translation_frames::FWD_FRAME_0:
             return translate_triplet((urange)[n * 3], (urange)[n * 3 + 1], (urange)[n * 3 + 2]);
             break;
         case translation_frames::REV_FRAME_0:
             return translate_triplet(complement((urange)[(urange).size() - n * 3 - 1]), complement((urange)[(urange).size() - n * 3 - 2]), complement((urange)[(urange).size() - n * 3 - 3]));
             break;
         case translation_frames::FWD_FRAME_1:
             return translate_triplet((urange)[n * 3 + 1], (urange)[n * 3 + 2], (urange)[n * 3 + 3]);
             break;
         case translation_frames::REV_FRAME_1:
             return translate_triplet(complement((urange)[(urange).size() - n * 3 - 2]), complement((urange)[(urange).size() - n * 3 - 3]), complement((urange)[(urange).size() - n * 3 - 4]));
             break;
         case translation_frames::FWD_FRAME_2:
             return translate_triplet((urange)[n * 3 + 2], (urange)[n * 3 + 3], (urange)[n * 3 + 4]);
             break;
         case translation_frames::REV_FRAME_2:
             return translate_triplet(complement((urange)[(urange).size() - n * 3 - 3]), complement((urange)[(urange).size() - n * 3 - 4]), complement((urange)[(urange).size() - n * 3 - 5]));
             break;
         default:
             throw std::invalid_argument(multiple_frame_error.c_str());
             break;
     }
    }

    //!\overload
    const_reference operator[](size_type const n) const
    {
        assert(n < size());
        switch (tf)
        {
            case translation_frames::FWD_FRAME_0:
                return translate_triplet((urange)[n * 3], (urange)[n * 3 + 1], (urange)[n * 3 + 2]);
                break;
            case translation_frames::REV_FRAME_0:
                return translate_triplet(complement((urange)[(urange).size() - n * 3 - 1]), complement((urange)[(urange).size() - n * 3 - 2]), complement((urange)[(urange).size() - n * 3 - 3]));
                break;
            case translation_frames::FWD_FRAME_1:
                return translate_triplet((urange)[n * 3 + 1], (urange)[n * 3 + 2], (urange)[n * 3 + 3]);
                break;
            case translation_frames::REV_FRAME_1:
                return translate_triplet(complement((urange)[(urange).size() - n * 3 - 2]), complement((urange)[(urange).size() - n * 3 - 3]), complement((urange)[(urange).size() - n * 3 - 4]));
                break;
            case translation_frames::FWD_FRAME_2:
                return translate_triplet((urange)[n * 3 + 2], (urange)[n * 3 + 3], (urange)[n * 3 + 4]);
                break;
            case translation_frames::REV_FRAME_2:
                return translate_triplet(complement((urange)[(urange).size() - n * 3 - 3]), complement((urange)[(urange).size() - n * 3 - 4]), complement((urange)[(urange).size() - n * 3 - 5]));
                break;
            default:
                throw std::invalid_argument(multiple_frame_error.c_str());
                break;
        }
    }
    //!\}

    //!\brief Implicit conversion to container types.
    template <RandomAccessContainer container_type>
    operator container_type()
    //!\cond
        requires std::is_same_v<aa27, value_type_t<container_type>>
    //!\endcond
    {
        container_type ret;
        ret.resize(size());
        std::copy(cbegin(), cend(), ret.begin());
        return ret;
    }

    //!\overload
    template <RandomAccessContainer container_type>
    operator container_type() const
    //!\cond
        requires std::is_same_v<aa27, value_type_t<container_type>>
    //!\endcond
    {
        container_type ret;
        ret.resize(size());
        std::copy(cbegin(), cend(), ret.begin());
        return ret;
    }
};

//!\brief Class template argument deduction for view_translate_single.
template <typename urng_t>
view_translate_single(urng_t &&, translation_frames const) -> view_translate_single<std::ranges::all_view<urng_t>>;


//!\brief Class template argument deduction for view_translate_single with default translation_frames.
template <typename urng_t>
view_translate_single(urng_t &&) -> view_translate_single<std::ranges::all_view<urng_t>>;

} // namespace seqan3::detail


namespace seqan3::detail
{

//!\brief Type trait that declares any seqan3::gap_decorator_anchor_set to be **NOT a view**.
template <typename type>
constexpr int enable_view<view_translate_single<type>> = 1;

template <typename type>
constexpr int enable_view<view_translate_single<type> const> = 1;

}

// ============================================================================
//  translate_single (adaptor object)
// ============================================================================

namespace seqan3::view
{

/*!\name Alphabet related views
 * \{
 */

/*!\brief A view that translates nucleotide into aminoacid alphabet for one of the six frames.
 * \tparam urng_t The type of the range being processed.
 * \param[in] urange The range being processed.
 * \param[in] tf A value of seqan3::translation_frames that indicates the desired frames.
 * \returns A range containing frames with aminoacid sequence. See below for the properties of the returned range.
 * \ingroup view
 *
 * \details
 *
 * This view can be used to translate nucleotide sequences into aminoacid sequences (see translation_frames for possible combination of frames).
 *
 * **Header**
 * ```cpp
 *      #include <seqan3/range/view/translate.hpp>
 * ```
 *
 * ### View properties
 *
 * | range concepts and reference_t  | `urng_t` (underlying range type)      | `rrng_t` (returned range type)                     |
 * |---------------------------------|:-------------------------------------:|:--------------------------------------------------:|
 * | std::ranges::InputRange         | *required*                            | *preserved*                                        |
 * | std::ranges::ForwardRange       | *required*                            | *preserved*                                        |
 * | std::ranges::BidirectionalRange | *required*                            | *preserved*                                        |
 * | std::ranges::RandomAccessRange  | *required*                            | *preserved*                                        |
 * | std::ranges::ContiguousRange    |                                       | *lost*                                             |
 * |                                 |                                       |                                                    |
 * | std::ranges::ViewableRange      | *required*                            | *guaranteed*                                       |
 * | std::ranges::View               |                                       | *guaranteed*                                       |
 * | std::ranges::SizedRange         | *required*                            | *preserved*                                        |
 * | std::ranges::CommonRange        |                                       | *guaranteed*                                       |
 * | std::ranges::OutputRange        |                                       | *lost*                                             |
 * | seqan3::ConstIterableRange      | *required*                            | *preserved*                                        |
 * |                                 |                                       |                                                    |
 * | seqan3::reference_t             | seqan3::NucleotideAlphabet            | seqan3::aa27                                       |
 *
 * * `urng_t` is the type of the range modified by this view (input).
 * * `rrng_type` is the type of the range returned by this view.
 * * for more details, see \ref view.
 *
 * ### Example
 *
 * Operating on a range of seqan3::dna5:
 * \snippet test/snippet/range/view/translate.cpp dna5
 * \hideinitializer
 */
inline constexpr auto translate_single = deep{detail::translate_fn<true>{}};

} // seqan3::view

// ============================================================================
//  view_translate (range definition)
// ============================================================================

namespace seqan3::detail
{

/*!\brief The return type of seqan3::view::translate.
 * \implements std::ranges::View
 * \implements std::ranges::SizedRange
 * \implements std::ranges::RandomAccessRange
 * \tparam urng_t The type of the range being translated.
 * \param[in] tf Translation frames to be used.
 * \ingroup view
 */
template <std::ranges::View urng_t>
//!\cond
    requires std::ranges::SizedRange<urng_t> &&
             std::ranges::RandomAccessRange<urng_t> &&
             NucleotideAlphabet<reference_t<urng_t>>
//!\endcond
class view_translate : public ranges::view_base
{
private:
    //!\brief The data members of view_translate_single.
    urng_t urange;
    //!\brief The frames that should be used for translation.
    translation_frames tf;
    //!\brief The selected frames corresponding to the frames required.
    std::vector<translation_frames> selected_frames{};

public:
    /*!\name Member types
     * \{
     */
    //!\brief The reference_type.
    using reference         = view_translate_single<urng_t>;
    //!\brief The const_reference type.
    using const_reference   = reference;
    //!\brief The value_type (which equals the reference_type with any references removed).
    using value_type        = reference;
    //!\brief The size_type.
    using size_type         = size_type_t<urng_t>;
    //!\brief A signed integer type, usually std::ptrdiff_t.
    using difference_type   = difference_type_t<urng_t>;
    //!\brief The iterator type of this view (a random access iterator).
    using iterator          = detail::random_access_iterator<view_translate>;
    //!\brief The const iterator type of this view (same as iterator, because it's a view).
    using const_iterator    = detail::random_access_iterator<view_translate const>;
    //!\}

protected:
    /*!\name Compatibility
     * \brief Static constexpr variables that emulate/encapsulate seqan3::Compatible (which doesn't work for types during their definition).
     * \{
     */
    //!\cond
    // unfortunately we cannot specialise the variable template so we have to add an auxiliary here
    template <typename t>
        requires (dimension_v<t> == dimension_v<value_type> + 1) &&
                 std::is_same_v<remove_cvref_t<innermost_value_type_t<value_type>>,
                                remove_cvref_t<innermost_value_type_t<t>>>
    static constexpr bool is_compatible_this_aux = true;
    //!\endcond
    //!\}

public:

    /*!\name Constructors, destructor and assignment
     * \{
     */
    view_translate()                                                 noexcept = default; //!< Defaulted.
    constexpr view_translate(view_translate const & rhs)             noexcept = default; //!< Defaulted.
    constexpr view_translate(view_translate && rhs)                  noexcept = default; //!< Defaulted.
    constexpr view_translate & operator=(view_translate const & rhs) noexcept = default; //!< Defaulted.
    constexpr view_translate & operator=(view_translate && rhs)      noexcept = default; //!< Defaulted.
    ~view_translate()                                                noexcept = default; //!< Defaulted.

    view_translate(urng_t _urange, translation_frames const _tf = translation_frames::SIX_FRAME)
        : urange{std::move(_urange)}, tf{_tf}
    {
        if ((_tf & translation_frames::FWD_FRAME_0) == translation_frames::FWD_FRAME_0)
            selected_frames.push_back(translation_frames::FWD_FRAME_0);
        if ((_tf & translation_frames::FWD_FRAME_1) == translation_frames::FWD_FRAME_1)
            selected_frames.push_back(translation_frames::FWD_FRAME_1);
        if ((_tf & translation_frames::FWD_FRAME_2) == translation_frames::FWD_FRAME_2)
            selected_frames.push_back(translation_frames::FWD_FRAME_2);
        if ((_tf & translation_frames::REV_FRAME_0) == translation_frames::REV_FRAME_0)
            selected_frames.push_back(translation_frames::REV_FRAME_0);
        if ((_tf & translation_frames::REV_FRAME_1) == translation_frames::REV_FRAME_1)
            selected_frames.push_back(translation_frames::REV_FRAME_1);
        if ((_tf & translation_frames::REV_FRAME_2) == translation_frames::REV_FRAME_2)
            selected_frames.push_back(translation_frames::REV_FRAME_2);
    }

    /*!\brief Construct from another range.
     * \param[in] _urange The underlying range (of ranges).
     * \param[in] _tf The frames that should be used for translation.
     */
    template<std::ranges::ViewableRange rng_t>
    //!\cond
        requires !std::Same<std::remove_reference_t<rng_t>, view_translate> &&
                 !std::Same<std::remove_reference_t<rng_t>, urng_t> &&
                 std::Constructible<urng_t, std::ranges::all_view<rng_t>>
    //!\endcond
    view_translate(rng_t && _urange, translation_frames const _tf = translation_frames::SIX_FRAME)
        : view_translate{std::view::all(std::forward<rng_t>(_urange)), _tf}
    {}
    //!\}

    /*!\name Iterators
     * \{
     */
    /*!\brief Returns an iterator to the first element of the container.
     * \returns Iterator to the first element.
     *
     * If the container is empty, the returned iterator will be equal to end().
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    iterator begin() noexcept
    {
        return {*this, 0};
    }

    //!\overload
    const_iterator begin() const noexcept
    {
        return {*this, 0};
    }

    //!\overload
    const_iterator cbegin() const noexcept
    {
        return begin();
    }

    /*!\brief Returns an iterator to the element following the last element of the container.
     * \returns Iterator to the first element.
     *
     * This element acts as a placeholder; attempting to dereference it results in undefined behaviour.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    iterator end() noexcept
    {
        return {*this, size()};
    }

    //!\overload
    const_iterator end() const noexcept
    {
        return {*this, size()};
    }

    //!\overload
    const_iterator cend() const noexcept
    {
        return end();
    }
    //!\}

    /*!\brief Returns the number of elements in the view.
     * \returns The number of elements in the container.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    size_type size() noexcept
    {
        return (size_type) selected_frames.size();
    }

    //!\overload
    size_type size() const noexcept
    {
        return (size_type) selected_frames.size();
    }

    /*!\name Element access
     * \{
     */
    /*!\brief Return the n-th element.
     * \param[in] n The element to retrieve.
     *
     * Accessing an element behind the last causes undefined behaviour. In debug mode an assertion checks the size of
     * the container.
     *
     * ### Exceptions
     *
     * Strong exception guarantee (never modifies data).
     *
     * ### Complexity
     *
     * Constant.
     */
    reference operator[](size_type const n)
    {
        assert(n < size());
        return urange | view::translate_single(selected_frames[n]);
    }

    //!\overload
    const_reference operator[](size_type const n) const
    {
        assert(n < size());
        return urange | view::translate_single(selected_frames[n]);
    }
    //!\}

    //!\brief Implicit conversion to container types.
    template <RandomAccessContainer container_type>
    operator container_type()
    //!\cond
        requires is_compatible_this_aux<container_type>
    //!\endcond
    {
        container_type ret;
        ret.resize(size());
        for (size_type i = 0; i < size(); i++)
            ret[i] = static_cast<value_type_t<container_type>>(operator[](i));
        return ret;
    }

    //!\overload
    template <RandomAccessContainer container_type>
    operator container_type() const
    //!\cond
        requires is_compatible_this_aux<container_type>
    //!\endcond
    {
        container_type ret;
        ret.resize(size());
        for (size_type i = 0; i < size(); i++)
            ret[i] = static_cast<value_type_t<container_type>>(operator[](i));
        return ret;
    }
};

//!\brief Class template argument deduction for view_translate.
template <typename urng_t>
//!\cond
    requires std::ranges::SizedRange<urng_t> &&
             std::ranges::RandomAccessRange<urng_t> &&
             NucleotideAlphabet<reference_t<urng_t>>
//!\endcond
view_translate(urng_t &&, translation_frames const = translation_frames{}) -> view_translate<std::ranges::all_view<urng_t>>;

} // namespace seqan3::detail

// ============================================================================
//  translate (adaptor object)
// ============================================================================

namespace seqan3::view
{

/*!\name Alphabet related views
 * \{
 */

/*!\brief A view that translates nucleotide into aminoacid alphabet with 1, 2, 3 or 6 frames.
 * \tparam urng_t The type of the range being processed.
 * \param[in] urange The range being processed.
 * \param[in] tf A value of seqan3::tanslation_frames that indicates the desired frames.
 * \returns A range of ranges containing frames with aminoacid sequence. See below for the properties of the returned range.
 * \ingroup view
 *
 * \details
 *
 * This view can be used to translate nucleotide sequences into aminoacid sequences (see translation_frames for possible combination of frames).
 *
 * **Header**
 * ```cpp
 *      #include <seqan3/range/view/translate.hpp>
 * ```
 *
 * ### View properties
 *
 * | range concepts and reference_t  | `urng_t` (underlying range type)      | `rrng_t` (returned range type)                     |
 * |---------------------------------|:-------------------------------------:|:--------------------------------------------------:|
 * | std::ranges::InputRange         | *required*                            | *preserved*                                        |
 * | std::ranges::ForwardRange       | *required*                            | *preserved*                                        |
 * | std::ranges::BidirectionalRange | *required*                            | *preserved*                                        |
 * | std::ranges::RandomAccessRange  | *required*                            | *preserved*                                        |
 * | std::ranges::ContiguousRange    |                                       | *lost*                                             |
 * |                                 |                                       |                                                    |
 * | std::ranges::ViewableRange      | *required*                            | *guaranteed*                                       |
 * | std::ranges::View               |                                       | *guaranteed*                                       |
 * | std::ranges::SizedRange         | *required*                            | *preserved*                                        |
 * | std::ranges::CommonRange        |                                       | *guaranteed*                                       |
 * | std::ranges::OutputRange        |                                       | *lost*                                             |
 * | seqan3::ConstIterableRange      | *required*                            | *preserved*                                        |
 * |                                 |                                       |                                                    |
 * | seqan3::reference_t             | seqan3::NucleotideAlphabet            | std::ranges::View && std::ranges::RandomAccessRange && std::ranges::SizedRange |
 *
 * * `urng_t` is the type of the range modified by this view (input).
 * * `rrng_type` is the type of the range returned by this view.
 * * for more details, see \ref view.
 *
 * ### Example
 *
 * Operating on a range of seqan3::dna5:
 * \snippet test/snippet/range/view/translate.cpp usage
 * \hideinitializer
 */
inline constexpr auto translate = deep{detail::translate_fn<false>{}};
//!\}

} // namespace seqan3::view

namespace seqan3::detail
{
template <bool single>
template <std::ranges::Range urng_t>
constexpr auto translate_fn<single>::operator()(urng_t && urange, translation_frames const tf) const
{
    static_assert(std::ranges::ViewableRange<urng_t>,
        "The range parameter to view::translate[_single] cannot be a temporary of a non-view range.");
    static_assert(std::ranges::SizedRange<urng_t>,
        "The range parameter to view::translate[_single] must model std::ranges::SizedRange.");
    static_assert(std::ranges::RandomAccessRange<urng_t>,
        "The range parameter to view::translate[_single] must model std::ranges::RandomAccessRange.");
    static_assert(NucleotideAlphabet<reference_t<urng_t>>,
        "The range parameter to view::translate[_single] must be over elements of seqan3::NucleotideAlphabet.");

    if constexpr (single)
    {
        static_assert(std::ranges::View<decltype(detail::view_translate_single{std::forward<urng_t>(urange), tf})>);
//         static_assert(std::ranges::ViewableRange<decltype(detail::view_translate_single{std::forward<urng_t>(urange), tf})>);
// //         static_assert(std::ranges::ViewableRange<detail::view_translate_single<std::ranges::all_view<urng_t>>>);
        return detail::view_translate_single{std::forward<urng_t>(urange), tf};
    }
    else
        return detail::view_translate{std::forward<urng_t>(urange), tf};
}

} // namespace seqan3::detail
