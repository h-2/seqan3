// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 * \brief Provides seqan3::eager_split_view
 */

#pragma once

#include <seqan3/std/concepts>
#include <seqan3/std/iterator>
#include <seqan3/std/ranges>

#include <seqan3/range/views/detail.hpp>

//-----------------------------------------------------------------------------
// Implementation of single pass input view.
//-----------------------------------------------------------------------------

namespace seqan3::detail
{
/*!\brief An input_iterator over the associated range.
 * \implements std::input_iterator
 * \ingroup views
 * \tparam view_type The type of the associated type.
 *
 * This iterator reduces every iterator type of the associated view to an single pass input iterator.
 */

class eager_split_iterator
{
private:
    std::string_view urange;
    char             delimiter{};
    std::string_view subrange;
    char const *     uend = nullptr;
    char const *     subend = nullptr;
    bool             skip_quotes = false;
    bool             in_quote = false;
public:

    /*!\name Associated types
     * \{
     */
    //!\brief Difference type.
    using difference_type       = std::ptrdiff_t;
    //!\brief Value type.
    using value_type            = std::string_view;
    //!\brief Pointer type.
    using pointer               = std::string_view const *;
    //!\brief Reference type.
    using reference             = std::string_view;
    //!\brief Iterator category.
    using iterator_category     = std::input_iterator_tag;
    //!\brief Iterator category.
    using iterator_concept      = std::forward_iterator_tag; // could be bi

    //!\}

    /*!\name Construction, destruction and assignment
     * \{
     */
    //!\brief Default construction.
    eager_split_iterator() = default;
    //!\brief Copy construction.
    constexpr eager_split_iterator(eager_split_iterator const & rhs) = default;
    //!\brief Move construction.
    constexpr eager_split_iterator(eager_split_iterator && rhs) = default;
    //!\brief Copy assignment.
    constexpr eager_split_iterator & operator=(eager_split_iterator const & rhs) = default;
    //!\brief Move assignment.
    constexpr eager_split_iterator & operator=(eager_split_iterator && rhs) = default;
    //!\brief Destruction.
    ~eager_split_iterator() = default;

    //!\brief Constructing from the underlying seqan3::eager_split_view.
    constexpr eager_split_iterator(std::string_view _urange, char const _delimiter, bool skip_quotes_) noexcept :
        urange{_urange},
        delimiter{_delimiter},
        uend{_urange.data() + _urange.size()},
        subend{_urange.data()},
        skip_quotes{skip_quotes_}
    {
        ++(*this);
    }
    //!\}

    /*!\name Access operations
     * \{
     */
    //!\brief Dereferences the cached iterator.
    constexpr reference operator*() const noexcept
    {
        return subrange;
    }

    //!\brief Returns pointer to the pointed-to object.
    pointer operator->() const noexcept
    {
        return &subrange;
    }
    //!\}

    /*!\name Iterator operations
     * \{
     */
    //!\brief Pre-increment.
    constexpr eager_split_iterator & operator++() noexcept
    {
        auto subbegin = subend;
        if (skip_quotes)
        {
            while ((subend < uend) && (in_quote || (*subend != delimiter)))
            {
                in_quote ^= (*subend == '\"');

                ++subend;
            }
            in_quote = false;
        }
        else
        {
            while ((subend < uend) && (*subend != delimiter))
                ++subend;
        }
        subrange = std::string_view{subbegin, static_cast<size_t>(subend - subbegin)};
        // move behind delimiter
        ++subend;
        return *this;
    }

    //!\brief Post-increment.
    constexpr auto operator++(int) noexcept
    {
        std::string_view tmp{**this};
        ++(*this);
        return tmp;
    }
    //!\}

    /*!\name Comparison operators
     * \{
     */
    //!\brief Compares for equality with sentinel.
    friend constexpr bool operator==(std::default_sentinel_t const &, eager_split_iterator const & rhs) noexcept
    {
        return rhs.subend > rhs.uend + 1; // subend is == uend + 1 on last element and becomes larger on next increment
    }

    //!\copydoc operator==
    friend constexpr bool operator==(eager_split_iterator const & lhs, std::default_sentinel_t const &) noexcept
    {
        return std::default_sentinel == lhs;
    }

    //!\brief Compares for equality with sentinel.
    friend constexpr bool operator!=(std::default_sentinel_t const &, eager_split_iterator const & rhs) noexcept
    {
        return !(rhs == std::default_sentinel);
    }

    //!\copydoc operator==
    friend constexpr bool operator!=(eager_split_iterator const & lhs, std::default_sentinel_t const &) noexcept
    {
        return !(lhs == std::default_sentinel);
    }

   //!\brief Compares for equality with self.
    friend constexpr bool operator==(eager_split_iterator const & lhs, eager_split_iterator const & rhs) noexcept
    {
        return std::tuple(lhs.urange.data(), lhs.urange.size(), lhs.delimiter, lhs.subrange.data(), lhs.subrange.size())
            == std::tuple(rhs.urange.data(), rhs.urange.size(), rhs.delimiter, rhs.subrange.data(), rhs.subrange.size());
    }

    //!\copydoc operator==
    friend constexpr bool operator!=(eager_split_iterator const & lhs, eager_split_iterator const & rhs) noexcept
    {
        return !(lhs == rhs);
    }
    //!\}
};


/*!\brief Adds eager_split behavior to the underlying range.
 * \tparam urng_t The underlying range type.
 * \implements std::ranges::input_range
 * \ingroup views
 */
//![view_def]
class eager_split_view : public std::ranges::view_interface<eager_split_view>
{
//![view_def]
private:

    std::string_view urange;
    char             delimiter;
    bool             skip_quotes;

    /*!\name Member types
     * \{
     */
    //!\brief Iterator type.
    using iterator = eager_split_iterator;
    //!\brief The sentinel type.
    using sentinel = std::default_sentinel_t;
    //\}

public:
    /*!\name Constructor, destructor, and assignment.
     * \{
     * \brief All standard functions are explicitly set to default.
     */
    //!\brief Default default-constructor.
    constexpr eager_split_view() = default;
    //!\brief Default copy-constructor.
    constexpr eager_split_view(eager_split_view const &) = default;
    //!\brief Default move-constructor.
    constexpr eager_split_view(eager_split_view &&) = default;
    //!\brief Default copy-assignment.
    constexpr eager_split_view & operator=(eager_split_view const &) = default;
    //!\brief Default move-assignment
    constexpr eager_split_view & operator=(eager_split_view &&) = default;
    //!\brief Default destructor.
    ~eager_split_view() = default;

    //!\brief Construction from the underlying view.
    explicit constexpr eager_split_view(std::string_view urng_, char delimiter_, bool skip_quotes_ = false) :
        urange{urng_}, delimiter{delimiter_}, skip_quotes{skip_quotes_}
    {}

    //!\}

    /*!\name Iterators
     * \{
     */
    /*!\brief Returns an iterator to the current begin of the underlying range.
     *
     * \details
     *
     * Subsequent calls to begin will result in different positions if the iterator was incremented
     * between the calls.
     */
    constexpr iterator begin() const
    {
        return {urange, delimiter, skip_quotes};
    }

    //!\brief Returns a sentinel.
    constexpr sentinel end() const
    {
        return std::default_sentinel;
    }
    //!\}
};

// ============================================================================
//  eager_split_fn (adaptor definition)
// ============================================================================

/*!\brief views::eager_split's range adaptor object type (non-closure).
 */
struct eager_split_fn
{
    //!\brief Store the argument and return a range adaptor closure object.
    constexpr auto operator()(char const delimiter, bool skip_quotes = false) const noexcept
    {
        return detail::adaptor_from_functor{*this, delimiter, skip_quotes};
    }

    /*!\brief Call the view's constructor with the underlying view as argument.
     * \param[in] urange    The input range to process. Must be a std::string_view.
     * \param[in] delimiter The delimiter to split on.
     * \returns A range of string views
     */
    constexpr auto operator()(std::string_view urange, char const delimiter, bool skip_quotes = false) const noexcept
    {
        return eager_split_view{urange, delimiter, skip_quotes};
    }
};


} // seqan3::detail

//-----------------------------------------------------------------------------
// View shortcut for functor.
//-----------------------------------------------------------------------------

//![adaptor_def]
namespace seqan3::views
{
/*!\name General purpose views
 * \{
 */

/*!\brief               A view adapter that decays most of the range properties and adds single pass behavior.
 * \tparam urng_t       The type of the range being processed. See below for requirements.
 * \param[in] urange    The range being processed.
 * \returns             A range with single pass input behavior. See below for the properties of the returned range.
 * \ingroup views
 *
 * \details
 *
 * \header_file{seqan3/utility/views/eager_split.hpp}
 *
 * TODO all of this
 *
 * ### View properties
 *
 *
 * | Concepts and traits              | `urng_t` (underlying range type)      | `rrng_t` (returned range type)                     |
 * |----------------------------------|:-------------------------------------:|:--------------------------------------------------:|
 * | std::ranges::input_range         | *required*                            | *preserved*                                        |
 * | std::ranges::forward_range       |                                       | *lost*                                             |
 * | std::ranges::bidirectional_range |                                       | *lost*                                             |
 * | std::ranges::random_access_range |                                       | *lost*                                             |
 * | std::ranges::contiguous_range    |                                       | *lost*                                             |
 * |                                  |                                       |                                                    |
 * | std::ranges::viewable_range      | *required*                            | *guaranteed*                                       |
 * | std::ranges::view                |                                       | *guaranteed*                                       |
 * | std::ranges::sized_range         |                                       | *lost*                                             |
 * | std::ranges::common_range        |                                       | *lost*                                             |
 * | std::ranges::output_range        |                                       | *preserved*                                        |
 * | seqan3::const_iterable_range     |                                       | *lost*                                             |
 * |                                  |                                       |                                                    |
 * | std::ranges::range_reference_t   |                                       | std::ranges::range_reference_t<urng_t>             |
 *
 * See the \link views views submodule documentation \endlink for detailed descriptions of the view properties.
 *
 * ### Thread safety
 *
 * Concurrent access to this view, e.g. while iterating over it, is not thread-safe and must be protected externally.
 *
 * ### Example
 *
 * \include test/snippet/utility/views/eager_split.cpp
 * \hideinitializer
 */
inline constexpr auto eager_split = detail::eager_split_fn{};

//!\}
} // namespace seqan3::views
//![adaptor_def]
