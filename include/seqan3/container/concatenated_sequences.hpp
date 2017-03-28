// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2017, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2017, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ============================================================================
// Author: Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
// ============================================================================

#pragma once

#include <vector>
#include <type_traits>

#include <range/v3/view/drop.hpp>
#include <range/v3/view/take.hpp>

#include "concepts.hpp"

namespace seqan3
{
/*!\file concatenated_sequences.hpp
 * \ingroup container
 * \brief Contains \ref concatenated_sequences
 */

/*! \brief Container that stores sequences concatenated internally.
 * \tparam inner_type The type of sequences that will be stored. Must satisfy \ref random_access_sequence_concept
 * \todo add complexity to function documentation
 */
template <typename inner_type,
          typename delimiters_type = std::vector<typename inner_type::size_type>>
    requires random_access_sequence_concept<inner_type>
    //TODO add requirements on delimiters_type
class concatenated_sequences
{
public:
    using value_type = ranges::any_random_access_view<typename inner_type::reference>;
    using reference = value_type;
    using const_reference = value_type const;
//     using iterator = detail::ra_iterator<concatenated_sequences>; //TODO must satisfy forward_iterator_concept and convertible to const_interator
//     using const_iterator = detail::ra_const_iterator<concatenated_sequences>; //TODO must satisfy forward_iterator_concept
    using difference_type = size_t;
    using size_type = size_t; // TODO must be the same as iterator_traits::difference_type for iterator and const_iterator


    /* rule of six */
    //!\name Rule of Six constructors and destructor
    //!\{
    concatenated_sequences()
    {
        clear();
    }
    concatenated_sequences(concatenated_sequences const &) = default;
    concatenated_sequences(concatenated_sequences &&) = default;
    concatenated_sequences & operator=(concatenated_sequences const &) = default;
    concatenated_sequences & operator=(concatenated_sequences &&) = default;
    ~concatenated_sequences() = default;
    //!\}



        template <typename type>
        requires input_range_concept<std::decay_t<type>> //&&
//                  std::is_same_v<std::decay_t<type>::value_type, inner_type>
    concatenated_sequences & operator=(type const & in)
    {
        clear();

        // benchmark between using join and insert
        data = in | ranges::view::join;

        for (auto const & val : in)
        {
//             data.insert(data.end(), val.begin(), val.end());
            delimiters.push_back(delimiters.back() + val.size());
        }

        return *this;
    }


//     iterator begin() const
//     {
//         return iterator{*this};
//     }
//
//     const_iterator cbegin() const
//     {
//         return const_iterator{*this};
//     }
//
//     iterator end() const
//     {
//         return iterator{*this, true};
//     }
//
//     const_iterator cend() const
//     {
//         return const_iterator{*this, true};
//     }



    /*!\name Element access
     * \{
     */
    /*!\brief Return the i-th element as a view.
     * \param i The element to retrieve.
     * \throws std::out_of_range If you access an element behind the last.
     * \returns A ranges::any_random_access_view on the underlying concatenated sequences that acts as a proxy for the element.
     */
    reference at(size_type const i)
    {
        //TODO add SEQAN_UNLIKELY
        if (i >= size())
            throw std::out_of_range{"Trying to access element behind the last in concatenated_sequences."};
        return (*this)[i];
    }

    /*!\brief Return the i-th element as a view.
     * \param i The element to retrieve.
     * \throws std::out_of_range If you access an element behind the last.
     * \returns A ranges::any_random_access_view on the underlying concatenated sequences that acts as a proxy for the element.
     */
    const_reference at(size_type const i) const
    {
        //TODO add SEQAN_UNLIKELY
        if (i >= size())
            throw std::out_of_range{"Trying to access element behind the last in concatenated_sequences."};
        return (*this)[i];
    }

    /*!\brief Return the i-th element as a view.
     * \param i The element to retrieve.
     * \returns A ranges::any_random_access_view on the underlying concatenated sequences that acts as a proxy for the element.
     *
     * Accessing an element behind the last causes undefined behaviour. In debug mode an assertion checks the size of the container.
     */
    reference operator[](size_type const i)
    {
        assert(i < size());
        return data | ranges::view::drop(delimiters[i]) | ranges::view::take(delimiters[i+1] - delimiters[i]);
    }

    /*!\brief Return the i-th element as a view.
     * \param i The element to retrieve.
     * \returns A ranges::any_random_access_view on the underlying concatenated sequences that acts as a proxy for the element.
     *
     * Accessing an element behind the last causes undefined behaviour. In debug mode an assertion checks the size of the container.
     */
    const_reference operator[](size_type const i) const
    {
        assert(i < size());
        return data | ranges::view::drop(delimiters[i]) | ranges::view::take(delimiters[i+1] - delimiters[i]);
    }

    /*!\brief Return the first element as a view. Calling front on an empty container is undefined.
     * \returns A ranges::any_random_access_view on the underlying concatenated sequences that acts as a proxy for the element.
     *
     * Calling front on an empty container is undefined. In debug mode an assertion checks the size of the container.
     */
    reference front()
    {
        assert(size() > 0);
        return (*this)[0];
    }

    /*!\brief Return the first element as a view.
     * \returns A ranges::any_random_access_view on the underlying concatenated sequences that acts as a proxy for the element.
     *
     * Calling front on an empty container is undefined. In debug mode an assertion checks the size of the container.
     */
    const_reference front() const
    {
        assert(size() > 0);
        return (*this)[0];
    }

    /*!\brief Return the last element as a view.
     * \returns A ranges::any_random_access_view on the underlying concatenated sequences that acts as a proxy for the element.
     *
     * Calling back on an empty container is undefined. In debug mode an assertion checks the size of the container.
     */
    reference back()
    {
        assert(size() > 0);
        return (*this)[size()-1];
    }

    /*!\brief Return the last element as a view.
     * \returns A ranges::any_random_access_view on the underlying concatenated sequences that acts as a proxy for the element.
     *
     * Calling back on an empty container is undefined. In debug mode an assertion checks the size of the container.
     */
    const_reference back() const
    {
        assert(size() > 0);
        return (*this)[size()-1];
    }
    //!\}

    /*!\name Capacity
     * \{
     */
    /*!\brief Checks whether the container is empty.
     * \returns `true` if the container is empty, `false` otherwise.
     */
    bool empty() const noexcept
    {
        return size() == 0;
    }

    /*!\brief Returns the number of elements in the container, i.e. std::distance(begin(), end()).
     * \returns The number of elements in the container.
     */
    size_type size() const noexcept
    {
        return delimiters.size() - 1;
    }

    /*!\brief Returns the maximum number of elements the container is able to hold due to system or library
     * implementation limitations, i.e. std::distance(begin(), end()) for the largest container.
     * \returns The number of elements in the container.
     *
     * This value typically reflects the theoretical limit on the size of the container. At runtime, the size
     * of the container may be limited to a value smaller than max_size() by the amount of RAM available.
     */
    size_type max_size() const noexcept
    {
        return delimiters.max_size() - 1;
    }

    /* How do we deal with these? Not clear if inner_type has them
    void reserve()
    {}

    void capacity()
    {}

    void shrink_to_fit()
    {}
    */
    //!\}

    //!\name Modifiers
    //!\{
    /*!\brief Removes all elements from the container.
     * \returns The number of elements in the container.
     */
    void clear() noexcept
    {
        data.clear();
        delimiters.clear();
        delimiters.push_back(0);
    }

    /*!\brief Inserts value before position in the container.
     * \param pos Iterator before which the content will be inserted. `pos` may be the end() iterator.
     * \param value Element value to insert.
     * \returns Iterator pointing to the inserted value.
     *
     * Causes reallocation if the new size() is greater than the old capacity(). If the new size() is greater
     * than capacity(), all iterators and references are invalidated. Otherwise, only the iterators and
     * references before the insertion point remain valid. The past-the-end iterator is also invalidated.
     *
     * \todo exception safety and complexity
     */
    iterator insert(const_iterator pos, value_type const & value)
    {
        auto pos_as_num = std::distance(pos, cbegin());
        data.insert(data.cbegin() + delimiters[pos_as_num], begin(value), end(value));
        delimiters.insert(delimiters.cbegin() + pos_as_num,
                          *(delimiters.cbegin() + pos_as_num));
        // TODO parallel execution policy or vectorization?
        std::for_each(delimiters.cbegin() + pos_as_num + 1,
                      delimiters.cend(),
                      [&value] (auto & d) { d += value });
        return begin() + pos_as_num;
    }
    // no specialization for temporaries, since we have to copy anyway

    /*!\brief Inserts count copies of value before position in the container.
     * \param pos Iterator before which the content will be inserted. `pos` may be the end() iterator.
     * \param count Number of copies.
     * \param value Element value to insert.
     * \returns Iterator pointing to the first element inserted, or pos if `count==0`.
     *
     * Causes reallocation if the new size() is greater than the old capacity(). If the new size() is greater
     * than capacity(), all iterators and references are invalidated. Otherwise, only the iterators and
     * references before the insertion point remain valid. The past-the-end iterator is also invalidated.
     *
     * \todo exception safety and complexity
     */
    iterator insert(const_iterator pos, size_type const count, value_type const & value)
    {
        // TODO SEQAN_UNLIKELY
        if (count == 0)
            return;
        auto pos_as_num = std::distance(pos, cbegin());
        auto repeated = value | ranges::view::repeat_n(count);
        data.insert(data.cbegin() + delimiters[pos_as_num], begin(repeated), end(repeated));
        delimiters.insert(delimiters.cbegin() + pos_as_num,
                          count,
                          *(delimiters.cbegin() + pos_as_num));
        // TODO parallel execution policy or vectorization?
        std::for_each(delimiters.cbegin() + pos_as_num + 1,
                      delimiters.cbegin() + pos_as_num + count,
                      [&value, factor = 1] (auto & d) mutable { d += value * factor++; });
        // TODO parallel execution policy or vectorization?
        std::for_each(delimiters.cbegin() + pos_as_num + count,
                      delimiters.cend(),
                      [&value] (auto & d) { d += value; });
        return begin() + pos_as_num;
    }

    /*!\brief Inserts elements from range `[first, last)` before position in the container.
     * \tparam iterator_type The type of the iterator, must satisfy input_iterator_concept.
     * \param pos Iterator before which the content will be inserted. `pos` may be the end() iterator.
     * \param first Begin of range to insert.
     * \param last Behind the end of range to insert.
     * \returns Iterator pointing to the first element inserted, or pos if `first==last`.
     *
     * The behavior is undefined if first and last are iterators into `*this`.
     *
     * Causes reallocation if the new size() is greater than the old capacity(). If the new size() is greater
     * than capacity(), all iterators and references are invalidated. Otherwise, only the iterators and
     * references before the insertion point remain valid. The past-the-end iterator is also invalidated.
     *
     * \todo exception safety and complexity
     */
    template <typename iterator_type>
    //  requires input_iterator_concept<iterator_type>
    iterator insert(const_iterator pos, iterator_type first, iterator_type last)
    {
        // TODO SEQAN_UNLIKELY
        if (first == last)
            return;
        //TODO
    }

    /*!\brief Inserts elements from initializer list before position in the container.
     * \param pos Iterator before which the content will be inserted. `pos` may be the end() iterator.
     * \param ilist Initializer list with values to insert.
     * \returns Iterator pointing to the first element inserted, or pos if `ilist` is empty.
     *
     * Causes reallocation if the new size() is greater than the old capacity(). If the new size() is greater
     * than capacity(), all iterators and references are invalidated. Otherwise, only the iterators and
     * references before the insertion point remain valid. The past-the-end iterator is also invalidated.
     *
     * \todo exception safety and complexity
     */
    iterator insert(const_iterator pos, std::initializer_list<value_type> ilist)
    {
        //TODO
    }

    // TODO emplace, erase, push_back, emplace_back, pop_back, resize

    //TODO document
    void swap(concatenated_sequences & rhs)
    {
        std::swap(data, rhs.data);
        std::swap(delimiters, rhs.delimiters);
    }

    void swap(concatenated_sequences && rhs)
    {
        std::swap(data, rhs.data);
        std::swap(delimiters, rhs.delimiters);
    }
    //!\}
protected:
    inner_type data;
    delimiters_type delimiters;
};

} // namespace seqan3


