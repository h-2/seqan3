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

#include <any>
#include <type_traits>
#include <vector>

#include <range/v3/view/drop.hpp>
#include <range/v3/view/take.hpp>

#include "concepts.hpp"

namespace seqan3
{
/*!\file concatenated_sequences.hpp
 * \ingroup container
 * \brief Contains concatenated_sequences
 */

/*!\brief Container that stores sequences concatenated internally.
 * \tparam inner_type The type of sequences that will be stored. Must satisfy seqan3::random_access_sequence_concept
 * \todo add complexity to function documentation
 */
template <typename inner_type,
          typename data_delimiters_type = std::vector<typename inner_type::size_type>>
    requires random_access_sequence_concept<inner_type>
    //TODO add requirements on data_delimiters_type
class concatenated_sequences
{
protected:
     //!\privatesection
    std::decay_t<inner_type> data_values;
    data_delimiters_type data_delimiters;
public:
    //!\publicsection
    using value_type = std::decay_t<inner_type>;
    using reference = decltype(data_values | ranges::view::slice(0, 1));
    using const_reference = decltype(data_values | ranges::view::slice(0, 1) | ranges::view::const_);
    using iterator = std::any;//detail::ra_iterator<concatenated_sequences>;
    using const_iterator = std::any;//detail::ra_const_iterator<concatenated_sequences>;
    using difference_type = size_t;
    using size_type = size_t; // TODO must be the same as iterator_traits::difference_type for iterator and const_iterator


    /* rule of six */
    //!\name Constructors, destructor and assignment operators
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


    template <typename type>
        requires input_range_concept<std::decay_t<type>> &&
                 std::is_same_v<std::decay_t<typename type::value_type>, value_type>
    concatenated_sequences & operator=(type const & in)
    {
        clear();

        // benchmark between using join and insert
        data_values = in | ranges::view::join;

        for (auto const & val : in)
        {
//             data_values.insert(data_values.end(), val.begin(), val.end());
            data_delimiters.push_back(data_delimiters.back() + val.size());
        }
        std::cout << data_values << '\n';
        for (auto && i : data_delimiters)
            std::cout << i << ", ";
        std::cout << '\n';
        return *this;
    }
    //!\}

    /*!\name Iterators
     * \{
     */
    /*!\brief Returns an iterator to the first element of the container.
     * \returns Iterator to the first element.
     *
     * If the container is empty, the returned iterator will be equal to end().
     */
    iterator begin() const noexcept
    {
        return iterator{*this};
    }

    /*!\brief Returns an iterator to the first element of the container.
     * \returns Iterator to the first element.
     *
     * If the container is empty, the returned iterator will be equal to end().
     */
    const_iterator cbegin() const noexcept
    {
        return const_iterator{*this};
    }

    /*!\brief Returns an iterator to the element following the last element of the container.
     * \returns Iterator to the first element.
     *
     * This element acts as a placeholder; attempting to access it results in undefined behavior.
     */
    iterator end() const noexcept
    {
        return iterator{*this, true};
    }

    /*!\brief Returns an iterator to the element following the last element of the container.
     * \returns Iterator to the first element.
     *
     * This element acts as a placeholder; attempting to access it results in undefined behavior.
     */
    const_iterator cend() const noexcept
    {
        return const_iterator{*this, true};
    }
    //\}

    /*!\name Element access
     * \{
     */
    /*!\brief Return the i-th element as a view.
     * \param i The element to retrieve.
     * \throws std::out_of_range If you access an element behind the last.
     * \returns A ranges::view on the underlying concatenated sequences that acts as a proxy for the element.
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
     * \returns A ranges::view on the underlying concatenated sequences that acts as a proxy for the element.
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
     * \returns A ranges::view on the underlying concatenated sequences that acts as a proxy for the element.
     *
     * Accessing an element behind the last causes undefined behaviour. In debug mode an assertion checks the size of
     * the container.
     */
    reference operator[](size_type const i)
    {
        assert(i < size());
        return data_values | ranges::view::slice(data_delimiters[i], data_delimiters[i+1]);
    }

    /*!\brief Return the i-th element as a view.
     * \param i The element to retrieve.
     * \returns A ranges::view on the underlying concatenated sequences that acts as a proxy for the element.
     *
     * Accessing an element behind the last causes undefined behaviour. In debug mode an assertion checks the size of
     * the container.
     */
    const_reference operator[](size_type const i) const
    {
        assert(i < size());
        return data_values | ranges::view::slice(data_delimiters[i], data_delimiters[i+1])
                           | ranges::view::const_;
    }

    /*!\brief Return the first element as a view. Calling front on an empty container is undefined.
     * \returns A ranges::view on the underlying concatenated sequences that acts as a proxy for the element.
     *
     * Calling front on an empty container is undefined. In debug mode an assertion checks the size of the container.
     */
    reference front()
    {
        assert(size() > 0);
        return (*this)[0];
    }

    /*!\brief Return the first element as a view.
     * \returns A ranges::view on the underlying concatenated sequences that acts as a proxy for the element.
     *
     * Calling front on an empty container is undefined. In debug mode an assertion checks the size of the container.
     */
    const_reference front() const
    {
        assert(size() > 0);
        return (*this)[0];
    }

    /*!\brief Return the last element as a view.
     * \returns A ranges::view on the underlying concatenated sequences that acts as a proxy for the element.
     *
     * Calling back on an empty container is undefined. In debug mode an assertion checks the size of the container.
     */
    reference back()
    {
        assert(size() > 0);
        return (*this)[size()-1];
    }

    /*!\brief Return the last element as a view.
     * \returns A ranges::view on the underlying concatenated sequences that acts as a proxy for the element.
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
        return data_delimiters.size() - 1;
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
        return data_delimiters.max_size() - 1;
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
        data_values.clear();
        data_delimiters.clear();
        data_delimiters.push_back(0);
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
        data_values.insert(data_values.cbegin() + data_delimiters[pos_as_num], begin(value), end(value));
        data_delimiters.insert(data_delimiters.cbegin() + pos_as_num,
                          *(data_delimiters.cbegin() + pos_as_num));
        // TODO parallel execution policy or vectorization?
        std::for_each(data_delimiters.cbegin() + pos_as_num + 1,
                      data_delimiters.cend(),
                      [&value] (auto & d) { d += value; });
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
            return pos;
        auto pos_as_num = std::distance(pos, cbegin());
        auto repeated = ranges::view::repeat_n(value, count);
        data_values.insert(data_values.cbegin() + data_delimiters[pos_as_num], begin(repeated), end(repeated));
        data_delimiters.insert(data_delimiters.cbegin() + pos_as_num,
                          count,
                          *(data_delimiters.cbegin() + pos_as_num));
        // TODO parallel execution policy or vectorization?
        std::for_each(data_delimiters.cbegin() + pos_as_num + 1,
                      data_delimiters.cbegin() + pos_as_num + count,
                      [&value, factor = 1] (auto & d) mutable { d += value * factor++; });
        // TODO parallel execution policy or vectorization?
        std::for_each(data_delimiters.cbegin() + pos_as_num + count,
                      data_delimiters.cend(),
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
            return pos;
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
        std::swap(data_values, rhs.data_values);
        std::swap(data_delimiters, rhs.data_delimiters);
    }

    void swap(concatenated_sequences && rhs)
    {
        std::swap(data_values, rhs.data_values);
        std::swap(data_delimiters, rhs.data_delimiters);
    }
    //!\}

};

} // namespace seqan3


