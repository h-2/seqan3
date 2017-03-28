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

template <typename inner_type,
          typename delimiters_type = std::vector<typename inner_type::size_type>>
    requires random_access_sequence_concept<inner_type>
    //TODO add requirements on delimiters_type
class concatenated_sequences
{
public:

    using value_type = ranges::any_random_access_view<typename inner_type::reference>;
    using reference = value_type;
    using const_reference = const value_type;
    using iterator = detail::ra_iterator<concatenated_sequences>; //TODO must satisfy forward_iterator_concept and convertible to const_interator
    using const_iterator = detail::ra_const_iterator<concatenated_sequences>; //TODO must satisfy forward_iterator_concept
    using difference_type = size_t;
    using size_type = size_t; // TODO must be the same as iterator_traits::difference_type for iterator and const_iterator


    /* rule of six */
    concatenated_sequences()
    {
        clear();
    }
    concatenated_sequences(concatenated_sequences const &) = default;
    concatenated_sequences(concatenated_sequences &&) = default;
    concatenated_sequences & operator=(concatenated_sequences const &) = default;
    concatenated_sequences & operator=(concatenated_sequences &&) = default;
    ~concatenated_sequences() = default;

    iterator begin() const
    {
        return iterator{*this, false};
    }

    const_iterator cbegin() const
    {
        return const_iterator{*this, false};
    }

    iterator end() const
    {
        return iterator{*this, true};
    }

    const_iterator cbegin() const
    {
        return const_iterator{*this, true};
    }

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

    size_type max_size() const
    {
        return delimiters.max_size() - 1;
    }

    bool empty() const
    {
        return size() == 0;
    }

    //TODO all the stuff from higher container concepts

    template <typename type>
        requires input_range_concept<std::decay_t<type>> //&&
//                  std::is_same_v<std::decay_t<type>::value_type, inner_type>
    concatenated_sequences & operator=(type const & in)
    {
        clear();

        // this would be nice, but doesn't work right now
        // data = in | ranges::view::concat;

        for (auto const & val : in)
        {
            data.insert(data.end(), val.begin(), val.end());
            delimiters.push_back(delimiters.back() + val.size());
        }

        return *this;
    }

    size_t size()
    {
        return delimiters.size() - 1;
    }

    void clear()
    {
        data.clear();
        delimiters.clear();
        delimiters.push_back(0);
    }

    reference operator[](size_type i)
    {
        assert(i < size());
        return data | ranges::view::drop(delimiters[i]) | ranges::view::take(delimiters[i+1] - delimiters[i]);
    }

protected:
    inner_type data;
    delimiters_type delimiters;
};

} // namespace seqan3


