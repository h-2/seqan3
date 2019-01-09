// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2018, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
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

#include <gtest/gtest.h>

#include <seqan3/core/metafunction/function.hpp>

constexpr
int    constexpr_nonvoid_free_fun(int i) { return i; }
int nonconstexpr_nonvoid_free_fun(int i) { return i; }

constexpr
int constexpr_nonvoid_free_fun_const_ref(int const & i) { return i; }
int nonconstexpr_nonvoid_free_fun_const_ref(int const & i) { return i; }

constexpr
void    constexpr_void_free_fun(int) { return; }
void nonconstexpr_void_free_fun(int) { return; }

struct constexpr_nonvoid_member_t
{
    int constexpr get_i(int i)
    {
        return i;
    }
};

struct constexpr_void_member_t
{
    void constexpr get_i(int)
    {}
};

struct nonconstexpr_nonvoid_member_t
{
    int get_i(int i)
    {
        return i;
    }
};

struct nonconstexpr_void_member_t
{
    void get_i(int)
    {}
};

TEST(metafunction, is_constexpr_invocable)
{
    int i = 32;
    int constexpr j = 42;

    EXPECT_TRUE((SEQAN3_IS_CONSTEXPR(constexpr_nonvoid_free_fun(3))));
    EXPECT_TRUE((SEQAN3_IS_CONSTEXPR(constexpr_nonvoid_free_fun(j))));
    EXPECT_TRUE((!SEQAN3_IS_CONSTEXPR(constexpr_nonvoid_free_fun(i))));
    EXPECT_TRUE((!SEQAN3_IS_CONSTEXPR(nonconstexpr_nonvoid_free_fun(3))));

    EXPECT_TRUE((SEQAN3_IS_CONSTEXPR(constexpr_nonvoid_free_fun_const_ref(static_cast<int const &>(3)))));
    EXPECT_TRUE((SEQAN3_IS_CONSTEXPR(constexpr_nonvoid_free_fun_const_ref(j))));
    EXPECT_TRUE((!SEQAN3_IS_CONSTEXPR(constexpr_nonvoid_free_fun_const_ref(i))));
    EXPECT_TRUE((!SEQAN3_IS_CONSTEXPR(nonconstexpr_nonvoid_free_fun_const_ref(static_cast<int const &>(3)))));

    EXPECT_TRUE((SEQAN3_IS_CONSTEXPR(constexpr_void_free_fun(3))));
    EXPECT_TRUE((SEQAN3_IS_CONSTEXPR(constexpr_void_free_fun(j))));
    EXPECT_TRUE((!SEQAN3_IS_CONSTEXPR(constexpr_void_free_fun(i))));
    EXPECT_TRUE((!SEQAN3_IS_CONSTEXPR(nonconstexpr_void_free_fun(3))));

    EXPECT_TRUE((SEQAN3_IS_CONSTEXPR(constexpr_nonvoid_member_t{}.get_i(3))));
    EXPECT_TRUE((SEQAN3_IS_CONSTEXPR(constexpr_void_member_t{}.get_i(3))));
    EXPECT_TRUE((!SEQAN3_IS_CONSTEXPR(nonconstexpr_nonvoid_member_t{}.get_i(3))));
    EXPECT_TRUE((!SEQAN3_IS_CONSTEXPR(nonconstexpr_void_member_t{}.get_i(3))));
}
