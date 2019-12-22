// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <deque>
#include <list>
#include <vector>

#include <benchmark/benchmark.h>

#include <seqan3/alphabet/all.hpp>
#include <seqan3/range/container/all.hpp>
#include <seqan3/test/seqan2.hpp>

template <typename t>
using sdsl_int_vec = sdsl::int_vector<sizeof(t) * 8>;

template <typename t>
using small_vec = seqan3::small_vector<t, 10'000>;

// ============================================================================
//  push_back
// ============================================================================

template <template <typename> typename container_t, typename alphabet_t>
void push_back(benchmark::State & state)
{
    alphabet_t a{};

    for (auto _ : state)
    {
        container_t<alphabet_t> c;
        for (size_t i = 0; i < 10'000; ++i)
            c.push_back(a);
        benchmark::DoNotOptimize(a = c.back());
    }

    state.counters["sizeof"] = sizeof(alphabet_t);
    if constexpr (seqan3::alphabet<alphabet_t>)
        state.counters["alph_size"] = seqan3::alphabet_size<alphabet_t>;
}

BENCHMARK_TEMPLATE(push_back, std::vector, char);
BENCHMARK_TEMPLATE(push_back, std::vector, uint8_t);
BENCHMARK_TEMPLATE(push_back, std::vector, uint16_t);
BENCHMARK_TEMPLATE(push_back, std::vector, uint32_t);
BENCHMARK_TEMPLATE(push_back, std::vector, uint64_t);
BENCHMARK_TEMPLATE(push_back, std::vector, seqan3::gap);
BENCHMARK_TEMPLATE(push_back, std::vector, seqan3::dna4);
BENCHMARK_TEMPLATE(push_back, std::vector, seqan3::gapped<seqan3::dna4>);
BENCHMARK_TEMPLATE(push_back, std::vector, seqan3::dna15);
BENCHMARK_TEMPLATE(push_back, std::vector, seqan3::aa27);
BENCHMARK_TEMPLATE(push_back, std::vector, seqan3::alphabet_variant<char, seqan3::dna4>);

BENCHMARK_TEMPLATE(push_back, std::deque, char);
BENCHMARK_TEMPLATE(push_back, std::deque, uint8_t);
BENCHMARK_TEMPLATE(push_back, std::deque, uint16_t);
BENCHMARK_TEMPLATE(push_back, std::deque, uint32_t);
BENCHMARK_TEMPLATE(push_back, std::deque, uint64_t);
BENCHMARK_TEMPLATE(push_back, std::deque, seqan3::gap);
BENCHMARK_TEMPLATE(push_back, std::deque, seqan3::dna4);
BENCHMARK_TEMPLATE(push_back, std::deque, seqan3::gapped<seqan3::dna4>);
BENCHMARK_TEMPLATE(push_back, std::deque, seqan3::dna15);
BENCHMARK_TEMPLATE(push_back, std::deque, seqan3::aa27);
BENCHMARK_TEMPLATE(push_back, std::deque, seqan3::alphabet_variant<char, seqan3::dna4>);

BENCHMARK_TEMPLATE(push_back, std::list, char);
BENCHMARK_TEMPLATE(push_back, std::list, uint8_t);
BENCHMARK_TEMPLATE(push_back, std::list, uint16_t);
BENCHMARK_TEMPLATE(push_back, std::list, uint32_t);
BENCHMARK_TEMPLATE(push_back, std::list, uint64_t);
BENCHMARK_TEMPLATE(push_back, std::list, seqan3::gap);
BENCHMARK_TEMPLATE(push_back, std::list, seqan3::dna4);
BENCHMARK_TEMPLATE(push_back, std::list, seqan3::gapped<seqan3::dna4>);
BENCHMARK_TEMPLATE(push_back, std::list, seqan3::dna15);
BENCHMARK_TEMPLATE(push_back, std::list, seqan3::aa27);
BENCHMARK_TEMPLATE(push_back, std::list, seqan3::alphabet_variant<char, seqan3::dna4>);

BENCHMARK_TEMPLATE(push_back, sdsl_int_vec, uint8_t);
BENCHMARK_TEMPLATE(push_back, sdsl_int_vec, uint16_t);
BENCHMARK_TEMPLATE(push_back, sdsl_int_vec, uint32_t);
BENCHMARK_TEMPLATE(push_back, sdsl_int_vec, uint64_t);

BENCHMARK_TEMPLATE(push_back, seqan3::bitcompressed_vector, char);
BENCHMARK_TEMPLATE(push_back, seqan3::bitcompressed_vector, uint32_t);
BENCHMARK_TEMPLATE(push_back, seqan3::bitcompressed_vector, seqan3::gap);
BENCHMARK_TEMPLATE(push_back, seqan3::bitcompressed_vector, seqan3::dna4);
BENCHMARK_TEMPLATE(push_back, seqan3::bitcompressed_vector, seqan3::gapped<seqan3::dna4>);
BENCHMARK_TEMPLATE(push_back, seqan3::bitcompressed_vector, seqan3::dna15);
BENCHMARK_TEMPLATE(push_back, seqan3::bitcompressed_vector, seqan3::aa27);
BENCHMARK_TEMPLATE(push_back, seqan3::bitcompressed_vector, seqan3::alphabet_variant<char, seqan3::dna4>);

BENCHMARK_TEMPLATE(push_back, small_vec, char);
BENCHMARK_TEMPLATE(push_back, small_vec, uint32_t);
BENCHMARK_TEMPLATE(push_back, small_vec, seqan3::gap);
BENCHMARK_TEMPLATE(push_back, small_vec, seqan3::dna4);
BENCHMARK_TEMPLATE(push_back, small_vec, seqan3::gapped<seqan3::dna4>);
BENCHMARK_TEMPLATE(push_back, small_vec, seqan3::dna15);
BENCHMARK_TEMPLATE(push_back, small_vec, seqan3::aa27);
BENCHMARK_TEMPLATE(push_back, small_vec, seqan3::alphabet_variant<char, seqan3::dna4>);

// ============================================================================
//  push_back SeqAn2
// ============================================================================

#if SEQAN3_HAS_SEQAN2

#include <seqan/sequence.h>

template <template <typename...> typename container_t, typename spec_t, typename alphabet_t>
void push_back2(benchmark::State & state)
{
    alphabet_t a{};

    for (auto _ : state)
    {
        container_t<alphabet_t, spec_t> c;
        for (size_t i = 0; i < 10'000; ++i)
            seqan::appendValue(c, a);
        a = seqan::back(c);
    }

    state.counters["sizeof"] = sizeof(alphabet_t);
    state.counters["alph_size"] = seqan::ValueSize<alphabet_t>::VALUE;
}

BENCHMARK_TEMPLATE(push_back, std::vector, seqan::Dna);
BENCHMARK_TEMPLATE(push_back, std::vector, seqan::Dna5);
BENCHMARK_TEMPLATE(push_back, std::vector, seqan::Iupac);
BENCHMARK_TEMPLATE(push_back, std::vector, seqan::AminoAcid);
BENCHMARK_TEMPLATE(push_back, std::vector, seqan::Dna5Q);

BENCHMARK_TEMPLATE(push_back2, seqan::String, seqan::Alloc<>, char);
BENCHMARK_TEMPLATE(push_back2, seqan::String, seqan::Alloc<>, uint8_t);
BENCHMARK_TEMPLATE(push_back2, seqan::String, seqan::Alloc<>, uint16_t);
BENCHMARK_TEMPLATE(push_back2, seqan::String, seqan::Alloc<>, uint32_t);
BENCHMARK_TEMPLATE(push_back2, seqan::String, seqan::Alloc<>, uint64_t);

BENCHMARK_TEMPLATE(push_back2, seqan::String, seqan::Alloc<>, seqan3::gap);
BENCHMARK_TEMPLATE(push_back2, seqan::String, seqan::Alloc<>, seqan3::dna4);
BENCHMARK_TEMPLATE(push_back2, seqan::String, seqan::Alloc<>, seqan3::gapped<seqan3::dna4>);
BENCHMARK_TEMPLATE(push_back2, seqan::String, seqan::Alloc<>, seqan3::dna15);
BENCHMARK_TEMPLATE(push_back2, seqan::String, seqan::Alloc<>, seqan3::aa27);
BENCHMARK_TEMPLATE(push_back2, seqan::String, seqan::Alloc<>, seqan3::alphabet_variant<char, seqan3::dna4>);

BENCHMARK_TEMPLATE(push_back2, seqan::String, seqan::Alloc<>, seqan::Dna);
BENCHMARK_TEMPLATE(push_back2, seqan::String, seqan::Alloc<>, seqan::Dna5);
BENCHMARK_TEMPLATE(push_back2, seqan::String, seqan::Alloc<>, seqan::Iupac);
BENCHMARK_TEMPLATE(push_back2, seqan::String, seqan::Alloc<>, seqan::AminoAcid);
BENCHMARK_TEMPLATE(push_back2, seqan::String, seqan::Alloc<>, seqan::Dna5Q);

BENCHMARK_TEMPLATE(push_back2, seqan::String, seqan::Packed<>, seqan::Dna);
BENCHMARK_TEMPLATE(push_back2, seqan::String, seqan::Packed<>, seqan::Dna5);
BENCHMARK_TEMPLATE(push_back2, seqan::String, seqan::Packed<>, seqan::Iupac);
BENCHMARK_TEMPLATE(push_back2, seqan::String, seqan::Packed<>, seqan::AminoAcid);
// BENCHMARK_TEMPLATE(push_back2, seqan::String, seqan::Packed<>, seqan::Dna5Q); // broken in SeqAn2

BENCHMARK_TEMPLATE(push_back2, seqan::String, seqan::Array<10'000>, seqan::Dna);
BENCHMARK_TEMPLATE(push_back2, seqan::String, seqan::Array<10'000>, seqan::Dna5);
BENCHMARK_TEMPLATE(push_back2, seqan::String, seqan::Array<10'000>, seqan::Iupac);
BENCHMARK_TEMPLATE(push_back2, seqan::String, seqan::Array<10'000>, seqan::AminoAcid);
BENCHMARK_TEMPLATE(push_back2, seqan::String, seqan::Array<10'000>, seqan::Dna5Q);

#endif

// ============================================================================
//  run
// ============================================================================

BENCHMARK_MAIN();
