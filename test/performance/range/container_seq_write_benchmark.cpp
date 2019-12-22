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
#include <seqan3/test/performance/sequence_generator.hpp>
#include <seqan3/test/seqan2.hpp>

template <typename t>
using sdsl_int_vec = sdsl::int_vector<sizeof(t) * 8>;

template <typename t>
using small_vec = seqan3::small_vector<t, 10'000>;

// ============================================================================
//  sequential_write
// ============================================================================

template <template <typename> typename container_t, typename alphabet_t>
void sequential_write(benchmark::State & state)
{
    auto cont_rando = seqan3::test::generate_sequence<alphabet_t>(10000, 0, 0);
    container_t<alphabet_t> source(cont_rando.begin(), cont_rando.end());

    alphabet_t a{};
    for (auto _ : state)
    {
        for (auto && c : source)
            c = a;
    }

    state.counters["sizeof"] = sizeof(alphabet_t);
    if constexpr (seqan3::alphabet<alphabet_t>)
        state.counters["alph_size"] = seqan3::alphabet_size<alphabet_t>;
}

BENCHMARK_TEMPLATE(sequential_write, std::vector, char);
BENCHMARK_TEMPLATE(sequential_write, std::vector, uint8_t);
BENCHMARK_TEMPLATE(sequential_write, std::vector, uint16_t);
BENCHMARK_TEMPLATE(sequential_write, std::vector, uint32_t);
BENCHMARK_TEMPLATE(sequential_write, std::vector, uint64_t);
BENCHMARK_TEMPLATE(sequential_write, std::vector, seqan3::gap);
BENCHMARK_TEMPLATE(sequential_write, std::vector, seqan3::dna4);
BENCHMARK_TEMPLATE(sequential_write, std::vector, seqan3::gapped<seqan3::dna4>);
BENCHMARK_TEMPLATE(sequential_write, std::vector, seqan3::dna15);
BENCHMARK_TEMPLATE(sequential_write, std::vector, seqan3::aa27);
BENCHMARK_TEMPLATE(sequential_write, std::vector, seqan3::alphabet_variant<char, seqan3::dna4>);

BENCHMARK_TEMPLATE(sequential_write, std::deque, char);
BENCHMARK_TEMPLATE(sequential_write, std::deque, uint8_t);
BENCHMARK_TEMPLATE(sequential_write, std::deque, uint16_t);
BENCHMARK_TEMPLATE(sequential_write, std::deque, uint32_t);
BENCHMARK_TEMPLATE(sequential_write, std::deque, uint64_t);
BENCHMARK_TEMPLATE(sequential_write, std::deque, seqan3::gap);
BENCHMARK_TEMPLATE(sequential_write, std::deque, seqan3::dna4);
BENCHMARK_TEMPLATE(sequential_write, std::deque, seqan3::gapped<seqan3::dna4>);
BENCHMARK_TEMPLATE(sequential_write, std::deque, seqan3::dna15);
BENCHMARK_TEMPLATE(sequential_write, std::deque, seqan3::aa27);
BENCHMARK_TEMPLATE(sequential_write, std::deque, seqan3::alphabet_variant<char, seqan3::dna4>);

BENCHMARK_TEMPLATE(sequential_write, std::list, char);
BENCHMARK_TEMPLATE(sequential_write, std::list, uint8_t);
BENCHMARK_TEMPLATE(sequential_write, std::list, uint16_t);
BENCHMARK_TEMPLATE(sequential_write, std::list, uint32_t);
BENCHMARK_TEMPLATE(sequential_write, std::list, uint64_t);
BENCHMARK_TEMPLATE(sequential_write, std::list, seqan3::gap);
BENCHMARK_TEMPLATE(sequential_write, std::list, seqan3::dna4);
BENCHMARK_TEMPLATE(sequential_write, std::list, seqan3::gapped<seqan3::dna4>);
BENCHMARK_TEMPLATE(sequential_write, std::list, seqan3::dna15);
BENCHMARK_TEMPLATE(sequential_write, std::list, seqan3::aa27);
BENCHMARK_TEMPLATE(sequential_write, std::list, seqan3::alphabet_variant<char, seqan3::dna4>);

BENCHMARK_TEMPLATE(sequential_write, sdsl_int_vec, uint8_t);
BENCHMARK_TEMPLATE(sequential_write, sdsl_int_vec, uint16_t);
BENCHMARK_TEMPLATE(sequential_write, sdsl_int_vec, uint32_t);
BENCHMARK_TEMPLATE(sequential_write, sdsl_int_vec, uint64_t);

BENCHMARK_TEMPLATE(sequential_write, seqan3::bitcompressed_vector, char);
BENCHMARK_TEMPLATE(sequential_write, seqan3::bitcompressed_vector, uint32_t);
BENCHMARK_TEMPLATE(sequential_write, seqan3::bitcompressed_vector, seqan3::gap);
BENCHMARK_TEMPLATE(sequential_write, seqan3::bitcompressed_vector, seqan3::dna4);
BENCHMARK_TEMPLATE(sequential_write, seqan3::bitcompressed_vector, seqan3::gapped<seqan3::dna4>);
BENCHMARK_TEMPLATE(sequential_write, seqan3::bitcompressed_vector, seqan3::dna15);
BENCHMARK_TEMPLATE(sequential_write, seqan3::bitcompressed_vector, seqan3::aa27);
BENCHMARK_TEMPLATE(sequential_write, seqan3::bitcompressed_vector, seqan3::alphabet_variant<char, seqan3::dna4>);

BENCHMARK_TEMPLATE(sequential_write, small_vec, char);
BENCHMARK_TEMPLATE(sequential_write, small_vec, uint32_t);
BENCHMARK_TEMPLATE(sequential_write, small_vec, seqan3::gap);
BENCHMARK_TEMPLATE(sequential_write, small_vec, seqan3::dna4);
BENCHMARK_TEMPLATE(sequential_write, small_vec, seqan3::gapped<seqan3::dna4>);
BENCHMARK_TEMPLATE(sequential_write, small_vec, seqan3::dna15);
BENCHMARK_TEMPLATE(sequential_write, small_vec, seqan3::aa27);
BENCHMARK_TEMPLATE(sequential_write, small_vec, seqan3::alphabet_variant<char, seqan3::dna4>);

// ============================================================================
//  SeqAn2: sequential_write
// ============================================================================

#if SEQAN3_HAS_SEQAN2

#include <seqan/sequence.h>

template <template <typename...> typename container_t, typename spec_t, typename alphabet_t>
void sequential_write2(benchmark::State & state)
{
    auto cont_rando = seqan3::test::generate_sequence_seqan2<alphabet_t>(10'000, 0, 0);
    container_t<alphabet_t, spec_t> source = cont_rando;

    alphabet_t a{};
    for (auto _ : state)
    {
        for (auto && c : source)
            c = a;
    }

    state.counters["sizeof"] = sizeof(alphabet_t);
    state.counters["alph_size"] = seqan::ValueSize<alphabet_t>::VALUE;
}

BENCHMARK_TEMPLATE(sequential_write2, seqan::String, seqan::Alloc<>, char);
BENCHMARK_TEMPLATE(sequential_write2, seqan::String, seqan::Alloc<>, uint8_t);
BENCHMARK_TEMPLATE(sequential_write2, seqan::String, seqan::Alloc<>, uint16_t);
BENCHMARK_TEMPLATE(sequential_write2, seqan::String, seqan::Alloc<>, uint32_t);
BENCHMARK_TEMPLATE(sequential_write2, seqan::String, seqan::Alloc<>, uint64_t);

BENCHMARK_TEMPLATE(sequential_write2, seqan::String, seqan::Alloc<>, seqan::Dna);
BENCHMARK_TEMPLATE(sequential_write2, seqan::String, seqan::Alloc<>, seqan::Dna5);
BENCHMARK_TEMPLATE(sequential_write2, seqan::String, seqan::Alloc<>, seqan::Iupac);
BENCHMARK_TEMPLATE(sequential_write2, seqan::String, seqan::Alloc<>, seqan::AminoAcid);
BENCHMARK_TEMPLATE(sequential_write2, seqan::String, seqan::Alloc<>, seqan::Dna5Q);

BENCHMARK_TEMPLATE(sequential_write2, seqan::String, seqan::Packed<>, seqan::Dna);
BENCHMARK_TEMPLATE(sequential_write2, seqan::String, seqan::Packed<>, seqan::Dna5);
BENCHMARK_TEMPLATE(sequential_write2, seqan::String, seqan::Packed<>, seqan::Iupac);
BENCHMARK_TEMPLATE(sequential_write2, seqan::String, seqan::Packed<>, seqan::AminoAcid);
// BENCHMARK_TEMPLATE(sequential_write2, seqan::String, seqan::Packed<>, seqan::Dna5Q); // broken in SeqAn2

BENCHMARK_TEMPLATE(sequential_write2, seqan::String, seqan::Array<10'000>, seqan::Dna);
BENCHMARK_TEMPLATE(sequential_write2, seqan::String, seqan::Array<10'000>, seqan::Dna5);
BENCHMARK_TEMPLATE(sequential_write2, seqan::String, seqan::Array<10'000>, seqan::Iupac);
BENCHMARK_TEMPLATE(sequential_write2, seqan::String, seqan::Array<10'000>, seqan::AminoAcid);
BENCHMARK_TEMPLATE(sequential_write2, seqan::String, seqan::Array<10'000>, seqan::Dna5Q);
#endif // SEQAN3_HAS_SEQAN2

// ============================================================================
//  run
// ============================================================================

BENCHMARK_MAIN();
