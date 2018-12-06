// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
//
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
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
// ==========================================================================
#if __has_include(<seqan/seq_io.h>)
    #include <seqan/seq_io.h>
#endif

#include <algorithm>
#include <cctype>
#include <cstring>

#include <benchmark/benchmark.h>

#include <sstream>

#include <range/v3/algorithm/equal.hpp>
#include <range/v3/view/transform.hpp>

#include <seqan3/alphabet/quality/all.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/input_format_concept.hpp>
#include <seqan3/io/sequence_file/output.hpp>
#include <seqan3/io/sequence_file/output_format_concept.hpp>
#include <seqan3/io/sequence_file/format_fasta.hpp>
#include <seqan3/range/view/convert.hpp>

using namespace seqan3;

static void write3(benchmark::State& state)
{
    std::stringstream ostream; //{"> seq\nACTAGACTAGCTACGATCAGCTACGATCAGCTACGA\n"};
    sequence_file_format_fasta format;
    sequence_file_output_options options;
    std::string id{"seq"};
    dna5_vector seq{"ACGTAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"_dna5};

    for (auto _ : state)
    {
        format.write(ostream, options, seq, id, std::ignore);
    }
}

BENCHMARK(write3);

#if __has_include(<seqan/seq_io.h>)

    static void write2(benchmark::State& state)
    {
        std::stringstream outStream;
        seqan::CharString meta = "seq";
        seqan::Dna5String seq = "ACGTAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";

        for (auto _ : state)
        {
            seqan::writeRecord(outStream, meta, seq, seqan::Fasta());
        }
    }

    BENCHMARK(write2);
#endif

static void read3(benchmark::State& state)
{
    std::string dummy_file{};
    for (size_t idx = 0; idx < 10000000; idx++)
        dummy_file += ">seq\nACTAGACTAGCTACGATCAGCTACGATCAGCTACGA\n";
    std::stringstream istream{dummy_file};
    sequence_file_format_fasta format;
    sequence_file_input_options<dna5, false> options;
    std::string id;
    dna5_vector seq;
    for (auto _ : state)
    {
        format.read(istream, options, seq, id, std::ignore);
	id.clear();
	seq.clear();
    }
}
BENCHMARK(read3);

static void read3_from_file(benchmark::State& state)
{
    // Create fasta file
    std::string filename{"Test_Format_Fasta_Benchmark.fasta"};
   /* sequence_file_output fout{filename};
    for (size_t idx = 0; idx < 10000000; idx++)
    {
        std::string id{ "seq" };
        dna5_vector seq{ "ACTAGACTAGCTACGATCAGCTACGATCAGCTACGA"_dna5 };
        fout.push_back(std::tie(seq, id));   // as a tuple
    }*/

    sequence_file_input fin{filename};
    auto file_it = fin.begin();
    std::string id;
    dna5_vector seq;

    for (auto _ : state)
    {
        [[maybe_unused]] volatile auto record = *file_it;
        ++file_it;
    }
}

BENCHMARK(read3_from_file);


#if __has_include(<seqan/seq_io.h>)

#include <fstream>
/*
std::fstream createFastAFile(seqan::CharString &tempFilename)
{
    using namespace seqan;

    tempFilename = SEQAN_TEMP_FILENAME();
    std::fstream myfile;
    myfile.open(toCString(tempFilename));
    SEQAN_ASSERT(file->is_open());
    for (size_t idx = 0; idx < 10000000; idx++)
        myfile << ">seq\nACTAGACTAGCTACGATCAGCTACGATCAGCTACGA\n";

    return myfile;
}*/

static void read2(benchmark::State& state)
    {
        seqan::CharString meta;
        seqan::Dna5String seq;
        std::string dummy_file{};
        for (size_t idx = 0; idx < 10000000; idx++)
            dummy_file += ">seq\nACTAGACTAGCTACGATCAGCTACGATCAGCTACGA\n";
        std::istringstream istream{dummy_file};

    	std::istreambuf_iterator<char> iter(istream);


        for (auto _ : state)
        {
            seqan::readRecord(meta,seq, iter, seqan::Fasta());
        }

    }
    BENCHMARK(read2);

    static void read2_from_file(benchmark::State& state)
   {
       seqan::CharString filename;
       //std::fstream *file = createFastAFile("Test_Format_Fasta_Benchmark.fasta");
       seqan::SeqFileIn reader("Test_Format_Fasta_Benchmark.fasta");
       seqan::CharString meta;
       seqan::Dna5String seq;
       for (auto _ : state)
       {
           seqan::readRecord(meta, seq, reader, seqan::Fasta());
       }

       //file.close();
   }

    BENCHMARK(read2_from_file);
#endif

BENCHMARK_MAIN();
