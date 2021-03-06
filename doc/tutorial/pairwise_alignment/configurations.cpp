//! [include]
#include <seqan3/alignment/configuration/all.hpp>
//! [include]

//! [include_aligned_ends]
#include <seqan3/alignment/configuration/align_config_aligned_ends.hpp>
//! [include_aligned_ends]

//! [include_scoring_scheme]
#include <seqan3/alignment/scoring/aminoacid_scoring_scheme.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
//! [include_scoring_scheme]

//! [include_gap_scheme]
#include <seqan3/alignment/scoring/gap_scheme.hpp>
//! [include_gap_scheme]

//! [include_result]
#include <seqan3/alignment/configuration/align_config_result.hpp>
//! [include_result]

//! [include_band]
#include <seqan3/alignment/band/static_band.hpp>
#include <seqan3/alignment/configuration/align_config_band.hpp>
//! [include_band]

//! [include_edit]
#include <seqan3/alignment/configuration/align_config_edit.hpp>
//! [include_edit]

#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/alphabet/aminoacid/all.hpp>

int main()
{
{
//! [aligned_ends]
using namespace seqan3;

front_end_first fef{std::true_type{}};
back_end_first bef{std::false_type{}};
front_end_second fes{true};
back_end_second bes{false};

auto cfg_1 = align_cfg::aligned_ends{end_gaps{fef, bef, fes, bes}};
auto cfg_2 = align_cfg::aligned_ends{end_gaps{fef, fes}};
//! [aligned_ends]
(void) cfg_1;
(void) cfg_2;
}

{
//! [scoring_scheme]
using namespace seqan3;

// Define a simple scoring scheme with match and mismatch cost and get the score.
nucleotide_scoring_scheme nc_scheme{match_score{4}, mismatch_score{-5}};
auto sc_nc = nc_scheme.score('A'_dna4, 'C'_dna4); // sc_nc == -5.

// Define a amino acid similarity matrix and get the score.
aminoacid_scoring_scheme aa_scheme{};
aa_scheme.set_similarity_matrix(aminoacid_similarity_matrix::BLOSUM30);
auto sc_aa = aa_scheme.score('M'_aa27, 'K'_aa27); // sc_aa == 2.
//! [scoring_scheme]
(void) sc_nc;
(void) sc_aa;
}

{
//! [gap_scheme]
using namespace seqan3;

// Define a gap scheme with custom gap scores.
gap_scheme g{gap_score{-1}, gap_open_score{-10}};

auto gap = g.get_gap_score();  // gap == -1
auto gap_open = g.get_gap_open_score(); // gap_open == -10
//! [gap_scheme]
(void) gap;
(void) gap_open;
}

{
//! [result]
using namespace seqan3;

// Configure the alignment to only compute the score.
auto cfg = align_cfg::result{with_score};
//! [result]
(void) cfg;
}

{
//! [band]
using namespace seqan3;

// Configure a banded alignment.
auto cfg = align_cfg::band{static_band{lower_bound{-4}, upper_bound{4}}};
//! [band]
(void) cfg;
}

{
//! [edit]
using namespace seqan3;

// Configure an edit distance alignment.
auto cfg = align_cfg::edit;
//! [edit]
(void) cfg;
}
}
