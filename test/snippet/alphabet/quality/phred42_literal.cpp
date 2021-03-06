#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/core/debug_stream.hpp>

using namespace seqan3;

int main()
{
    // directly assign to a std::vector<phred42> using a string literal
    std::vector<phred42> qual_vec = "###!"_phred42;

    // This is the same as a sequence of char literals:
    std::vector<phred42> qual_vec2 = {'#'_phred42, '#'_phred42, '#'_phred42, '!'_phred42};

    debug_stream << ranges::equal(qual_vec, qual_vec2) << std::endl; // prints 1 (true)
}
