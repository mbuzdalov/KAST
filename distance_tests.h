#include <ranges>
#include <vector>
#include <iostream>

#include <seqan3/alphabet/all.hpp>

#include "distance.h"
#include "utils.h"

using namespace seqan3::literals;

template<typename T, typename R>
void append_range(std::vector<T> &destination, R const &range) {
    destination.insert(destination.end(), range.begin(), range.end());
}

/*
   Prep the query and reference sequences and perform counts
*/
void prep(std::vector<unsigned> &qrycounts, std::vector<unsigned> &refcounts, unsigned k)
{
   std::vector<seqan3::dna5> qryseq { "AGGCAGCGTACGAACCTACTGGAGTTGCGGTATGGGACCAGGCGACCTCTGATGCAGAGATACAGGAGCGCCGCGCCGGGTCTTCCTTGTAGAAGTCCTG"_dna5 };
   std::vector<seqan3::dna5> refseq { "CGGAGACCTCCGTGGACGGGGAAGTCCTGCGCGGGTCAGACGTACGCCCCGATTAGTTGCCCGGACGCCCGGTTGGCAGAAGTGACGGCGACTGCCCTCA"_dna5 };

   // TODO: reverse-complementing onto itself should possibly be doable without copying,
   // and possibly be contained in utils
   std::vector<seqan3::dna5> qryseqrc = qryseq;
   append_range(qryseq, "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"_dna5); // this should probably the same size as options.klen
   append_range(qryseq, qryseqrc | std::views::reverse | seqan3::views::complement);

   std::vector<seqan3::dna5> refseqrc = refseq;
   append_range(refseq, "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"_dna5); // this should probably the same size as options.klen
   append_range(refseq, refseqrc | std::views::reverse | seqan3::views::complement);

   count_kmers(qrycounts, qryseq, k);
   count_kmers(refcounts, refseq, k);
}

void prep(std::vector<unsigned> &qrycounts, std::vector<unsigned> &refcounts,
          std::vector<double> &qrymarkov, std::vector<double> &refmarkov,
          unsigned k, unsigned markov_order)
{
   std::vector<seqan3::dna5> qryseq { "AGGCAGCGTACGAACCTACTGGAGTTGCGGTATGGGACCAGGCGACCTCTGATGCAGAGATACAGGAGCGCCGCGCCGGGTCTTCCTTGTAGAAGTCCTG"_dna5 };
   std::vector<seqan3::dna5> refseq { "CGGAGACCTCCGTGGACGGGGAAGTCCTGCGCGGGTCAGACGTACGCCCCGATTAGTTGCCCGGACGCCCGGTTGGCAGAAGTGACGGCGACTGCCCTCA"_dna5 };

   // TODO: reverse-complementing onto itself should possibly be doable without copying
   // and possibly be contained in utils
   std::vector<seqan3::dna5> qryseqrc = qryseq;
   append_range(qryseq, "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"_dna5); // this should probably the same size as options.klen
   append_range(qryseq, qryseqrc | std::views::reverse | seqan3::views::complement);

   std::vector<seqan3::dna5> refseqrc = refseq;
   append_range(refseq, "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"_dna5); // this should probably the same size as options.klen
   append_range(refseq, refseqrc | std::views::reverse | seqan3::views::complement);

   count_kmers(qrycounts, qryseq, k);
   count_kmers(refcounts, refseq, k);

   markov(qrymarkov, qryseq, k, markov_order);
   markov(refmarkov, refseq, k, markov_order);
}

void prep_aa(std::vector<unsigned> &qrycounts, std::vector<unsigned> &refcounts,
             unsigned k)
{
   std::vector<seqan3::aa20> qryseq { "MTEITAAMVKELRESTGAGMMDCKNALSETNGDFDKAVQLLREKGLGKAAKKADRLAAEGLVSVKVSDDFTIAAMRPSYLSYEDLDMTFVENEYKALVAELEKENEERRRLKDPNKPEHKIPQFASRKQLSDAILKEAEEKIKEELKAQGKPEKIWDNIIPGKMNSFIADNSQLDSKLTLMGQFYVMDDKKTVEQVIAEKEKEFGGKIKIVEFICFEVGEGLEKKTEDFAAEVAAQL"_aa20 };
   std::vector<seqan3::aa20> refseq { "SATVSEINSETDFVAKNDQFIALTKDTTAHIQSNSLQSVEELHSSTINGVKFEEYLKSQIATIGENLVVRRFATLKAGANGVVNGYIHTNGRVGVVIAAACDSAEVASKSRDLLRQICMH"_aa20 };

   count_kmers(qrycounts, qryseq, k);
   count_kmers(refcounts, refseq, k);
}

void prep_aa(std::vector<unsigned> &qrycounts, std::vector<unsigned> &refcounts,
             std::vector<double> &qrymarkov, std::vector<double> &refmarkov,
             unsigned k, unsigned markov_order)
{
   std::vector<seqan3::aa20> qryseq { "MTEITAAMVKELRESTGAGMMDCKNALSETNGDFDKAVQLLREKGLGKAAKKADRLAAEGLVSVKVSDDFTIAAMRPSYLSYEDLDMTFVENEYKALVAELEKENEERRRLKDPNKPEHKIPQFASRKQLSDAILKEAEEKIKEELKAQGKPEKIWDNIIPGKMNSFIADNSQLDSKLTLMGQFYVMDDKKTVEQVIAEKEKEFGGKIKIVEFICFEVGEGLEKKTEDFAAEVAAQL"_aa20 };
   std::vector<seqan3::aa20> refseq { "SATVSEINSETDFVAKNDQFIALTKDTTAHIQSNSLQSVEELHSSTINGVKFEEYLKSQIATIGENLVVRRFATLKAGANGVVNGYIHTNGRVGVVIAAACDSAEVASKSRDLLRQICMH"_aa20 };

   count_kmers(qrycounts, qryseq, k);
   count_kmers(refcounts, refseq, k);

   markov(qrymarkov, qryseq, k, markov_order);
   markov(refmarkov, refseq, k, markov_order);
}

void prep_raa(std::vector<unsigned> &qrycounts, std::vector<unsigned> &refcounts,
              unsigned k)
{
   std::vector<seqan3::aa10murphy> qryseq { "MTEITAAMVKELRESTGAGMMDCKNALSETNGDFDKAVQLLREKGLGKAAKKADRLAAEGLVSVKVSDDFTIAAMRPSYLSYEDLDMTFVENEYKALVAELEKENEERRRLKDPNKPEHKIPQFASRKQLSDAILKEAEEKIKEELKAQGKPEKIWDNIIPGKMNSFIADNSQLDSKLTLMGQFYVMDDKKTVEQVIAEKEKEFGGKIKIVEFICFEVGEGLEKKTEDFAAEVAAQL"_aa10murphy };
   std::vector<seqan3::aa10murphy> refseq { "SATVSEINSETDFVAKNDQFIALTKDTTAHIQSNSLQSVEELHSSTINGVKFEEYLKSQIATIGENLVVRRFATLKAGANGVVNGYIHTNGRVGVVIAAACDSAEVASKSRDLLRQICMH"_aa10murphy };

   count_kmers(qrycounts, qryseq, k);
   count_kmers(refcounts, refseq, k);
}

void prep_raa(std::vector<unsigned> &qrycounts, std::vector<unsigned> &refcounts,
              std::vector<double> &qrymarkov, std::vector<double> &refmarkov,
              unsigned k, unsigned markov_order)
{
   std::vector<seqan3::aa10murphy> qryseq { "MTEITAAMVKELRESTGAGMMDCKNALSETNGDFDKAVQLLREKGLGKAAKKADRLAAEGLVSVKVSDDFTIAAMRPSYLSYEDLDMTFVENEYKALVAELEKENEERRRLKDPNKPEHKIPQFASRKQLSDAILKEAEEKIKEELKAQGKPEKIWDNIIPGKMNSFIADNSQLDSKLTLMGQFYVMDDKKTVEQVIAEKEKEFGGKIKIVEFICFEVGEGLEKKTEDFAAEVAAQL"_aa10murphy };
   std::vector<seqan3::aa10murphy> refseq { "SATVSEINSETDFVAKNDQFIALTKDTTAHIQSNSLQSVEELHSSTINGVKFEEYLKSQIATIGENLVVRRFATLKAGANGVVNGYIHTNGRVGVVIAAACDSAEVASKSRDLLRQICMH"_aa10murphy };

   count_kmers(qrycounts, qryseq, k);
   count_kmers(refcounts, refseq, k);

   markov(qrymarkov, qryseq, k, markov_order);
   markov(refmarkov, refseq, k, markov_order);
}

/////////////////////////////////////////////////////////

#define FAIL_IF_NOT_WITHIN(found, expected, eps) { \
    double found_val = found; \
    if (found_val < expected - eps || found_val > expected + eps) { \
        seqan3::debug_stream << "Test " << __func__ << "/" << #found << " failed: expected " << expected << " found " << found_val << std::endl; \
        exit(1); \
    } else { \
        seqan3::debug_stream << "Test " << __func__ << "/" << #found << " OK" << std::endl; \
    } \
}

void tests_prep_dna_3() {
    std::vector<unsigned> qrycounts, refcounts;
    std::vector<double> qrymarkov, refmarkov;
    prep(qrycounts, refcounts, qrymarkov, refmarkov, 3, 0);

    FAIL_IF_NOT_WITHIN(d2(refcounts, qrycounts), 0.10619, 0.0001);
    FAIL_IF_NOT_WITHIN(euler(refcounts, qrycounts), 0.10306, 0.0001);
    FAIL_IF_NOT_WITHIN(manhattan(refcounts, qrycounts), 0.63265, 0.0001);
    FAIL_IF_NOT_WITHIN(bray_curtis_distance(refcounts, qrycounts), 0.316326530612, 0.0001);
    FAIL_IF_NOT_WITHIN(normalised_google_distance(refcounts, qrycounts), 0.316326530612, 0.0001);
    FAIL_IF_NOT_WITHIN(chebyshev(refcounts, qrycounts), 0.03061, 0.0001);

    FAIL_IF_NOT_WITHIN(d2s(refcounts, qrycounts, refmarkov, qrymarkov), 0.43246389442301436, 0.0001);
    FAIL_IF_NOT_WITHIN(d2star(refcounts, qrycounts, refmarkov, qrymarkov), 0.4027100011247771, 0.0001);
}

void tests_prep_dna_5() {
    std::vector<unsigned> qrycounts, refcounts;
    std::vector<double> qrymarkov, refmarkov;
    prep(qrycounts, refcounts, qrymarkov, refmarkov, 5, 0);

    FAIL_IF_NOT_WITHIN(d2(refcounts, qrycounts), 0.3537, 0.0001);
    FAIL_IF_NOT_WITHIN(euler(refcounts, qrycounts), 0.09317, 0.0001);
    FAIL_IF_NOT_WITHIN(manhattan(refcounts, qrycounts), 1.45833, 0.0001);
    FAIL_IF_NOT_WITHIN(bray_curtis_distance(refcounts, qrycounts), 0.729166666667, 0.0001);
    FAIL_IF_NOT_WITHIN(normalised_google_distance(refcounts, qrycounts), 0.729166666667, 0.0001);
    FAIL_IF_NOT_WITHIN(chebyshev(refcounts, qrycounts), 0.01042, 0.0001);

    FAIL_IF_NOT_WITHIN(d2s(refcounts, qrycounts, refmarkov, qrymarkov), 0.3808187760303444, 0.0001);
    FAIL_IF_NOT_WITHIN(d2star(refcounts, qrycounts, refmarkov, qrymarkov), 0.4333069320392635, 0.0001);
}

void tests_prep_dna_7() {
    std::vector<unsigned> qrycounts, refcounts;
    std::vector<double> qrymarkov, refmarkov;
    prep(qrycounts, refcounts, qrymarkov, refmarkov, 7, 0);

    FAIL_IF_NOT_WITHIN(d2(refcounts, qrycounts), 0.47872, 0.0001);
    FAIL_IF_NOT_WITHIN(euler(refcounts, qrycounts), 0.10092, 0.0001);
    //FAIL_IF_NOT_WITHIN(manhattan(refcounts, qrycounts), ???, 0.0001);
    FAIL_IF_NOT_WITHIN(bray_curtis_distance(refcounts, qrycounts), 0.957446808511, 0.0001);
    FAIL_IF_NOT_WITHIN(normalised_google_distance(refcounts, qrycounts), 0.957446808511, 0.0001);
    FAIL_IF_NOT_WITHIN(chebyshev(refcounts, qrycounts), 0.00532, 0.0001);

    FAIL_IF_NOT_WITHIN(d2s(refcounts, qrycounts, refmarkov, qrymarkov), 0.3111486571852624, 0.0001);
    FAIL_IF_NOT_WITHIN(d2star(refcounts, qrycounts, refmarkov, qrymarkov), 0.4811460538701716, 0.0001);
}

void tests_prep_dna_9() {
    std::vector<unsigned> qrycounts, refcounts;
    std::vector<double> qrymarkov, refmarkov;
    prep(qrycounts, refcounts, qrymarkov, refmarkov, 9, 0);

    FAIL_IF_NOT_WITHIN(d2(refcounts, qrycounts), 0.49457, 0.0001);
    FAIL_IF_NOT_WITHIN(euler(refcounts, qrycounts), 0.10369, 0.0001);
    //FAIL_IF_NOT_WITHIN(manhattan(refcounts, qrycounts), ???, 0.0001);
    FAIL_IF_NOT_WITHIN(bray_curtis_distance(refcounts, qrycounts), 0.989130434783, 0.0001);
    FAIL_IF_NOT_WITHIN(normalised_google_distance(refcounts, qrycounts), 0.989130434783, 0.0001);
    FAIL_IF_NOT_WITHIN(chebyshev(refcounts, qrycounts), 0.00543, 0.0001);

    FAIL_IF_NOT_WITHIN(d2s(refcounts, qrycounts, refmarkov, qrymarkov), 0.3105504265761406, 0.0001);
    FAIL_IF_NOT_WITHIN(d2star(refcounts, qrycounts, refmarkov, qrymarkov), 0.4938153388000316, 0.0001);
}

/////////////////////////////////////////////////////////

void tests_prep_aa_3() {
    std::vector<unsigned> qrycounts, refcounts;
    std::vector<double> qrymarkov, refmarkov;
    prep_aa(qrycounts, refcounts, qrymarkov, refmarkov, 3, 0);

    FAIL_IF_NOT_WITHIN(d2(refcounts, qrycounts), 0.47198, 0.0001);
    FAIL_IF_NOT_WITHIN(euler(refcounts, qrycounts), 0.11299, 0.0001);
    FAIL_IF_NOT_WITHIN(manhattan(refcounts, qrycounts), 1.91489, 0.0001);
    FAIL_IF_NOT_WITHIN(bray_curtis_distance(refcounts, qrycounts), 0.943343, 0.0001);
    FAIL_IF_NOT_WITHIN(normalised_google_distance(refcounts, qrycounts), 0.957447, 0.0001);
    FAIL_IF_NOT_WITHIN(chebyshev(refcounts, qrycounts), 0.0169492, 0.0001);

    FAIL_IF_NOT_WITHIN(d2s(refcounts, qrycounts, refmarkov, qrymarkov), 0.377069, 0.0001);
    FAIL_IF_NOT_WITHIN(d2star(refcounts, qrycounts, refmarkov, qrymarkov), 0.497807, 0.0001);
}

void tests_prep_aa_4() {
    std::vector<unsigned> qrycounts, refcounts;
    std::vector<double> qrymarkov, refmarkov;
    prep_aa(qrycounts, refcounts, qrymarkov, refmarkov, 4, 0);

    //FAIL_IF_NOT_WITHIN(d2(refcounts, qrycounts), ???, 0.0001);
    //FAIL_IF_NOT_WITHIN(euler(refcounts, qrycounts), ???, 0.0001);
    //FAIL_IF_NOT_WITHIN(manhattan(refcounts, qrycounts), ???, 0.0001);
    //FAIL_IF_NOT_WITHIN(bray_curtis_distance(refcounts, qrycounts), ???, 0.0001);
    //FAIL_IF_NOT_WITHIN(normalised_google_distance(refcounts, qrycounts), ???, 0.0001);
    //FAIL_IF_NOT_WITHIN(chebyshev(refcounts, qrycounts), ???, 0.0001);

    FAIL_IF_NOT_WITHIN(d2s(refcounts, qrycounts, refmarkov, qrymarkov), 0.379908, 0.0001);
    FAIL_IF_NOT_WITHIN(d2star(refcounts, qrycounts, refmarkov, qrymarkov), 0.500178, 0.0001);
}

void tests_prep_aa_5() {
    std::vector<unsigned> qrycounts, refcounts;
    std::vector<double> qrymarkov, refmarkov;
    prep_aa(qrycounts, refcounts, qrymarkov, refmarkov, 5, 0);

    FAIL_IF_NOT_WITHIN(d2(refcounts, qrycounts), 0.5, 0.0001);
    FAIL_IF_NOT_WITHIN(euler(refcounts, qrycounts), 0.113633, 0.0001);
    FAIL_IF_NOT_WITHIN(manhattan(refcounts, qrycounts), 2, 0.0001);
    FAIL_IF_NOT_WITHIN(bray_curtis_distance(refcounts, qrycounts), 1, 0.0001);
    FAIL_IF_NOT_WITHIN(normalised_google_distance(refcounts, qrycounts), 1, 0.0001);
    FAIL_IF_NOT_WITHIN(chebyshev(refcounts, qrycounts), 0.00862069, 0.0001);

    FAIL_IF_NOT_WITHIN(d2s(refcounts, qrycounts, refmarkov, qrymarkov), 0.391587, 0.0001);
    FAIL_IF_NOT_WITHIN(d2star(refcounts, qrycounts, refmarkov, qrymarkov), 0.500023, 0.0001);
}

void tests_prep_aa_6() {
    std::vector<unsigned> qrycounts, refcounts;
    prep_aa(qrycounts, refcounts, 6);

    FAIL_IF_NOT_WITHIN(d2(refcounts, qrycounts), 0.5, 0.0001);
    //FAIL_IF_NOT_WITHIN(euler(refcounts, qrycounts), ???, 0.0001);
    //FAIL_IF_NOT_WITHIN(manhattan(refcounts, qrycounts), ???, 0.0001);
    //FAIL_IF_NOT_WITHIN(bray_curtis_distance(refcounts, qrycounts), ???, 0.0001);
    //FAIL_IF_NOT_WITHIN(normalised_google_distance(refcounts, qrycounts), ???, 0.0001);
    //FAIL_IF_NOT_WITHIN(chebyshev(refcounts, qrycounts), ???, 0.0001);
}

/////////////////////////////////////////////////////////

void tests_prep_raa_3() {
    std::vector<unsigned> qrycounts, refcounts;
    std::vector<double> qrymarkov, refmarkov;
    prep_raa(qrycounts, refcounts, qrymarkov, refmarkov, 3, 0);

    FAIL_IF_NOT_WITHIN(d2(refcounts, qrycounts), 0.301479, 0.0001);
    FAIL_IF_NOT_WITHIN(euler(refcounts, qrycounts), 0.114879, 0.0001);
    FAIL_IF_NOT_WITHIN(manhattan(refcounts, qrycounts), 1.42243, 0.0001);
    FAIL_IF_NOT_WITHIN(bray_curtis_distance(refcounts, qrycounts), 0.72238, 0.0001);
    FAIL_IF_NOT_WITHIN(normalised_google_distance(refcounts, qrycounts), 0.791489, 0.0001);
    FAIL_IF_NOT_WITHIN(chebyshev(refcounts, qrycounts), 0.0254237, 0.0001);

    FAIL_IF_NOT_WITHIN(d2s(refcounts, qrycounts, refmarkov, qrymarkov), 0.449371, 0.0001);
    FAIL_IF_NOT_WITHIN(d2star(refcounts, qrycounts, refmarkov, qrymarkov), 0.49648, 0.0001);
}

void tests_prep_raa_4() {
    std::vector<unsigned> qrycounts, refcounts;
    std::vector<double> qrymarkov, refmarkov;
    prep_raa(qrycounts, refcounts, qrymarkov, refmarkov, 4, 0);

    FAIL_IF_NOT_WITHIN(d2(refcounts, qrycounts), 0.445397, 0.0001);
    FAIL_IF_NOT_WITHIN(euler(refcounts, qrycounts), 0.112743, 0.0001);
    FAIL_IF_NOT_WITHIN(manhattan(refcounts, qrycounts), 1.83761, 0.0001);
    FAIL_IF_NOT_WITHIN(bray_curtis_distance(refcounts, qrycounts), 0.903134, 0.0001);
    FAIL_IF_NOT_WITHIN(normalised_google_distance(refcounts, qrycounts), 0.92735, 0.0001);
    FAIL_IF_NOT_WITHIN(chebyshev(refcounts, qrycounts), 0.017094, 0.0001);

    FAIL_IF_NOT_WITHIN(d2s(refcounts, qrycounts, refmarkov, qrymarkov), 0.391654, 0.0001);
    FAIL_IF_NOT_WITHIN(d2star(refcounts, qrycounts, refmarkov, qrymarkov), 0.502855, 0.0001);
}

void tests_prep_raa_5() {
    std::vector<unsigned> qrycounts, refcounts;
    std::vector<double> qrymarkov, refmarkov;
    prep_raa(qrycounts, refcounts, qrymarkov, refmarkov, 5, 0);

    FAIL_IF_NOT_WITHIN(d2(refcounts, qrycounts), 0.484987, 0.0001);
    FAIL_IF_NOT_WITHIN(euler(refcounts, qrycounts), 0.112819, 0.0001);
    FAIL_IF_NOT_WITHIN(manhattan(refcounts, qrycounts), 1.95708, 0.0001);
    FAIL_IF_NOT_WITHIN(bray_curtis_distance(refcounts, qrycounts), 0.971347, 0.0001);
    FAIL_IF_NOT_WITHIN(normalised_google_distance(refcounts, qrycounts), 0.978541, 0.0001);
    FAIL_IF_NOT_WITHIN(chebyshev(refcounts, qrycounts), 0.0172414, 0.0001);

    FAIL_IF_NOT_WITHIN(d2s(refcounts, qrycounts, refmarkov, qrymarkov), 0.367139, 0.0001);
    FAIL_IF_NOT_WITHIN(d2star(refcounts, qrycounts, refmarkov, qrymarkov), 0.499936, 0.0001);
}

void tests_prep_raa_6() {
    std::vector<unsigned> qrycounts, refcounts;
    prep_raa(qrycounts, refcounts, 6);

    FAIL_IF_NOT_WITHIN(d2(refcounts, qrycounts), 0.5, 0.0001);
    FAIL_IF_NOT_WITHIN(euler(refcounts, qrycounts), 0.114044, 0.0001);
    FAIL_IF_NOT_WITHIN(manhattan(refcounts, qrycounts), 2.0, 0.0001);
    FAIL_IF_NOT_WITHIN(bray_curtis_distance(refcounts, qrycounts), 1, 0.0001);
    FAIL_IF_NOT_WITHIN(normalised_google_distance(refcounts, qrycounts), 1, 0.0001);
    FAIL_IF_NOT_WITHIN(chebyshev(refcounts, qrycounts), 0.00869565, 0.0001);
}
