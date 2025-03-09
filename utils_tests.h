#ifndef __KAST_UTILS_TESTS_H__
#define __KAST_UTILS_TESTS_H__

#include <algorithm>
#include <seqan3/alphabet/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include "utils.h"

using namespace seqan3::literals;

#define FAIL_IF_NOT_EQUAL(expected, found) { \
    auto found_val = found; \
    auto expected_val = expected; \
    if (found_val != expected_val) { \
        seqan3::debug_stream << "Test " << __func__ << "/" << #found << ": expected " << expected_val << " found " << found_val << std::endl; \
        exit(1); \
    } else { \
        seqan3::debug_stream << "Test " << __func__ << "/" << #found << " OK" << std::endl; \
    } \
}

void dna5_masked_single_bit() {
    auto seq = "AGGCAGCGTACGAACCTACTGGAGTTGCGGTATGGGACCAGGCGACCTCTGATGCAGAGATACAGGAGCGCCGCGCCGGGTCTTCCTTGTAGAAGTCCTG"_dna5;
    std::vector<unsigned> counts;
    std::vector<std::string> masks { "10000" };

    count_kmers(counts, seq, 5, 1, masks);

    FAIL_IF_NOT_EQUAL(4, counts.size());
    FAIL_IF_NOT_EQUAL(21, counts[rank_of_char_as<seqan3::dna4>('A')]);
    FAIL_IF_NOT_EQUAL(24, counts[rank_of_char_as<seqan3::dna4>('C')]);
    FAIL_IF_NOT_EQUAL(33, counts[rank_of_char_as<seqan3::dna4>('G')]);
    FAIL_IF_NOT_EQUAL(18, counts[rank_of_char_as<seqan3::dna4>('T')]);
}

void aa_masked_single_bit() {
    auto seq = "MVLTIYPDELVQIVS"_aa27;
    std::vector<unsigned> counts;
    std::vector<std::string> masks { "10000" };

    count_kmers(counts, seq, 5, 1, masks);

    FAIL_IF_NOT_EQUAL(27, counts.size());
    FAIL_IF_NOT_EQUAL(1, counts[rank_of_char_as<seqan3::aa27>('M')]);
    FAIL_IF_NOT_EQUAL(2, counts[rank_of_char_as<seqan3::aa27>('V')]);
    FAIL_IF_NOT_EQUAL(2, counts[rank_of_char_as<seqan3::aa27>('L')]);
    FAIL_IF_NOT_EQUAL(1, counts[rank_of_char_as<seqan3::aa27>('T')]);
    FAIL_IF_NOT_EQUAL(1, counts[rank_of_char_as<seqan3::aa27>('I')]);
    FAIL_IF_NOT_EQUAL(1, counts[rank_of_char_as<seqan3::aa27>('Y')]);
    FAIL_IF_NOT_EQUAL(1, counts[rank_of_char_as<seqan3::aa27>('P')]);
    FAIL_IF_NOT_EQUAL(1, counts[rank_of_char_as<seqan3::aa27>('D')]);
    FAIL_IF_NOT_EQUAL(1, counts[rank_of_char_as<seqan3::aa27>('E')]);
    FAIL_IF_NOT_EQUAL(27 - 9, std::count(counts.begin(), counts.end(), 0));
}

void raa_masked_single_bit() {
    auto seq = "MVLTIYPDELVQIVS"_aa10murphy;
    std::vector<unsigned> counts;
    std::vector<std::string> masks { "10000" };

    count_kmers(counts, seq, 5, 1, masks);

    FAIL_IF_NOT_EQUAL(10, counts.size());
    FAIL_IF_NOT_EQUAL(0, counts[rank_of_char_as<seqan3::aa10murphy>('A')]);
    FAIL_IF_NOT_EQUAL(2, counts[rank_of_char_as<seqan3::aa10murphy>('B')]);
    FAIL_IF_NOT_EQUAL(0, counts[rank_of_char_as<seqan3::aa10murphy>('C')]);
    FAIL_IF_NOT_EQUAL(1, counts[rank_of_char_as<seqan3::aa10murphy>('F')]);
    FAIL_IF_NOT_EQUAL(0, counts[rank_of_char_as<seqan3::aa10murphy>('G')]);
    FAIL_IF_NOT_EQUAL(0, counts[rank_of_char_as<seqan3::aa10murphy>('H')]);
    FAIL_IF_NOT_EQUAL(6, counts[rank_of_char_as<seqan3::aa10murphy>('I')]);
    FAIL_IF_NOT_EQUAL(0, counts[rank_of_char_as<seqan3::aa10murphy>('K')]);
    FAIL_IF_NOT_EQUAL(1, counts[rank_of_char_as<seqan3::aa10murphy>('P')]);
    FAIL_IF_NOT_EQUAL(1, counts[rank_of_char_as<seqan3::aa10murphy>('S')]);
}

#endif // __KAST_UTILS_TESTS_H__
