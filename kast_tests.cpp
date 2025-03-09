#include "distance_tests.h"
#include "utils_tests.h"

int main() {
    split_by_space_tests();

    dna5_masked_single_bit();
    aa_masked_single_bit();
    raa_masked_single_bit();

    tests_prep_dna_3();
    tests_prep_dna_5();
    tests_prep_dna_7();
    tests_prep_dna_9();
    tests_prep_aa_3();
    tests_prep_aa_4();
    tests_prep_aa_5();
    tests_prep_aa_6();
    tests_prep_raa_3();
    tests_prep_raa_4();
    tests_prep_raa_5();
    tests_prep_raa_6();

    seqan3::debug_stream << "All tests OK" << std::endl;
    return 0;
}
