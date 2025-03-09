#include <iostream>

#include "distance_tests.h"
//#include "utils_tests.h"

int main() {
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
    return 0;
}

/*
SEQAN_BEGIN_TESTSUITE(KAST_tests)
{
   // Call roundUp() tests
   SEQAN_CALL_TEST(d2_dna);
   SEQAN_CALL_TEST(euler_dna);
   SEQAN_CALL_TEST(manhattan_dna);
   SEQAN_CALL_TEST(bc_dna);
   SEQAN_CALL_TEST(ngd_dna);
   SEQAN_CALL_TEST(chebyshev_dna);
   SEQAN_CALL_TEST(d2s_dna);
   SEQAN_CALL_TEST(d2star_dna);

   SEQAN_CALL_TEST(d2_aa);
   SEQAN_CALL_TEST(euler_aa);
   SEQAN_CALL_TEST(manhattan_aa);
   SEQAN_CALL_TEST(bc_aa);
   SEQAN_CALL_TEST(ngd_aa);
   SEQAN_CALL_TEST(chebyshev_aa);
   SEQAN_CALL_TEST(d2s_aa);
   SEQAN_CALL_TEST(d2star_aa);

   SEQAN_CALL_TEST(d2_raa);
   SEQAN_CALL_TEST(euler_raa);
   SEQAN_CALL_TEST(manhattan_raa);
   SEQAN_CALL_TEST(bc_raa);
   SEQAN_CALL_TEST(ngd_raa);
   SEQAN_CALL_TEST(chebyshev_raa);
   SEQAN_CALL_TEST(d2s_raa);
   SEQAN_CALL_TEST(d2star_raa);

   //SEQAN_CALL_TEST(count_dna);

   SEQAN_CALL_TEST(mask_count_dna);
   SEQAN_CALL_TEST(mask_count_raa);
   SEQAN_CALL_TEST(mask_count_aa);

   //SEQAN_CALL_TEST(rounding_test_2);
   //SEQAN_CALL_TEST(rounding_test_3);

   // Call checkSorted() tests
   //SEQAN_CALL_TEST(check_sorted_test_1);
   //SEQAN_CALL_TEST(check_sorted_test_2);
}
SEQAN_END_TESTSUITE

*/
