#include <CUnit/CUnit.h>
#include <CUnit/Basic.h>

#include "../helpers_int.h"
#include "test_frobenius_int.h"
#include "test_miller_rabin_int.h"

#include "../helpers.h"
#include "test_miller_rabin_long.h"

int main()
{
	CU_pSuite suite;

	if (CUE_NOMEMORY == CU_initialize_registry())
		die("could not initialize CUnit registry: out of memory\n");

	init_int();

	suite = CU_add_suite("Miller-Rabin (int)", NULL, NULL);
	CU_ADD_TEST(suite, test_miller_rabin_powm_int);
	//CU_ADD_TEST(suite, test_miller_rabin_primes_int);

	suite = CU_add_suite("Miller-Rabin (GMP)", init, cleanup);
	CU_ADD_TEST(suite, test_miller_rabin_some_numbers);
	CU_ADD_TEST(suite, test_miller_rabin_primes);
	CU_ADD_TEST(suite, test_miller_rabin_composites);
	CU_ADD_TEST(suite, test_miller_rabin_composites2);
	CU_ADD_TEST(suite, test_miller_rabin_both);

	suite = CU_add_suite("Frobenius (int)", NULL, NULL);
	CU_ADD_TEST(suite, test_frobenius_int_int_sqrt);
	CU_ADD_TEST(suite, test_frobenius_int_gcd);
	CU_ADD_TEST(suite, test_frobenius_int_is_square);
	CU_ADD_TEST(suite, test_frobenius_int_jacobi);
	CU_ADD_TEST(suite, test_frobenius_int_split);
	CU_ADD_TEST(suite, test_frobenius_mult_mod_int);
	CU_ADD_TEST(suite, test_frobenius_powm_mod_int);
	CU_ADD_TEST(suite, test_frobenius_int_get_random_int);
	CU_ADD_TEST(suite, test_frobenius_squares_int);
	CU_ADD_TEST(suite, test_frobenius_trial_division_int);
	CU_ADD_TEST(suite, test_frobenius_int_rqft_small_primes);
	CU_ADD_TEST(suite, test_frobenius_int_rqft_small_composites);
	CU_ADD_TEST(suite, test_frobenius_primelist_int);
	CU_ADD_TEST(suite, test_frobenius_problematic_primes_int);
	CU_ADD_TEST(suite, test_frobenius_larger_primes_int);

	CU_basic_set_mode(CU_BRM_VERBOSE);
	CU_basic_run_tests();

	CU_cleanup_registry();
}
