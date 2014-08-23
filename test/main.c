#include <CUnit/CUnit.h>
#include <CUnit/Basic.h>

#include "../helpers_int.h"
#include "test_frobenius_int.h"
#include "test_miller_rabin_int.h"

#include "../helpers.h"
#include "test_frobenius_long.h"
#include "test_miller_rabin_long.h"

int main()
{
	CU_pSuite suite;

	if (CUE_NOMEMORY == CU_initialize_registry())
		die("could not initialize CUnit registry: out of memory\n");

//	suite = CU_add_suite("Miller-Rabin (int)", init_int, NULL);
//	CU_ADD_TEST(suite, test_mr_powm_int);
//	CU_ADD_TEST(suite, test_mr_primes_int);

//	suite = CU_add_suite("Miller-Rabin (GMP)", init, cleanup);
//	CU_ADD_TEST(suite, test_mr_some_numbers);
//	CU_ADD_TEST(suite, test_mr_primes);
//	CU_ADD_TEST(suite, test_mr_composites);
//	CU_ADD_TEST(suite, test_mr_composites2);
//	CU_ADD_TEST(suite, test_mr_both);

//	suite = CU_add_suite("Frobenius (int)", init_int, NULL);
//	CU_ADD_TEST(suite, frob_int_int_sqrt);
//	CU_ADD_TEST(suite, frob_int_gcd);
//	CU_ADD_TEST(suite, frob_int_is_square);
//	CU_ADD_TEST(suite, frob_int_jacobi);
//	CU_ADD_TEST(suite, frob_int_split);
//	CU_ADD_TEST(suite, frob_mult_mod_int);
//	CU_ADD_TEST(suite, frob_powm_mod_int);
//	CU_ADD_TEST(suite, frob_int_get_random_int);
//	CU_ADD_TEST(suite, frob_squares_int);
//	CU_ADD_TEST(suite, frob_trial_division_int);
//	CU_ADD_TEST(suite, frob_int_rqft_small_primes);
//	CU_ADD_TEST(suite, frob_int_rqft_small_composites);
//	CU_ADD_TEST(suite, frob_primelist_int);
//	CU_ADD_TEST(suite, frob_problematic_primes_int);
//	CU_ADD_TEST(suite, frob_larger_primes_int);

	suite = CU_add_suite("Frobenius (GMP)", init, cleanup);
//	CU_ADD_TEST(suite, frob_mult_x);
//	CU_ADD_TEST(suite, frob_sigma_basics);
//	CU_ADD_TEST(suite, frob_sigma_short_integer);
//	CU_ADD_TEST(suite, frob_sigma_power);
	CU_ADD_TEST(suite, frob_power_basics);
//	CU_ADD_TEST(suite, frob_power_x_lucas);
	CU_ADD_TEST(suite, frob_inverse);
	CU_ADD_TEST(suite, frob_fast_algorithm1);
	CU_ADD_TEST(suite, frob_fast_algorithm2);
	CU_ADD_TEST(suite, frob_fast_algorithm3);
	CU_ADD_TEST(suite, frob_composites);
	CU_ADD_TEST(suite, frob_split);
	CU_ADD_TEST(suite, frob_mult_mod);
	CU_ADD_TEST(suite, frob_powm_mod);
	CU_ADD_TEST(suite, frob_squares);
	CU_ADD_TEST(suite, frob_trial_division);
	CU_ADD_TEST(suite, frob_rqft_small_primes);
	CU_ADD_TEST(suite, frob_larger_primes);
//	CU_ADD_TEST(suite, frob_primelist);

	CU_basic_set_mode(CU_BRM_VERBOSE);
	CU_basic_run_tests();

	CU_cleanup_registry();
}
