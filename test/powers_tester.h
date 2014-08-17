unsigned long n_;
mpz_t MODULUS, POLY(f), POLY(foo), POLY(bar), POLY(x), baz;
mpz_inits(MODULUS, POLY(foo), POLY(f), POLY(bar), POLY(x), baz, NULL);

// Initialize POLY(x)
n_ = 0x7ffffffff;
mpz_set_ui(n, n_);
mpz_urandomm(n, r_state, n);
mpz_add(n, n, n);
mpz_add_ui(n, n, 0x70000001);
n_ = mpz_get_ui(n);
mpz_set_ui(x_x, 1);

for (unsigned long c_ = 1; c_ < 10; c_++) {
	if (jacobi(n_ - c_, n_) != 1)
		continue;
	mpz_set_ui(c, c_);

	for (unsigned long b_ = 1; b_ < 10; b_++) {
		if (jacobi(b_*b_+4*c_, n_) != -1)
			continue;
		mpz_set_ui(b, b_);

		mpz_urandomm(f_x, r_state, n);
		mpz_urandomm(f_1, r_state, n);

		RUN_TEST
	}
}

mpz_clears(MODULUS, POLY(f), POLY(foo), POLY(bar), POLY(x), baz, NULL);

#undef RUN_TEST
