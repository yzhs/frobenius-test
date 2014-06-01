#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

#include <gmp.h>

#include "helpers.h"
#include "frobenius.h"

int main()
{
       mpz_t tmp, tmp0, tmp1, tmp2;

       mpz_inits(tmp, tmp0, tmp1, tmp2, NULL);
       init();

       mpz_set_ui(tmp, 1215239);
       assert(RQFT(tmp, 1) != composite);
       mpz_set_ui(tmp, 1215237);
       assert(RQFT(tmp, 1) == composite);
       mpz_set_str(tmp, "2147483659", 10);
       assert(RQFT(tmp, 1) == probably_prime);
       mpz_set_str(tmp, "32317006071311007300714876688669951960444102669715484032130"\
                       "34542752465513886789089319720141152291346368871796092189801949411955"\
                       "91504909210950881523864482831206308773673009960917501977503896521067"\
                       "96057638384067568276792218642619756161838094338476170470581645852036"\
                       "30504288757589154106580860755239912393038552191433338966834242068497"\
                       "47865645694948561760353263220580778056593310261927084603141502585928"\
                       "64177116725943603718461857357598351152301645904403697613233287231227"\
                       "12568471082020972515710172693132346967854258065669793504599726835299"\
                       "86382155251663894373355436021354332296046453184786049521481935558536"\
                       "11059596231637", 10);
       assert(RQFT(tmp, 1) == probably_prime);

       cleanup();
       mpz_clears(tmp, tmp0, tmp1, tmp2, NULL);
}
