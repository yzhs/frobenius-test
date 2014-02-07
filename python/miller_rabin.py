import unittest
import sys
if sys.version > '3':
	long = int

from random import randint

def split(n):
	d = n - 1
	s = 0
	while d % 2 == 0:
		s += 1
		d //= 2
	return (s, d)

# This function checks whether a given number n is a prime or not, using the
# Miller-Rabin primality test.  This is a probabilistic test which randomly
# chooses an integer a as a base and checks whether n satisfies a certain
# property (which depends on b).  If it does, n is a prime for at least three
# out of four of the possible values of a, if it does not, it is certainly not
# prime.
# The implementation is taken from the following pseudo code found on
# http://en.wikipedia.org/wiki/Miller-Rabin_primality_test:
# 
#   Input: n > 3, an odd integer to be tested for primality;
#   Input: k, a parameter that determines the accuracy of the test
#   Output: composite if n is composite, otherwise probably prime
#   write n - 1 as 2^s*d with d odd by factoring powers of 2 from n - 1
#   LOOP: repeat k times:
#      pick a random integer a in the range 2, n - 2
#      x <- a^d mod n
#      if x = 1 or x = n - 1 then do next LOOP
#      for r = 1 .. s
#         x <- x^2 mod n
#         if x = 1 then return composite
#         if x = n - 1 then do next LOOP
#      return composite
#   return probably prime
# 
# The function returns true if it found no evidence, that n might be composite
# and false if it found a counter example.

def miller_rabin(n, k):
	(s, d) = split(n)
	nm1 = n - 1
	for _ in range(k):
		a = randint(2, n-2)
		x = pow(a, d, n)
		if x == 1 or x == nm1:
			continue
		for r in range(1, s+1):
			x = pow(x, 2, n)
			if x == 1:
				return False
			if x == nm1:
				break
		if x != nm1:
			return False
	return True

class TestMillerRabin(unittest.TestCase):
	def setUp(self):
		self.n = 1000
		self.large_primes = []
		primes = self.large_primes
		with open("primelist.txt") as f:
			data = f.read().split("\n")[:-1]
			for p in data:
				primes.append(int(p))

	def test_primes(self):
		for i in range(self.n):
			n = randint(0, len(self.large_primes))
			self.assertTrue(miller_rabin(self.large_primes[n], 1))

	def test_composites(self):
		counter = 0
		for i in range(self.n):
			n = 2*(randint(self.large_primes[1]+1, self.large_primes[2]-1)//2) + 1
			if miller_rabin(n, 1):
				counter += 1
		self.assertLess(counter, self.n//4)

	def test_some_stuff(self):
		self.assertTrue(miller_rabin(7, 1))
		self.assertFalse(miller_rabin(611879**2*611957**4, 1))
		self.assertFalse(miller_rabin(1235790412356789098765432827498274392743929834792843282734279348239482349**9, 1))

if __name__ == "__main__":
	unittest.main()

# vim: set ts=4 sw=4 noet :
