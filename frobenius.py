import array
from math import ceil, floor, log, sqrt
import operator
from random import randint

from miller_rabin import split
from polynomial import *

# {{{ helper functions
def sieve(n):
	prime_list = []
	ar = array.array('b', [0,1]*((n+2) // 2))
	ar[1] = 0
	ar[2] = 1
	l = ceil(sqrt(n))
	for i in range(3, l, 2):
		if ar[i] == 0:
			continue
		for j in range(i*i, n+1, i):
			ar[j] = 0
	for i in range(2, n+1):
		if ar[i] == 1:
			prime_list.append(i)
	return prime_list

prime_list = sieve(50000)

def primes(n):
	for p in prime_list:
		if p > n:
			break
		yield p

def even(x): return x % 2 == 0
def odd(x): return x % 2 == 1

def gcd(a, b):
	while b != 0:
		(a, b) = (b, a % b)
	return a

def polygcd(a, b):
	while b.degree() == 0:
		(a, b) = (b, a % b)
	return a

def jacobi(x, y):
	if y < 3 or even(y):
		print("jacobi(x,y): y must be an odd integer >= 3, found {}".format(y))
		return 0
	res = 1
	while True:
		x = x % y
		if x == 0:
			return 0
		while even(x):
			x //= 2
			m8 = y % 8
			if m8 == 3 or m8 == 5:
				res = -res
		if x == 1:
			return res
		(x, y) = (y, x)
		if x % 4 == 3 and y % 4 == 3:
			res = -res

def pow_mod(b, e, m):
	if isinstance(m, tuple) or isinstance(m, list):
		return pow_mod(b, e, m[-1]) % m[0]
	if isinstance(b, Polynomial):
		res = poly(b.n, (1,))
	else:
		res = 1
	while e != 0:
		if odd(e):
			res = (res * b) % m
		b = (b * b) % m
		e //= 2
	return res
# }}}


def debug(s):
	if False:
		print(s)

def isqrt(n):
	try:
		r = floor(sqrt(n))
	except OverflowError:
		r = 2 * 10**floor(log(n)/log(100))
	if r*r > n:
		r = n // r
	while r < n//r-1:
		r = (r + n//r)//2
	return r

def is_square(n):
	r = isqrt(n)
	return r*r == n

# {{{ QFT
# The Quadratic Frobenius Test (QFT) with parameters (b,c) consists of the
# following.
def QFT(n, b, c, B):
	# Suppose n>1 is odd, (b^2+4c over n)=-1, and (-c over n)=1.
	#assert(n > 1), "{} <= 1".format(n)
	#assert(jacobi(b**2+4*c, n) == -1), "b={}, c={}, n={}".format(b, c, n)
	#assert(jacobi(-c, n) == 1), "c={}, n={}".format(c, n)

	# (1) Test n for divisibility by primes less than or equal to min{B,
	# sqrt(n)}.  If it is divisible by one of these primes, declare n to be
	# composite and stop.
#	for p in primes(min(B, isqrt(n))):
#		if n % p == 0:
#			#debug("found a prime divisor p={} of n={}".format(p, n))
#			return False # composite

	# (2) Test whether sqrt(n) in ℤ.  If it is, declare n to be composite and
	# stop.
#	if floor(sqrt(n))**2 == n:
#		debug("n={}={}^2 is a perfect square and therefore composite".format(n, floor(sqrt(n))))
#		return False # composite

	x = get_x(n)
	m = x**2-x*b-c
	#debug("doing computations modulo ({}, {})".format(n, m))

	# (3) Compute x^((n+1)/2) mod (n, x^2-bx-c).  If x^((n+1)/2) not in ℤ/nℤ,
	# declare n to be composite and stop.
	foo = pow_mod(x, (n + 1) // 2, m)
	if foo.degree() > 0:
		debug("n={} is composite by the criterion checked in step (3): {}".format(n, foo))
		return False # composite

	# (4) Compute x^(n+1) mod (n, x^2-bx-c).  If x^(n+1) not congruent -c,
	# declare n to be composite and stop.
	foo = pow_mod(x, n + 1, m)
	if foo != poly(n, (-c,)):
		debug("n={} is composite by the criterion checked in step (4)".format(n))
		return False # composite

	# (5) Let n^2-1=2^r*s, where s is odd.  If x^s not congruent 1 mod (n,
	# x^2-bx-c), and x^(2^j*s) not congruent -1 mod (n, x^2-bx-c) for all
	# 0≤j≤r-2, declare n to be composite and stop.
	# If n is not declared composite in Steps 1—5, declare n to be a probable
	# prime.
	(r, s) = split(n**2)
	assert(2**r*s == n**2-1), "2^{}*{} == {}*{} != {} = {}^2-1}".format(r, s, 2**r, s, n**2-1, n)
	foo = pow_mod(x, s, m)
	if foo == poly(n, (1,)):
		#debug("as asserted by step (5a), n={} is probably prime".format(n))
		return True # probably prime
	else:
		pass#debug("x^{} == {}".format(s, foo))
	for j in range(0, r-2 + 1):
		bar = pow_mod(x, 2**j * s, m)
		if bar == poly(n, (-1,)):
			#debug("as asserted by step (5b), n={} is probably prime".format(n))
			return True # probably prime
		else:
			pass#debug("x^(2^{}*{}) == {}".format(j, s, foo))
	return False # composite

def RQFT(n, B):
	assert(n > 1)
	assert(odd(n))

	if n in prime_list:
		return True
	fail = False
	b = c = 0
	j1 = j2 = 0

	for p in primes(min(B, isqrt(n))):
		if n % p == 0:
			#debug("found a prime divisor p={} of n={}".format(p, n))
			return False # composite

	# (2) Test whether sqrt(n) in ℤ.  If it is, declare n to be composite and
	# stop.
	if is_square(n):
		debug("n={}={}^2 is a perfect square and therefore composite".format(n, floor(sqrt(n))))
		return False # composite

	for _ in range(B):
		fail = False
		# (1) Choose pairs (b,c) at random with 1≤b,c≤n until one is found with
		# (b^2+4c over n)=-1 and (-c over n)=1, or with gcd(b^2+4c,n), gcd(b,n) or
		# gcd(c,n) a nontrivial divisor of n.  However, if the latter case occurs
		# before the former, declar n to be composite and stop.  If after B pairs
		# are tested, none is found satisfying the above conditions, declare n to be
		# a probable prime and stop.
		b = randint(1, n)
		c = randint(1, n)
		j1 = jacobi(b*b+4*c, n)
		j2 = jacobi(-c, n)
		if j1 == -1 and j2 == 1:
			if gcd(b**2+4*c, n) not in [1, n] \
			or gcd(b, n) not in [1, n] \
			or gcd(c, n) not in [1, n]:
				return False
			break

	if jacobi(b*b+4*c, n) != -1 or jacobi(-c, n) != 1:
		print("n = {}, b = {}, c = {}, (b^2+4c/n) = {}, (-c/n) = {}".format(n, b, c, jacobi(b**2+4*c, n), jacobi(-c, n)))
		return None


	# (2) Perform the QFT with parameters (b,c).
	return QFT(n, b, c, B)
# }}}

B = 50000

class TestFrobenius(unittest.TestCase):
	def setUp(self):
		self.n = 1000
		self.large_primes = []
		self.composites = []
		primes = self.large_primes
		with open("primelist.txt") as f:
			data = f.read().split("\n")[:-1]
			for p in data:
				primes.append(int(p))
		with open("composites.txt") as f:
			data = f.read().split("\n")[:-1]
			for c in data:
				self.composites.append(int(c))

	def test_primes(self):
		for i in range(self.n):
			n = randint(0, len(self.large_primes))
			self.assertTrue(RQFT(self.large_primes[n], B) != False)

	def test_squarefree_composites(self):
		for n in self.composites:
			self.assertFalse(RQFT(n, B))

	def test_composites(self):
		counter = errors = 0
		for _ in range(self.n):
			i = randint(0, len(self.large_primes)-1)
			if self.large_primes[i] + 2 == self.large_primes[i+1]:
				i+=1
			n = randint(self.large_primes[i]+1, self.large_primes[i+1]-1)
			while n % 2 == 0:
				n = randint(self.large_primes[i]+1, self.large_primes[i+1]-1)
			# [n] is composite by construction
			self.assertFalse(RQFT(n, B))

	def test_some_stuff(self):
		self.assertTrue(RQFT(7, B))
		self.assertFalse(RQFT(611879**2*611957**4, B))
		self.assertFalse(RQFT(1235790412356789098765432827498274392743929834792843282734279348239482349**9, B))

if __name__ == "__main__":
	unittest.main()

# vim: set ts=4 sw=4 noet :
