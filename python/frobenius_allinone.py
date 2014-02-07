import array
from math import ceil, floor, log, sqrt
import operator
from random import randint, seed
import operator
import unittest
from small_primes import prime_list

from functools import reduce

def gcd(a, b):
	while b != 0:
		(a, b) = (b, a % b)
	return a

def split(n):
	d = n - 1
	s = 0
	while d % 2 == 0:
		s += 1
		d //= 2
	return (s, d)

def mult_mod(n, f, g, b, c):
	if f[0] == 0:
		x = f[1]
		return [(x*y) % n for y in g]
	# (αx + β) (γx + δ) mod (x² - bx - c) = (αγb + αδ + βγ) x + αγc + βδ
	return [(f[0]*g[0]*b + f[0]*g[1] + f[1]*g[0]) % n, (f[0]*g[0]*c + f[1]*g[1]) % n]

def pow_mod(n, base, exp, b, c):
	"""efficiently compute b^e % m
	base: [int]
	exp: int
	b, c: int
	"""
	res = [0, 1]
	while exp != 0:
		if exp % 2 == 1:
			res = mult_mod(n, res, base, b, c)
		base = mult_mod(n, base, base, b, c)
		exp //= 2
	return res

def jacobi(x, y):
	"""efficiently compute the jacobi symbol (x/y)"""
	#assert y > 3 and y % 2 == 1,("jacobi(x,y): y must be an odd integer >= 3, found {}".format(y))
	res = 1
	while True:
		x = x % y
		if x == 0:
			return 0
		while x % 2 == 0:
			x //= 2
			m8 = y % 8
			if m8 == 3 or m8 == 5:
				res = -res
		if x == 1:
			return res
		(x, y) = (y, x)
		if x % 4 == 3 and y % 4 == 3:
			res = -res

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

# taken from http://facthacks.cr.yp.to/product.html
def producttree(X):
	result = [X]
	while len(X) > 1:
		X = [reduce(lambda a, b: a * b, X[i*2:(i+1)*2]) for i in range((len(X)+1)//2)]
		result.append(X)
	return result

# dito
def remaindersusingproducttree(n,T):
	result = [n]
	for t in reversed(T):
		result = [result[i//2] % t[i] for i in range(len(t))]
	return result

def remainders(n,X):
	return remaindersusingproducttree(n,producttree(X))

foo = producttree(prime_list)[-1][0]

def trial_division(n, B):
	#return gcd(n, foo) not in [1, n]
	i = 0
	b = min(B, isqrt(n), len(prime_list))
	while prime_list[i] <= b:
		if n % prime_list[i] == 0:
			return True # found a divisor
		i+=1
	return False

def QFT(n, b, c, B, use_rqft=False):
	if not use_rqft:
		if trial_division(n, B):
			return False # composite

		if is_square(n):
			return False # composite

	x = [1, 0]

	foo = pow_mod(n, x, (n + 1) // 2, b, c)
	if foo[0] != 0:
		return False # composite

	if (foo[1]*foo[1]) % n != n-c:
		return False # composite

	(r, s) = split(n**2)
	foo = pow_mod(n, x, s, b, c)
	if foo[0] == 0 and foo[1] == 1:
		return True # probably prime
	for j in range(0, r-2 + 1):
		if foo[0] == 0 and foo[1] == n-1:
			return True # probably prime
		foo = mult_mod(n, foo, foo, b, c)
	return False # composite


def RQFT(n, B):
	#assert(n > 1)
	#assert(n % 2 == 1)

	b = c = 0
	j1 = j2 = 0

	if trial_division(n, B):
		return False # composite

	if is_square(n):
		return False # composite

	for _ in range(B):
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

	return QFT(n, b, c, B, True)

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
			self.assertTrue(RQFT(self.large_primes[n], B))

	def test_composites(self):
		counter = 0
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

	def test_squarefree_composites(self):
		for n in self.composites:
			self.assertFalse(RQFT(n, B))

if __name__ == "__main__":
	seed()
	unittest.main()

# vim: set ts=4 sw=4 noet :
