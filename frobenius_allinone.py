import array
from math import ceil, floor, log, sqrt
import operator
from random import randint, seed
import operator
import unittest
from small_primes import prime_list

def even(x):
	return x % 2 == 0
def odd(x):
	return x % 2 == 1

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

def Polynomial(n, coeff):
	coeff = [0] * max(3-len(coeff), 0) + [x % n for x in coeff]
	while len(coeff) > 1 and coeff[0] == 0:
		coeff.pop(0)
	return coeff

def isInt(p):
	return len(p) == 1 
    #return p[0] == 0 and p[1] == 0

def mult_mod(n, a, b, m):
	if isInt(a) == 1:
		x = a[0]
		return Polynomial(n, [x*y for y in b])
	# (αx + β) (γx + δ) mod (x² - bx - c) = (αγb + αδ + βγ) x + αγc + βδ
	return Polynomial(n, [-a[0]*b[0]*m[1] + a[0]*b[1] + a[1]*b[0], -a[0]*b[0]*m[2] + a[1]*b[1]])

def pow_mod(n, b, e, m):
	"""efficiently compute b^e % m
	b: Polynomial
	e: int
	m: Polynomial
	"""
	res = Polynomial(n, [1])
	while e != 0:
		if odd(e):
			res = mult_mod(n, res, b, m)
		b = mult_mod(n, b, b, m)
		e //= 2
	return res

def jacobi(x, y):
	"""efficiently compute the jacobi symbol (x/y)"""
	assert y > 3 and odd(y),("jacobi(x,y): y must be an odd integer >= 3, found {}".format(y))
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

def QFT(n, b, c, B, use_rqft=False):
	if not use_rqft:
		trial_division(n, B)

		if is_square(n):
			return False # composite

	x = Polynomial(n, [1, 0])
	m = Polynomial(n, [1, -b, -c])

	foo = pow_mod(n, x, (n + 1) // 2, m)
	if not isInt(foo):
		return False # composite

	foo = mult_mod(n, foo, foo, m)
	#foo = pow_mod(n, x, n + 1, m)
	if foo != [-c % n]:
		return False # composite

	(r, s) = split(n**2)
	foo = pow_mod(n, x, s, m)
	if foo == [1]:
		return True # probably prime
	for j in range(0, r-2 + 1):
		#foo = pow_mod(n, x, 2**j * s, m)
		if foo == [n-1]:
			return True # probably prime
		foo = mult_mod(n, foo, foo, m)
	return False # composite

def trial_division(n, B):
	i = 0
	b = min(B, isqrt(n), len(prime_list))
	while prime_list[i] <= b:
		if n % prime_list[i] == 0:
			return False # composite
		i+=1


def RQFT(n, B):
	assert(n > 1)
	assert(odd(n))

	if n in prime_list:
		return True
	fail = False
	b = c = 0
	j1 = j2 = 0

	trial_division(n, B)

	if is_square(n):
		return False # composite

	for _ in range(B):
		fail = False
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
			self.assertTrue(RQFT(self.large_primes[n], B) != False)

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

# vim: setlocal ts=4 sw=4 noet :
