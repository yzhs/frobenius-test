import operator

import unittest

def xgcd(a, b):
	if b == 0:
		return (1, 0)
	else:
		(q, r) = divmod(a, b)
		(s, t) = xgcd(b, r)
		return (t, s-q*t)

def inverse_mod(m, x):
	return xgcd(x, m)[0]

def sign(n):
	if n < 0:
		return -1
	elif n > 0:
		return 1
	else:
		return 0

def zip_with(op, left, right):
	l = max(len(left), len(right))
	result = [0]*l
	for i in range(l):
		if i < len(left):
			result[i] = left[i]
		if i < len(right):
			result[i] = op(result[i], right[i])
	return result

# {{{ Polynomials modulo n
class Polynomial:
	"""A univariate polynomial represented as a list of coefficients."""
	def __init__(self, n, coeff):
		"""The list of coefficients is ordered by decreasing degree of
		the variable's power."""
		self.n = n
		self.coeff = list([x % n for x in coeff])
		while len(self.coeff) > 1 and self.coeff[0] == 0:
			self.coeff.pop(0)

	def degree(self):
		return len(self.coeff) - 1
	def componentwise(self, other, op):
		if isinstance(other, Polynomial):
			assert self.n == other.n
			return Polynomial(self.n, map(lambda x: x%self.n, zip_with(op, self.coeff[::-1], other.coeff[::-1])[::-1]))
		elif isinstance(other, int):
			tmp = self.coeff[:]
			tmp[-1] = op(tmp[-1], other)
			return Polynomial(self.n, tmp)
		else:
			raise Exception("Second operand has wrong type: found {} instead of \
					Polynomial or int".format(str(type(other))))

	def __add__(self, other):
		return self.componentwise(other, operator.add)
	def __sub__(self, other):
		return self.componentwise(other, operator.sub)
	def __mul__(self, other):
		if isinstance(other, int):
			return self * Polynomial(self.n, (other,))
		assert self.n == other.n
		coeff = [0]*(1 + self.degree() + other.degree())
		for i in range(len(self.coeff)):
			for j in range(len(other.coeff)):
				coeff[i+j] += self.coeff[i]*other.coeff[j]
		return Polynomial(self.n, coeff)

	# In ZZ/pZZ truediv and floordiv are exactly the same.
	def __truediv__(self, other): return self.__divmod__(other)[0]
	def __floordiv__(self, other): return self.__divmod__(other)[0]
	def __mod__(self, other): return self.__divmod__(other)[1]
	def __divmod__(self, other):
		n = self.n
		if isinstance(other, int):
			other = other % n
			assert other != 0
			other = Polynomial(self.n, (other,))
		assert n == other.n
		"""calculate (self / other, self % other) efficiently"""
		dividend, divisor = self.coeff[:], other.coeff
		assert [] not in [dividend, divisor], "one argument is not a valid polynomial (empty list of coefficients)"
		if len(divisor) == 1:
			return (Polynomial(n, [x * inverse_mod(n, divisor[0]) for x in dividend]),
					Polynomial(n, (0,)))
		deg1 = self.degree()
		deg2 = other.degree()
		if deg2 > deg1:
			return (Polynomial(n, (0,)), self)
		res = []
		while len(dividend) >= len(divisor):
			c = dividend[0]*inverse_mod(n, divisor[0]) % n
			res.append(c)
			for j in range(len(divisor)):
				dividend[j] = (dividend[j] - c*divisor[j]) % n
			dividend.pop(0)
		return (Polynomial(n, res), Polynomial(n, dividend))
	def __pow__(base, exp):
		if not isinstance(exp, int):
			raise NotImplemented("Polynomial ** x where x is of type {} is not implemented".format(type(exp)))
		b = Polynomial(base.n, base.coeff)
		res = poly(base.n, (1,))
		while exp != 0:
			if exp % 2 == 1:
				res = res * b
			b = b * b
			exp //= 2
		return res

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

	def __neg__(self):
		return Polynomial(self.n, map(operator.neg, self.coeff))

	def __eq__(self, other):
		if isinstance(other, int):
			return self.degree() == 0 and other == self.coeff[0]
		assert self.n == other.n
		return self.coeff == other.coeff
	def __ne__(self, other):
		if isinstance(other, int):
			return self.degree() != 0 or other != self.coeff[0]
		assert self.n == other.n
		return self.coeff != other.coeff

	def __radd__(self, other): return self + other
	def __rsub__(self, other): return (-self) + other
	def __rmul__(self, other): return self * other
	def __rfloordiv__(self, other): return self // Polynomial(self.n, (other,))
	def __rtruediv__(self, other): return self // Polynomial(self.n, (other,))
	def __rmod__(self, other): return self % Polynomial(self.n, (other,))

	def __str__(self):
		res = ""
		i = len(self.coeff)-1
		j = 0
		if self.coeff[0] < 0:
			res += "-"
		while j < len(self.coeff):
			c = self.coeff[j]
			if i == 0 or j == len(self.coeff)-1:
				res += " + "
				res += str(c)
			else:
				if j != 0:
					res += " + "
				if c != 1:
					res += str(abs(c))
				res += "x"
				if i != 1:
					res += "^" + str(i)
				while j < len(self.coeff) - 1 and self.coeff[j+1] == 0:
					j += 1
					i -= 1
			i -= 1
			j += 1
		return res + " mod " + str(self.n)
	def __repr__(self):
		return "Polynomial({}, {})".format(self.n, repr(self.coeff))

def poly(n, coeff):
	return Polynomial(n, coeff)

def get_x(m):
	return poly(m, (1,0))
# }}}

class TestPolynomial(unittest.TestCase):
	def setUp(self):
		self.x = get_x(13)
		x = self.x
		self.p = x**4 - 4*x**3 - 2*x + 5
		self.q = x**2 - 7*x + 1

	def test_notation(self):
		x, p, q = self.x, self.p, self.q
		self.assertEqual(p, poly(13, (1,-4,0,-2,5)))
		self.assertEqual(q, poly(13, (1,-7,1)))
	def test_str(self):
		x, p, q = self.x, self.p, self.q
		self.assertEqual(str(p).replace(" ", ""), "x^4+9x^3+11x+5mod13")
		self.assertEqual(str(q).replace(" ", ""), "x^2+6x+1mod13")
	def test_plus(self):
		x, p, q = self.x, self.p, self.q
		self.assertEqual(p+q, x**4 - 4*x**3 + x**2 - 9*x + 6)
	def test_minus(self):
		x, p, q = self.x, self.p, self.q
		self.assertEqual(p-q, x**4 - 4*x**3 - x**2 + 5*x + 4)
	def test_times(self):
		x, p, q = self.x, self.p, self.q
		self.assertEqual(p*q, x**6 - 11*x**5 + 29*x**4 - 6*x**3 + 19*x**2 -37*x + 5)
	def test_div(self):
		x, p, q = self.x, self.p, self.q
		self.assertEqual(p//q, x**2 + 3*x + 20)
	def test_mod(self):
		x, p, q = self.x, self.p, self.q
		self.assertEqual(p%q, x*135 - 15)
	def foo(self):
		x, p, q = self.x, self.p, self.q
		self.assertEqual(gcd(p, q), None) 

# vim: set ts=4 sw=4 noet :
