from functools import cache

import sage.libs.lrcalc.lrcalc as lrcalc
from sage.combinat.sf.sfa import is_SymmetricFunction

schur = SymmetricFunctions(QQ).schur()

@cache
def product_schurs(mu1, mu2):
	return schur(mu1) * schur(mu2)

class PolynomialCoefficientSchur:
	
	def __init__(self, coefficient_ring):
		self.__coefficients = {}
		if not sage.rings.polynomial.polynomial_ring.is_PolynomialRing(coefficient_ring):
			raise TypeError("coefficient_ring should be a univariate polynomial ring")
		self.__coefficient_ring = coefficient_ring
		return
	
	@classmethod
	def from_symmetric_function(cls, symmetric_function, coefficient_ring):
		if not is_SymmetricFunction(symmetric_function):
			raise TypeError("symmetric_function should be a symmetric function")
		result = cls(coefficient_ring)
		degree = symmetric_function.degree()
		for d in range(degree, -1, -1): # d = degree, degree-1, ..., 0
			for mu in Partitions(d).list():
				coefficient = symmetric_function.scalar(schur(mu))
				if coefficient == 0: continue
				result.set_coefficient(mu, coefficient)
				symmetric_function -= coefficient * schur(mu)
				if symmetric_function == 0:
					break
			if symmetric_function == 0:
				break
		return result
	
	def coefficient_ring(self):
		return self.__coefficient_ring
	
	def schurs(self):
		return self.__coefficients.keys()
	
	def coefficients(self):
		return self.__coefficients.values()
	
	def coefficient(self, mu):
		return self.__coefficients[Partition(mu)] if Partition(mu) in self.__coefficients else 0
	
	def set_coefficient(self, mu, coefficient):
		if not coefficient in self.coefficient_ring():
			raise TypeError(f"coefficient should be in {self.coefficient_ring()}")
		self.__coefficients[Partition(mu)] = coefficient
		return
	
	def evaluate(self, num):
		if not num in QQ:
			raise TypeError(f"num should be in QQ")
		result = 0
		for mu in self.schurs():
			f = self.coefficient(mu)
			result += f(num) * schur(mu)
		return result
	
	def degree(self):
		return max(map(sum, self.schurs()))
	
	def __repr__(self):
		if not bool(self.schurs()):
			return "0"
		result = []
		schurs = list(self.schurs())
		schurs.sort(key = lambda mu: tuple([-sum(mu)] + list(mu)))
		schurs.reverse()
		for mu in schurs:
			coefficient = self.coefficient(mu)
			if coefficient in QQ:
				if coefficient > 0:
					if coefficient == 1:
						result.append(f"s{mu}")
					else:
						result.append(f"{coefficient}*s{mu}")
				else:
					if len(result) == 0:
						if coefficient == -1:
							result.append(f"-s{mu}")
						else:
							result.append(f"{coefficient}*s{mu}")
					else:
						result[-1] = "-"
						if coefficient == -1:
							result.append(f"s{mu}")
						else:
							result.append(f"{-coefficient}*s{mu}")
			else:
				result.append(f"({self.coefficient(mu)})*s{mu}")
			result.append("+")
		result.pop()
		return " ".join(result)
	
	def __eq__(self, other):
		if other == 0:
			return not bool(self.schurs())
		if not isinstance(other, PolynomialCoefficientSchur):
			return NotImplemented
		if self.coefficient_ring() != other.coefficient_ring():
			raise TypeError("the coefficient rings of the two PolynomialCoefficientSchur objects should be the same")
		if self.schurs() != other.schurs():
			return False
		for mu in self.schurs():
			if self.coefficient(mu) != other.coefficient(mu):
				return False
		return True
	
	def __pos__(self):
		return self
	
	def __neg__(self):
		result = PolynomialCoefficientSchur(self.coefficient_ring())
		for mu in self.schurs():
			result.set_coefficient(mu, - self.coefficient(mu))
		return result
	
	def __add__(self, other):
		if is_SymmetricFunction(other):
			other = PolynomialCoefficientSchur.from_symmetric_function(other, self.coefficient_ring())
		if not isinstance(other, PolynomialCoefficientSchur):
			return NotImplemented
		if self.coefficient_ring() != other.coefficient_ring():
			raise TypeError("the coefficient rings of the two PolynomialCoefficientSchur objects should be the same")
		result = PolynomialCoefficientSchur(self.coefficient_ring())
		for mu in self.schurs():
			if mu in other.schurs():
				coefficient = self.coefficient(mu) + other.coefficient(mu)
				if coefficient != 0:
					result.set_coefficient(mu, coefficient)
			else:
				result.set_coefficient(mu, self.coefficient(mu))
		for mu in other.schurs() - self.schurs():
			result.set_coefficient(mu, other.coefficient(mu))
		return result
	
	def __sub__(self, other):
		return self + (- other)
	
	def scalar_multiplication(self, scalar):
		if not scalar in self.coefficient_ring():
			raise TypeError(f"scalar should be in {self.coefficient_ring()}")
		result = PolynomialCoefficientSchur(self.coefficient_ring())
		for mu in self.schurs():
			result.set_coefficient(mu, self.coefficient(mu) * scalar)
		return result
	
	def __mul__(self, other):
		if is_SymmetricFunction(other):
			other = PolynomialCoefficientSchur.from_symmetric_function(other, self.coefficient_ring())
		if not isinstance(other, PolynomialCoefficientSchur):
			return NotImplemented
		if self.coefficient_ring() != other.coefficient_ring():
			raise TypeError("the coefficient rings of the two PolynomialCoefficientSchur objects should be the same")
		result = PolynomialCoefficientSchur(self.coefficient_ring())
		for mu1 in self.schurs():
			for mu2 in other.schurs():
				coefficient = self.coefficient(mu1) * other.coefficient(mu2)
				product = PolynomialCoefficientSchur.from_symmetric_function(product_schurs(mu1, mu2), self.coefficient_ring())
				result += (product.scalar_multiplication(coefficient))
		return result
	
	def __rmul__(self, other):
		return self.__mul__(other)
	
	def __pow__(self, other):
		if not other in ZZ or other < 0:
			return NotImplemented
		if other == 0:
			return schur([])
		if other == 1:
			return self
		else:
			return (self ** (other >> 1)) * (self ** ((other + 1) >> 1))









