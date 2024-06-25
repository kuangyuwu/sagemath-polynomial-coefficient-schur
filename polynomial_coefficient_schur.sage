from sage.combinat.sf.sfa import is_SymmetricFunction
from sage.rings.polynomial.polynomial_ring import is_PolynomialRing

schur = SymmetricFunctions(QQ).schur()

class PolynomialCoefficientSchur:

	def __init__(self, sym_func=None, coeff_ring=None, coeff_dict=None):
		
		if coeff_ring is not None and not is_PolynomialRing(coeff_ring):
			raise TypeError("Invalid coeff_ring: not a univariate polynomial ring")
		
		self._coeff_ring = coeff_ring
		self._coeff_dict = {}

		if sym_func is not None:
			if coeff_dict is not None:
				raise ValueError("__init__ does not accept both sym_func and coeff_dict")
			self._from_sym_func(sym_func)
			
		elif coeff_dict is not None:
			self._set_coeff(coeff_dict)
		
		return
	
	def _from_sym_func(self, sym_func):
		if not is_SymmetricFunction(sym_func):
			raise TypeError("Invalid sym_func: not a symmetric function")
		degree = sym_func.degree()
		for d in range(degree, -1, -1): # d = degree, degree-1, ..., 0
			for part in Partitions(d).list():
				coeff = sym_func.scalar(schur(part))
				if coeff == 0: continue
				self[part] = coeff
				sym_func -= coeff * schur(part)
				if sym_func == 0:
					break
			if sym_func == 0:
				break
		return
	
	def _set_coeff(self, coeff_dict):
		for key in coeff_dict:
			self[key] = coeff_dict[key]
		return

	def __setitem__(self, key, value):
		part = Partition(key)
		if value == 0:
			if part in self._coeff_dict:
				del self._coeff_dict[part]
			return
		if self._coeff_ring is None:
			ring = value.parent()
			if is_PolynomialRing(ring):
				self._coeff_ring = ring
			elif value not in QQ:
				raise TypeError("Invalid value: not a univariate polynomial over QQ")
			else:
				value = QQ(value)
		else:
			if not value in self._coeff_ring:
				raise TypeError(f"Invalid value: not in {self._coeff_ring}")
		self._coeff_dict[part] = value
		return

	def __getitem__(self, key):
		part = Partition(key)
		return self._coeff_dict[part] if part in self._coeff_dict else 0
	
	def __repr__(self):
		if not self._coeff_dict:
			return "0"
		result = []
		parts = list(self._coeff_dict.keys())
		parts.sort(key = lambda mu: tuple([-sum(mu)] + list(mu)))
		parts.reverse()
		for part in parts:
			coeff = self._coeff_dict[part]
			if coeff in QQ:
				if coeff > 0:
					if coeff == 1:
						result.append(f"s{part}")
					else:
						result.append(f"{coeff}*s{part}")
				else:
					if len(result) == 0:
						if coeff == -1:
							result.append(f"-s{part}")
						else:
							result.append(f"{coeff}*s{part}")
					else:
						result[-1] = "-"
						if coeff == -1:
							result.append(f"s{part}")
						else:
							result.append(f"{-coeff}*s{part}")
			else:
				result.append(f"({coeff})*s{part}")
			result.append("+")
		result.pop()
		return " ".join(result)

	def coefficient_ring(self):
		return self._coeff_ring
	
	def schurs(self):
		return self._coeff_dict.keys()
	
	def coefficients(self):
		return self._coeff_dict.values()
	
	def coefficient(self, mu):
		return self._coeff_dict[Partition(mu)] if Partition(mu) in self._coeff_dict else 0
	
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









