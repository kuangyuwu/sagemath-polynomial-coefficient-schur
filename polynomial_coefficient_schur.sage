from sage.combinat.sf.sfa import is_SymmetricFunction
from sage.rings.polynomial.polynomial_ring import is_PolynomialRing

schur = SymmetricFunctions(QQ).schur()

class PolynomialCoefficientSchur:


	def __init__(self, sym_func=None, coeff_ring=None, coeff_dict=None):

		r"""
		Parameters:
		symm_func (optional): the (constant-coefficient) symmetric function to be casted to polynomial-coefficient Schur
		coeff_ring (optional): the polynomial ring that the coefficients live in
		coeff_dict (optional): the dictionary representing the coefficient of the polynomial-coefficient Schur, with each key being a partition and the corresponding value being the coefficient of the Schur indexed by the partition

		Raises:
		ValueError: if both sym_func and coeff_dict are provided
		TypeError: if the coeff_ring is not a univariate polynomial ring over QQ
		           if the sym_func is not QQ-coefficient
				   if a key in coeff_dict is not a partition
				   if a value in coeff_dict is not in coeff_ring (if not None)
				   if two values in coeff_dict are not in the same univariate polynomial rings
		"""
		
		if coeff_ring is not None and (not is_PolynomialRing(coeff_ring) or coeff_ring.base_ring() != QQ):
			raise TypeError("Invalid coeff_ring: not a univariate polynomial ring over QQ")
		
		self._coeff_ring = coeff_ring
		self._coeff_dict = {}

		if sym_func is not None:
			if coeff_dict is not None:
				raise ValueError("__init__ does not accept both sym_func and coeff_dict")
			self._from_sym_func(sym_func)
			
		elif coeff_dict is not None:
			self._set_coeff(coeff_dict)
		
		return
	

	def degree(self):

		r"""
		returns the degree of the polynomial-coefficient Schur
		the degree is the highest sum of all partitions corresponding to Schur functions with non-zero coefficients

		Returns:
		int: the degree of the polynomial-coefficient Schur
		"""

		return max(map(sum, self._coeff_dict.keys()))
	

	def evaluate(self, num):

		r"""
		substitute the number for the variable in the coefficients and return the result
		(if coeff_ring is None (implying all coefficients are constant), the method returns self as a constant-coefficient symmetric function)

		Parameters:
		num: a rational number to be substituted

		Returns:
		symmetric function (in Schur bases)
		"""

		if num not in QQ:
			raise TypeError(f"Invalid num: not in QQ")
		result = 0
		if self._coeff_ring is None:
			for part in self._coeff_dict:
				f = self._coeff_dict[part]
				result += f * schur(part)
		else:
			for part in self._coeff_dict:
				f = self._coeff_dict[part]
				result += f(num) * schur(part)
		return result
		

	def coeff_ring(self):

		r"""
		returns the coefficient ring (or None if not set)

		Returns:
		univariate polynomial ring over QQ (or None)
		"""

		return self._coeff_ring
	

	def smul(self, scalar):

		r"""
		multiply self by the scalar and return the result (scalar multiplication)

		Parameters:
		scalar: a polynomial or a rational number

		Returns:
		polynomial-coefficient-schur

		Raises:
		TypeError: if the scalar is not in the coefficient ring (if already set)
		           if the scalar is not a univariate polynomial over QQ
		"""

		if self._coeff_ring is None:
			ring = scalar.parent()
			if is_PolynomialRing(ring):
				self._coeff_ring = ring
			elif scalar not in QQ:
				raise TypeError("Invalid scalar: not a univariate polynomial over QQ")
			else:
				scalar = QQ(scalar)
		elif scalar not in self._coeff_ring:
			raise TypeError(f"Invalid scalar: not in {self._coeff_ring}")
		if self == 0:
			return PolynomialCoefficientSchur()
		result = PolynomialCoefficientSchur()
		for part in self._coeff_dict:
			result[part] = self[part] * scalar
		return result
	
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
			if is_PolynomialRing(ring) and ring.base_ring() == QQ:
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
	

	
	def __eq__(self, other):
		if other == 0:
			return not self._coeff_dict
		if is_SymmetricFunction(other):
			return self == PolynomialCoefficientSchur(other)
		if not isinstance(other, PolynomialCoefficientSchur):
			return NotImplemented
		if self._coeff_dict.keys() != other._coeff_dict.keys():
			return False
		for part in self._coeff_dict:
			if self._coeff_dict[part] != other._coeff_dict[part]:
				return False
		return True
	
	def __pos__(self):
		return self
	
	def __neg__(self):
		result = PolynomialCoefficientSchur()
		for part in self._coeff_dict:
			result[part] = - self[part]
		return result
	
	def __add__(self, other):
		if is_SymmetricFunction(other):
			return self + PolynomialCoefficientSchur(other)
		if not isinstance(other, PolynomialCoefficientSchur):
			return NotImplemented
		if self._coeff_ring is not None and other._coeff_ring is not None and self._coeff_ring != other._coeff_ring:
			raise ValueError("Invalid operands: different coefficient rings")
		if self == 0:
			return other
		if other == 0:
			return self
		result = PolynomialCoefficientSchur()
		for part in self._coeff_dict:
			result[part] = self[part] + other[part]
		for part in other._coeff_dict:
			if result[part] == 0:
				result[part] = self[part] + other[part]
		return result
	
	def __sub__(self, other):
		return self + (- other)
	
	def __mul__(self, other):
		if is_SymmetricFunction(other):
			return self * PolynomialCoefficientSchur(other)
		if not isinstance(other, PolynomialCoefficientSchur):
			return NotImplemented
		if self._coeff_ring is not None and other._coeff_ring is not None and self._coeff_ring != other._coeff_ring:
			raise ValueError("Invalid operands: different coefficient rings")
		if self == 0 or other == 0:
			return PolynomialCoefficientSchur()
		result = PolynomialCoefficientSchur()
		for part1 in self._coeff_dict:
			for part2 in other._coeff_dict:
				coeff = self[part1] * other[part2]
				product = PolynomialCoefficientSchur(schur[part1] * schur[part2])
				result += product.smul(coeff)
		return result
	
	def __rmul__(self, other):
		return self.__mul__(other)
	
	def __pow__(self, other):
		if other not in ZZ or other < 0:
			return NotImplemented
		if other == 0:
			return schur([])
		if other == 1:
			return self
		else:
			return self * (self ** (other - 1))









