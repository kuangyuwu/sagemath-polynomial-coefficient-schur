from sage.combinat.sf.sfa import is_SymmetricFunction
from sage.rings.polynomial.polynomial_ring import is_PolynomialRing

schur = SymmetricFunctions(QQ).schur()

class PolynomialCoefficientSchur:


	def __init__(self, sym_func=None, coeff_ring=None, coeff_dict=None):
		r"""
		Initialize a PolynomialCoefficientSchur object.

		Parameters:
		sym_func (symmetric function, optional): The constant-coefficient symmetric function to be cast to polynomial-coefficient Schur.
		coeff_ring (PolynomialRing, optional): The polynomial ring that the coefficients belong to.
		coeff_dict (dict, optional): Dictionary representing the coefficients of the polynomial-coefficient Schur. Each key is a partition and the corresponding value is the coefficient of the Schur function indexed by that partition.

		Raises:
		ValueError: If both sym_func and coeff_dict are provided.
		TypeError: If coeff_ring is not a univariate polynomial ring over QQ.
				   If sym_func is not QQ-coefficient.
				   If a key in coeff_dict is not a partition.
				   If a value in coeff_dict is not in coeff_ring (if not None).
				   If two values in coeff_dict are not in the same univariate polynomial ring.
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
		Return the degree of the polynomial-coefficient Schur function.
		Return -1 if self == 0.

		The degree is the highest sum of all partitions corresponding to Schur functions with non-zero coefficients.

		Returns:
		int: The degree of the polynomial-coefficient Schur function.
		"""
		if self == 0:
			return -1
		return max(map(sum, self._coeff_dict.keys()))
	

	def evaluate(self, num):
		r"""
		Substitute a number for the variable in the coefficients and return the result.
		If coeff_ring is None (implying all coefficients are constant), the method returns self as a constant-coefficient symmetric function.

		Parameters:
		num (Rational): A rational number to be substituted.

		Returns:
		SymmetricFunction: The resulting symmetric function in Schur basis.

		Raises:
		TypeError: If num is not a rational number.
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

	def __call__(self, num):
		return self.evaluate(num)	

	def coeff_ring(self):
		r"""
		Return the coefficient ring (or None if not set).

		Returns:
		PolynomialRing or None: The coefficient ring.
		"""
		return self._coeff_ring
	

	def smul(self, scalar):
		r"""
		Multiply self by the scalar and return the result (scalar multiplication)

		Parameters:
		scalar (Polynomial or Rational): a polynomial or a rational number

		Returns:
		PolynomialCoefficientSchur: The resulting polynomial-coefficient Schur function.

		Raises:
		TypeError: If the scalar is not in the coefficient ring (if already set)
		           If the scalar is not a univariate polynomial over QQ
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
		r"""
		Convert a symmetric function to a polynomial-coefficient Schur function.

		Parameters:
		sym_func (SymmetricFunction): The symmetric function to convert.

		Raises:
		TypeError: If sym_func is not a symmetric function.
		"""
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
		r"""
		Set the coefficients of the polynomial-coefficient Schur function from a dictionary.

		Parameters:
		coeff_dict (dict): Dictionary representing the coefficients of the polynomial-coefficient Schur.
		"""
		for key in coeff_dict:
			self[key] = coeff_dict[key]
		return

	def __setitem__(self, key, value):
		r"""
		Set the coefficient of a particular Schur function indexed by a partition.

		Parameters:
		key (Partition): The partition indexing the Schur function.
		value (Polynomial or Rational): The coefficient to set.

		Raises:
		TypeError: If value is not a valid coefficient in the polynomial ring.
		"""
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
		r"""
		Get the coefficient of a particular Schur function indexed by a partition.

		Parameters:
		key (Partition): The partition indexing the Schur function.

		Returns:
		Polynomial or rational number: The coefficient of the Schur function.
		"""
		part = Partition(key)
		return self._coeff_dict[part] if part in self._coeff_dict else 0
	
	def __repr__(self):
		r"""
		Return a string representation of the polynomial-coefficient Schur function.

		Returns:
		str: The string representation of the polynomial-coefficient Schur function.
		"""
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
		r"""
		Check if two polynomial-coefficient Schur functions are equal.

		Parameters:
		other (PolynomialCoefficientSchur or SymmetricFunction): The other function to compare.

		Returns:
		bool: True if the functions are equal, False otherwise.

		Raises:
		NotImplemented Error: If the other function is not PolynomialCoefficientSchur.
		"""
		if other == 0:
			return not self._coeff_dict
		if is_SymmetricFunction(other):
			return self == PolynomialCoefficientSchur(other)
		if not isinstance(other, PolynomialCoefficientSchur):
			raise NotImplementedError
		if self._coeff_dict.keys() != other._coeff_dict.keys():
			return False
		for part in self._coeff_dict:
			if self._coeff_dict[part] != other._coeff_dict[part]:
				return False
		return True
	
	def __pos__(self):
		r"""
		Return the positive of the polynomial-coefficient Schur function (self).

		Returns:
		PolynomialCoefficientSchur: The positive of the function.
		"""
		return self
	
	def __neg__(self):
		r"""
		Return the negation of the polynomial-coefficient Schur function.

		Returns:
		PolynomialCoefficientSchur: The negation of the function.
		"""
		result = PolynomialCoefficientSchur()
		for part in self._coeff_dict:
			result[part] = - self[part]
		return result
	
	def __add__(self, other):
		r"""
		Add two polynomial-coefficient Schur functions.

		Parameters:
		other (PolynomialCoefficientSchur or a symmetric function): The other function to add.

		Returns:
		PolynomialCoefficientSchur: The sum of the two functions.

		Raises:
		NotImplemented Error: If the other function is not PolynomialCoefficientSchur or a symmetric function.
		ValueError: If the coefficient rings of the two functions are different.
		"""
		if is_SymmetricFunction(other):
			return self + PolynomialCoefficientSchur(other)
		if not isinstance(other, PolynomialCoefficientSchur):
			raise NotImplementedError
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
		r"""
		Subtract one polynomial-coefficient Schur function from another.

		Parameters:
		other (PolynomialCoefficientSchur or a symmetric function): The other function to subtract.

		Returns:
		PolynomialCoefficientSchur: The result of the subtraction.

		Raises:
		NotImplemented Error: If the other function is not PolynomialCoefficientSchur or a symmetric function.
		ValueError: If the coefficient rings of the two functions are different.
		"""
		return self + (- other)
	
	def __mul__(self, other):
		r"""
		Multiply two polynomial-coefficient Schur functions.

		Parameters:
		other (PolynomialCoefficientSchur or a symmetric function): The other function to multiply.

		Returns:
		PolynomialCoefficientSchur: The product of the two functions.

		Raises:
		NotImplemented Error: If the other function is not PolynomialCoefficientSchur or a symmetric function.
		ValueError: If the coefficient rings of the two functions are different.
		"""
		if is_SymmetricFunction(other):
			return self * PolynomialCoefficientSchur(other)
		if not isinstance(other, PolynomialCoefficientSchur):
			raise NotImplementedError("Use the method smul() for scalar multiplication")
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
		r"""
		Raise the polynomial-coefficient Schur function to a power.

		Parameters:
		other (int): The exponent to raise to.

		Returns:
		PolynomialCoefficientSchur: The result of the exponentiation.

		Raises:
		NotImplemented Error: If the exponent is not a non-negative integer.
		"""
		if other not in ZZ or other < 0:
			raise NotImplementedError
		if other == 0:
			return schur([])
		if other == 1:
			return self
		else:
			return self * (self ** (other - 1))









