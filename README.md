# Polynomial-coefficient Schur symmetric functions (Sagemath)

This is a SageMath implementation of polynomial-coefficient Schur functions, a generalization of [symmetric functions](https://doc.sagemath.org/html/en/reference/combinat/sage/combinat/sf/sfa.html) in the Schur bases.

The code provides a class, `PolynomialCoefficientSchur`, which allows for operations such as addition, subtraction, multiplication, and evaluation of polynomial-coefficient Schur functions.

## Installation

To use this code, you need to have SageMath installed. You can download SageMath from its [official website](https://www.sagemath.org/).

Once SageMath is installed, you can save the code locally and load it in a SageMath environment. (Replace the `path/to` with the directory where you saved the file.)

```python
load("path/to/polynomial_coefficient_schur.sage")
```

## Usage

### Class: `PolynomialCoefficientSchur`

This class represents polynomial-coefficient Schur functions.

#### Initialization `PolynomialCoefficientSchur(sym_func=None, coeff_ring=None, coeff_dict=None)`

You can initialize a `PolynomialCoefficientSchur` object in three ways:

1. Zero function:

   ```python
   F_1 = PolynomialCoefficientSchur()
   print(F_1) # Output: 0
   ```

2. With a constant-coefficient symmetric function:

   ```python
   schur = SymmetricFunctions(QQ).schur()
   sym_func = 2*schur([2,1]) - schur([1,1,1])
   F_2 = PolynomialCoefficientSchur(sym_func=sym_func)
   print(F_2) # Output: 2*s[2, 1] - s[1, 1, 1]
   ```

3. With a dictionary of coefficients:
   ```python
   R.<n> = QQ[] # Declare R as the Univariate Polynomial Ring in n over Rational Field
   coeff_dict = {Partition([2,1]): 2*n - 1, Partition([3]): n^2}
   F_3 = PolynomialCoefficientSchur(coeff_dict=coeff_dict, coeff_ring=R)
   print(F_3) # Output: (n^2)*s[3] + (2*n - 1)*s[2, 1]
   ```

##### Coefficient Ring `coeff_ring`

In 1. and 2., the coefficient ring of `F_1` and `F_2` is set to the default value `None`.
Alternatively, the coefficient ring may be specified at initialization.

```python
R.<n> = QQ[]
F_1 = PolynomialCoefficientSchur(coeff_ring=R)
F_2 = PolynomialCoefficientSchur(sym_func=sym_func, coeff_ring=R)
```

In 3. `coeff_ring=R` can be omitted, in which case the coefficient ring is still automatically set to be the Univariate Polynomial Ring in n over Rational Field.

```python
F_3 = PolynomialCoefficientSchur(coeff_dict=coeff_dict)
```

#### Operator Overloading

The class supports various operators for convenience:

- Set and get coefficients `[]` :

  ```python
  R.<n> = QQ[] # Declare R as the Univariate Polynomial Ring in n over Rational Field
  F = PolynomialCoefficientSchur()
  print(F) # Output: 0
  F[[2,1,1]] = 2*n^3 + n # Set the coefficient of s[2,1,1] to be 2*n^3 + n
  print(F) # Output: (2*n^3 + n)*s[2, 1, 1]
  print(F[[2,1,1]]) # Print the coefficient of s[2,1,1] Output: 2*n^3 + n
  ```

- Equality `==` :

  ```python
  F_1 = PolynomialCoefficientSchur()
  print(F_1) # Output: 0

  R.<n> = QQ[] # Declare R as the Univariate Polynomial Ring in n over Rational Field
  F_2 = PolynomialCoefficientSchur()
  F_2[[2]] = n + 1
  F_2[[1,1]] = 2*n
  print(F_2) # Output: (n + 1)*s[2] + (2*n)*s[1, 1]

  coeff_dict = {Partition([2]): n + 1, Partition([1,1]): 2*n}
  F_3 = PolynomialCoefficientSchur(coeff_dict=coeff_dict)
  print(F_3) # Output: (n + 1)*s[2] + (2*n)*s[1, 1]

  print(F_1 == F_2) # Output: False
  print(F_2 == F_3) # Output: True
  print(F_2 != F_3) # Output: False
  ```

- Addition `+`, subtraction `-` :

  ```python
  R.<n> = QQ[] # Declare R as the Univariate Polynomial Ring in n over Rational Field
  F_1 = PolynomialCoefficientSchur(coeff_dict={Partition([3]): n^2 + 3*n + 2, Partition([2,1]): n})
  print(F_1) # Output: (n^2 + 3*n + 2)*s[3] + (n)*s[2, 1]

  F_2 = PolynomialCoefficientSchur(coeff_dict={Partition([3]): 2*n^2 + 1, Partition([1,1,1]): 2*n + 3})
  print(F_2) # Output: (2*n^2 + 1)*s[3] + (2*n + 3)*s[1, 1, 1]

  print(F_1 + F_2) # Output: (3*n^2 + 3*n + 3)*s[3] + (n)*s[2, 1] + (2*n + 3)*s[1, 1, 1]
  print(F_1 - F_2) # Output: (-n^2 + 3*n + 1)*s[3] + (n)*s[2, 1] + (-2*n - 3)*s[1, 1, 1]
  ```

- Multiplication `*`, exponentiation: `^` :

  ```python
  R.<n> = QQ[] # Declare R as the Univariate Polynomial Ring in n over Rational Field
  F_1 = PolynomialCoefficientSchur(coeff_dict={Partition([2]): n})
  print(F_1) # Output: (n)*s[2]
  F_2 = PolynomialCoefficientSchur(coeff_dict={Partition([1,1]): n - 2})
  print(F_2) # Output: (n - 2)*s[1, 1]

  print(F_1 * F_2) # Output: (n^2 - 2*n)*s[3, 1] + (n^2 - 2*n)*s[2, 1, 1]
  print(F_1 ^ 2) # Output: (n^2)*s[4] + (n^2)*s[3, 1] + (n^2)*s[2, 2]
  ```

- Evaluation `()`: This is equivalent to the method `.evaluate(num)`

  ```python
  R.<n> = QQ[] # Declare R as the Univariate Polynomial Ring in n over Rational Field
  F = PolynomialCoefficientSchur(coeff_dict={Partition([3]): n - 2, Partition([2,1]): n^3})
  print(F) # Output: (n - 2)*s[3] + (n^3)*s[2, 1]

  print(F(-1)) # Output: -s[2, 1] - 3*s[3]
  print(F(2)) # Output: 8*s[2, 1]
  ```

#### Methods

- `coeff_ring()`: Returns the coefficient ring (or None if not set).

  ```python
  F_1 = PolynomialCoefficientSchur()
  print(F_1.coeff_ring()) # Output: None

  R.<n> = QQ[] # Declare R as the Univariate Polynomial Ring in n over Rational Field
  F_2 = PolynomialCoefficientSchur(coeff_ring=R)
  print(F_2.coeff_ring()) # Output: Univariate Polynomial Ring in n over Rational Field

  S.<x> = QQ[] # Declare S as the Univariate Polynomial Ring in x over Rational Field
  F_3 = PolynomialCoefficientSchur(coeff_dict={Partition([2,1]): 5*x})
  print(F_3.coeff_ring()) # Output: Univariate Polynomial Ring in x over Rational Field
  ```

- `degree()`: Returns the degree of the polynomial-coefficient Schur function. (Returns -1 if self == 0.)

  ```python
  R.<n> = QQ[] # Declare R as the Univariate Polynomial Ring in n over Rational Field
  F = PolynomialCoefficientSchur(coeff_dict={Partition([3]): n + 1, Partition([2,1,1]): -n^2})
  print(F) # Output: (n + 1)*s[3] + (-n^2)*s[2, 1, 1]

  print(F.degree()) # Output: 4
  # The highest-degree schur with nonzero coefficient is s[2,1,1], which has degree 2 + 1 + 1 = 4
  ```

- `evaluate(num)`: Substitutes a number for the variable in the coefficients and returns the result.

  ```python
  R.<n> = QQ[] # Declare R as the Univariate Polynomial Ring in n over Rational Field
  F = PolynomialCoefficientSchur(coeff_dict={Partition([3]): n - 2, Partition([2,1]): n^3})
  print(F) # Output: (n - 2)*s[3] + (n^3)*s[2, 1]

  print(F.evaluate(-1)) # Output: -s[2, 1] - 3*s[3]
  print(F.evaluate(2)) # Output: 8*s[2, 1]
  ```

- `smul(scalar)`: Multiplies the polynomial-coefficient Schur function by a scalar and returns the result.

  ```python
  R.<n> = QQ[] # Declare R as the Univariate Polynomial Ring in n over Rational Field
  F = PolynomialCoefficientSchur(coeff_dict={Partition([2]): n, Partition([1,1]): -2})
  print(F) # Output: (n)*s[2] - 2*s[1, 1]

  print(F.smul(-3)) # Output: (-3*n)*s[2] + 6*s[1, 1]
  print(F.smul(n + 1)) # Output: (n^2 + n)*s[2] + (-2*n - 2)*s[1, 1]
  ```
