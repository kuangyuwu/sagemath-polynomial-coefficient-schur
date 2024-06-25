PATH = ""
load(f"{PATH}/polynomial_coefficient_schur.sage")
R.<n> = PolynomialRing(QQ)
S.<k> = PolynomialRing(QQ)
schur = SymmetricFunctions(QQ).schur()

a = PolynomialCoefficientSchur()
print(a, a._coeff_dict, a._coeff_ring)
            # 0 {} None
print(a[[2]])
            # 0
try:
    a[[2]] = sqrt(3)
except TypeError as err:
    print(f"TypeError: \"{err}\"")
            # TypeError: "Invalid value: not a univariate polynomial over QQ"
a[[2]] = 3*n^2 + 2
print(a[[2]])
            # 3*n^2 + 2
print(a._coeff_ring)
            # Univariate Polynomial Ring in n over Rational Field
try:
    a[[2]] = k
except TypeError as err:
    print(f"TypeError: \"{err}\"")
            # TypeError: "Invalid value: not in Univariate Polynomial Ring in n over Rational Field"

a[[2]] += 4*n
print(a[[2]])
            # 3*n^2 + 4*n + 2
a[[2]] -= 2*n + 3
print(a[[2]])
            # 3*n^2 + 2*n - 1
a[[2]] *= (n + 1)
print(a[[2]])
            # 3*n^3 + 5*n^2 + n - 1
a[[2]] /= (n + 1)
print(a[[2]])
            # 3*n^2 + 2*n - 1
try:
    a[[2]] /= n^2
except TypeError as err:
    print(f"TypeError: \"{err}\"")
            # TypeError: "Invalid value: not in Univariate Polynomial Ring in n over Rational Field"
print(a[[2]])
            # 3*n^2 + 2*n - 1

a[[2,1]] = n^3
print(a._coeff_dict)
            # {[2]: 3*n^2 + 2*n - 1, [2, 1]: n^3}
a[[2,1]] = 0
print(a._coeff_dict)
            # {[2]: 3*n^2 + 2*n - 1}

try:
    print(a[[1,2]])
except ValueError as err:
    print(f"ValueError: \"{err}\"")
            # ValueError: "[1, 2] is not an element of Partitions"
try:
    a[[1,2]] = 3
except ValueError as err:
    print(f"ValueError: \"{err}\"")
            # ValueError: "[1, 2] is not an element of Partitions"




b = PolynomialCoefficientSchur(coeff_ring=R)
print(b, b._coeff_dict, b._coeff_ring)
            # 0 {} Univariate Polynomial Ring in n over Rational Field
b[[1,1]] = 4*n
try:
    a[[2]] = sqrt(3)
except TypeError as err:
    print(f"TypeError: \"{err}\"")
            # TypeError: "Invalid value: not in Univariate Polynomial Ring in n over Rational Field"
try:
    a[[2]] = k
except TypeError as err:
    print(f"TypeError: \"{err}\"")
            # TypeError: "Invalid value: not in Univariate Polynomial Ring in n over Rational Field"
print(b)
            # (4*n)*s[1, 1]

c = PolynomialCoefficientSchur(3 * schur([1,1]) - 2 * schur([]))
print(c, c._coeff_dict, c._coeff_ring)
            # -2*s[] + 3*s[1, 1] {[1, 1]: 3, []: -2} None
c2 = PolynomialCoefficientSchur(3 * schur([1,1]) - 2 * schur([]), R)
print(c2, c2._coeff_dict, c2._coeff_ring)
            # -2*s[] + 3*s[1, 1] {[1, 1]: 3, []: -2} Univariate Polynomial Ring in n over Rational Field
try:
    c3 = PolynomialCoefficientSchur(3 * schur([1,1]) - 2 * schur([]), coeff_dict={Partition([1, 1]): 3, Partition([]): -2})
except ValueError as err:
    print(f"ValueError: \"{err}\"")
            # ValueError: "__init__ does not accept both sym_func and coeff_dict"


d = PolynomialCoefficientSchur(coeff_dict={Partition([1, 1]): 1, Partition([]): -1, Partition([3,2,1]): -k^3, Partition([2,1,1]): 2/3})
print(d, d._coeff_dict, d._coeff_ring)
            # -s[] + s[1, 1] + (-k^3)*s[3, 2, 1] {[1, 1]: 1, []: -1, [3, 2, 1]: -k^3} Univariate Polynomial Ring in k over Rational Field