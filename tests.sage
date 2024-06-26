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
    c4 = PolynomialCoefficientSchur(58)
except TypeError as err:
    print(f"TypeError: \"{err}\"")
            # TypeError: "Invalid sym_func: not a symmetric function"
try:
    c3 = PolynomialCoefficientSchur(3 * schur([1,1]) - 2 * schur([]), coeff_dict={Partition([1, 1]): 3, Partition([]): -2})
except ValueError as err:
    print(f"ValueError: \"{err}\"")
            # ValueError: "__init__ does not accept both sym_func and coeff_dict"


d = PolynomialCoefficientSchur(coeff_dict={Partition([1, 1]): 1, Partition([]): -1, Partition([3,2,1]): -k^3, Partition([2,1,1]): 2/3})
print(d, d._coeff_dict, d._coeff_ring)
            # -s[] + s[1, 1] + (-k^3)*s[3, 2, 1] {[1, 1]: 1, []: -1, [3, 2, 1]: -k^3} Univariate Polynomial Ring in k over Rational Field





print(a.degree(), b.degree(), c.degree(), d.degree())
            # 2 2 2 6
a2 = a.evaluate(1/2)
print(a2, type(a2))
            # 3/4*s[2] <class 'sage.combinat.sf.schur.SymmetricFunctionAlgebra_schur_with_category.element_class'>
try:
    a2 = a.evaluate(sqrt(2))
except TypeError as err:
    print(f"TypeError: \"{err}\"")
            # TypeError: "Invalid num: not in QQ"
e = PolynomialCoefficientSchur(coeff_dict={Partition([2, 1]): 2*n - 4, Partition([2,2,1]): -n^4})
print(e.degree(), e.evaluate(2))
            # 5 -16*s[2, 2, 1]





e2 = PolynomialCoefficientSchur()
print(e2 == 0)
            # True
e2[[2,1]] = 2*n - 4
e2[[2,2,1]] = -n^4
print(e == e2)
            # True
c4 = 3 * schur([1,1]) - 2 * schur([])
print(c == c4)
            # True
print(d == e, d != e)
            # False True
f1 = PolynomialCoefficientSchur(coeff_dict={Partition([3]): n})
f2 = PolynomialCoefficientSchur(coeff_dict={Partition([3]): k})
print(f1 == f2, f1 != f2)
            # False True


g1 = PolynomialCoefficientSchur(coeff_dict={Partition([3]): n, Partition([2,1]): 3, Partition([1,1,1]): 2*n^2 + 3, Partition([1,1]): 2})
g2 = PolynomialCoefficientSchur(coeff_dict={Partition([3]): -n, Partition([2,1]): 3, Partition([1,1,1]): n^2 + 4*n, Partition([3,1]): n - 4})
g3 = PolynomialCoefficientSchur((1/2) * schur([3]) - 3 * schur([2,2])) 
print(+g1)
            # 2*s[1, 1] + (n)*s[3] + 3*s[2, 1] + (2*n^2 + 3)*s[1, 1, 1]
print(-g1)
            # -2*s[1, 1] + (-n)*s[3] - 3*s[2, 1] + (-2*n^2 - 3)*s[1, 1, 1]
print(g1 + g2)
            # 2*s[1, 1] + 6*s[2, 1] + (3*n^2 + 4*n + 3)*s[1, 1, 1] + (n - 4)*s[3, 1]
print(g1 - g2)
            # 2*s[1, 1] + (2*n)*s[3] + (n^2 - 4*n + 3)*s[1, 1, 1] + (-n + 4)*s[3, 1]
print(g1 + g3)
            # 2*s[1, 1] + (n + 1/2)*s[3] + 3*s[2, 1] + (2*n^2 + 3)*s[1, 1, 1] - 3*s[2, 2]
print(g1 - g3)
            # 2*s[1, 1] + (n - 1/2)*s[3] + 3*s[2, 1] + (2*n^2 + 3)*s[1, 1, 1] + 3*s[2, 2]
g4 = PolynomialCoefficientSchur(coeff_dict={Partition([1]): k})
try:
    print(g1 + g4)
except ValueError as err:
    print(f"ValueError: \"{err}\"")
            # ValueError: "Invalid operands: different coefficient rings"




print(g1.smul(-3/2))
            # -3*s[1, 1] + (-3/2*n)*s[3] - 9/2*s[2, 1] + (-3*n^2 - 9/2)*s[1, 1, 1]
print(g1.smul(2 * n - 1))
            # (4*n - 2)*s[1, 1] + (2*n^2 - n)*s[3] + (6*n - 3)*s[2, 1] + (4*n^3 - 2*n^2 + 6*n - 3)*s[1, 1, 1] 
h1 = PolynomialCoefficientSchur(coeff_dict={Partition([2]): n - 1, Partition([1,1]): 1/2})
h2 = 2 * schur([2]) - (1/4) * schur([1,1])
print(g1 * h1)
            # (2*n - 2)*s[3, 1] + s[2, 2] + (2*n - 1)*s[2, 1, 1] + s[1, 1, 1, 1] + (n^2 - n)*s[5] + (n^2 + 5/2*n - 3)*s[4, 1] + (n^2 + 2*n - 3/2)*s[3, 2] + (2*n^3 - 2*n^2 + 13/2*n - 9/2)*s[3, 1, 1] + (n^2 + 3*n)*s[2, 2, 1] + (2*n^3 - n^2 + 3*n)*s[2, 1, 1, 1] + (n^2 + 3/2)*s[1, 1, 1, 1, 1]
print(h1 * g1)
            # (2*n - 2)*s[3, 1] + s[2, 2] + (2*n - 1)*s[2, 1, 1] + s[1, 1, 1, 1] + (n^2 - n)*s[5] + (n^2 + 5/2*n - 3)*s[4, 1] + (n^2 + 2*n - 3/2)*s[3, 2] + (2*n^3 - 2*n^2 + 13/2*n - 9/2)*s[3, 1, 1] + (n^2 + 3*n)*s[2, 2, 1] + (2*n^3 - n^2 + 3*n)*s[2, 1, 1, 1] + (n^2 + 3/2)*s[1, 1, 1, 1, 1]
print(g1 * h2)
            # 4*s[3, 1] - 1/2*s[2, 2] + 7/2*s[2, 1, 1] - 1/2*s[1, 1, 1, 1] + (2*n)*s[5] + (7/4*n + 6)*s[4, 1] + (2*n + 21/4)*s[3, 2] + (4*n^2 - 1/4*n + 45/4)*s[3, 1, 1] + (-1/2*n^2 + 9/2)*s[2, 2, 1] + (7/2*n^2 + 9/2)*s[2, 1, 1, 1] + (-1/2*n^2 - 3/4)*s[1, 1, 1, 1, 1]
print(h2 * g1)
            # 4*s[3, 1] - 1/2*s[2, 2] + 7/2*s[2, 1, 1] - 1/2*s[1, 1, 1, 1] + (2*n)*s[5] + (7/4*n + 6)*s[4, 1] + (2*n + 21/4)*s[3, 2] + (4*n^2 - 1/4*n + 45/4)*s[3, 1, 1] + (-1/2*n^2 + 9/2)*s[2, 2, 1] + (7/2*n^2 + 9/2)*s[2, 1, 1, 1] + (-1/2*n^2 - 3/4)*s[1, 1, 1, 1, 1]
print(h1 ^ 2, h1 * h1)
            # (n^2 - 2*n + 1)*s[4] + (n^2 - n)*s[3, 1] + (n^2 - 2*n + 5/4)*s[2, 2] + (n - 3/4)*s[2, 1, 1] + 1/4*s[1, 1, 1, 1] (n^2 - 2*n + 1)*s[4] + (n^2 - n)*s[3, 1] + (n^2 - 2*n + 5/4)*s[2, 2] + (n - 3/4)*s[2, 1, 1] + 1/4*s[1, 1, 1, 1]
print(h1 ^ 3, h1 * h1 * h1)
            # (n^3 - 3*n^2 + 3*n - 1)*s[6] + (2*n^3 - 9/2*n^2 + 3*n - 1/2)*s[5, 1] + (3*n^3 - 15/2*n^2 + 27/4*n - 9/4)*s[4, 2] + (n^3 - 9/4*n + 5/4)*s[4, 1, 1] + (n^3 - 3/2*n^2 + 5/8)*s[3, 3] + (2*n^3 - 3*n^2 + 3/2*n - 1/4)*s[3, 2, 1] + (3/2*n^2 - 3/2*n + 1/8)*s[3, 1, 1, 1] + (n^3 - 3*n^2 + 15/4*n - 13/8)*s[2, 2, 2] + (3/2*n^2 - 9/4*n + 9/8)*s[2, 2, 1, 1] + (3/4*n - 1/2)*s[2, 1, 1, 1, 1] + 1/8*s[1, 1, 1, 1, 1, 1] (n^3 - 3*n^2 + 3*n - 1)*s[6] + (2*n^3 - 9/2*n^2 + 3*n - 1/2)*s[5, 1] + (3*n^3 - 15/2*n^2 + 27/4*n - 9/4)*s[4, 2] + (n^3 - 9/4*n + 5/4)*s[4, 1, 1] + (n^3 - 3/2*n^2 + 5/8)*s[3, 3] + (2*n^3 - 3*n^2 + 3/2*n - 1/4)*s[3, 2, 1] + (3/2*n^2 - 3/2*n + 1/8)*s[3, 1, 1, 1] + (n^3 - 3*n^2 + 15/4*n - 13/8)*s[2, 2, 2] + (3/2*n^2 - 9/4*n + 9/8)*s[2, 2, 1, 1] + (3/4*n - 1/2)*s[2, 1, 1, 1, 1] + 1/8*s[1, 1, 1, 1, 1, 1]

zero = PolynomialCoefficientSchur()
print(zero == 0)
            # True
print(+zero, type(+zero))
            # 0 <class '__main__.PolynomialCoefficientSchur'>
print(-zero, type(-zero))
            # 0 <class '__main__.PolynomialCoefficientSchur'>
print(zero + h1)
            # (n - 1)*s[2] + 1/2*s[1, 1]
print(h1 + zero)
            # (n - 1)*s[2] + 1/2*s[1, 1]
print(zero - h1)
            # (-n + 1)*s[2] - 1/2*s[1, 1]
print(h1 - zero)
            # (n - 1)*s[2] + 1/2*s[1, 1]
print(h1.smul(0), type(h1.smul(0)))
            # 0 <class '__main__.PolynomialCoefficientSchur'>
print(zero.smul(-5/7), type(zero.smul(-5/7)))
            # 0 <class '__main__.PolynomialCoefficientSchur'>
print(zero.smul(n^8), type(zero.smul(n^8)))
            # 0 <class '__main__.PolynomialCoefficientSchur'>
print(zero * h1, type(zero * h1))
            # 0 <class '__main__.PolynomialCoefficientSchur'>
print(h1 * zero, type(h1 * zero))
            # 0 <class '__main__.PolynomialCoefficientSchur'>
print(zero ^ 11, type(zero ^ 11))
            # 0 <class '__main__.PolynomialCoefficientSchur'>