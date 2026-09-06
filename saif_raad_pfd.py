import sympy as sp
from math import gcd
import itertools

class SaifPFD:
    """
    Python Library implementing Saif Raad's Unified Non-Recursive Framework 
    for High-Order Partial Fraction Decomposition (PFD).
    
    Author: Saif Raad Abdul-Mawla Lazim
    Reference: Saif Raad (2026), "A Unified Non-Recursive Framework for 
    High-Order Partial Fraction Decomposition and Structural Challenges".
    Compatible with Python 3.8+.
    """
    def __init__(self):
        self.x = sp.Symbol('x')

    # =========================================================================
    # Section 1: Monomial & Binomial Denominators
    # =========================================================================
    def formula_1_1(self, a, b, m, n, x=None):
        """Theorem 1.1: Single-step expansion for monomial-binomial denominator."""
        if x is None:
            x = self.x
        a, b = sp.S(a), sp.S(b)
        return 1 / (a * x**n) - (b * x**(m - n)) / (a * (b * x**m + a))

    def formula_1_2_series(self, a, b, m, n, k=None, x=None):
        """Theorem 1.2: Iterative expansion up to k-th step."""
        if x is None:
            x = self.x
        a, b = sp.S(a), sp.S(b)
        if k is None:
            k = int(sp.ceiling(n / m))

        sum_terms = 0
        for i in range(1, k + 1):
            num = (b**(i - 1)) * ((-1)**(i - 1))
            den = (a**i) * (x**(n - m * (i - 1)))
            sum_terms += num / den

        remainder = ((-1)**k) * (b**k * x**(k * m - n)) / ((a**k) * (b * x**m + a))
        return sum_terms + remainder

    # =========================================================================
    # Section 2: General Polynomial / Multinomial Denominators
    # =========================================================================
    def formula_2_1(self, a, b, c, m, k, n, x=None):
        """Theorem 2.1: Decomposition for 1 / (x^n * (c*x^k + b*x^m + a))"""
        if x is None:
            x = self.x
        a, b, c = sp.S(a), sp.S(b), sp.S(c)
        return 1 / (a * x**n) - (c * x**(k - n) + b * x**(m - n)) / (a * (c * x**k + b * x**m + a))

    def formula_2_2(self, coeffs, powers, a_last, n, x=None):
        """Theorem 2.2: Expansion for general multinomial denominators."""
        if x is None:
            x = self.x
        a_last = sp.S(a_last)
        coeffs = [sp.S(c) for c in coeffs]
        k = len(coeffs)
        term1 = 1 / (a_last * x**n)
        num_sum = sum(coeffs[i] * x**(powers[i] - n) for i in range(k))
        den_poly = sum(coeffs[i] * x**powers[i] for i in range(k)) + a_last
        return term1 - num_sum / (a_last * den_poly)

    def formula_2_3(self, coeffs, powers, a_last, n, r=None, x=None):
        """Theorem 2.3: General r-th iterative expansion."""
        if x is None:
            x = self.x
        a_last = sp.S(a_last)
        coeffs = [sp.S(c) for c in coeffs]
        k = len(coeffs)

        if r is None:
            r = 1
            min_m = min(powers)
            while r * min_m < n:
                r += 1

        main_sum = 0
        for s in range(1, r + 1):
            s_sum = 0
            for indices in itertools.product(range(k), repeat=s):
                prod_a = sp.Mul(*[coeffs[idx] for idx in indices])
                sum_m = sum(powers[idx] for idx in indices)
                term = (prod_a * (x**(sum_m - n))) / (a_last**s * x**n)
                s_sum += term
            main_sum += ((-1)**(s + 1)) * s_sum

        rem_sum = 0
        for indices in itertools.product(range(k), repeat=r):
            prod_a = sp.Mul(*[coeffs[idx] for idx in indices])
            sum_m = sum(powers[idx] for idx in indices)
            rem_sum += prod_a * x**(sum_m - n)

        den_poly = sum(coeffs[i] * x**powers[i] for i in range(k)) + a_last
        remainder = ((-1)**r) * rem_sum / ((a_last**r) * den_poly)

        return main_sum + remainder

    # =========================================================================
    # Section 3: Dual Binomial Denominators
    # =========================================================================
    def formula_3_1(self, a, b, c, d, m, n, x=None):
        """Theorem 3.1: Decoupling dual-binomial fields without polynomial expansion."""
        if x is None:
            x = self.x
        a, b, c, d = sp.S(a), sp.S(b), sp.S(c), sp.S(d)
        r = gcd(int(n), int(m))
        k = int(m // r)
        h = int(n // r)

        s = {}
        for i in range(1, h + 1):
            for s_val in range(-1, k + 1):
                cond = i * k - s_val * h
                if 0 < cond <= h:
                    s[i] = s_val
                    break
        s[0] = -1

        det = (b**h) * (c**k) + (a**h) * (d**k) * ((-1)**(k + h - 1))

        num_A = 0
        for i in range(1, h + 1):
            s_i = s[i]
            coeff = ((-1)**(s_i + i)) * (a**(i - 1)) * (d**(s_i + 1)) * (b**(h - i)) * (c**(k - (s_i + 1)))
            num_A += coeff * (x**(n - (i * k - s_i * h) * r))

        num_B = 0
        for i in range(1, h + 1):
            s_prev = s[i - 1]
            s_curr = s[i]
            for j in range(s_prev + 1, s_curr + 1):
                coeff = ((-1)**(j + i - 1)) * (a**(i - 1)) * (d**j) * (b**(h + 1 - i)) * (c**(k - (j + 1)))
                num_B += coeff * (x**(m - (i * k - j * h) * r))

        return num_A / (det * (d * x**n + c)) + num_B / (det * (b * x**m + a))

    # =========================================================================
    # Section 4: Arbitrary Powers & Higher-Order Expansions
    # =========================================================================
    def formula_4_1(self, a, b, m, n, h_power, max_terms=5, x=None):
        """Theorem 4.1: Non-recursive expansion for high-multiplicity power h."""
        if x is None:
            x = self.x
        a, b = sp.S(a), sp.S(b)
        k = int(sp.ceiling(n / m))

        finite_sum = 0
        for i in range(1, k + 1):
            num = ((-1)**(i - 1)) * sp.factorial(i - 1) * (b**(i - (h_power + 1)))
            den = sp.factorial(h_power) * sp.factorial(i - 1 + h_power) * (a**(i + h_power)) * (x**(n - m * (i - (h_power + 1))))
            finite_sum += num / den

        infinite_sum = 0
        for j in range(max_terms):
            num = ((-1)**(j + k)) * sp.factorial(j + k) * (b**(j + k - h_power)) * (x**(m * (j + k - h_power) - n))
            den = sp.factorial(h_power) * sp.factorial(j + k - h_power) * (a**(j + k + 1))
            infinite_sum += num / den

        return finite_sum + infinite_sum