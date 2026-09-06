# Saif_Raad_pfd

A Unified Non-Recursive Framework for High-Order Partial Fraction Decomposition (PFD).

## Author
**Saif R. Lazim**  
*Department of Mathematics, College of Science, University of Basrah*

## Overview
`Saif_Raad_pfd` is a Python library designed to compute exact, symbolic non-recursive partial fraction decompositions based on the mathematical research paper:
> *"A Unified Non-Recursive Framework for High-Order Partial Fraction Decomposition and Structural Challenges"* (Saif R. Lazim, 2026).

It relies on exact rational representation via `SymPy` to prevent floating-point inaccuracies.

## Features
- Non-recursive partial fraction expansion for high-order denominators.
- Decoupling of dual-binomial denominator structures without polynomial expansion.
- General multi-index summation expansions (r-th iterative expansions).
- Symbolic precision with exact fraction evaluation.
- Fully compatible with Python 3.8+.

## Installation & Usage

Ensure SymPy is installed:
pip install sympy

Import and use the library in Python:

```python
from saif_raad_pfd import SaifPFD

pfd = SaifPFD()

# Example: Theorem 2.3 General Expansion
result = pfd.formula_2_3(coeffs=[2, 3], powers=[2, 4], a_last=5, n=3, r=2)
print(result)