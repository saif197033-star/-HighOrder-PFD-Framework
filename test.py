from saif_raad_pfd import SaifPFD

pfd = SaifPFD()

print("=== Theorem 1.2 Test ===")
print(pfd.formula_1_2_series(a=5, b=6, m=5, n=18))

print("\n=== Theorem 2.1 Test ===")
print(pfd.formula_2_1(a=2, b=3, c=4, m=2, k=4, n=5))

print("\n=== Theorem 3.1 Test ===")
print(pfd.formula_3_1(a=1, b=2, c=3, d=4, m=3, n=2))

print("\n=== Theorem 4.1 Test ===")
print(pfd.formula_4_1(a=2, b=3, m=2, n=5, h_power=2))

input("\nPress Enter to exit...")
