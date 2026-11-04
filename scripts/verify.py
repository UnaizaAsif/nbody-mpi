import sys
a = [float(x) for x in open(sys.argv[1]).read().split()]
b = [float(x) for x in open(sys.argv[2]).read().split()]
assert len(a) == len(b), "different sizes"
d = max(abs(x - y) for x, y in zip(a, b))
tol = 1e-9
print(f"max difference = {d:.3e}  ->", "PASS" if d < tol else "FAIL", f"(tolerance {tol})")
