#!/usr/bin/env python
# coding: utf-8

# # Finite field element class

# In[1]:


MOD = 31

class FieldElement:
    def __init__(self, value):
        """Initialize a field element"""
        self.value = value % MOD

    def __add__(self, other):
        """Add two field elements."""
        return FieldElement((self.value + other.value) % MOD)

    def __sub__(self, other):
        """Subtract two field elements."""
        return FieldElement((self.value - other.value) % MOD)

    def __mul__(self, other):
        """Multiply two field elements."""
        return FieldElement((self.value * other.value) % MOD)

    def __truediv__(self, other):
        """Divide two field elements (self * other^-1)."""
        return self * other.inv()

    def __neg__(self):
        """Negate a field element."""
        return FieldElement(-self.value % MOD)

    def __repr__(self):
        """String representation of the field element."""
        return f"{self.value}"

    def __eq__(self, other):
        """Check equality of two field elements."""
        return self.value == other.value


    def __pow__(self, exponent):
        """Compute exponentiation"""
        if exponent < 0:
            raise ValueError("Negative exponents are not supported")
        result = FieldElement(1)
        base = self
        while exponent > 0:
            if exponent % 2 == 1:
                result *= base
            base = base.square()
            exponent //= 2
        return result

    def inv(self):
        """Compute the multiplicative inverse using Fermat's little theorem."""
        if self.value == 0:
            raise ZeroDivisionError("Cannot invert zero")
        return FieldElement(pow(self.value, MOD - 2, MOD))

    def square(self):
        """Compute the square of the field element."""
        return self * self

    def double(self):
        """Double the field element."""
        return self + self

    @classmethod
    def one(cls):
        """Return the multiplicative identity (1)."""
        return cls(1)

    @classmethod
    def zero(cls):
        """Return the additive identity (0)."""
        return cls(0)


# In[2]:


a = FieldElement(5)
b = FieldElement(10)

print(f"a + b: {a + b}")  # 15
print(f"a - b: {a - b}")  # 26 (5 - 10 ≡ 26 mod 31)
print(f"a * b: {a * b}")  # 19 (50 ≡ 19 mod 31)
print(f"a / b: {a / b}")  # 16 (5 * 10^{-1} ≡ 16 mod 31)

print(f"a * a.inv(): {a * a.inv()}")  # 1

print(f"a ** 2: {a ** 2}")       # 25
print(f"a.double(): {a.double()}")  # 10


# # Circle point class for cirlce group math

# In[3]:


class CirclePoint:
    """A point (x, y) on the circle where x^2 + y^2 = 1, Circle Curve"""
    def __init__(self, x, y):
        """Initialize a point with coordinates x and y, ensuring it lies on the circle."""
        if (x.square() + y.square()).value != 1:
            raise ValueError("Point does not lie on the circle: x^2 + y^2 ≠ 1")
        self.x = x
        self.y = y

    def __repr__(self):
        """String representation of the circle point."""
        return f"Point({self.x}, {self.y})"

    @classmethod
    def identity(cls):
        """Return the identity element (1, 0) of the circle group."""
        return cls(FieldElement.one(), FieldElement.zero())

    def add(self, other):
        """Perform group operation: (x1,y1)・(x2,y2) = (x1*x2 - y1*y2, x1*y2 + x2*y1)."""
        nx = self.x * other.x - self.y * other.y
        ny = self.x * other.y + other.x * self.y
        return CirclePoint(nx, ny)

    def double(self):
        """Apply squaring map: π(x,y) = (2*x^2 - 1, 2*x*y), since x^2 + y^2 = 1."""
        xx = self.x.square().double() - FieldElement.one()
        yy = self.x.double() * self.y
        return CirclePoint(xx, yy)

    def power(self, n):
        """Compute the nth power of the point."""
        if n < 0:
            raise ValueError("Negative exponents are not supported")
        out = CirclePoint.identity()
        base = CirclePoint(self.x, self.y)
        tmp = n
        while tmp > 0:
            if tmp & 1:
                out = out.add(base)
            base = base.double()
            tmp //= 2
        return out

    def inverse(self):
        """Return the inverse element (x, -y)."""
        return CirclePoint(self.x, -self.y)

    def __eq__(self, other):
        """Check equality of two circle points."""
        if not isinstance(other, CirclePoint):
            return False
        return self.x == other.x and self.y == other.y

    def __hash__(self):
        """Hash function for use in sets."""
        return hash((self.x.value, self.y.value))

def find_generator():
    """Find a generator for the circle group where x^2 + y^2 = 1)."""
    for xx in range(1, MOD):
        for yy in range(1, MOD):
            if (xx * xx + yy * yy) % MOD == 1:
                return CirclePoint(FieldElement(xx), FieldElement(yy))
    return CirclePoint.identity()

def get_nth_generator(logn):
    """Get the generator for subgroup of size 2^logn."""
    if (MOD + 1) % (1 << logn) != 0:
        raise ValueError(f"Cannot generate subgroup of size 2^{logn}")
    step = (MOD + 1) // (1 << logn)
    g = find_generator()
    return g.power(step)


# In[4]:


p1 = CirclePoint(FieldElement(13), FieldElement(7))
p2 = CirclePoint(FieldElement(30), FieldElement(0))

# Identity
identity = CirclePoint.identity()
print(f"Identity: {identity}")

# group operation
p3 = p1.add(p2)
print(f"p1・p2 = {p3}")

# squaring map
p1_double = p1.double()
print(f"π(p1) = {p1_double}")

# Inverse
p1_inv = p1.inverse()
print(f"p1's inverse = {p1_inv}")
p1_plus_inv = p1.add(p1_inv)
print(f"p1・p1_inv = {p1_plus_inv}")

# power
p1_power_31 = p1.power(31)
p1_power_32 = p1.power(32)
print(f"p1^{31} = {p1_power_31}")
print(f"p1^{31} = p1^{-1} : {p1_power_31 == p1.inverse()}")
print(f"p1^{32} = {p1_power_32}")


# In[5]:


p1 = CirclePoint(FieldElement(1), FieldElement(0))
p2 = CirclePoint(FieldElement(30), FieldElement(0))
q = CirclePoint(FieldElement(2), FieldElement(11))

a = q.add(p1)
print(a)
b = q.add(p2)
print(b)


# # Generate subgroup
# 

# In[6]:


def generate_subgroup(k: int) -> list[CirclePoint]:
    """Generate a subgroup of size 2^k using the generator."""
    g_k = get_nth_generator(k)
    p = CirclePoint.identity()
    return [ (p := p.add(g_k)) if i else p
             for i in range(1 << k) ]


def print_all_subgroups(max_k: int = 5) -> None:
    """Print all subgroups up to size 2^max_k."""
    print("Subgroups of C(F_p):")
    for k in range(max_k + 1):
        print(f"Size {1 << k}: ", generate_subgroup(k))


print_all_subgroups()

G3 = generate_subgroup(3)
print("G3: size 8 subgroup")
for pt in G3:
    print(pt)


# # Generate twin coset

# In[7]:


def generate_twin_coset(Q, G_n_minus_1):
    """Generate twin coset: D = Q・G_(n-1) ∪ Q^{-1}・G_(n-1)."""
    Q_inv = Q.inverse()
    coset_Q = [Q.add(g) for g in G_n_minus_1]
    coset_Q_inv = [Q_inv.add(g) for g in G_n_minus_1]
    return coset_Q + coset_Q_inv

G2 = generate_subgroup(2)
Q = CirclePoint(FieldElement(13), FieldElement(7))
twin_cosets = generate_twin_coset(Q, G2)
print(twin_cosets)


# # Generate standard position coset D = Q・G_n

# In[8]:


def generate_standard_position_coset(Q, G_n):
    """Generate standard position coset D = Q・G_n."""
    return [Q.add(g) for g in G_n]

G3 = generate_subgroup(3)
Q = CirclePoint(FieldElement(13), FieldElement(7))
D = generate_standard_position_coset(Q, G3)
print("Standard Position Coset D (size 8):")
for pt in D:
    print(pt)


# # Standard position cose and twin coset

# In[9]:


D_std = generate_standard_position_coset(Q, G3)
D_twin = generate_twin_coset(Q, G2)
print("Q・G3:", D_std)
print("(Q・G2) ∪ (Q^-1・G2):", D_twin)
if set(D_std) == set(D_twin):
    print("Equal: standard position coset Q·G3 = Q·G2 ∪ Q^-1·G2")
else:
    print("Not Equal")


# # Decompose standard position coset into twin cosets

# In[10]:


def decompose_to_twin_cosets(D, n, m):
    """Decompose standard position coset D (size 2^m) into twin cosets of size 2^n."""
    assert len(D) == 1 << m
    G_n_minus_1 = generate_subgroup(n - 1)
    twin_cosets = []
    for k in range((1 << m) // (1 << n)):
        Q_k = Q.power(4 * k + 1)
        twin_coset_k = generate_twin_coset(Q_k, G_n_minus_1)
        twin_cosets.append(twin_coset_k)
    return twin_cosets

print("Standard Positon Coset :", D)
twin_cosets = decompose_to_twin_cosets(D, 2, 3)
for i, tc in enumerate(twin_cosets):
    print(f"Twin-Coset {i} (size 4):")
    for pt in tc:
        print(pt)


# # Halve standard position coset using squaring map

# In[11]:


def halve_standard_position_coset_by_squaring_map(standard_coset):
    """Apply squaring map to reduce coset size from 2^n to 2^(n-1)."""
    new_points = [pt.double() for pt in standard_coset]
    unique_dict = {}
    for pt in new_points:
        key = (pt.x.value, pt.y.value)
        unique_dict[key] = pt
    unique = list(unique_dict.values())
    unique.sort(key=lambda p: (p.x.value, p.y.value))
    return unique

D_halved = halve_standard_position_coset_by_squaring_map(D)
print("Result of halving Standard Position Coset D by square mapping (size):", len(D_halved))
for pt in D_halved:
    print(pt)


# # Circle FFT and IFFT

# In[12]:


def circle_fft(evaluations, circle_domain):
    """Transform evaluations on circle domain to polynomial coefficients using Circle FFT."""
    if len(evaluations) != len(circle_domain) or len(evaluations) & (len(evaluations) - 1):
        raise ValueError("evaluations and circle_domain must have equal power-of-2 length")
    arr = [FieldElement(e) for e in evaluations]
    coefficients = _fft_quotient_map_y(arr, circle_domain)
    return [c.value for c in coefficients]

def _fft_quotient_map_y(values, circle_domain):
    """FFT layer using quotient map φ_y for y-component decomposition."""
    n = len(values)
    if n == 1:
        return [values[0]]
    half = n // 2
    # Split into left and right halves
    left = values[:half]
    right = list(reversed(values[half:]))
    y_factors = [pt.y for pt in circle_domain[:half]]
    two = FieldElement(2)
    # Compute even and odd parts with respect to y
    # Even part: f_0(x) = (f(x,y)+f(x,-y))/2
    f_even_y = [(left[i] + right[i]) / two for i in range(half)]
    # Odd part: f_1(x) = (f(x,y)-f(x,-y))/2*y
    f_odd_y = [(left[i] - right[i]) / (y_factors[i] * two) for i in range(half)]
    x_domain = [pt.x for pt in circle_domain]
    half_x_domain = x_domain[:half]
    # Recurse on x-component
    coeffs_even = _fft_squaring_map_x(f_even_y, half_x_domain)
    coeffs_odd = _fft_squaring_map_x(f_odd_y, half_x_domain)
    return [coeff for pair in zip(coeffs_even, coeffs_odd) for coeff in pair]

def _fft_squaring_map_x(values, x_domain):
    """FFT layer using squaring map π: x → (2*x^2 - 1) for x-component decomposition."""
    n = len(values)
    if n == 1:
        return [values[0]]
    half = n // 2
    # Split into left and right halves
    left = values[:half]
    right = list(reversed(values[half:]))
    x_factors = x_domain[:half]
    two = FieldElement(2)
    # Compute even and odd parts with respect to x
    # Even part: f_even(π(x)) = (f_0(x)+f_0(-x))/2
    f_even_x = [(left[i] + right[i]) / two for i in range(half)]
    # Odd part: f_odd(π(x)) = (f_0(x)-f_0(-x))/2*x
    f_odd_x = [(left[i] - right[i]) / (x_factors[i] * two) for i in range(half)]
    # Compute next x-domain using squaring map
    # π(x) = 2x^2 - 1
    next_x_domain = [(x.square().double() - FieldElement.one()) for x in x_domain[:half]]
    # Recurse on the next domain
    coeffs_even = _fft_squaring_map_x(f_even_x, next_x_domain)
    coeffs_odd = _fft_squaring_map_x(f_odd_x, next_x_domain)
    return [coeff for pair in zip(coeffs_even, coeffs_odd) for coeff in pair]

def circle_ifft(coefficients, circle_domain):
    """Reconstruct evaluations from Circle FFT coefficients using Inverse Circle FFT."""
    if len(coefficients) != len(circle_domain) or len(coefficients) & (len(coefficients) - 1):
        raise ValueError("coefficients and circle_domain must have equal power-of-2 length")
    arr = [FieldElement(c) for c in coefficients]
    evaluations = _ifft_quotient_map_y(arr, circle_domain)
    return [e.value for e in evaluations]

def _ifft_quotient_map_y(coefficients, circle_domain):
    """IFFT layer using quotient map φ_y for y-component reconstruction."""
    n = len(coefficients)
    if n == 1:
        return [coefficients[0]]
    half = n // 2
    # Split coefficients into even and odd indices
    coeffs_even = coefficients[::2]
    coeffs_odd = coefficients[1::2]
    x_domain = [pt.x for pt in circle_domain]
    y_factors = [pt.y for pt in circle_domain]
    half_x_domain = x_domain[:half]
    # Recurse on x-component
    f_even_y = _ifft_squaring_map_x(coeffs_even, half_x_domain)
    f_odd_y = _ifft_squaring_map_x(coeffs_odd, half_x_domain)
    # Reconstruct left and right halves
    evaluations_left = [f_even_y[i] + (y_factors[i] * f_odd_y[i]) for i in range(half)]
    evaluations_right = [f_even_y[i] - (y_factors[i] * f_odd_y[i]) for i in range(half)]
    evaluations_right.reverse()
    return evaluations_left + evaluations_right

def _ifft_squaring_map_x(coefficients, x_domain):
    """IFFT layer using squaring map π for x-component reconstruction."""
    n = len(coefficients)
    if n == 1:
        return [coefficients[0]]
    half = n // 2
    # Split coefficients into even and odd indices
    coeffs_even = coefficients[::2]
    coeffs_odd = coefficients[1::2]
    x_factors = x_domain[:half]
    # Compute next x-domain using squaring map
    next_x_domain = [(x.square().double() - FieldElement.one()) for x in x_domain[:half]]
    # Recurse on the next domain
    f_even_x = _ifft_squaring_map_x(coeffs_even, next_x_domain)
    f_odd_x = _ifft_squaring_map_x(coeffs_odd, next_x_domain)
    # Reconstruct left and right halves
    evaluations_left = [f_even_x[i] + (x_factors[i] * f_odd_x[i]) for i in range(half)]
    evaluations_right = [f_even_x[i] - (x_factors[i] * f_odd_x[i]) for i in range(half)]
    evaluations_right.reverse()
    return evaluations_left + evaluations_right

Q = CirclePoint(FieldElement(24), FieldElement(18))
Domain = generate_standard_position_coset(Q, G3)
evals = [1, 2, 4, 8, 16, 1, 2, 4]
print("Initial evaluations:", evals)
coefs = circle_fft(evals, Domain)
print("FFT coefficients:", coefs)
final_evals = circle_ifft(coefs, Domain)
print("IFFT reconstructed evaluations:", final_evals)
print("Matches original:", final_evals == evals)


# # Verify FFT with basis functions

# In[13]:


def check_with_basis(coefs, circle_domain):
    """Re-evaluate coefficients using basis functions."""
    arrC = [FieldElement(c) for c in coefs]
    out = []
    for pt in circle_domain:
        acc = FieldElement(0)
        for i, cval in enumerate(arrC):
            acc += cval * basis_func(pt, i)
        out.append(acc.value)
    return out

def basis_func(pt, idx):
    """Compute basis function for index idx at point."""
    b = FieldElement(1)
    bitpos = 0
    tmp = idx
    while tmp > 0:
        if bitpos == 0:
            if tmp & 1:
                b = b * pt.y
        else:
            if tmp & 1:
                b = b * repeated_vn(pt.x, bitpos)
        tmp >>= 1
        bitpos += 1
    return b

def repeated_vn(x, k):
    """Repeatedly apply v_n(x) = 2*x^2 - 1, k-1 times."""
    tmp = x
    for _ in range(k - 1):
        tmp = tmp.square().double() - FieldElement.one()
    return tmp

reeval = check_with_basis(coefs, Domain)
print("Re-evaluated with basis:", reeval)
print("Matches original:", reeval == evals)


# # Benchmarking FFT

import time

Q = CirclePoint(FieldElement(24), FieldElement(18))
Domain = generate_standard_position_coset(Q, G3)

evals = [1, 2, 4, 8, 16, 1, 2, 4]

# warm up
for _ in range(100):
    circle_fft(evals, Domain)

repeat = 10000

start = time.perf_counter()

for _ in range(repeat):
    coefs = circle_fft(evals, Domain)

end = time.perf_counter()

total_time = end - start
avg_time = total_time / repeat

print(f"FFT size       : {len(evals)}")
print(f"Iterations     : {repeat}")
print(f"Total time     : {total_time:.6f} s")
print(f"Average time   : {avg_time * 1e6:.3f} us")