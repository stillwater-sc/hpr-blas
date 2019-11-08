# applications/polynomial

Polynomial evaluation examples

# Challange Problem

For the more typical case where we are using rounded real numbers for everything, the solution is described in many publications by Ulrich Kulisch. The most recent one is in his Computer Arithmetic and Validity, 2nd edition, Section 9.6.2 "Accurate evaluation of polynomials". I reference this in the posits4.pdf document that's on posithub.org. A challenge problem from Kulisch:

P[t_] := 8118t^4 -11482t^3 +t^2 +5741t-2030

evaluated at a t value near 1/√2, where it has an exact root. The polynomial is cast into the form of a lower-triangular system:

x_1 =    8118
x_2 = x1 * t - 11482
x_3 = x2 * t + 1
x_4 = x3 * t + 5741
x_5 = x4 * t - 2030

where x_5 is the desired polynomial value. Solving the system, refining the answer using a residual computed with the quire, and perhaps iterating, nails the correct value within 0.5 ULP.
