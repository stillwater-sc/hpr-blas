# examples/numerical

Numerical Analysis Examples

# Newton-Raphson

A Newton-Raphson algorithm using valids.

It is a demonstration of the benefits of custom posit/valid configurations and the error-free linear algebra.

In numerical analysis, Newton's method (also known as the Newton–Raphson method), named after Isaac Newton and Joseph Raphson, is a method for finding successively better approximations to the roots (or zeroes) of a real-valued function. It is one example of a root-finding algorithm.

{\displaystyle x:f(x)=0\,.} x:f(x)=0\,.

The Newton–Raphson method in one variable is implemented as follows:

The method starts with a function f defined over the real numbers x, the function's derivative f ′, and an initial guess x0 for a root of the function f. If the function satisfies the assumptions made in the derivation of the formula and the initial guess is close, then a better approximation x1 is

{\displaystyle x_{1}=x_{0}-{\frac {f(x_{0})}{f'(x_{0})}}\,.} x_{1}=x_{0}-{\frac {f(x_{0})}{f'(x_{0})}}\,.
Geometrically, (x1, 0) is the intersection of the x-axis and the tangent of the graph of f at (x0, f (x0)).

The process is repeated as

{\displaystyle x_{n+1}=x_{n}-{\frac {f(x_{n})}{f'(x_{n})}}\,} x_{n+1}=x_{n}-{\frac {f(x_{n})}{f'(x_{n})}}\,
until a sufficiently accurate value is reached.
