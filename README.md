# Saddle-free Newton.

Assume f is a C2-smooth function from U \in Rn to R. Given that min(f) on U is greater that $-\infinity$ find $argmin(f).$

This is a common task that is usually solved by Gradient Descent, and by Newton method if f is convex.
However, if f is not convex and has saddle points they are attractors for Gradient descent, which very often stops it from converging.
Moreover, if n grows saddle points become a main obstacle for Gradient descent, because probability of a function to have a local
minimum drops exponentially quick.

In 2014 Bengio et al. introduced and proved a novel method, that they called "Saddle-free Newton". It is designed to overcome
saddle points for an arbitrary function. 

Though this method is very powerful, a public implementation is lacking. This repository is an attempt to fill this gap.
