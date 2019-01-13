# Saddle-free Newton.

Assume f is a C2-smooth function from U \in Rn to R. Given that $argmin(f) \in \int(U)$ find $argmin(f).$

This is a common task that is usually solved by Gradient descent, and by Newton method if $f$ is convex.
However, if $f$ is not convex and has saddle points they are attractors for Gradient descent, which very often stops it from converging.
Moreover, if $n$ grows, saddle points become a main obstacle for Gradient descent, because probability of a function to have a local
minimum drops exponentially quick.

In 2014 Dauphin et al. introduced and proved a novel method, that they called "Saddle-free Newton". It is designed to overcome
saddle points for an arbitrary function. 

Though this method is proven powerful (See figure 4 in chapter 7 in https://arxiv.org/abs/1406.2572), a public implementation is lacking. This package is an attempt to fill this gap.
