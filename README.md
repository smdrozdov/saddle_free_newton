# Saddle-free Newton.

Assume f is a C2-smooth function from $U \in \mathbb{R}^n$ to $\mathbb{R}$. Given that $\argmin(f) \in \int(U)$ find $\argmin(f).$

This is a common task usually solved by Gradient descent, and by Newton method if $f$ is convex.
However, if $f$ is not convex and has saddle points, these saddle points are attractors for Gradient descent process, which very often stops method from converging to the real minimum.
Moreover, if $n$ grows, saddle points become the only blocker for Gradient descent. The probability of getting blocked in the local minimum drops as exp(-n).

In 2014 Dauphin et al. introduced and proved a novel method, that they called "Saddle-free Newton" [1]. It is designed to overcome saddle points of arbitrary function. 

Though this method is proven powerful (See figure 4 in chapter 7 in [1]), a public implementation is lacking. This package is an attempt to fill this gap.


References:
[1] Yann N. Dauphin et al. (2014) Identifying and attacking the saddle point problem in high-dimensional non-convex optimization, https://arxiv.org/abs/1406.2572
