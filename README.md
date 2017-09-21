# 3DEPG
3D extended phase graph header library

Inspired by the work of Aaron Hess and Matthew Robson (http://onlinelibrary.wiley.com/doi/10.1002/mrm.26213), this small C++ header library enables extended phase graph simulations with spoiling in three dimensions. Given the exponential explosion in spin states that 3D spoiling can produce, this library employs a few tricks to try to keep things efficient:

1. It uses a specific layout of the spin states in memory, along with the MKL, to make excitation and relaxation operations very fast.

2. It uses a sparse representation of the 3D lattice of spin states, storing only occupied states and creating new entries only when needed.

I make no promises this is really the fastest/most-compact/smartest way to do it, but its what I found to be practical. If you have suggestions for how to make it more efficient, please make a pull request!
