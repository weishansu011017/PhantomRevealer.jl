# PhantomRevealer

A Julia interface for analyzing dump files from the Phantom Smoothed Particle Hydrodynamics code.

Most of the analysis is based on SPH interpolation. Check out Price2010 for further information.

## Installation

### 1. Install Julia

Julia is required before installing this package. You can check this website: [Julia Downloads](https://julialang.org/downloads/) for further information.

After installing Julia, type 

```bash
julia
```
to start the julia REPL. 

### 2. Install this package
In the Julia REPL, activate the Pkg module to install the package:

```julia
using Pkg
```

Next, install this package directly from the Git repository:
```julia
Pkg.add(url="https://github.com/weishansu011017/PhantomRevealer.git")
```

## Usage

(Not done yet, you can check out the files under `example` to learn it.)

## References

Price, Daniel J. "Smoothed Particle Hydrodynamics and Magnetohydrodynamics" (doi: 10.1016/j.jcp.2010.12.011)

Phantom homepage: [Phantom SPH](https://phantomsph.github.io)

Sarracen documentation: [Sarracen Documentation](https://sarracen.readthedocs.io/en/latest/)