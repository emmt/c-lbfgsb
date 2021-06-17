# Wrappers for the L-BFGS-B algorithm

This repository provides thin C and [Yorick](http://yorick.github.com/)
wrappers for using L-BFGS-B algorithm by Ciyou Zhu, Richard Byrd, Jorge Nocedal
and Jose Luis Morales.

L-BFGS-B is a numerical method to minimize a multi-variate differentiable
objective function possibly under separble bound constraints.  The user is
required to provide the bounds and an initial solution and to compute the
objective function and its gradient. L-BFGS-B is a quasi-Newton method with low
memory requirements (*"L"* is for *"Limited memory"*) and which can optionally
take into account separable bound constraints (the final "*B*") on the
variables.  To determine efficient search directions, L-BFGS-B approximates the
Hessian of the objective function by a a limited memory version of the
Broyden-Fletcher-Goldfarb-Shanno model (*"BFGS"* for short).


## Installation

### To install the library

Go to [`src`](./src) directory.  Edit `Makefile` to match your settings.
Then build and install the library and header file:

```sh
make
make install PREFIX=...
```

where `...` denotes a directory where to install files.

To test the software:

```sh
make check
```


### To install the Yorick plug-in

Follow instructions in [yorick/README.md](yorick/README.md) file, it is not
needed to install the library first.


## License(s)

The FORTRAN code of L-BFGS-B (version 3.0) by Ciyou Zhu, Richard Byrd, Jorge
Nocedal and Jose Luis Morales is in directory [`lbfgsb-3.0`](./lbfgsb-3.0].
This code has been released under the [“*New BSD
License*”](./lbfgsb-3.0/License.txt) (aka “*Modified BSD License*” or
“*3-clause license*”) and is freely available
[hure](http://users.iems.northwestern.edu/~nocedal/lbfgsb.html).

The C and Yorick parts (in directories [`src`](./src) and [`yorick`](./yorick))
are released under the [Simplified BSD 3-Clause License](./LICENSE.md).
