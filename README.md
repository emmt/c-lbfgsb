# A C wrapper for L-BFGS-B algorithm

This repository provides thin wrappers for using L-BFGS-B algorithm in C and
Yorick.


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


## License(s)

The FORTRAN code of L-BFGS-B (version 3.0) by Ciyou Zhu, Richard Byrd, Jorge
Nocedal and Jose Luis Morales is in directory [`lbfgsb-3.0`](./lbfgsb-3.0].
This code has been released under the [“*New BSD
License*”](./lbfgsb-3.0/License.txt) (aka “*Modified BSD License*” or
“*3-clause license*”) and is freely available
[hure](http://users.iems.northwestern.edu/~nocedal/lbfgsb.html).

The C and Yorick parts (in directories [`src`](./src) and [`yorick`](./yorick))
are released under the [Simplified BSD 3-Clause License](./LICENSE.md).
