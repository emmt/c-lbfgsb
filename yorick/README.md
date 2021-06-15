Yorick plugin for L-BFGS-B algorithm
====================================


Usage
-----

Minimizing a multi-variate function under box constraints with L-BFGS-B
algorithm is best described by examples.  See functions `lbfgsb-test1`,
`lbfgsb-test2`, and `lbfgsb-test3` in [`lbfgsb-tests.i`](./lbfgsb-tests.i).


Installation
------------

In short, building and installing the plug-in can be as quick as:
```sh
cd $BUILD_DIR
$SRC_DIR/configure
make
make install
```
where `$BUILD_DIR` is the build directory (at your convenience) and `$SRC_DIR`
is the source directory of the plug-in code.  The build and source directories
can be the same in which case, call `./configure` to configure for building.
Call `$SRC_DIR/configure -h` for a description of available options.

If the plug-in has been properly installed, it is sufficient to use any of its
functions to automatically load the plug-in.  You may force the loading of the
plug-in by something like:
```cpp
#include "lbfgsb.i"
```
or
```cpp
require, "lbfgsb.i";
```
in your code.

More detailled installation explanations are given below.

1. You must have Yorick installed on your machine.

2. Unpack the plug-in code somewhere.

3. Configure for compilation.  There are two possibilities:

   For an **in-place build**, go to the source directory of the plug-in code
   and run the configuration script:
   ```sh
   cd SRC_DIR
   ./configure
   ```
   To see the configuration options, type:
   ```sh
   ./configure --help
   ```

   To compile in a **different build directory**, say `$BUILD_DIR`, create the
   build directory, go to the build directory, and run the configuration
   script:
   ```sh
   mkdir -p $BUILD_DIR
   cd $BUILD_DIR
   $SRC_DIR/configure
   ```
   where `$SRC_DIR` is the path to the source directory of the plug-in code.
   To see the configuration options, type:
   ```sh
   $SRC_DIR/configure --help
   ```

4. Compile the code:
   ```sh
   make clean
   make
   ```

5. Install the plug-in in Yorick directories:
   ```sh
   make install
   ```
