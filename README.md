# ijvectorstates

Quad-state triangular bitmap matrix used to cache contribution states of a set of M N-dimensional integer vectors representing electron selections.

In short, a program uses a series of M integer vectors, each of dimension N; any vectors differing by 3 or more elements lead to a zero Hamiltonian matrix element.  Checking two integer vectors I and J repeatedly has been identified as a possible optimization point for the code.

For M vectors, there are M!/(2!\*(M-2)!) = M\*(M+1)/2 pair-wise combinations.  Comparing a vector to itself implicitly yields zero differences -- we don't need to cache that, so we can eliminate M of the combinations, leaving M\*(M-1)/2.  (That's also the dimension of the upper (or lower) triangle of a matrix minus the diagonal.)

A Fortran LOGICAL defaults to being the same size as the INTEGER type, but can be altered to a single byte using `LOGICAL*1` as the type.  This module is designed to use memory as efficiently as possible, representing each of the M\*(M-1)/2 states as 2 bits (valued [0,3]).  Thus, it packs the state of 4 elements into each byte.

## Using an IJVectorStates

After declaring a variable (let's call it `IJS`) it must be initialized:

```
Integer, Parameter              :: M = 50
Integer                         :: IJ(N, M), I_idx, J_idx, S
Type(IJVectorStatesType)        :: IJS

Call IJVectorStatesInit(IJS, M)
```

Where M is the number of integer vectors.  The instance starts off with each element set to `IJVectorStateUnchecked` to indicate it has not been initialized; at this point, any quick queries of the cached state will return `IJVectorStateTrue` (1) for like-index and `IJVectorStateUnchecked` (3) for all other indice combinations:

```
Write(*,*) IJVectorStatesCheckQuick(IJS, 1, 1)
```

### Initialization

The IJVectorStates object can be used in two ways:

1. Pass all vector combinations to `IJVectorStatesCheck()` to populate the cache, then use `IJVectorStatesCheckQuick()` throughout the rest of the code.
2. Use `IJVectorStatesCheck()` everywhere in the program to populate the cache only as needed.

Scheme 1 will usually be the most performant and also reduces the number of arguments passed each time a state check is performed (3 arguments for `IJVectorStatesCheckQuick()` versus 6 for `IJVectorStatesCheck()`):

```
Do I_idx = 1, M - 1
    Do J_idx = I_Idx + 1, M
        IJVectorStatesCheck(IJS, N, IJ(:,I_idx), I_idx, IJ(:,J_idx), J_idx)
    End Do
End Do
```

### Checking cached state

To check the cached state only, the vectors do not need to be passed, just the vector indices:

```
S = IJVectorStatesCheckQuick(IJS, 1, 1)
```

## Example

The `ftest.F90` program is an example that uses the IJVectorStates module.  The included `Makefile` has compiler-specific includes (`Makefile.inc.gfortran` and `Makefile.inc.ifort`) that set compiler and flags; copy one of these as `Makefile.inc` or create a symbolic link named `Makefile.inc` to setup the build.  Then type `make`.

The `ftest` executable produced should output the following when run:

```
$ ./ftest
 Initial state array

1 3 3 3 3 3 3 3 3 3 3 3
3 1 3 3 3 3 3 3 3 3 3 3
3 3 1 3 3 3 3 3 3 3 3 3
3 3 3 1 3 3 3 3 3 3 3 3
3 3 3 3 1 3 3 3 3 3 3 3
3 3 3 3 3 1 3 3 3 3 3 3
3 3 3 3 3 3 1 3 3 3 3 3
3 3 3 3 3 3 3 1 3 3 3 3
3 3 3 3 3 3 3 3 1 3 3 3
3 3 3 3 3 3 3 3 3 1 3 3
3 3 3 3 3 3 3 3 3 3 1 3
3 3 3 3 3 3 3 3 3 3 3 1

 First pass with IJVectorStatesCheck()

1 1 0 1 1 0 0 0 0 0 0 0
1 1 0 0 0 0 0 0 0 0 0 0
0 0 1 0 0 0 0 0 0 0 0 0
1 0 0 1 0 0 0 0 0 0 0 0
1 0 0 0 1 0 0 0 0 0 0 0
0 0 0 0 0 1 0 0 0 0 0 0
0 0 0 0 0 0 1 1 0 1 1 0
0 0 0 0 0 0 1 1 0 0 1 0
0 0 0 0 0 0 0 0 1 0 0 0
0 0 0 0 0 0 1 0 0 1 0 0
0 0 0 0 0 0 1 1 0 0 1 0
0 0 0 0 0 0 0 0 0 0 0 1

 All cached values, hopefully (no 3s)

1 1 0 1 1 0 0 0 0 0 0 0
1 1 0 0 0 0 0 0 0 0 0 0
0 0 1 0 0 0 0 0 0 0 0 0
1 0 0 1 0 0 0 0 0 0 0 0
1 0 0 0 1 0 0 0 0 0 0 0
0 0 0 0 0 1 0 0 0 0 0 0
0 0 0 0 0 0 1 1 0 1 1 0
0 0 0 0 0 0 1 1 0 0 1 0
0 0 0 0 0 0 0 0 1 0 0 0
0 0 0 0 0 0 1 0 0 1 0 0
0 0 0 0 0 0 1 1 0 0 1 0
0 0 0 0 0 0 0 0 0 0 0 1
```

The matrix form being output is symmetrical, thank goodness!
