# Converging Periodic Boundary Conditions and Detection of Topological Gaps on Regular Hyperbolic Tessellations
Fabian R. Lux<sup>1</sup> and Emil Prodan<sup>1</sup><br />
<sup>1</sup>*Department of Physics and, Department of Mathematical Sciences, Yeshiva University, New York, NY 10016, USA*

DOI: [10.13140/RG.2.2.14897.04968](http://dx.doi.org/10.13140/RG.2.2.14897.04968)

This is a repository for code, related to our paper on the "Converging Periodic Boundary Conditions and Detection of Topological Gaps on Regular Hyperbolic Tessellations".

## Abstract

Tessellations of the hyperbolic spaces by regular polygons are becoming popular because they support discrete quantum and classical models displaying unique spectral and topological characteristics. Resolving the true bulk spectra and the thermodynamic response functions of these models requires converging periodic boundary conditions and our work delivers a practical solution for this open problem on generic {p,q}-tessellations. This enables us to identify the true spectral gaps of bulk Hamiltonians and, as an application, we construct all but one topological models that deliver the topological gaps predicted by the K-theory of the lattices. We demonstrate the emergence of the expected topological spectral flows whenever two such bulk models are deformed into each other and, additionally, we prove the emergence of topological channels whenever a soft physical interface is created between different topological classes of Hamiltonians.

## Introduction

For detailed information, we refer to our [paper](http://dx.doi.org/10.13140/RG.2.2.14897.04968). Very briefly, our goal is to study the group

$$\Delta_{\lbrace p, q \rbrace}^+ = \langle A, B | A^p = B^q= (AB)^2 \rangle,$$ 

which is known as the proper triangle group. For hyperbolic crystals this group acts as a group of symmetries and contains infinitely many group elements. In a hyperbolic $\Delta_{\lbrace p, q \rbrace}^+$ Cayley crystal, the lattice sites can be unambigously labelled by elements of the group. In other words, the action of the group leaves no lattice site fixed.

Since the group is infinite, some practial choices have to be made in order to simulate the physics of such a Cayley crystal on a computer. We explore two alternatives: open and periodic boundary conditions.

Both cases make use of a faithful representation

$$ \sigma \colon \Delta_{\lbrace p, q \rbrace}^+ \to \mathrm{GL}(3, \mathbb{Z}[\xi_{2pq}]),$$

where $\xi_n = 2 \cos(2\pi/n)$ is an algebraic number and $\mathbb{Z}[\xi_{n}]$ is the ring of integer polynomials in $\xi_{n}$.
There is a minimal, monic polynomial of degree $\varphi(n)/2$ which has $\xi_n$ as its root and where $\varphi(n)$ is Euler's totient function.
This means  $\xi_n^{\varphi(n)/2}$ can be expressed as a linear combination of lower order terms. We construct the minimal polynomial and realize $ \mathbb{Z}[\xi_{n}]$ as the vector space $\mathbb{Z}^{\varphi(n)/2}$. A rule of multiplication $\mathbb{Z}^{\varphi(n)/2}\times \mathbb{Z}^{\varphi(n)/2} \to \mathbb{Z}^{\varphi(n)/2}$ is inherited from $ \mathbb{Z}[\xi_{n}]$.

With open boundary conditions, we grow the group "outwards" by the alternating completion of A and B cycles up to a predefined order. We do so by using $\sigma$ to supply the group multiplication and use it as a means of distinguishes the group elements from each other. This is possible, since $\sigma$ is faithful.

For periodic boundary conditions, $\sigma(\Delta_{\lbrace p, q \rbrace}^+) \mod N$ is a finite group, where the modulo operation acts on each individual coefficient of $ \mathbb{Z}[\xi_{2pq}]$. We then still generate the group elements as before, but there is in principal no need for truncating the process, since it will terminate automatically as soon as no new group elements are being generated anymore.

This repository contains a code which implements this procedure. The group generation algorithm is very efficient and is capable of generating thousands of group elements in a fraction of a second using a modern CPU.

## About the code

### Overview

The code base is structured in three segments, each with a different purpose and programming language:

1. `./src/` and `./include`: C++ group generation algorithms (gives list of group elements, and the right-regular representation), exact diagonilization of very large matrices in Fortran.
2. `./py`: Python scripts to set up different Hamiltonians from the right-regular representation (bulk-spectra, projection operators, Y-junction)
3. `./wls`: WolframLanguage scripts for the creation of the ring-reduction and the group generators for given $p$ and $q$. Additional scripts for plotting on the hyperbolic disk such as the local density of states.

The docstrings of the respective source files and the `README.md` files contain more information.

### Installing the C++ and Fortran code

#### Requirements

The main code base relies on `C++17` and the libraries

- (optional) `OpenBLAS64`, http://www.openblas.net/ (with `long` datatype for array indexing instead of `int`)
- (optional) `Boost`, https://www.boost.org/ 
- (optional) `Eigen3`, https://eigen.tuxfamily.org/index.php?title=Main_Page

The environment variables `Boost_DIR`, `Eigen3_DIR` and `OpenBLAS64_DIR` need to point to the root directory of the respective libraries. 

#### Compilation

The code can be installed by 
```bash
mkdir build
cd build
cmake ..
make
```

To test if the code works as expected, one can run 

```bash
make test
```

Note that the test functionality has very little coverage of the code base and is very rudimentary. Only some imports are tested.

Finally one needs to set the environment variable `HYPERBOLIC_BUILD` with the absolute path of the build directory that was just created and the variable `HYPERBOLIC_DIR` to the directory which contains the file `group_specs.inp` and will store group-related data.

### How to use the code

#### The input file

Create the text file `group_specs.inp` in the directory specified by `HYPERBOLIC_DIR`. A valid example could be the following:

```
p 5
q 4
N 2
```

If $N>0$, the proper $p$-$q$-triangle group will generated modulo $N$. If $N<0$, the first $N$ cycles of the of the triangle group are completed by starting with a $p$-cycle.

There is a unique file name prefix corresponding to each set of $p$, $q$ and $N$. If $N>0$, outpot file names are of the form

`$HYPERBOLIC_DIR/<p>_<q>_modulo_<N>.*`

where `<x>` needs to be replaced with the actual value of `x` and `*` denotes a possible file name extension. If $N<0$, the corresponding string is

`$HYPERBOLIC_DIR/<p>_<q>_open_<-N>.*`

In order for the code to function, it requires information about the ring reduction and the group generators for given $p$ and $q$. For convenience those have already been computed for many cases of practical relevance.

#### Generating the group

The executable `generate_group` will read `group_specs.inp` and generate the group elements. The result is stored in a `.words` file. For example, a typical `.words` file could look like the following:

```
159 AABBAABBAAB
158 ABBAAABBAAAAB
153 AAABBAABBAA
152 AAABBAABBA
150 AABBBAABBBA
148 AABBBAABBA
...
```

where the integer describes the index of the group element, keeping track of the order in which group elements where discovered. The integer is then followed by a string which describes a word in the proper-triangle group. `A` is the generator of `p`-cycles, while `B` is the generator of `q`-cycles. Once created, the matrix representation of the group elements can be reconstructed from the information in this file.

#### Generating the right regular representation

The executable `generate_reg` will read `group_specs.inp` and the corresponding `.words` file which needs to exist beforehand. From this information it will construct the right-regular representation of `A`, `B` and their corresponding inverse operations `iA`, `iB`.  For each of these operators `X` an output file 

`$HYPERBOLIC_DIR/<p>_<q>_modulo_<N>.<X>`

or 

`$HYPERBOLIC_DIR/<p>_<q>_open_<-N>.<X>`

is created. Let $g_i$ and $g_j$ denote group elements with index $i$ and $j$ in the `.words` file.
A typical example outpout can then look like the following:

```
18 17
19 18
4 19
54 20
-1 21
...
```

where a row with non-negative entries  $i,j$ means that

$$ g_i = g_j X^{-1}$$

if $i=-1$, then $g_j X^{-1}$ is not contained in the `.words` file, which can happen for the case of open boundary conditions. Setting $i=-1$ in this case makes it possible to catch this exception.
