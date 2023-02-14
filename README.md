# Periodic Boundary Conditions on Regular Hyperbolic Tessellations for Fast-Converging Accurate Bulk Simulations
Fabian R. Lux<sup>1</sup> and Emil Prodan<sup>1</sup><br />
<sup>1</sup>*Department of Physics and, Department of Mathematical Sciences, Yeshiva University, New York, NY 10016, USA*

DOI: [XY](http://dx.doi.org/?)

This is a repository for code, related to our paper on the "Periodic Boundary Conditions on Regular Hyperbolic Tessellations for Fast-Converging Accurate Bulk Simulations".

## Abstract

TODO

## About the code

### Requirements

The main code base relies on `C++17` and the libraries

- `Boost`, https://www.boost.org/
- `Eigen3`, https://eigen.tuxfamily.org/index.php?title=Main_Page

The environment variables `Boost_DIR` and `Eigen3_DIR` need to point to the root directory of the respective libaries. For this one can modify and use the configure scripts `./configure.bat` on Windows and `./configure.sh` on Linux or Mac.

### Installing the code

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