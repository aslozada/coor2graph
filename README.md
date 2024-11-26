# coor2graph

This software supports the build system CMake.

**Compilers**
1. gfortran(<=14.2.1 20240805)
2. ifort(<=2021.13.0 20240602)

### CMake

The CMake build system requires both make and CMake to be installed, the latter has to be version 3.30.3 or newer.

Build `coor2graph` with CMake works with the following chain of commands:

```bash
cmake -B build -DBUILD_STATIC=ON
make -C build
```

coord2graph uses the `networkx` python library. Required dependencies:
* numpy
* pandas
* networkx
* matplotlib

#### Alternative
Create an enviroment
* In ubuntu
```conda create --name myenvironment python=3```
```conda activate myenviroment```

to install required packages (if necessary): 
```conda install <package-name>[=version]```




### How to use?

### Basic use
```coor2graph --input <file>.gro --rcut <#> --pair <sym1> <sym2> <distance>  --graph <prefix> --measure <networkx measure>```

#### Example
![lattice_model](https://github.com/user-attachments/assets/d167f95e-ea9a-4b8c-b741-99ee7b053892)


