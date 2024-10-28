# coor2graph

This software supports the build system CMake.

**Compilers**
1. gfortran(<=14.2.1 20240805)
2. ifort(<=2021.13.0 20240602)

### CMake

The CMake build system requires both make and CMake to be installed, the latter has to be version 3.30.3 or newer.

Build `coor2graph` with CMake works with the following chain od commands:

```bash
cmake -B build -DBUILD_STATIC=ON
make -C build
```

coord2graph call the `networkx` python library. Required depedencies:
* numpy
* pandas
* networkx
* matplotlib

### How to use?

### Basic use
```coor2graph --input <file>.gro --pure y --rcut # --cdistance #  --graph <prefix> --measure <function>```


