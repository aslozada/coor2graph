# coor2graph
This software build an adjacency matrix from coordinate file. The coordinate file can be obtained from a molecular dynamics or stochastic simulation. 
Using the adjacency matrix and the Python library Networkx, this program enables the analysis of graphs or complex networks resulting from the molecular simulation.
Periodic boundary effects can be evaluated.

* Options
- [x] pair
- [ ] triad
- [ ] 2D lattice
- [ ] 3D lattice
- [ ] Non-cubic pbc

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
* scipy (for large matrix)

#### Alternative
Create an enviroment
* In ubuntu
```conda create --name myenvironment python=3```
```conda activate myenvironment```

to install required packages (if necessary): 
```conda install <package-name>[=version]```

### GRO file format
```
MD of 9 ising spins, t= 0.0
    9
    1ISING  SP1    1   0.000   0.000   0.000  0.0000  0.0000  0.0000
    2ISING  SP1    2   0.100   0.000   0.000  0.0000  0.0000  0.0000
    3ISING  SP1    3   0.200   0.000   0.000  0.0000  0.0000  0.0000
    4ISING  SP1    4   0.000   0.100   0.000  0.0000  0.0000  0.0000
    5ISING  SP1    5   0.100   0.100   0.000  0.0000  0.0000  0.0000
    6ISING  SP1    6   0.200   0.100   0.000  0.0000  0.0000  0.0000
    7ISING  SP1    7   0.000   0.200   0.000  0.0000  0.0000  0.0000
    8ISING  SP1    8   0.100   0.200   0.000  0.0000  0.0000  0.0000
    9ISING  SP1    9   0.200   0.200   0.000  0.0000  0.0000  0.0000
   0.40000   0.40000   0.40000
```
Positions (x,y,z) and box(3) units in nm

### How to use?

### Basic use
```coor2graph --input <file>.gro --rcut <#> --pbc <y|n> --pair <sym1> <sym2> <distance>  --graph <prefix> --measure <networkx measure> --frequency <y|n> <#|-1>```

### Currently avaliable graph measures (TODO) (Update from networkx.org)

* Degree centrality
  
   degree | in_degree | out_degree
* Eigenvector centrality
  
   eigenvector | katz
* Closeness centrality
  
   closeness |incremental_closeness
* Current flow closeness
  
   flow_centrality | information
* Betweenness

  betweenness | betweeness_subset | edge_betweenness
* Current flow betweenness

  flow_betweenness edge_flow_betweenness
* Others

harmonic_centrality | dispersion | laplacian

#### Example (Ising model)
First-neighbor interaction

![lattice_model](https://github.com/user-attachments/assets/d167f95e-ea9a-4b8c-b741-99ee7b053892)

#### Closeness 
![prefix_1](https://github.com/user-attachments/assets/56d2b87c-0d38-4e54-bbb4-07c52ec8ce49)



#### Example (3-sites water model)
![graphs](https://github.com/user-attachments/assets/c34049fb-dfe1-4a94-82e4-0a49088c6b3c)

#### Properties along a trajectory
![histogram](https://github.com/user-attachments/assets/bc0c0655-ea7d-4162-ba48-a46713f79d65)



### PBC effects on centrality degree
(with pbc conditions)
![wpbc_1](https://github.com/user-attachments/assets/a0b49eff-c0eb-4476-ae05-77bfb62f7a5d)
(without pbc conditions)
![woutpbc_1](https://github.com/user-attachments/assets/cd325594-ee55-47e8-9158-71a081fcd6ef)


.
