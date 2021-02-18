# Geometrically Non-Linear Hexahedral Finite Element
Total Lagrangian approach is used to account for the geometrically non-linear behaviour in the code provided here.

## Input Parameters
| Variable | Matrix Sizing | Format | Description |
| --- | --- | --- | --- |
| `E` | real scalar | 1e9 | Young's modulus in kPa |
| `v` | real scalar | 0 | Poisson's ratio |
| `coord` | 3n * 3 | [X<sub>1</sub> Y<sub>1</sub> Z<sub>1</sub>] | Coordinates of the nodes |
| `etopol` | n * 8 | [n<sub>1</sub> n<sub>2</sub> n<sub>3</sub> n<sub>4</sub> n<sub>5</sub> n<sub>6</sub> n<sub>7</sub> n<sub>8</sub>] | Element topology |
| `bc` | nbc * 2 | [i u<sub>i</sub>] | Boundary conditions: degrees of freedom, <math>i</math> having defined displacements of u<sub>i</sub> |
| `f` | 3n * 1 | [f<sub>1</sub> f<sub>2</sub> ... f<sub>3n</sub>]<sup>T</sup> | Total applied force at each degree of freedom in kN |

where `n` is the number of elements, `nbc` is the number of nodes that have Dirichlet boundary condition (displacement is known). Note that the element topology has to be written in the format shown below
<p align="center">
  <img src="../Linear/8_Hexahedral_Nodal_Numbering.png" width="350" title="hover text">
</p>

## Output Parameters
| Variable | Matrix Sizing | Format | Description |
| --- | --- | --- | --- |
| `d` | 3n * 1 | [d<sub>1</sub> d<sub>2</sub> ... d<sub>3n</sub>]<sup>T</sup> | Displacement at each degree of freedoms in m |

## Files
1. TLHexahedral.m - The main MATLAB files to run the Total-Lagrangian Geometrically Non-Linear Hexahedral FE analysis.
2. TLFE.m - Contains the functions to be ran at the element-level looping.
3. GNLcantilever_endload.m - A function file to provide a sample of input for the Geometrically Non-Linear Hexahedral FE analysis.

## Algorithm Process
<p align="center">
  <img src="Geometrically_Non-Linear_Hexahedral_Process.png" width="1654" title="hover text">
</p>
Note that the number in top left of each procedure box relates to TLHexahedral.m while the number in the top right of each procedure box relates to the line number in the TLHE.m.
