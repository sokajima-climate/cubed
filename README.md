# cubed
A Python implementation of the subroutine for calculating curvature vorticity and curvature.
https://sites.google.com/view/cubed/

### Logo
<img src="https://github.com/user-attachments/assets/5a5a4292-cb49-4b01-b888-4d3b6845b5bb" width="320">\

## Features
- Computes curvature vorticity and curvature from wind field data.
- Handles periodic boundary conditions for longitude.
- Provides a structured approach to numerical differentiation.

## Description of the methodology
(The following is based on Okajima et al. (2021). You may also refer to the publication.)
Vorticity can be decomposed locally into shear and curvature terms as follows (e.g., Holton 2004);\
```math
\zeta=-\frac{\partial V}{\partial n}+\frac{V}{R_s}
```
where V denotes scalar wind speed, n the direction perpendicular to the flow, and Rs the radius of curvature. The first and second terms of the RHS in the above equation represent shear vorticity and curvature vorticity, respectively, which can be calculated as\
```math
-\frac{\partial V}{\partial n}=-\frac{1}{V^2}(-u v u_x-v^2 v_x+u^2 u_y+u v v_y)
```
```math
\frac{V}{R_s}=\frac{1}{V^2}(-u v u_x+u^2 v_x-v^2 u_y+u v v_y)
```
where u and v denote local zonal and meridional wind velocities, respectively, and subscripts x and y partial derivatives in the zonal and meridional directions, respectively. 

The curvature, defined as\
```math
\kappa_2=\frac{1}{R_s}=\frac{1}{V^3}(-u v u_x+u^2 v_x-v^2 u_y+u v v_y)
```
can be derived from the definition of curvature of a two-dimensional curve implicitly represented by\
```math
\psi(x,y)=c
```
```math
\kappa_2=\frac{\begin{vmatrix}
  \psi_{xx} & \psi_{xy} & \psi_x \\
  \psi_{yx} & \psi_{yy} & \psi_y \\
  \psi_x & \psi_y & 0
\end{vmatrix}}{(\psi_x^2+\psi_y^2)^{3/2}}
```
The curvature or curvature vorticity enables us to circumvent the difficulties in determining areas of cyclonic and anticyclonic circulations, because it is free from shear vorticity and thus extracts vortex circulation with a certain radius.

The curvature-based methodology enables us to evaluate contributions from cyclonic and anticyclonic eddies separately to Eulerian statistics, by accumulating instantaneous contributions only at grid points where cyclonic or anticyclonic curvature is observed.

##### Examples:

- Eulerian eddy statistics (such as v't'850, v'v'300) and eddy feedback forcing of westerly wind acceleration: Okajima, S., H. Nakamura, Y. Kaspi: Cyclonic and anticyclonic contributions to atmospheric energetics, Scientific Reports, 11, 13202, 2021.
- Air-sea turbulent heat flux, precipitation, and associated hydrological cycle: Okajima, S., Nakamura, H., Spengler, T. (2024). Midlatitude oceanic fronts strengthen the hydrological cycle between cyclones and anticyclones. Geophysical Research Letters, 51, e2023GL106187.

- Here is a press release for the curvature-based methodology: https://www.u-tokyo.ac.jp/focus/en/press/z0508_00181.html

##### "Materialization" of cyclonic eddies based on CUBED\
<img src="https://github.com/user-attachments/assets/d52a92cc-1bdc-4b0e-9f27-85b4f79c712b" width="320">\
<img src="https://github.com/user-attachments/assets/14c2425a-4527-4650-9316-adcde1da3c20" width="480">\
You can find more detailed information about the "materialization" in the article in Tenki (in Japanese).
- 岡島悟 (2022): ３次元格子データのmaterialization―3Dプリンタによる”物質化”―，「天気」Vol.69, No.8.PDF

## files
### cuv_calc_xarray.py (Recommended)
Script to calculate curvature and curvature vorticity based on the curvature-based methodology CUBED, using xarray. 
### cuv_calc.py
Script to calculate curvature and curvature vorticity based on the curvature-based methodology CUBED, using numpy. 
### cuv_calc_numba.py
Script to calculate curvature and curvature vorticity based on the curvature-based methodology CUBED, using numpy and jit (numba).
(** I have not check the functionality of cuv_calc_numba.py yet)

## Requirements
- Python 3.x
- NumPy
- Xarray (optional)

## Reference
- Okajima, S., Nakamura, H. & Kaspi, Y. Cyclonic and anticyclonic contributions to atmospheric energetics. Sci Rep 11, 13202 (2021). https://doi.org/10.1038/s41598-021-92548-7
- Okajima, S., Nakamura, H., & Kaspi, Y. (2024). Anticyclonic suppression of the North Pacific transient eddy activity in midwinter. Geophysical Research Letters, 51, e2023GL106932. https://doi.org/10.1029/2023GL106932
- Okajima, S., Nakamura, H., & Kaspi, Y. (2024). Reply to comment by Chang on “anticyclonic suppression of the North Pacific transient eddy activity in midwinter”. Geophysical Research Letters, 51, e2024GL111599. https://doi.org/10.1029/2024GL111599

## Installation
Clone this repository and install dependencies if necessary:
```bash
git clone https://github.com/sokajima-climate/cubed.git
cd cubed
pip install numpy
