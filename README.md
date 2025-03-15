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
Vorticity can be decomposed locally into shear and curvature terms as follows (e.g., Holton 2004)
```math
\zeta=-\frac{\partial V}{\partial n}+\frac{V}{R_s}
```
where V denotes scalar wind speed, n the direction perpendicular to the flow, and Rs the radius of curvature. The first and second terms of the RHS in the above equation represent shear vorticity and curvature vorticity, respectively, which can be calculated as
```math
-\frac{\partial V}{\partial n}=-\frac{1}{V^2}(-u v u_x-v^2 v_x+u^2 u_y+u v v_y)
```
```math
\frac{V}{R_s}=\frac{1}{V^2}(-u v u_x+u^2 v_x-v^2 u_y+u v v_y)
```
where u and v denote local zonal and meridional wind velocities, respectively, and subscripts x and y partial derivatives in the zonal and meridional directions, respectively. 

The curvature, defined as
```math
\kappa_2=\frac{1}{R_s}=\frac{1}{V^3}(-u v u_x+u^2 v_x-v^2 u_y+u v v_y)
```
can be derived from the definition of curvature of a two-dimensional curve implicitly represented by
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

## Q&A:
##### Is it possible to calculate the curvature of grid points directly based on the height field? 
Yes, it is possible to calculate the curvature of grid points based on the height field, by assuming geostrophic equilibrium. Under the assumption of geostrophic equilibrium, the stream function is proportional to the geopotential; ψ = gz/f. Using this relationship, you can calculate the curvature from the geopotential height field.
(The equation of the curvature from the stream function appears to be the same as the definition of curvature of an implicitly represented two-dimensional curve in the above equation.)
##### Did you filter the curvature while calculating the curvature from horizontal winds? 
No, we did not filter the curvature. This is because we wanted to obtain cyclonic & anticyclonic domains that correspond directly to those that appear on weather charts. (Those systems are identified based on unfiltered SLP, for example) If you apply filtering, we will obtain cyclonic/anticyclonic regions that are artificially too symmetric because of filtering. This will make it unable to investigate cyclone-anticyclone asymmetry.
##### I am still not sure which threshold to adopt.
Your curvature threshold will depend on what you are to focus on. Our methodology evaluates the local curvature (or equivalently, curvature radius) of the flow, which represents the flow topology in a physically straightforward fashion.
If you want to partition all grid points into cyclonic and anticyclonic regions, a threshold of zero curvature would work well. If you are to obtain "transition zones" between cyclonic and anticyclonic regions, a non-zero curvature threshold would be effective. It will be according to the features you want to focus on - for transient eddies, the curvature radius of 2,500km or 3,000km would work as in our paper. Ultimately, it will correspond to what you regard as "cyclonic and anticyclonic regions", which will depend on purposes and contexts. 

## Examples:

- Eulerian eddy statistics (such as v't'850, v'v'300) and eddy feedback forcing of westerly wind acceleration: Okajima, S., H. Nakamura, Y. Kaspi: Cyclonic and anticyclonic contributions to atmospheric energetics, Scientific Reports, 11, 13202, 2021.
- Air-sea turbulent heat flux, precipitation, and associated hydrological cycle: Okajima, S., Nakamura, H., Spengler, T. (2024). Midlatitude oceanic fronts strengthen the hydrological cycle between cyclones and anticyclones. Geophysical Research Letters, 51, e2023GL106187.

- Here is a press release for the curvature-based methodology: https://www.u-tokyo.ac.jp/focus/en/press/z0508_00181.html

### "Materialization" of cyclonic eddies based on CUBED\
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
