# cubed
A Python implementation of the subroutine for calculating curvature vorticity and curvature.
https://sites.google.com/view/cubed/

## Logo
![cubed](https://github.com/user-attachments/assets/f98caa43-f8b0-4292-9534-e56b564f0706)

## Features
- Computes curvature vorticity and curvature from wind field data.
- Handles periodic boundary conditions for longitude.
- Provides a structured approach to numerical differentiation.

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
