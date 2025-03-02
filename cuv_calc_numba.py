def cuv_calc_numba(lon, lat, ut, vt):
    """
    Python version of the cuv_calc Fortran subroutine with transposed array indexing.
    
    Parameters:
    -----------
    lon, lat : numpy.ndarray
        1D arrays containing longitude and latitude values
    ut, vt : numpy.ndarray
        2D arrays of size (nlat, nlon) containing u and v components
        
    Returns:
    --------
    tuple of numpy.ndarray
        (cvort, cuv, tp) arrays of size (nlat, nlon)
    """
    # lat = ut["lat"].values
    # lon = ut["lon"].values
    nlat = len(lat)
    nlon = len(lon)

    # Initialize output arrays
    cvort = np.zeros((nlat, nlon), dtype=np.float64)
    cuv = np.zeros((nlat, nlon), dtype=np.float64)
    tp = np.zeros((nlat, nlon), dtype=np.float64)
    
    # Initialize intermediate arrays
    dux = np.zeros((nlat, nlon), dtype=np.float64)
    duy = np.zeros((nlat, nlon), dtype=np.float64)
    dvx = np.zeros((nlat, nlon), dtype=np.float64)
    dvy = np.zeros((nlat, nlon), dtype=np.float64)
    
    # Constants
    RE = 6.371e6
    pi = math.pi
    dist = 400.0
    ddeg = 2.0 * pi * RE / 360.0
    
    # Compute derivatives
    _compute_derivatives(nlon, nlat, lon, lat, ut, vt, dux, duy, dvx, dvy, tp, pi, ddeg)
    
    # Compute curvature vortex and curvature
    _compute_curvature(nlon, nlat, ut, vt, dux, duy, dvx, dvy, cvort, cuv)
       
    # Note: Smoothing operations are commented out in original code
    # Implement smoothing functions if needed
    
    return cvort, cuv, tp

@jit(nopython=True, parallel=True)
def _compute_derivatives(nlon, nlat, lon, lat, ut, vt, dux, duy, dvx, dvy, tp, pi, ddeg):
    """
    Compute derivatives of ut and vt with transposed indexing
    """
    # Use numba's parallel range for parallelization (equivalent to OpenMP in Fortran)
    for ilat in prange(1, nlat-1):
        # Compute y-derivatives for all longitudes
        for ilon in range(nlon):
            duy[ilat, ilon] = (ut[ilat+1, ilon] - ut[ilat-1, ilon]) / (lat[ilat+1] - lat[ilat-1]) / ddeg
            dvy[ilat, ilon] = (vt[ilat+1, ilon] - vt[ilat-1, ilon]) / (lat[ilat+1] - lat[ilat-1]) / ddeg
        
        # Compute x-derivatives for interior longitudes
        for ilon in range(1, nlon-1):
            cos_lat = math.cos(lat[ilat] * pi / 180.0)
            dux[ilat, ilon] = (ut[ilat, ilon+1] - ut[ilat, ilon-1]) / (lon[ilon+1] - lon[ilon-1]) / ddeg / cos_lat
            dvx[ilat, ilon] = (vt[ilat, ilon+1] - vt[ilat, ilon-1]) / (lon[ilon+1] - lon[ilon-1]) / ddeg / cos_lat
            tp[ilat, ilon] = ut[ilat, ilon] + ut[ilat, ilon-1] + ut[ilat, ilon+1] + ut[ilat-1, ilon] + ut[ilat+1, ilon]
        
        # Handle boundary conditions (wrapping around)
        cos_lat = math.cos(lat[ilat] * pi / 180.0)
        
        # At ilon = 0
        dux[ilat, 0] = (ut[ilat, 1] - ut[ilat, nlon-1]) / ((lon[1] - lon[nlon-1]) + 360.0) / ddeg / cos_lat
        dvx[ilat, 0] = (vt[ilat, 1] - vt[ilat, nlon-1]) / ((lon[1] - lon[nlon-1]) + 360.0) / ddeg / cos_lat
        
        # At ilon = nlon-1
        dux[ilat, nlon-1] = (ut[ilat, 0] - ut[ilat, nlon-2]) / ((lon[0] - lon[nlon-2]) + 360.0) / ddeg / cos_lat
        dvx[ilat, nlon-1] = (vt[ilat, 0] - vt[ilat, nlon-2]) / ((lon[0] - lon[nlon-2]) + 360.0) / ddeg / cos_lat
        
        # dtopo at boundaries
        tp[ilat, 0] = ut[ilat, 0] + ut[ilat, nlon-1] + ut[ilat, 1] + ut[ilat-1, 0] + ut[ilat+1, 0]
        tp[ilat, nlon-1] = ut[ilat, nlon-1] + ut[ilat, nlon-2] + ut[ilat, 0] + ut[ilat-1, nlon-1] + ut[ilat+1, nlon-1]

@jit(nopython=True, parallel=True)
def _compute_curvature(nlon, nlat, ut, vt, dux, duy, dvx, dvy, ddd1, ddd2):
    """
    Compute curvature values with transposed indexing
    """
    for ilat in prange(1, nlat-1):
        for ilon in range(nlon):
            uu = ut[ilat, ilon]**2
            vv = vt[ilat, ilon]**2
            uv = ut[ilat, ilon] * vt[ilat, ilon]
            
            # Avoid division by zero
            denominator = uu + vv
            if denominator > 0:
                numerator = (-uv * dux[ilat, ilon] + uu * dvx[ilat, ilon] - 
                             vv * duy[ilat, ilon] + uv * dvy[ilat, ilon])
                
                ddd1[ilat, ilon] = numerator / denominator
                ddd2[ilat, ilon] = numerator / (denominator**1.5)
            # If denominator is zero, ddd1 and ddd2 remain zero as initialized

# Example of how to use the function:
"""
import numpy as np

# Define input parameters
nlon, nlat = 360, 180
lon = np.linspace(0, 359, nlon)
lat = np.linspace(-89.5, 89.5, nlat)
undef = -9999.0
ut = np.random.rand(nlat, nlon)  # Example u-component data
vt = np.random.rand(nlat, nlon)  # Example v-component data

# Run the function
cvort, cuv, tp = cuv_calc_numba(lon, lat, ut, vt)
"""
