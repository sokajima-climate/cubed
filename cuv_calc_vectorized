def cuv_calc_vectorized(lon, lat, ut, vt):
    """
    Vectorized version of cuv_calc using NumPy operations.
    This can be faster than the loop version but may use more memory.
    
    Parameters:
    lon -- longitude array
    lat -- latitude array
    ut -- u component array (2D: lat x lon)
    vt -- v component array (2D: lat x lon)
    
    Returns:
    cvort, cuv, tp -- curvature vorticity, curvature, and topology arrays
    """
    # Get dimensions
    nlat = len(lat)
    nlon = len(lon)
    
    # Initialize output arrays
    cvort = np.zeros((nlat, nlon), dtype=np.float64)
    cuv = np.zeros((nlat, nlon), dtype=np.float64)
    tp = np.zeros((nlat, nlon), dtype=np.float64)
    
    # Constants
    RE = 6.371e6
    pi = np.pi
    ddeg = 2.0 * pi * RE / 360.0
    
    # Create arrays for derivatives
    dux = np.zeros((nlat, nlon))
    duy = np.zeros((nlat, nlon))
    dvx = np.zeros((nlat, nlon))
    dvy = np.zeros((nlat, nlon))
    
    # Interior points for y-derivatives
    duy[1:-1, :] = (ut[2:, :] - ut[:-2, :]) / np.tile((lat[2:] - lat[:-2])[:, np.newaxis], (1, nlon)) / ddeg
    dvy[1:-1, :] = (vt[2:, :] - vt[:-2, :]) / np.tile((lat[2:] - lat[:-2])[:, np.newaxis], (1, nlon)) / ddeg
    
    # Compute cosine of latitude for interior points
    cos_lat = np.cos(lat[1:-1] * pi / 180.0)[:, np.newaxis]
    
    # Interior points for x-derivatives (fixing the broadcasting issue)
    for ilat in range(1, nlat-1):
        cos_lat_val = math.cos(lat[ilat] * pi / 180.0)
        
        # Interior points
        dux[ilat, 1:-1] = (ut[ilat, 2:] - ut[ilat, :-2]) / (lon[2:] - lon[:-2]) / ddeg / cos_lat_val
        dvx[ilat, 1:-1] = (vt[ilat, 2:] - vt[ilat, :-2]) / (lon[2:] - lon[:-2]) / ddeg / cos_lat_val
        
        # Wrap-around for longitude at boundaries
        dux[ilat, 0] = (ut[ilat, 1] - ut[ilat, nlon-1]) / ((lon[1] - lon[nlon-1]) + 360.0) / ddeg / cos_lat_val
        dvx[ilat, 0] = (vt[ilat, 1] - vt[ilat, nlon-1]) / ((lon[1] - lon[nlon-1]) + 360.0) / ddeg / cos_lat_val
        
        dux[ilat, nlon-1] = (ut[ilat, 0] - ut[ilat, nlon-2]) / ((lon[0] - lon[nlon-2]) + 360.0) / ddeg / cos_lat_val
        dvx[ilat, nlon-1] = (vt[ilat, 0] - vt[ilat, nlon-2]) / ((lon[0] - lon[nlon-2]) + 360.0) / ddeg / cos_lat_val
        
    # Compute topology (this can be vectorized for interior points, but boundaries need special handling)
    for ilat in range(1, nlat-1):
        for ilon in range(1, nlon-1):
            tp[ilat, ilon] = ut[ilat, ilon] + ut[ilat, ilon-1] + ut[ilat, ilon+1] + ut[ilat-1, ilon] + ut[ilat+1, ilon]
        
        # Boundary conditions
        tp[ilat, 0] = ut[ilat, 0] + ut[ilat, nlon-1] + ut[ilat, 1] + ut[ilat-1, 0] + ut[ilat+1, 0]
        tp[ilat, nlon-1] = ut[ilat, nlon-1] + ut[ilat, nlon-2] + ut[ilat, 0] + ut[ilat-1, nlon-1] + ut[ilat+1, nlon-1]
    
    # Compute curvature values
    uu = ut[1:-1, :]**2
    vv = vt[1:-1, :]**2
    uv = ut[1:-1, :] * vt[1:-1, :]
    
    # Compute numerator for all points
    numerator = (-uv * dux[1:-1, :] + uu * dvx[1:-1, :] - 
                 vv * duy[1:-1, :] + uv * dvy[1:-1, :])
    
    # Compute denominator and avoid division by zero
    denominator = uu + vv
    mask = denominator > 0
    
    # Apply formula where denominator is non-zero
    cvort[1:-1, :] = np.where(mask, numerator / denominator, 0)
    cuv[1:-1, :] = np.where(mask, numerator / denominator**1.5, 0)
    
    return cvort, cuv, tp

    
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
cvort, cuv, tp = cuv_calc_vectorized(lon, lat, ut, vt)
"""
