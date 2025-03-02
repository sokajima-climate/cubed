def cuv_calc_xarray_vectorized(ut, vt, lon=None, lat=None):
    """
    Vectorized version of cuv_calc for xarray DataArrays.
    This handles dimension mismatches carefully.
    
    Parameters and returns same as cuv_calc_xarray.
    """
    # Get dimensions from input data
    if lat is None:
        lat = ut.lat.values
    if lon is None:
        lon = ut.lon.values
    
    nlat = len(lat)
    nlon = len(lon)
    
    # Convert to numpy arrays for calculations
    ut_np = ut.values
    vt_np = vt.values
    
    # Create output numpy arrays
    ddd1_np = np.zeros_like(ut_np)
    ddd2_np = np.zeros_like(ut_np)
    dtopo_np = np.zeros_like(ut_np)
    
    # Temporary arrays for derivatives
    dux_np = np.zeros_like(ut_np)
    duy_np = np.zeros_like(ut_np)
    dvx_np = np.zeros_like(ut_np)
    dvy_np = np.zeros_like(ut_np)
    
    # Constants
    RE = 6.371e6
    pi = np.pi
    ddeg = 2.0 * pi * RE / 360.0
    
    # Pre-compute values needed for derivatives
    lat_diff = np.zeros(nlat)
    for i in range(1, nlat-1):
        lat_diff[i] = lat[i+1] - lat[i-1]
    
    lon_diff = np.zeros(nlon)
    for i in range(1, nlon-1):
        lon_diff[i] = lon[i+1] - lon[i-1]
    
    cos_lat_values = np.cos(np.radians(lat))
    
    # Compute derivatives - carefully to avoid dimension mismatch
    for ilat in range(1, nlat-1):
        duy_np[ilat, :] = (ut_np[ilat+1, :] - ut_np[ilat-1, :]) / lat_diff[ilat] / ddeg
        dvy_np[ilat, :] = (vt_np[ilat+1, :] - vt_np[ilat-1, :]) / lat_diff[ilat] / ddeg
        
        cos_lat = cos_lat_values[ilat]
        for ilon in range(1, nlon-1):
            dux_np[ilat, ilon] = (ut_np[ilat, ilon+1] - ut_np[ilat, ilon-1]) / lon_diff[ilon] / ddeg / cos_lat
            dvx_np[ilat, ilon] = (vt_np[ilat, ilon+1] - vt_np[ilat, ilon-1]) / lon_diff[ilon] / ddeg / cos_lat
            
        # Handle boundary conditions with special care for longitude wrap-around
        dux_np[ilat, 0] = (ut_np[ilat, 1] - ut_np[ilat, nlon-1]) / ((lon[1] - lon[nlon-1]) + 360.0) / ddeg / cos_lat
        dvx_np[ilat, 0] = (vt_np[ilat, 1] - vt_np[ilat, nlon-1]) / ((lon[1] - lon[nlon-1]) + 360.0) / ddeg / cos_lat
        
        dux_np[ilat, nlon-1] = (ut_np[ilat, 0] - ut_np[ilat, nlon-2]) / ((lon[0] - lon[nlon-2]) + 360.0) / ddeg / cos_lat
        dvx_np[ilat, nlon-1] = (vt_np[ilat, 0] - vt_np[ilat, nlon-2]) / ((lon[0] - lon[nlon-2]) + 360.0) / ddeg / cos_lat
    
    # Compute dtopo
    for ilat in range(1, nlat-1):
        for ilon in range(1, nlon-1):
            dtopo_np[ilat, ilon] = ut_np[ilat, ilon] + ut_np[ilat, ilon-1] + ut_np[ilat, ilon+1] + ut_np[ilat-1, ilon] + ut_np[ilat+1, ilon]
        
        # Boundary conditions
        dtopo_np[ilat, 0] = ut_np[ilat, 0] + ut_np[ilat, nlon-1] + ut_np[ilat, 1] + ut_np[ilat-1, 0] + ut_np[ilat+1, 0]
        dtopo_np[ilat, nlon-1] = ut_np[ilat, nlon-1] + ut_np[ilat, nlon-2] + ut_np[ilat, 0] + ut_np[ilat-1, nlon-1] + ut_np[ilat+1, nlon-1]
    
    # Vectorized computation of curvature
    for ilat in range(1, nlat-1):
        for ilon in range(nlon):
            uu = ut_np[ilat, ilon]**2
            vv = vt_np[ilat, ilon]**2
            uv = ut_np[ilat, ilon] * vt_np[ilat, ilon]
            
            denominator = uu + vv
            if denominator > 0:
                numerator = (-uv * dux_np[ilat, ilon] + uu * dvx_np[ilat, ilon] - 
                            vv * duy_np[ilat, ilon] + uv * dvy_np[ilat, ilon])
                
                ddd1_np[ilat, ilon] = numerator / denominator
                ddd2_np[ilat, ilon] = numerator / (denominator**1.5)
    
    # Convert back to xarray DataArrays
    cvort = xr.DataArray(
        ddd1_np,
        coords=ut.coords,
        dims=ut.dims,
        name='curvature_vorticity'
    )
    
    cuv = xr.DataArray(
        ddd2_np,
        coords=ut.coords,
        dims=ut.dims,
        name='curvature'
    )

    tp = xr.DataArray(
        dtopo_np,
        coords=ut.coords,
        dims=ut.dims,
        name='topography'
    )
    
    return cvort, cuv, tp

# Usage example
"""
# import xarray as xr
# cvort, cuv, tp = cuv_calc_xarray_vectorized(u, v)
"""
