#!/usr/bin/env python3
from netCDF4 import Dataset
import numpy as np

def gaussian2d(sigma, x0, y0, dx, dy, nx, ny):
    """
    Returns a 2D Gaussian of the form exp(-1*i((x-x0)^2 + (y-y0)^2)/(2*sigma^2))
    """
    g2d = np.zeros((nx, ny), dtype=np.float32)
    for i in range(nx):
        for j in range(ny):
            g2d[i,j] = np.exp(-1.0*(((i*dx-x0)/(2.0*sigma))**2 + 
                                    ((j*dy-y0)/(2.0*sigma))**2))
    return g2d 

if __name__ == "__main__":
    sigma = 0.1
    dx = 0.1
    dy = 0.1 
    nx = 100
    ny = 100
    x0 = 5.0
    y0 = 5.0
    ncfile = Dataset('testpy.nc', mode='w', format='NETCDF3_CLASSIC')
    nx_dim = ncfile.createDimension("nx", 101)
    ny_dim = ncfile.createDimension("ny", 101)
    h = ncfile.createVariable('h', np.float32, ('nx','ny'))
    h[:,:] = gaussian2d(sigma, x0, y0, dx, dy, nx, ny)
    ncfile.close()
