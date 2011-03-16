'''
View GIS Elevation data for the ocean floor. And the land too, I guess... Why would you want to do that anyway?

'''
import scipy.io as sio
import compbio.config as config
import compbio.utils.colors as mycolors
from matplotlib import pyplot as plt
from numpy import * 
def view(lat_range = [9.,15.], lon_range = [140., 150.] ):
  '''View GIS Eleveation data for the ocean floor.

kwargs:
 lat_range  [9,13]
 lon_range  [140,144]
'''
  nc = sio.netcdf_file(config.dataPath('random/ocean_topo.nc'))
  xvals = nc.variables['x'].data
  yvals = nc.variables['y'].data
  idxs = [nonzero(logical_and(greater(xvals,lon_range[0]),
                              less(xvals, lon_range[1])))[0],
          nonzero(logical_and(greater(yvals,lat_range[0]),
                              less(yvals, lat_range[1])))[0]]
  cmap = mycolors.blackbody()                            
  arr = array(nc.variables['z'].data)
  subimg = arr[:,idxs[0]][idxs[1],:]
  plt.imshow(subimg[::-1], cmap = cmap, interpolation = 'nearest')
  
