def DigitalResults(footprint, tseb_pt_1, tseb_pt_2, dir_out, 
                   n_rn, n_h, n_le, n_g, n_t_et, pixel_size):
'''
parameters:
footprint: directory of the footprint image.
tseb_pt_1: directory of the TSEB result.
tseb_pt_2: directory of the TSEB ancillary result.
dir_out: directory of the outputs from this scripts. they are transform results and they can be deleted.
n_rn, n_h, n_le, n_g: the layer number of the net radiation, sensible heat flux, latent heat flux, and soil surface heat flux.
n_t_et: the layer number of the ratio between canopy latent heat flux and total latent heat flux.
pixel_size: the pixel size (e.g., 3.6 meter by 3.6 meter)

return:
Net radiation, sensible heat flux, latent heat flux, soil surface heat flux, and canopy latent heat flux within the footprint gained from TSEB model.
'''
    import arcpy
    import numpy as np
    import pandas as pd
    
    cellsize = str(pixel_size)+" "+str(pixel_size)

    arcpy.Resample_management(in_raster=footprint, 
                              out_raster=dir_out+"\\footprint_resample.tif", 
                              cell_size=cellsize, 
                              resampling_type="BILINEAR")
    Grid_Describe = arcpy.Describe(tseb_pt_1)
    Grid_Extent = Grid_Describe.extent
    extent = "{} {} {} {}".format(Grid_Extent.XMin, Grid_Extent.YMin, Grid_Extent.XMax, Grid_Extent.YMax)
    
    arcpy.Clip_management(in_raster=dir_out+"\\footprint_resample.tif", 
                          rectangle=extent, 
                          out_raster=dir_out+"\\footprint_clip.tif", 
                          in_template_dataset=tseb_pt_1, 
                          nodata_value="0.000000e+00", 
                          clipping_geometry="NONE", maintain_clipping_extent="MAINTAIN_EXTENT")

    raster_footprint = arcpy.RasterToNumPyArray(dir_out+"\\footprint_clip.tif", nodata_to_value=-9999)
    raster_footprint[raster_footprint>1] = 0
    raster_footprint[raster_footprint<0] = 0
    raster_tseb = arcpy.RasterToNumPyArray(tseb_pt_1, nodata_to_value=-9999)
    raster_tseb_ancillary = arcpy.RasterToNumPyArray(tseb_pt_2, nodata_to_value=-9999)

    out_rn = raster_tseb[n_rn,:,:]*raster_footprint
    out_rn[out_rn>1500] = 0
    out_rn[out_rn<-100] = 0
    out_h = raster_tseb[n_h,:,:]*raster_footprint
    out_h[out_h>1500] = 0
    out_h[out_h<-100] = 0
    out_le = raster_tseb[n_le,:,:]*raster_footprint
    out_le[out_le>1500] = 0
    out_le[out_le<-100] = 0
    out_g = raster_tseb[n_g,:,:]*raster_footprint
    out_g[out_g>1500] = 0
    out_g[out_g<-100] = 0
    out_t = out_le*raster_tseb_ancillary[n_t_et,:,:]
    out_t[out_t>1500] = 0
    out_t[out_t<-100] = 0

    out_rn = np.nansum(out_rn)
    out_h = np.nansum(out_h)
    out_le = np.nansum(out_le)
    out_g = np.nansum(out_g)
    out_t = np.nansum(out_t)
    print("Net radiation:",round(out_rn,3))
    print("Sensible heat flux:",round(out_h,3))
    print("Latent heat flux:",round(out_le,3))
    print("Soil surface heat flux:",round(out_g,3))
    print("Canopy latent heat flux:",round(out_t,3))
    return(out_rn, out_h, out_le, out_g, out_t)