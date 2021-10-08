def DigitalResults(footprint, tseb_r_1, tseb_r_2, dir_out, 
                   lai_image,
                   n_rn, n_h, n_le, n_g, n_t_et, pixel_size,
                   upper_boundary, lower_boundary, delete_tmp_files="Yes"):
    '''
    parameters:
    footprint: directory of the footprint image.
    tseb_r_1: directory of the TSEB result: multiple-layer image.
    tseb_r_2: directory of the TSEB ancillary result: multiple-layer image.
    dir_out: directory of the outputs from this scripts. they are transform results and they can be deleted.
    lai_image: directory of the LAI image.
    n_rn, n_h, n_le, n_g: the layer number of the net radiation, sensible heat flux, latent heat flux, and soil surface heat flux.
    n_t_et: the layer number of the ratio between canopy latent heat flux and total latent heat flux.
    pixel_size: the pixel size (e.g., 3.6 meter by 3.6 meter).
    upper_boundary: the upper threshold for all fluxes at one pixel which does not make sense, e.g., 10,000 W/m2 for LE at one pixel.
    lower_boundary: the lower threshold for all fluxes at one pixel which does not make sense, e.g., -1,500 W/m2 for LE at one pixel.
    delete_tmp_files: Default is "Yes", and this means the temporary (middle products) files will be deleted at the end. Any other input (string) results in saving the temporary files.

    return:
    Net radiation, sensible heat flux, latent heat flux, soil surface heat flux, canopy latent heat flux, and LAI within the footprint area.
    '''
    import arcpy
    import numpy as np
    import pandas as pd
    import os
    
    cellsize = str(pixel_size)+" "+str(pixel_size)

    arcpy.Resample_management(in_raster=footprint, 
                              out_raster=dir_out+"\\footprint_resample.tif", 
                              cell_size=cellsize, 
                              resampling_type="BILINEAR")
    Grid_Describe = arcpy.Describe(tseb_r_1)
    Grid_Extent = Grid_Describe.extent
    extent = "{} {} {} {}".format(Grid_Extent.XMin, Grid_Extent.YMin, Grid_Extent.XMax, Grid_Extent.YMax)
    
    arcpy.Clip_management(in_raster=dir_out+"\\footprint_resample.tif", 
                          rectangle=extent, 
                          out_raster=dir_out+"\\footprint_clip.tif", 
                          in_template_dataset=tseb_r_1, 
                          nodata_value="0.000000e+00", 
                          clipping_geometry="NONE", maintain_clipping_extent="MAINTAIN_EXTENT")

    raster_footprint = arcpy.RasterToNumPyArray(dir_out+"\\footprint_clip.tif", nodata_to_value=-9999)
    raster_footprint[raster_footprint>1] = np.nan
    raster_footprint[raster_footprint<0] = np.nan
    raster_tseb = arcpy.RasterToNumPyArray(tseb_r_1, nodata_to_value=-9999)
    raster_tseb_ancillary = arcpy.RasterToNumPyArray(tseb_r_2, nodata_to_value=np.nan)
    raster_lai = arcpy.RasterToNumPyArray(lai_image, nodata_to_value=np.nan)

    raster_rn = raster_tseb[n_rn,:,:]
    raster_rn[raster_rn>upper_boundary] = np.nan
    raster_rn[raster_rn<lower_boundary] = np.nan
    out_rn = raster_rn*raster_footprint
    
    raster_h = raster_tseb[n_h,:,:]
    raster_h[raster_h>upper_boundary] = np.nan
    raster_h[raster_h<lower_boundary] = np.nan
    out_h = raster_h*raster_footprint
    
    raster_le = raster_tseb[n_le,:,:]
    raster_le[raster_le>upper_boundary] = np.nan
    raster_le[raster_le<lower_boundary] = np.nan
    out_le = raster_le*raster_footprint
    
    raster_g = raster_tseb[n_g,:,:]
    raster_g[raster_g>upper_boundary] = np.nan
    raster_g[raster_g<lower_boundary] = np.nan
    out_g = raster_g*raster_footprint
    
    raster_t = raster_tseb_ancillary[n_t_et,:,:]
    raster_t[raster_t>1] = np.nan
    raster_t[raster_t<0] = np.nan
    out_t = out_le*raster_t
    
    raster_tet = raster_footprint*0+1
    out_tet = raster_tet*raster_tseb_ancillary[n_t_et,:,:]
    out_tet = np.nanmean(out_tet)    
    
    raster_lai[raster_lai>5] = np.nan
    raster_lai[raster_lai<0] = np.nan
    out_lai = raster_lai*(raster_footprint*0+1)

    out_rn = np.nansum(out_rn)
    out_h = np.nansum(out_h)
    out_le = np.nansum(out_le)
    out_g = np.nansum(out_g)
    out_t = np.nansum(out_t)
    out_lai = np.nanmean(out_lai)
    print("Rn - Net radiation:",round(out_rn,3))
    print("H - Sensible heat flux:",round(out_h,3))
    print("LE - Latent heat flux:",round(out_le,3))
    print("G - Soil surface heat flux:",round(out_g,3))
    print("T - Canopy latent heat flux:",round(out_t,3))
    print("ET partitioning:",round(out_tet,3))
    print("LAI:",round(out_lai,3))
    
    if delete_tmp_files == "Yes":
        os.remove(dir_out+"\\footprint_resample.tif")
        os.remove(dir_out+"\\footprint_clip.tif")
    else:
        print("Temporary files are saved in the output folder.")
        
    return(out_rn, out_h, out_le, out_g, out_t, out_tet, out_lai)
