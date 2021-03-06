# TSEB_Components
Digital number of important components from the TSEB model: Rn, H, LE, G, LE (Canopy), and T/ET.
Digital number of LAI within the footprint area.

## Brief introduction
The main idea is get the digital value of the LAI and the main components from the TSEB results within the footprint. Net radiation, sensible heat flux, latent heat flux, soil heat flux, canopy latent heat flux, ET partitioning, and LAI.

## Parameters in this script
- `footprint`: footprint image.
- `tseb_pt_1`: a main multiple-image result from the TSEB model.
- `tseb_pt_2`: the ancillary image result from the TSEB model.
- `dir_out`: a folder directory the temporary results are saved.
- `lai_image`: the LAI image in order to get the LAI value within the footprint area. This value is averaged within the footprint area.
- `n_rn`: the number layer of the net radiation.
- `n_h`: the number layer of the net radiation.
- `n_le`: the number layer of the latent heat flux.
- `n_g`: the number layer of the soil surface heat flux.
- `n_t_et`: the number layer of the ratio between canopy latent heat flux and total latent heat flux.
- `pixel_size`: the pixel size of the image, e.g., 3.6 meter by 3.6 meter.
- `upper_boundary`: the upper threshold for all fluxes at one pixel which does not make sense, e.g., 10,000 W/m2 for LE at one pixel.
- `lower_boundary`: the lower threshold for all fluxes at one pixel which does not make sense, e.g., -1,500 W/m2 for LE at one pixel.
- `delete_tmp_files`: the default string is "Yes", and this means the temporary (middle products) files will be deleted at the end. Any other input (string) results in saving the temporary files

## Returned digital values
- `out_rn`: Net radiation within the footprint area.
- `out_h`: Sensible heat flux within the footprint area.
- `out_le`: Latent heat flux within the footprint area.
- `out_g`: Ground heat flux within the footprint area.
- `out_t`: Canopy latent heat flux within the footprint area.
- `out_tet`: ET partitioning within the footprint area.
- `out_lai`: Averaged LAI value within the footprint area.

## How to cite:
Please cite the paper below when you are using this script for your paper work.<br>
[Evapotranspiration partitioning assessment using a machine-learning-based leaf area index and the two-source energy balance model with sUAV information](https://www.researchgate.net/publication/350820947_Evapotranspiration_partitioning_assessment_using_a_machine-learning-based_leaf_area_index_and_the_two-source_energy_balance_model_with_sUAV_information)

## Contact
[Rui Gao](https://www.researchgate.net/profile/Rui-Gao-55): Rui.Gao@usu.edu | Rui.Ray.Gao@gmail.com
