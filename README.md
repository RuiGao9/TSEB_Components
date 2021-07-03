# TSEB_Components
Digital number of important components from the TSEB model: Rn, H, LE, G, and LE (Canopy)

## Brief introduction
The main idea is get the digital number of the main components from the TSEB results within the footprint. Net radiation, sensible heat flux, latent heat flux, soil heat flux, and canopy latent heat flux.

## Parameters in this script
- `footprint`: footprint image.
- `tseb_pt_1`: a main multiple-image result from the TSEB model.
- `tseb_pt_2`: the ancillary image result from the TSEB model.
- `dir_out`: a folder directory the temporary results are saved.
- `n_rn`: the number layer of the net radiation.
- `n_h`: the number layer of the net radiation.
- `n_le`: the number layer of the latent heat flux.
- `n_g`: the number layer of the soil surface heat flux.
- `n_t_et`: the number layer of the ratio between canopy latent heat flux and total latent heat flux.
- `pixel_size`: the pixel size of the image, e.g., 3.6 meter by 3.6 meter.

## Contact
[Rui Gao](https://www.researchgate.net/profile/Rui-Gao-55)<br>
Rui.Gao@usu.edu
