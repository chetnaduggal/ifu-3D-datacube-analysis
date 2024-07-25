# Integral-field (3D) datacube analysis for a quasar-host merger system

![SINFONI H-band datacube](https://github.com/chetnaduggal/IFU-datacube-analysis/assets/67710398/85aeb7fe-adb3-45ae-8e2d-e0d86ea1ad63)

The code and results in this repository are part of an ongoing study of the host galaxy of quasar 3C 297. The observations were obtained (Cycle 97 program 097.B-0452(A); PI G. Tremblay) from the now-decomissioned SINFONI instrument at the ESO Very Large Telescope (VLT), Chile. 
The data are now in public domain and the raw files are available [here](http://archive.eso.org/wdb/wdb/eso/eso_archive_main/query?prog_id=097.B-0452(A)&max_rows_returned=10000).    

----------------------------

### Data

The H-band (1.45−1.85 μm) data were reduced using the standard ESO-SINFONI pipeline. After sky subtraction and wavelength calibration, the final co-added data cube has a spatial scale of 0.125′′ per pix. The field of view is 8′′ × 8′′. The angular resolution is ∼0.4′′ as sampled by the point-spread function (PSF) from the standard star observations. At the target redshift, this translates to spatial resolution of ∼3 kpc. 

The galaxy has a highly perturbed morphology, with a prominent arc towards north of the core, which is cospatial with the jet axis (green contours from radio overlay). This suggests its a site of jet-ISM interaction from the expanding radio lobes. 

A quick look at the line profiles over the field of view showed that the _broad line emission is limited to the nuclear region_-- which were found to be due to kpc-scale outflowing (ionized) gas. The northern arc-like feature only shows narrow line emission and is likely a site of active star formation.  

Thus, some of the analyses are divided spatially into the "nuclear" and "northern arc" regions.

### Codes

- `pipelined_FOV_fitting_&_maps.ipynb` and `pipelined_FOV_fitting_&_maps.py` |  spectrum modelling
- `3c297_broad&narrow_line_maps.ipynb` |
- `spaxelwise_Nuclear7x7_NorthernArc.ipynb` |  Script to calculate ...
- `3c297_overlays.ipynb` |
- `WHAN_diagnostic_spaxelwise.ipynb` |
- `telluricOH_fitting.ipynb` |
- `generate_error_cube.ipynb` |


----------------------------

![](https://github.com/chetnaduggal/ifu-3D-datacube-analysis/blob/f7bdc12496aa005c6ef48d8a0624e7f3d6bae30b/int-regions-spectra.pdf)






