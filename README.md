# Integral-field (3D) datacube analysis for a quasar-host merger system

![movie](https://github.com/chetnaduggal/IFU-datacube-analysis/assets/67710398/85aeb7fe-adb3-45ae-8e2d-e0d86ea1ad63)

The codes and results in this repository are part of an ongoing study of the host galaxy of quasar 3C 297. The observations were obtained (Cycle 97 program 097.B-0452A; PI G. Tremblay) from the now-decomissioned SINFONI instrument at the ESO Very Large Telescope (VLT), Chile. 
The data are now in public domain and the raw files are available [here](http://archive.eso.org/wdb/wdb/eso/eso_archive_main/query?prog_id=097.B-0452(A)&max_rows_returned=10000).    

----------------------------

### Data

The H-band (1.45−1.85 μm) data were reduced using the standard ESO-SINFONI pipeline. After sky subtraction and wavelength calibration, the final co-added data cube has a spatial scale of 0.125′′ per pix. The field of view is 8′′ × 8′′. The angular resolution is ∼0.4′′ as sampled by the point-spread function (PSF) from the standard star observations. At the target redshift, this translates to spatial resolution of ∼3 kpc. 

The galaxy has a highly perturbed morphology, with a prominent arc towards north of the core. A quick look at the line profiles over the field of view (FOV) showed that the _broad line emission is limited to the nuclear region_-- which also shows blue-to-red shift roughly around the twin AGN (blue '+' marks). The northern arc-like feature only shows narrow line emission and bright UV emission in the HST-SINFONI overlay.

Thus, some of the analyses are divided spatially into the "nuclear" and "northern arc" regions.

### Codes

- `pipelined_FOV_fitting_&_maps` |  The main code-- automated spectral modelling in the FOV of the 30x34 mini-cube centered on the core. The line flux, velocity and velocity dispersion information is also recorded for each spaxel for 2D map construction.

The .py file shows the code at a glance, while the .ipynb file with the fitting outputs for 1020 spaxels is too large for Github to display.   

- `3c297_broad&narrow_line_maps.ipynb` |  Creating line flux, velocity and velocity dispersion maps for all the detected lines in 3c297 spectrum (Halpha, [N II] and [S II]). 
- `spaxelwise_Nuclear7x7_NorthernArc.ipynb` |  This script does spaxel-wise line fitting on the nuclear 7x7 spaxel region and the northern arc. Each output shows the resulting fit values and plots the model+residuals on the observed spectrum. 
- `3c297_overlays.ipynb` |   Overlays for comparing .....
- `WHAN_diagnostic_spaxelwise.ipynb` |  Script to perform WHAN analysis ...
- `telluricOH_fitting.ipynb` |  Computing instrumental line broadening from telluric OH lines in the pre-background subtraction datacube.
- `generate_error_cube.ipynb` |  This script creates an "error cube" form the noise in the observed data, for fittting uncertainty computations.
- `3x3central_leastsq_coupledfit.ipynb` |  A simple script that illustrates the 6-component Gaussian fitting used for spectral modelling in our analysis. Integrating the spectra of the central 3x3 spaxel region, the modelling shows a clear broad component in Halpha.


----------------------------

![image](https://github.com/chetnaduggal/ifu-3D-datacube-analysis/blob/dbb2273ca01a6ae0a13d48cb7c630dbadd3ec4de/int-regions-spectra.png)






