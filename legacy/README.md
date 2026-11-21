# Original Fortran and MATLAB (2000–2005)

`fortran/` is PLANT1D as committed to github.com/mjpuma/SoilWaterFlow in 2018: block-centered finite differences for the 1-D unsaturated flow equation with transpiration and evaporation, by M. A. Celia, A. J. Guswa and M. J. Puma. The `.IN` files are its example (Burkea sandy clay, 400 days). The compiled Mac binary and output files from that repository are not kept here; the outputs are in `tests/data/fortran_reference` (gzipped).

`matlab/` holds the original pre- and post-processing scripts (`presoil*.m` write the input files and generate storms; `rootprofile.m` and `betaweight.m` compute beta-distributed root weights).

The Python port in `src/soilwaterflow` reproduces this code bit for bit; see the main README for the build command and the comparison.
