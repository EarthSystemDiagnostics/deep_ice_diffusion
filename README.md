# deep_ice_diffusion
Code for paper "Novel approach to estimate the water isotope diffusion length in deep ice cores and application to MIS 19 in the EPICA Dome C ice core"

Must run scripts in a specific order:
1. Run Bayesian_conventional.R and conventional_bias_simulated.R (order here does not matter)
2. Run MIS19_diffusion_length.R (this is the main script)
3. Run interpolation_error.R (this is for an appendix figure)

The first 2 scripts produce data which is called by the main script. All figures are plotted in the main script except one appendix figure, Fig. A1, which is plotted in the final script.
