# data_viewer

An interactive Matplotlib tool for inspecting galaxies one at a time. In a single window it shows a galaxy's spectra, its photometry, its best-fit model and where it sits among the full sample.

I built it in 2017, during my PhD at the Laboratoire d'Astrophysique de Marseille, for my own analysis of galaxies from the VIMOS Ultra-Deep Survey (VUDS).

## Used for

[Khusanova et al. 2020, *UV and Lyα luminosity functions of galaxies and star formation rate density at the end of HI reionization from the VIMOS UltraDeep Survey (VUDS)*](https://doi.org/10.1051/0004-6361/201935400), A&A 634, A97 ([arXiv:1903.01884](https://arxiv.org/abs/1903.01884)).

The paper started from 111 candidate high-redshift (z > 5) galaxies. A single emission line or continuum break can be mistaken for a feature of a much closer galaxy, so each candidate had to be checked individually before it could be counted at z > 5. I used data_viewer to look at the selected candidates one by one and to check whether the 1D and 2D spectra, the photometry and the best-fit SED at the spectroscopic redshift agree. The final sample has 52 galaxies.

## What it does

The window brings together everything known about one galaxy, and every panel updates when you select a different galaxy.

```
┌──────────────────┬─────────────────────────────────────────────┐
│ SFR vs. M*       │ 2D spectrum, full range                     │
│ click a galaxy   │ 1D spectrum + photometry, zoom window shaded│
├──────────────────┼──────────────────────────┬──────────────────┤
│ sSFR vs. M*      │ 2D spectrum, zoomed      │ sliders:         │
│                  │ 1D spectrum, zoomed      │ Vmin, Vmax, z cut│
├──────────────────┴──────────────────────────┼──────────────────┤
│ SED: spectrum, photometry, best-fit template│ selected line    │
│ [Previous] [Next]       [Ext_law] [Flux/Mag]│ region           │
└─────────────────────────────────────────────┴──────────────────┘
```

- **Sample overview.** Star-formation rate and specific SFR are plotted against stellar mass for the whole sample, with 68% error bars. Clicking a point loads that galaxy, and the current galaxy is shown in red. Galaxies without a mass estimate are placed along the edge of the plot, so they can still be clicked.
- **Spectra.** The 2D and 1D spectra are shown over the full wavelength range and in a zoom window, with the broadband photometry converted to flux and overplotted. Clicking the full spectrum moves the zoom window. Two sliders set the contrast of the 2D images.
- **Line measurement.** Two clicks on the continuum of the zoomed spectrum measure the equivalent width of a line against a linear continuum between those points.
- **SED.** The observed spectrum is shown with the photometry and the best-fit Bruzual & Charlot (BC03) template. The template is redshifted, reddened with the galaxy's E(B−V) and normalized to its i-band magnitude through the filter curve. Buttons switch the extinction law (Calzetti, SMC) and the units (flux or AB magnitude).
- **Sample selection.** A slider sets the redshift cut and reloads the sample, which covers three VUDS fields: COSMOS, ECDFS and VVDS-02h.
- **Template generation.** For galaxies whose template does not exist yet, the tool writes a shell script that passes the fit parameters to BC03's `galaxevpl` and runs it.

`auxiliaryfunctions.py` holds the supporting code:
- catalog parsing;
- magnitude/flux conversions and synthetic photometry through filter curves;
- cosmology helpers (luminosity distance, age of the universe at a given redshift);
- a rest-frame line list;
- an implementation of IRAF's zscale display algorithm, adapted from STScI's numdisplay.

## Stack

Python 2.7 · Matplotlib (widgets, pick and click events) · NumPy · SciPy · Astropy (FITS) · BC03 stellar population models

## Status

This is an archived snapshot from 2017, kept as it was written. It targets Python 2.7 and the Matplotlib of that time. It reads the survey spectra, the SED-fitting catalogs and the BC03 models from local paths, and none of these are part of this repository.
