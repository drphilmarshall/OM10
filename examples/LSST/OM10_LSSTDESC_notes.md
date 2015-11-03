# Cosmology Forecasts from OM10 Time Delays with CosmoSIS

* LSST DESC Hack Day, October 30 2015, Argonne National Lab *

David Finley, Elise Jennings, Phil Marshall

First:
- PM and DF created [OM10 fork for DF](https://github.com/davidfinley/OM10)
- PM modified OM10 to allow use of astropy

- PM set up the `.ini` files so that `cosmosis` can be run from the OM10 examples folder:
```
bash
export COSMOSIS_SRC_DIR = /Users/pjm/work/stronglensing/LSST/DESC/Cosmology/cosmosis
source ${COSMOSIS_SRC_DIR}/config/setup-cosmosis
```

MOCK DATA:
- Used OM10 catalog to make 466 $(z_d, z_s)$ pairs at `OM10_LSSTDESC_zd_zs.txt`
- Lens system selection was as follows:
```python
maglim = 22.0   # Faintest image must be well detected
area = 18000.0  # LSST WFD survey area (sq deg)
IMSEP > 1.0     # Well-resolved (!) quasar images (arcesec)
APMAG_I < 21.0  # Bright lens galaxies for follow-up velocity dispersions
DELAY > 10.0    # Longish time delays (days)
```

- EJ used [Suyu et al (2010)](http://arxiv.org/pdf/0910.2773.pdf) to create CosmoSIS modules to take $(z_d, z_s)$ from above and append mock observations $(D_{\Delta t}, \sigma_{\Delta t})$. We assumed $\sigma_{\Delta t} = 0.05$ (5% precision per lens).
- DF and PM ran the modules on 466 lenses with:
```
    cosmosis OM10_LSSTDESC_generate.ini
```
This results in a new plain text data file, with 4 columns: `OM10_LSSTDESC_mock_data.txt`

- EJ created CosmoSIS modules for MOCK COSMOLOGICAL PARAMETER INFERENCE
- DF and PM used the second CosmoSIS module created by EJ in the `develop` branch of the CosmoSIS repository to take the mock data for 466 lenses and infer $H_0$, etc,
```
	cosmosis OM10_LSSTDESC_inference.ini
```
saving the output as `OM10_LSSTDESC_inferred_parameters.txt`. This takes a few minutes.

- The `postprocess` script can be used to make some standard-format plots of the MCMC-sampled posterior PDF, marginalized to 1 and 2D.
```
postprocess  --burn 5000 -o OM10_LSSTDESC_plots -p OM10_LSSTDESC_inferred OM10_LSSTDESC_inference.ini
```
This produces a lot of plots in a new folder, `OM10_LSSTDESC_plots`. To view them all (on a Mac), do
```
open OM10_LSSTDESC_plots/OM10_LSSTDESC_inferred_*.png
```
- Compressed 1D marginalized arameter inferences can be read from one of the files in the plots folder, eg,
```
more OM10_LSSTDESC_inferred_means.txt

# parameter mean std_dev
# OM10_LSSTDESC_inferred_parameters.txt
cosmological_parameters--omega_m   5.673062e-01   1.784710e-01
cosmological_parameters--h0   7.130534e-01   6.988178e-03
cosmological_parameters--omega_b   4.199131e-02   1.098027e-02
cosmological_parameters--omega_k   -2.133684e-05   1.641963e-04
cosmological_parameters--w   -9.135106e-01   2.050209e-01
cosmological_parameters--wa   -6.870567e-03   2.631039e-01
like   -inf   nan
```
The `omega_m` PDF is broad, but does seem to be shifted somewhat to higher values (0.57 +/- 0.18) - we should re-run with a different sample and check this. The Hubble constant `h` (0.713 +/- 0.007) is possibly also shifted low (but less significantly). Nice to see 1% precision even with freely varying curvature and Dark Energy parameters. Likewise, nice to see `omega_k` come out to be (-0.00002 +/- 0.00016). Its possible that this is a result of mismatched sampling distribution (when generating data) and likelihood; it could also just be a phase space effect, from the large uniform prior volume in the other parameters. First step should be to make a "corner plot" showing all 1 and 2-D marginalized PDFs.

- PM wrote two `cosmosis` issues about 1) [the samples with `-inf` log likelihood being accepted by `emcee` and then used by `postprocess`](https://bitbucket.org/joezuntz/cosmosis/issues/119/do-emcee-and-postprocess-deal-correctly), and 2) [the contents (header line? credible intervals?) of the `medians.txt` file.](https://bitbucket.org/joezuntz/cosmosis/issues/118/medians-file-contains-incorrect-header)
