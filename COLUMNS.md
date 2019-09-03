
| Parameter | Type   	    | Original Name  | Description  |
|---	    |---	        |---	         |---	        |
| LENSID    | long int      |  id            | Unique identifier for each row |
| FLAGTYPE  | int           |  type          | TODO |
| NIMG      | int           |  nim           | Number of lensed images |
| ZLENS     | double        |  zd            | Redshift of the lens |
| VELDISP   | double        |  sigma         | Velocity dispersion of the lens (km/s) |
| ELLIP     | double        |  epsd          | Ellipticity of the lens |
| PHIE      | double        |  phid          | Position angle of the lens (degrees East of North TODO is this right?) |
| GAMMA     | double        |  gammax        | Magnitude of external shear |
| PHIG      | double        |  phix          | Orientation of external shear (degrees East of North) |
| ZSRC      | double        |  zs            | Redshift of the source |
| XSRC      | double        |  xs            | Intrinsic X position of the source |
| YSRC      | double        |  ys            | Intrinsic Y position of the source |
| MAGI_IN   | double        |  ms            | i band apparent magnitude of source if there were no lensing (mags) |
| MAGI      | double        |  m3            | Apparent i magnitude for 3rd brightest image (mags) |
| IMSEP     | double        |  dtheta        | Maximum separation between any two images (arcsec) (TODO from footnote 1. Is this right?) |
| XIMG      | double arr(4) |  xi            | X position of each image |
| YIMG      | double arr(4) |  yi            | Y position of each image |
| MAG       | double arr(4) |  mui           | Apparent magnitude of each image (mags) |
| DELAY     | double arr(4) |  ti            | Time delay of each image (TODO units) |
| KAPPA     | double arr(4) |  kappa         | Convergence at each image location |
| FSTAR     | double arr(4) |  fstar         | Stellar mass fraction at each image location |
| Dd        | double        |  Dd            | Angular diameter distance from the observer to the lens (Mpc) |
| DDLUM     | double        |  DLd           | Luminosity distance from the observer to the lens (Mpc) |
| ABMAG_I   | double        |  absmd         | Absolute i band magnitude of the lens (mags) |
| APMAG_I   | double        |  md            | Apparent i band magnitude of the lens (mags) |
| KCORR     | double        |  Kcorrd        | K correction between observed i band and rest frame i band (APMAG_I - 5log10(DDLUM) - KCORR = ABMAG_I) |
| DS        | double        |  Ds            | Angular diameter distance to the source (Mpc) |
| DDS       | double        |  Dds           | Angular diameter distance from the lens to the source (Mpc) |
| SIGCRIT   | double        |  Sigcrit       | Lensing critical density: (c^2/4\pi G) * DS/(DDS * Dd) (TODO units) |
| DSLUM     | double        |  DLs           | Luminosity distance from the observer to the source (Mpc) |
| L_I       | double        |  -             | TODO |
| REFF      | double        |  -             | TODO |
| REFF_T    | double        |  -             | TODO |
