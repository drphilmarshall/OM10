# OM10

<a href='https://travis-ci.org/drphilmarshall/OM10'>
<img src='https://secure.travis-ci.org/drphilmarshall/OM10.png?branch=master'></a>

Tools for working with the Oguri &amp; Marshall (2010) mock catalog of strong gravitational lenses

### Getting started

Edit the following lines and add them to your `.login` file:

    setenv OM10_DIR ${WORK_DIR}/OM10
    setenv PYTHONPATH ${OM10_DIR}:${PYTHONPATH}

Now you can just:

    import om10

You'll need `pyfits` and `matplotlib` too.

### Example use

Read in the master FITS catalog and look up one system:    

    db = om10.DB(catalog="data/qso_mock.fits")

    id = 7176527
    lens = db.get_lens(id)

Make a plot of it:

    om10.plot_lens(lens)

Select a mock LSST sample:

    lenses = db.select_random(maglim=23.3,area=20000.0,IQ=0.7)

### License, Credits

This code is distributed under the MIT license, and is being developed by Marshall and Baumer (KIPAC), Liao (UCLA) and Agnello (UCLA).

If you use the OM10 mock lens catalog in your research, please cite [Oguri &amp; Marshall (2010)](http://adsabs.harvard.edu/abs/2010MNRAS.405.2579O). Here's the bibtex for you!

    @OM10,
       author = {{Oguri}, M. and {Marshall}, P.~J.},
        title = "{Gravitationally lensed quasars and supernovae in future wide-field optical imaging surveys}",
      journal = {\mnras},
    archivePrefix = "arXiv",
       eprint = {1001.2037},
     primaryClass = "astro-ph.CO",
     keywords = {gravitational lensing: strong, cosmological parameters, cosmology: theory},
         year = 2010,
        month = jul,
       volume = 405,
        pages = {2579-2593},
          doi = {10.1111/j.1365-2966.2010.16639.x},
       adsurl = {http://adsabs.harvard.edu/abs/2010MNRAS.405.2579O},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    }
