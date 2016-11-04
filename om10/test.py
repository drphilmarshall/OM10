# ======================================================================

import os,unittest
import om10

# ======================================================================

class TestDB(unittest.TestCase):
    """
    Tests the OM10 db class.

    Notes
    -----
    Execute these tests with:
        nosetests
    from anywhere in the module, provided you have run
        pip install nose

    This file is part of the OM10 project, distributed under the
    MIT License, by Phil Marshall (KIPAC).
    Please cite: Oguri & Marshall (2010), MNRAS, 405, 2579.
    """
    # ------------------------------------------------------------------

    def setUp(self):
        self.db = om10.DB(catalog=os.path.expandvars("data/qso_mock.fits"))
        return

    def tearDown(self):
        return

    def test_download(self):
        original = self.catalog
        self.db.download()
        self.assertEqual(db.catalog, original)
        # Should really test that the files are the same length too.

    def test_get_lens(self):
        # Get one lens:
        id = 7176527
        lens = self.db.get_lens(id)
        self.assertIsNotNone(lens)
        self.assertAlmostEqual(lens['ZLENS'][0],0.556)
        self.assertAlmostEqual(lens['ZSRC'][0],1.88)
        self.assertAlmostEqual(lens['MAG'][0][0],2.9873)
        return

# ======================================================================

if __name__ == '__main__':

    unittest.main()

# ======================================================================

# Get one lens:

#     id = 7176527
#     lens = db.get_lens(id)

#     if lens is not None:
#         print "Lens ",id," has zd,zs = ",lens.ZLENS[0],lens.ZSRC[0]
#         print "and has images with magnifications: ",lens.MAG[0]

# # Look up one system:
#
#     id = 7176527
#     lens = db.get_lens(id)
#
#     if lens is not None:
#         print "Lens ",id," has zd,zs = ",lens.ZLENS[0],lens.ZSRC[0]
#         print "and has images with magnifications: ",lens.MAG[0]

# # To make a mock catalog of SDSS lenses:
#
#     db.select_random(maglim=19.1,area=8000.0,IQ=1.4)
#     db.write_table("OM10_SQLS_mock_lensed_quasars.fits")

# # To select a mock catalog of LSST lenses:
#
#     db.select_random(maglim=23.3,area=20000.0,IQ=0.75)
#     print db.Nlenses," LSST lenses, with zd = ",db.sample.ZLENS

# # To make a mock catalog of KIDS lenses:
#
#     # db.select_random(maglim=22.9,area=1500.0,IQ=0.7,Nlens=1e7)
#     db.select_random(maglim=22.9,area=1500.0,IQ=0.7)
#     db.write_table("OM10_KiDS_mock_lensed_quasars.fits")

# # To make a mock catalog of PS1 lenses:
#
#     db.select_random(maglim=21.4,area=30000.0,IQ=1.0)
#     db.write_table("OM10_PS1_mock_lensed_quasars.fits")
#
# # and export them for plotting:
#
#     pars = ['ZLENS','ZSRC','APMAG_I','MAGI','IMSEP']
#     db.export_to_cpt(pars,"OM10_PS1_mock_lensed_quasars.cpt")

# # To make a mock catalog of LSST lenses:
#
#     db.select_random(maglim=21.5,area=20000.0,IQ=0.75)
#     print db.Nlenses," LSST lenses"
#
#     good = db.sample[numpy.where(db.sample.IMSEP > 1.0)]
#     print "Number with imsep > 1.0 arcsec = ",len(good)
#
#     bright = good[numpy.where(good.APMAG_I < 22.0)]
#     print "Number of these with md < 22 = ",len(bright)
#
#     lagged = bright[numpy.where(numpy.max(bright.DELAY,axis=1) > 10.0)]
#     print "Number of these with time delay > 10 days = ",len(lagged)
#
#     nearby = lagged[numpy.where((lagged.ZLENS > 0.1) * (lagged.ZLENS < 0.6))]
#     print "Number of these with 0.1 < zd < 0.6 = ",len(nearby)

# To make a mock catalog of DES time delay lenses:
#
#     db.select_random(maglim=18.0,area=5000.0,IQ=0.9)
#     db.write_table("OM10_DES_mock_time-delay_lensed_quasars.fits")

# and export them for plotting:
#
#     pars = ['ZLENS','APMAG_I','IMSEP']
#     db.export_to_cpt(pars,"OM10_DES_mock_lensed_quasars_lenses.cpt")
#     pars = ['ZSRC','MAGI','IMSEP']
#     db.export_to_cpt(pars,"OM10_DES_mock_lensed_quasars_sources.cpt")
#     pars = ['ZLENS','ZSRC','APMAG_I','MAGI','IMSEP']
#     db.export_to_cpt(pars,"OM10_DES_mock_lensed_quasars.cpt")

# # These files are designed to be plotted with CornerPlotter.py:
#
# CornerPlotter.py \
#   -o OM10_DES_mock_lensed_quasars_both.png \
#   OM10_DES_mock_lensed_quasars_sources.cpt,blue,shaded \
#   OM10_DES_mock_lensed_quasars_lenses.cpt,orange,shaded
#
# CornerPlotter.py \
#   -o OM10_DES_mock_lensed_quasars.png \
#   OM10_DES_mock_lensed_quasars.cpt,black,shaded
#
# CornerPlotter.py \
#   -o OM10_PS1-vs-DES_mock_lensed_quasars.png \
#   OM10_DES_mock_lensed_quasars.cpt,black,shaded \
#   OM10_PS1_mock_lensed_quasars.cpt,blue,outlines

# This script is part of the pappy module, available from
#   http://github.com/drphilmarshall/pappy


# # Read in LRGs from CFHTLS:
#
#     db.get_LRGs()
#
# # Associate LRGs with sample - this appends the CFHTLS magnitudes in all filters to each lens,
# # based on the i magnitude and redshift:
#
#     db.match_LRGs()

# # To select 10 lenses detectable with PS1 at each epoch:
#
#     db.select_random(maglim=21.4,area=30000.0,IQ=1.0,Nlens=10)
#     print db.Nlenses," representative PS1 3pi lenses, with zd = ", \
#       db.sample.ZLENS
#     # print "ugriz = ", \
#     #   db.sample.uMAG_LRG,db.sample.gMAG_LRG,db.sample.rMAG_LRG, \
#     #   db.sample.iMAG_LRG,db.sample.zMAG_LRG


# 10-sigma detection in a single epoch?
# surveys = PS1-3PI PS1-MDS DES-WL KIDS  HSC-WIDE HSC-DEEP LSST  SDSS-S82x100
# maglims = 21.4    23.3    23.6   22.9  24.9     25.3     23.3  21.3
# areas   = 30000   84      5000   1500  1500     30       20000 30000        # survey area in sq deg
# psfs    = 1.0     1.0     0.9    0.7   0.75     0.75     0.75  1.4          # PSF FWHM in arcsec
# Note that these numbers give lower yields that OM10, by about a factor of 2:
# this is just due to the single epoch requirement, in the stacked images we
# should be able to go much deeper.

# ======================================================================
