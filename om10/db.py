# ======================================================================

import sys,os,subprocess
import numpy as np
import os
from numpy import *
import math
from astropy.table import Table, hstack
from sklearn.neighbors import NearestNeighbors
from sklearn import preprocessing

# from astropy.table import Table
import astropy.io.fits as pyfits

import om10

vb = True

# ======================================================================

class DB(object):

    """
    NAME
        DB

    PURPOSE
        Read in an OM10 catalog and store it as a 'database', in the
        loosest sense of the word.

    COMMENTS

    INITIALISATION
        catalog            OM10 FITS lens catalog.

    METHODS
        get_lens(ID)       Extract one lens object
        select_random(Nlens=None,maglim=99.0,area=100000.0,IQ=0.0)
                           Extract a random sample of lenses that meet
                             some simple observational criteria.
        make_table()       Make a FITS table given raw text files.
        write_table(file)  Write out a new FITS table to file
        paint(lens/qso,catalog)
                           Find a matching real object from the supplied
                             catalog, and transfer its colors
        place(catalog)     Find a matching LRG and transfer its sky position.

    BUGS
      - LRG properties have not been added to table!

    AUTHORS
      This file is part of the OM10 project, distributed under the
      GPL v2, by Phil Marshall (KIPAC).
      Please cite: Oguri & Marshall (2010), MNRAS, 405, 2579.

    HISTORY
      2013-04-07  Started Marshall (Oxford)
    """
    # ------------------------------------------------------------------

    def __init__(self,catalog=None,generate=False):

        self.name = 'OM10 database'
        self.catalog = os.path.expandvars(catalog)

        # Make a FITS table from the supplied text catalogs, if required:
        if generate:
            self.make_table()
            self.catalog = os.path.expandvars("$OM10_DIR/data/qso_mock.fits")
            self.write_table(self.catalog)

        # If given a catalog, read it in:
        if self.catalog is not None:
            #self.lenses = pyfits.getdat(self.catalog)
            self.lenses = Table.read(self.catalog,format='fits')
            #print self.lenses

        # No down-sampling has been done yet, but methods operate on sample:
        self.sample = self.lenses.copy()

        # Count lenses:
        try: self.Nlenses = len(self.lenses['LENSID'])
        except: self.Nlenses = 0

        return

    # ------------------------------------------------------------------

    def write_table(self,catalog):
        try: os.remove(catalog)
        except OSError: pass

        if len(self.sample) == len(self.lenses):
            pyfits.writeto(catalog,self.lenses)
        else:
            pyfits.writeto(catalog,self.sample)

        if vb: print "om10.DB: wrote catalog of ",self.Nlenses," OM10 lenses to file at "+catalog
        return

    # ------------------------------------------------------------------

    def get_lens(self,ID):

        try: rec = self.lenses[self.lenses['LENSID'] == ID]
        except: rec = None

        return rec

    # ------------------------------------------------------------------

    def select_random(self,Nlens=None,maglim=99.0,area=100000.0,IQ=0.0):

        # Select all lenses that meet the rough observing criteria:

        try:
            sample = self.sample.copy()
            sample = sample[sample['MAGI'] < maglim]
            sample = sample[sample['IMSEP'] > 0.67*IQ]
        except:
            if vb: print "om10.DB: selection yields no lenses"
            return None

        # Compute expected number of lenses in survey:

        if Nlens is None:
            N = int(len(sample) * (area / 20000.0) * 0.2)
        else:
            N = Nlens
        if vb: print "om10.DB: selection yields ",N," lenses"
        if N > len(sample):
            print "om10.db: Warning: too few lenses in catalog, returning ",len(sample)," instead"
            N = len(sample)

        # Shuffle sample and return only this, or the required, number of systems:

        index = range(len(sample))
        np.random.shuffle(index)
        index = index[0:N]

        self.sample = sample[index]
        self.Nlenses = len(self.sample)

        return

    # ------------------------------------------------------------------

    def get_sky_positions(self,dmag=0.2,dz=0.2,input_cat='$OM10_DIR/data/CFHTLS_LRGs.txt'):

        LRGfile = os.path.expandvars(input_cat)
        try:
            d = np.loadtxt(LRGfile)
        except IOError:
            print "Cannot find LRG catalog!"
        if vb: print "om10.DB: read in LRG sky position data from ",LRGfile

        # Put LRG parameters in LRG structure:
        # RA DEC z mag_u mag_g mag_r mag_i mag_z

        self.LRGs = {}
        self.LRGs['RA']       = np.array(d[:, 0])
        self.LRGs['DEC']      = np.array(d[:, 1])
        self.LRGs['redshift'] = np.array(d[:, 2])
        self.LRGs['g-r']      = np.array(d[:, 4]) - np.array(d[:, 5])
        self.LRGs['r-i']      = np.array(d[:, 5]) - np.array(d[:, 6])
        self.LRGs['i-z']      = np.array(d[:, 6]) - np.array(d[:, 7])
        self.LRGs['mag_i']    = np.array(d[:, 6])
        features = np.array([self.LRGs['redshift'], self.LRGs['g-r'], self.LRGs['r-i'], self.LRGs['i-z'], self.LRGs['mag_i']]).transpose()
        self.LRGs['feature_scaler'] = preprocessing.StandardScaler().fit(features)
        scaled_features = self.LRGs['feature_scaler'].transform(features)
        self.LRGs['nbrFinder'] = NearestNeighbors(n_neighbors=1,algorithm='auto',metric='euclidean').fit(scaled_features)

        print "Mean LRG RA,DEC,z = ",np.average(self.LRGs['RA']),np.average(self.LRGs['DEC']),np.average(self.LRGs['redshift']),np.average(self.LRGs['mag_i']);
        print "Mean LRG i,(g-r) = ",np.average(self.LRGs['RA']),np.average(self.LRGs['DEC']),np.average(self.LRGs['redshift']),np.average(self.LRGs['mag_i']);

        if vb: print "om10.DB: number of LRGs stored = ",len(self.LRGs['redshift'])

        return

    # ------------------------------------------------------------------

    def assign_sky_positions(self,verbose=False):

        #try:
        #    tmp = self.sample.['MAGG_LENS'][0]
        #except :

        reallyverbose = verbose

        # Prepare new columns for LRG properties:
        self.sample['RA'] = 0.0
        self.sample['DEC'] = 0.0

        scaler = self.LRGs['feature_scaler']
        index_list = []

        for lens in self.sample:
            lens_features = np.array([lens['ZLENS'], lens['MAGG_LENS']-lens['MAGR_LENS'], \
            lens['MAGR_LENS']-lens['MAGI_LENS'], lens['MAGI_LENS']-lens['MAGZ_LENS'], lens['APMAG_I']])

            scaled_lens_features = scaler.transform(lens_features)
            distance, index = self.LRGs['nbrFinder'].kneighbors(scaled_lens_features)
            index_list.append(index)
            lens['RA'] = self.LRGs['RA'][index]
            lens['DEC'] = self.LRGs['DEC'][index]

            if reallyverbose:
                print "  Lens i,z: ",self.sample['APMAG_I'][k],self.sample['ZLENS'][k]
                print "  Lens RA,DEC: ",self.sample['RA'][k],self.sample['DEC'][k]

        return index_list

# ----------------------------------------------------------------------------

    def make_sim_input_catalog(self):
        n_obj = len(self.sample)
        n_tot_img = np.sum(self.sample['NIMG'])
        output_cols=['LENSID','RA','DEC','XIMG','YIMG','G','R','I','Z']

        sim_cat = Table(np.zeros((n_tot_img+n_obj,len(output_cols)), \
        dtype=('>i4', '>f8', '>f8', '>f8', '>f8', '>f8', '>f8', '>f8', '>f8')),names=output_cols)

        out_idx = 0
        for lens in self.sample:
            sim_cat[out_idx] = (lens['LENSID'],lens['RA'],lens['DEC'],0,0,lens['MAGG_LENS'], \
                                lens['MAGR_LENS'], lens['MAGI_LENS'], lens['MAGZ_LENS'])
            out_idx += 1
            mag_adjust = 2.5*np.log10(abs(lens['MAG'][lens['MAG'] != 0]))
            for img in np.arange(lens['NIMG']):
                sim_cat[out_idx] = (lens['LENSID'],lens['RA']+lens['XIMG'][img]/(np.cos(np.deg2rad(lens['DEC']))*3600.0),lens['DEC']+lens['YIMG'][img]/3600.0,\
                                    lens['XIMG'][img],lens['YIMG'][img],lens['MAGG_SRC']-mag_adjust[img], \
                                    lens['MAGR_SRC']-mag_adjust[img], lens['MAGI_SRC']-mag_adjust[img],\
                                    lens['MAGZ_SRC']-mag_adjust[img])
                out_idx += 1
        return sim_cat

# ----------------------------------------------------------------------------

    def paint(self,Nmax=None,verbose=False,lrg_input_cat='$OM10_DIR/data/LRGo.txt',qso_input_cat='$OM10_DIR/data/QSOo.txt'):
        ## read data from SDSS
        f=open(os.path.expandvars(lrg_input_cat),'r')
        lrg=loadtxt(f)
        f.close()
        #print lrg[0,0],lrg.shape
        g=open(os.path.expandvars(qso_input_cat),'r')
        qso=loadtxt(g)
        g.close()
        #print qso[0,0],qso.shape

        ###MY OWN REDSHIFT ONLY MATCHING HERE:

        lens_props = ['MAGG_LENS','MAGR_LENS','MAGI_LENS','MAGZ_LENS', \
        'MAGW1_LENS','MAGW2_LENS','MAGW3_LENS','MAGW4_LENS', 'SDSS_FLAG_LENS']

        src_props = ['MAGG_SRC','MAGR_SRC','MAGI_SRC','MAGZ_SRC', \
        'MAGW1_SRC','MAGW2_SRC','MAGW3_SRC','MAGW4_SRC', 'SDSS_FLAG_SRC']

        tmp_lens = Table(np.zeros((len(self.sample),len(lens_props)),dtype='f8'),names=lens_props)
        tmp_src = Table(np.zeros((len(self.sample),len(src_props)),dtype='f8'),names=src_props)

        if verbose: print 'setup done'

        lrg_sort = lrg[np.argsort(lrg[:,0]),:]
        qso_sort = qso[np.argsort(qso[:,0]),:]
        lens_count = 0

        for lens in self.sample:

            #paint lens
            ind = np.searchsorted(lrg_sort[:,0],lens['ZLENS'])
            if ind >= len(lrg_sort): ind = len(lrg_sort) - 1
            tmp_lens[lens_count] = lrg_sort[ind,6:] - lrg_sort[ind,8] + lens['APMAG_I'] #assign colors, not mags
            #paint source
            qso_ind = np.searchsorted(qso_sort[:,0],lens['ZSRC'])
            if qso_ind >= len(qso_sort): qso_ind = len(qso_sort) - 1
            tmp_src[lens_count] = qso_sort[qso_ind,1:] - qso_sort[qso_ind,3] + lens['MAGI']

            lens_count += 1

        self.sample = hstack([self.sample,tmp_lens,tmp_src])

        return

# ======================================================================

if __name__ == '__main__':

# Some examples!

# To make the FITS catalog from the master text catalogs:
#
#     db = om10.DB(generate=True)

# To read in an old FITS catalog:

    db = om10.DB(catalog=os.path.expandvars("$OM10_DIR/data/qso_mock.fits"))

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

# To make a mock catalog of LSST lenses:

#     db.select_random(maglim=21.5,area=20000.0,IQ=0.75)
#     print db.Nlenses," LSST lenses"
    db.select_random(maglim=23.3,area=18000.0,IQ=0.75)
    print db.Nlenses," LSST lenses"

    good = db.sample[np.where(db.sample.IMSEP > 1.0)]
    print "Number with imsep > 1.0 arcsec = ",len(good)

    bright = good[np.where(good.APMAG_I < 22.0)]
    print "Number of these with md < 22 = ",len(bright)

    lagged = bright[np.where(np.max(bright.DELAY,axis=1) > 10.0)]
    print "Number of these with time delay > 10 days = ",len(lagged)

    nearby = lagged[np.where((lagged.ZLENS > 0.1) * (lagged.ZLENS < 0.6))]
    print "Number of these with 0.1 < zd < 0.6 = ",len(nearby)

# Example outputs:

# Mag limit 21.5:
# 813  LSST lenses
# Number with imsep > 1.0 arcsec =  581
# Number of these with md < 22 =  523
# Number of these with time delay > 10 days =  505
# Number of these with 0.1 < zd < 0.6 =  254

# Mag limit 23.3:
# 2813  LSST lenses
# Number with imsep > 1.0 arcsec =  1911
# Number of these with md < 22 =  1614
# Number of these with time delay > 10 days =  1559
# Number of these with 0.1 < zd < 0.6 =  795

# To make a mock catalog of DES time delay lenses:
#
#     db.select_random(maglim=20.0,area=5000.0,IQ=0.9)
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


# Read in LRGs from CFHTLS:
    db.get_sky_positions()

# Associate LRGs with sample - this appends the CFHTLS magnitudes in all filters to each lens,
# based on the i magnitude and redshift:

    db.assign_sky_positions()

# How many got placed properly?

    good = db.lenses[np.where(db.lenses.RA > 0.0)]
    bad = db.lenses[np.where(db.lenses.RA < 0.0)]
    print "No. of OM10 lenses with matching LRG sky positions = ",len(good)
    print "  mean,min,max redshift = ",np.average(good.ZLENS),np.min(good.ZLENS),np.max(good.ZLENS)
    print "No. of OM10 lenses with no matching sky position = ",len(bad),np.min(bad.ZLENS),np.max(bad.ZLENS)
    print "  mean,min,max redshift = ",np.average(bad.ZLENS)

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
