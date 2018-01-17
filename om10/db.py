#  ======================================================================
from __future__ import print_function

import sys,os,subprocess
import numpy as np
import os, urllib
from numpy import *
import math
from astropy.table import Table, hstack
import astropy.io.fits as pyfits
from sklearn.neighbors import NearestNeighbors
from sklearn import preprocessing

from lenspop import population_functions, distances
from stellarpop import tools
import matplotlib.pyplot as plt
import scipy
import astropy.cosmology as cosmo
import time

import om10

# ======================================================================


class DB(object):
    """
    Reads in an OM10 catalog and stores it as a 'database', in the
    loosest sense of the word.

    Parameters
    ----------
    catalog : string
        OM10 FITS lens catalog name.

    Notes
    -----
      This file is part of the OM10 project, distributed under the
      MIT License by Phil Marshall (KIPAC).
      Please cite: Oguri & Marshall (2010), MNRAS, 405, 2579.
    """
    # ------------------------------------------------------------------

    def __init__(self, catalog=None, generate=False, vb=True):

        self.name = 'OM10 database'
        self.vb = vb

        if catalog is None:
            # Use the one that comes with the package:
            self.catalog = \
                os.path.join(os.path.dirname(os.path.realpath(__file__)),'../data/qso_mock.fits')
        else:
            self.catalog = os.path.expandvars(catalog)

        self.setup()

        return

    # ------------------------------------------------------------------

    def setup(self):
        """
        Read in the catalog and set up an initial (super) sample.
        """
        # Read in the catalog:
        self.lenses = Table.read(self.catalog, format='fits')
        t = Table(np.arange(len(self.lenses)).reshape(len(self.lenses), 1), names=('weight',), dtype=('f4',))
        for i in range(len(self.lenses)):
            t['weight']=1;
        self.lenses.add_columns(t.columns.values())
        if self.vb:
            print('OM10: Full db.lenses table contains {:d} systems'.format(len(self.lenses)))

        # No down-sampling has been done yet, but all methods operate
        # on a "sample" - so make a copy:
        self.sample = self.lenses.copy()
        self.Nlenses = len(self.sample)
        if self.vb:
            print('OM10: Initial db.sample contains {:d} systems'.format(self.Nlenses))

        return

    def reset(self):
        self.setup()
        return

    # ------------------------------------------------------------------

    def download(self):
        """
        Downloads a copy of the primary OM10 FITS table.

        Notes
        -----
        This could be useful, in case the one that came with the
        package gets deleted, or you just want a local copy. The new
        catalog will be placed in the current working directory, and the filename (stored in `db.catalog`) updated.
        """
        url = 'https://github.com/drphilmarshall/OM10/raw/master/data/qso_mock.fits'
        self.catalog = url.split('/')[-1]
        if self.vb: print('OM10: Looking for local catalog {:s}'.format(self.catalog))
        if not os.path.isfile(self.catalog):
            urllib.urlretrieve(url, self.catalog)
            if self.vb: print('OM10: Downloaded catalog: {:s}'.format(self.catalog))
        else:
            if self.vb: print('OM10: File already exists, no need to download.')
        return

    # ------------------------------------------------------------------

    def write_table(self,catalog):
        try: os.remove(catalog)
        except OSError: pass

        self.sample = np.array(self.sample)
        if len(self.sample) == len(self.lenses):
            pyfits.writeto(catalog,self.lenses)
        else:
            pyfits.writeto(catalog,self.sample)
        if self.vb: print('OM10: Wrote catalog of {:d} OM10 lenses to file at {:s}'.format(self.Nlenses, catalog))
        return

    # ------------------------------------------------------------------

    def get_lens(self,ID):

        try: rec = self.lenses[self.lenses['LENSID'] == ID]
        except: rec = None

        if self.vb:
            print('OM10: Extracted OM10 lens number {:d}:'.format(ID))
            print(rec)

        return rec

    # ------------------------------------------------------------------

    def select_random(self, Nlens=None, maglim=99.0, area=100000.0,
                      IQ=0.0):
        """
        Selects an appropriately-sized random sample of lenses that
        meet the rough observing criteria given.

        Parameters
        ----------
        Nlens : int, optional
            Specific desired number of lenses
        maglim : float
            10-sigma point source detection limit
        area : float
            Total survey area, in square degrees
        IQ : float
            Median survey image quality, in arcsec

        Notes
        -----
        If `Nlens` is not specified, it is calculated based on the OM10
        model. The full OM10 catalog contains 100,000 sq degrees worth
        of lenses. The detection criteria assumed are given in the OM10
        paper: we assume that the 3rd brightest quasar image must be
        brighter than the given `maglim`, and the image separation must
        be greater than 0.67 times the given `IQ`.
        """

        try:
            sample = self.sample.copy()
            sample = sample[sample['MAGI'] < maglim]
            sample = sample[sample['IMSEP'] > 0.67*IQ]
        except:
            if self.vb: print('OM10: Selection yields no lenses')
            return None

        # Compute expected number of lenses in survey:
        if Nlens is None:
            N = int(len(sample) * (area / 20000.0) * 0.2)
        else:
            N = Nlens
        if self.vb: print('OM10: selection yields {:d} lenses'.format(N))
        if N > len(sample):
            print('OM10: Warning: too few lenses in catalog, returning {:d} instead'.format(len(sample)))
            N = len(sample)

        # Shuffle sample and return only this, or the required, number of systems:
        index = range(len(sample))
        np.random.shuffle(index)
        index = index[0:N]

        self.sample = sample[index]

        return

    # ------------------------------------------------------------------

    def get_sky_positions(self, dmag=0.2, dz=0.2,
                          input_cat='$OM10_DIR/data/CFHTLS_LRGs.txt'):

        LRGfile = os.path.expandvars(input_cat)
        try:
            d = np.loadtxt(LRGfile)
        except IOError:
            print('Cannot find LRG catalog!')
        if self.vb: print('OM10: read in LRG sky position data from {:s}'.format(LRGfile))

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

        print('Mean LRG RA,DEC,z = ', np.average(self.LRGs['RA']), np.average(self.LRGs['DEC']), np.average(self.LRGs['redshift']), np.average(self.LRGs['mag_i']))
        print('Mean LRG i,(g-r) = ', np.average(self.LRGs['RA']), np.average(self.LRGs['DEC']), np.average(self.LRGs['redshift']), np.average(self.LRGs['mag_i']))

        if self.vb: print('om10.DB: number of LRGs stored = ', len(self.LRGs['redshift']))

        return

    # ------------------------------------------------------------------

    def assign_sky_positions(self, verbose=False):

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
                print('  Lens i,z: ', self.sample['APMAG_I'][k], self.sample['ZLENS'][k])
                print('  Lens RA,DEC: ', self.sample['RA'][k], self.sample['DEC'][k])

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

# The paint method became really long, so needed to decompose this part out
    def calculate_rest_frame_r_magnitude(self, sed, veldisp, redshift, cosmo):
        """
        Computes rest-frame r-band magnitude of a lens galaxy

        Parameters
        ----------
        sed : string
            Name of SED to use
        redshift : float
            Redshift of object
        d : float
            Distance modulus to object
        veldisp : float
            For lens galaxies, the velocity dispersion can be passed in to
            provide the absolute magnitude via the Fundamental Plane

        Returns
        -------
        RF_Rmag_app : float
            Reference r-band apparent magnitude
        offset_RF_abs : float
            Magnitude offset for converting absolute to apparent magnitude

        Notes
        -----
        We don't need this function when painting quasars, because the OM10
        catalog contains a reliable i-band apparent magnitude for each source.
        """
        lenspop_constructor = population_functions.LensPopulation_()
        # Reference Frame Absolute R magnitude
        RF_RMag_abs, _ = lenspop_constructor.EarlyTypeRelations(veldisp)
        RMag_abs = tools.ABFilterMagnitude(Rfilter, sed, redshift)
        distModulus = cosmo.distmod(redshift).value
        Rmag_app = RMag_abs + distModulus
        offset_abs_app = RMag_abs - Rmag_app
        offset_RF_abs = RF_RMag_abs - RMag_abs
        RF_Rmag_app = RF_RMag_abs - offset_abs_app
        return RF_Rmag_app, offset_RF_abs, distModulus

# ----------------------------------------------------------------------------

    def paint(self, Nmax=None, verbose=False,
              lrg_input_cat='$OM10_DIR/data/LRGo.txt',
              qso_input_cat='$OM10_DIR/data/QSOo.txt',
              synthetic=False):
        """
        Add new columns to the table, for the magnitudes in various filters.

        Parameters
        ----------
        synthetic : boolean
            Use `lenspop` to make synthetic magnitudes in various filters
        target : string
            Paint lenses ('lens') or sources ('source')
        lrg_input_cat : string
            Name of LRG catalog, if not using synthetic paint
        qso_input_cat : string
            Name of QSO catalog, if not using synthetic paint
        verbose : boolean
           print progress to stdout

        Notes
        -----
        * Synthetic painting is very slow, as we loop over each object.
        * The treatment of QSOs may be flawed: the offset calculation has not
          been tested.

        """

        if synthetic==False:
        # read data from SDSS
            f=open(os.path.expandvars(lrg_input_cat),'r')
            lrg=loadtxt(f)
            f.close()
            g=open(os.path.expandvars(qso_input_cat),'r')
            qso=loadtxt(g)
            g.close()

        ###MY OWN REDSHIFT ONLY MATCHING HERE:

            lens_props = ['MAGG_LENS','MAGR_LENS','MAGI_LENS','MAGZ_LENS', \
            'MAGW1_LENS','MAGW2_LENS','MAGW3_LENS','MAGW4_LENS', 'SDSS_FLAG_LENS']

            src_props = ['MAGG_SRC','MAGR_SRC','MAGI_SRC','MAGZ_SRC', \
            'MAGW1_SRC','MAGW2_SRC','MAGW3_SRC','MAGW4_SRC', 'SDSS_FLAG_SRC']

            tmp_lens = Table(np.zeros((len(self.sample),len(lens_props)),dtype='f8'),names=lens_props)
            tmp_src = Table(np.zeros((len(self.sample),len(src_props)),dtype='f8'),names=src_props)

            if verbose: print('setup done')

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


        if synthetic==True:
            lens_count = 0
            total = len(self.sample)
            Rfilter = tools.filterfromfile('r_SDSS')
            Ufilter = tools.filterfromfile('u_SDSS')
            # sort the Ufilter array
            Ufilterarg=np.sort(Ufilter[1])
            Ufilter = (Ufilter[0], Ufilterarg, 1)
            Gfilter = tools.filterfromfile('g_SDSS')
            Ifilter = tools.filterfromfile('i_SDSS')
            Zfilter = tools.filterfromfile('z_SDSS')
            self.Nlenses=len(self.sample)
            bands = ('r_SDSS_lens', 'g_SDSS_lens', 'i_SDSS_lens', 'z_SDSS_lens', 'u_SDSS_lens','r_SDSS_quasar', 'g_SDSS_quasar', 'i_SDSS_quasar', 'z_SDSS_quasar', 'u_SDSS_quasar')
            if verbose: print('OM10: computing synthetic magnitudes in the following bands: ', bands)
            # call a distance class constructor
            d = distances.Distance()
            # number of data in the table of calculated magnitude
            totalEntrees = self.Nlenses*10.0
            t = Table(np.arange(totalEntrees).reshape(self.Nlenses, 10),
                      names=bands)
            Lsed = tools.getSED('BC_Z=1.0_age=9.000gyr')
            Qsed = tools.getSED('agn')
            from astropy.cosmology import FlatLambdaCDM
            cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
            lenspop_constructor = population_functions.LensPopulation_()
            for lens in self.sample:
                # calculate the quasar magnitude
                redshift = lens['ZSRC']
                RF_Imag_app_q = lens['MAGI_IN']
                Qoffset = RF_Imag_app_q - tools.ABFilterMagnitude(Ifilter, Qsed, redshift)
                RF_Rmag_app_q = tools.ABFilterMagnitude(Rfilter, Qsed, redshift) + Qoffset
                RF_Gmag_app_q = tools.ABFilterMagnitude(Gfilter, Qsed, redshift) + Qoffset
                RF_Zmag_app_q = tools.ABFilterMagnitude(Zfilter, Qsed, redshift) + Qoffset
                if(redshift<3.9):
                    RF_Umag_app_q = tools.ABFilterMagnitude(Ufilter, Qsed, redshift) + Qoffset
                elif(redshift>=3.9):
                    RF_Umag_app_q = 99
                # calculate the lens magnitude
                veldisp = np.atleast_1d(lens['VELDISP'])
                redshift = lens['ZLENS']
                # Reference Frame Absolute R magnitude
                RF_RMag_abs, _ = lenspop_constructor.EarlyTypeRelations(veldisp)
                RMag_abs = tools.ABFilterMagnitude(Rfilter, Lsed, redshift)
                distMod = cosmo.distmod(redshift).value
                Rmag_app = RMag_abs + distMod
                offset_abs_app = RMag_abs - Rmag_app
                offset_RF_abs = RF_RMag_abs - RMag_abs
                RF_Rmag_app = RF_RMag_abs - offset_abs_app
                # Get filters and calculate magnitudes for each filter:
                RF_Umag_app = tools.ABFilterMagnitude(Ufilter, Lsed, redshift) + offset_RF_abs + distMod
                RF_Gmag_app = tools.ABFilterMagnitude(Gfilter, Lsed, redshift) + offset_RF_abs + distMod
                RF_Imag_app = tools.ABFilterMagnitude(Ifilter, Lsed, redshift) + offset_RF_abs + distMod
                RF_Zmag_app = tools.ABFilterMagnitude(Zfilter, Lsed, redshift) + offset_RF_abs + distMod
                t['u_SDSS_lens'][lens_count] = RF_Umag_app
                t['r_SDSS_lens'][lens_count] = RF_Rmag_app
                t['g_SDSS_lens'][lens_count] = RF_Gmag_app
                t['i_SDSS_lens'][lens_count] = RF_Imag_app
                t['z_SDSS_lens'][lens_count] = RF_Zmag_app
                t['u_SDSS_quasar'][lens_count] = RF_Umag_app_q
                t['r_SDSS_quasar'][lens_count] = RF_Rmag_app_q
                t['g_SDSS_quasar'][lens_count] = RF_Gmag_app_q
                t['i_SDSS_quasar'][lens_count] = RF_Imag_app_q
                t['z_SDSS_quasar'][lens_count] = RF_Zmag_app_q
                lens_count = lens_count + 1
                dot = np.mod(lens_count, total/np.min([79,total])) == 0
                if verbose and dot:
                    print('.', end="")

            # Update the sample by adding the table of calculated magnitude
    	    self.sample.add_columns(t.columns.values())
            self.lenses = self.sample.copy()

        return

# ----------------------------------------------------------------------------

    def gaussian_reweight(self, mean, stdev):
        """
        Add new columns to the table, for the magnitudes in various filters.

        Parameters
        ----------
        mean : float
            The mean of the parent gaussian distribution.
        stdev : float
            The standard deviation of the parent gaussian distribution.
        Returns
        ----------
        None
        Notes
        -----
        * This method adds one column named "weights" in the lens.sample.
        """
        # fit gaussian function
        # fit 2nd degree polinomial function and normalize
        # weight = gaussian(x)/original(x)
        import matplotlib.pyplot as plt
        plt.ioff()
        self.Nlenses=len(self.sample)
        t = Table(np.arange(self.Nlenses).reshape(self.Nlenses, 1), names=('weight',), dtype=('f4',))
        n, bins, patches = plt.hist(self.sample['ZLENS'], bins='auto', normed = True)
        bin_centers = bins[:-1] + 0.5 * (bins[1:] - bins[:-1])
        # calculate polynomial
        def bestFitHist(x, a, b, c, d, e):
            return a*x*x*x*x+b*x*x*x+c*x*x+d*x+e
        param, cov = scipy.optimize.curve_fit(bestFitHist, bin_centers, n)
        xNumbers = np.arange(0, 3, 0.05)
        yHist = bestFitHist(xNumbers, param[0], param[1], param[2], param[3], param[4])
        def gauss_function(x):
            return np.exp(-(x-mean)**2/(2*stdev**2))
        yGauss = gauss_function(xNumbers)
        for (lens, lenscount) in zip(self.sample, range(len(self.sample))):
            redshift = lens['ZLENS']
            histogram = bestFitHist(redshift, param[0], param[1], param[2], param[3], param[4])
            gauss = gauss_function(redshift)
            weight = gauss/histogram
            # Here, rejection sampling
            import random
            if weight<random.random():
                weight = 0
            t['weight'][lenscount] = weight
        self.sample.remove_columns(['weight'])
        self.sample.add_columns(t.columns.values())

        return

# ======================================================================

    def noissify_quasars(self, stdev):
        for lensCount in range(len(self.sample)):
            noise = np.random.normal(scale=stdev)
            self.sample[lensCount]['u_SDSS_quasar'] = self.sample[lensCount]['u_SDSS_quasar'] + noise
            noise = np.random.normal(scale=stdev)
            self.sample[lensCount]['r_SDSS_quasar'] = self.sample[lensCount]['r_SDSS_quasar'] + noise
            noise = np.random.normal(scale=stdev)
            self.sample[lensCount]['g_SDSS_quasar'] = self.sample[lensCount]['g_SDSS_quasar'] + noise
            noise = np.random.normal(scale=stdev)
            self.sample[lensCount]['z_SDSS_quasar'] = self.sample[lensCount]['z_SDSS_quasar'] + noise
            noise = np.random.normal(scale=stdev)
            self.sample[lensCount]['i_SDSS_quasar'] = self.sample[lensCount]['i_SDSS_quasar'] + noise

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
#         print("Lens ",id," has zd,zs = ",lens.ZLENS[0],lens.ZSRC[0])
#         print("and has images with magnifications: ",lens.MAG[0])

# # Look up one system:
#
#     id = 7176527
#     lens = db.get_lens(id)
#
#     if lens is not None:
#         print("Lens ",id," has zd,zs = ",lens.ZLENS[0],lens.ZSRC[0])
#         print("and has images with magnifications: ",lens.MAG[0])

# # To make a mock catalog of SDSS lenses:
#
#     db.select_random(maglim=19.1,area=8000.0,IQ=1.4)
#     db.write_table("OM10_SQLS_mock_lensed_quasars.fits")

# # To select a mock catalog of LSST lenses:
#
#     db.select_random(maglim=23.3,area=20000.0,IQ=0.75)
#     print(db.Nlenses," LSST lenses, with zd = ",db.sample.ZLENS)

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
#     print(db.Nlenses," LSST lenses")
    db.select_random(maglim=23.3,area=18000.0,IQ=0.75)
    print(db.Nlenses," LSST lenses")

    good = db.sample[np.where(db.sample.IMSEP > 1.0)]
    print("Number with imsep > 1.0 arcsec = ",len(good))

    bright = good[np.where(good.APMAG_I < 22.0)]
    print("Number of these with md < 22 = ",len(bright))

    lagged = bright[np.where(np.max(bright.DELAY,axis=1) > 10.0)]
    print("Number of these with time delay > 10 days = ",len(lagged))

    nearby = lagged[np.where((lagged.ZLENS > 0.1) * (lagged.ZLENS < 0.6))]
    print("Number of these with 0.1 < zd < 0.6 = ",len(nearby))

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
    print("No. of OM10 lenses with matching LRG sky positions = ",len(good))
    print("  mean,min,max redshift = ",np.average(good.ZLENS),np.min(good.ZLENS),np.max(good.ZLENS))
    print("No. of OM10 lenses with no matching sky position = ",len(bad),np.min(bad.ZLENS),np.max(bad.ZLENS))
    print("  mean,min,max redshift = ",np.average(bad.ZLENS))

# # To select 10 lenses detectable with PS1 at each epoch:
#
#     db.select_random(maglim=21.4,area=30000.0,IQ=1.0,Nlens=10)
#     print(db.Nlenses," representative PS1 3pi lenses, with zd = ", \
#       db.sample.ZLENS)
#     # print("ugriz = ", \
#     #   db.sample.uMAG_LRG,db.sample.gMAG_LRG,db.sample.rMAG_LRG, \
#     #   db.sample.iMAG_LRG,db.sample.zMAG_LRG)


# 10-sigma detection in a single epoch?
# surveys = PS1-3PI PS1-MDS DES-WL KIDS  HSC-WIDE HSC-DEEP LSST  SDSS-S82x100
# maglims = 21.4    23.3    23.6   22.9  24.9     25.3     23.3  21.3
# areas   = 30000   84      5000   1500  1500     30       20000 30000        # survey area in sq deg
# psfs    = 1.0     1.0     0.9    0.7   0.75     0.75     0.75  1.4          # PSF FWHM in arcsec
# Note that these numbers give lower yields that OM10, by about a factor of 2:
# this is just due to the single epoch requirement, in the stacked images we
# should be able to go much deeper.

# ======================================================================
