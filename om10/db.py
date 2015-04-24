# ======================================================================

import pyfits,sys,os,subprocess
import numpy as np
import os

# from astropy.table import Table

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
        self.catalog = catalog

        # Make a FITS table from the supplied text catalogs, if required:
        if generate:
            self.make_table()
            self.catalog = os.path.expandvars("$OM10_DIR/data/qso_mock.fits")
            self.write_table(self.catalog)

        # If given a catalog, read it in:
        if self.catalog is not None:
            self.lenses = pyfits.getdata(self.catalog)
            # self.lenses = Table.read(self.catalog,format='fits')
            # print self.lenses

        # No down-sampling has been done yet:
        self.sample = None

        # Count lenses:
        try: self.Nlenses = len(self.lenses.LENSID)
        except: self.Nlenses = 0

        return

    # ------------------------------------------------------------------

    def make_table(self):

        # Read in lens data from original text catalog:

        lensfile = os.path.expandvars("$OM10_DIR/data/qso_mock_result+lensphot+lensgeom.dat")
        d = np.loadtxt(lensfile)
        if vb: print "om10.DB: read in lens data from ",lensfile

        # Put lens parameters in columns:

        nim     = np.array(d[:, 0],dtype=np.int)
        zd      = np.array(d[:, 1])
        sigma   = np.array(d[:, 2])
        zs      = np.array(d[:, 3])
        ms      = np.array(d[:, 4])
        m3      = np.array(d[:, 5])
        dtheta  = np.array(d[:, 6])
        epsd    = np.array(d[:, 7])
        phid    = np.array(d[:, 8])
        gammax  = np.array(d[:, 9])
        phix    = np.array(d[:,10])
        xs      = np.array(d[:,11])
        ys      = np.array(d[:,12])
        type    = np.array(d[:,13],dtype=np.int)
        id      = np.array(d[:,14],dtype=np.int)
        Dd      = np.array(d[:,15])
        DLd     = np.array(d[:,16])
        absmd   = np.array(d[:,17])
        md      = np.array(d[:,18])
        Kcorrd  = np.array(d[:,19])
        Ds      = np.array(d[:,20])
        Dds     = np.array(d[:,21])
        Sigcrit = np.array(d[:,22])
        DLs     = np.array(d[:,23])

        size = len(id)

        if vb: print "om10.DB: stored lens parameters for",size,"lenses"

        if vb: print "om10.DB: pulling out image configurations..."

        xs  = np.zeros([size])
        ys  = np.zeros([size])
        xi  = np.zeros([size,4])
        yi  = np.zeros([size,4])
        ti  = np.zeros([size,4])
        mui = np.zeros([size,4])

        # Loop through image log file, line by line, filling up image
        # parameter arrays:

        logfile = os.path.expandvars("$OM10_DIR/data/qso_mock_log.dat")

        with open(logfile, "r") as file:

            count = 0
            for line in file:

                x = line.split()

                # print "line = ",x," length = ",len(x)

                if (len(x) == 5):
                    k = np.where(id == int(x[4]))
                    k = k[0]
                    xs[k] = float(x[1])
                    ys[k] = float(x[2])
                    # print "k = ",k," xs,ys = ",xs[k],ys[k]
                    i = 0
                    count += 1
                else:
                   xi[k,i]   = float(x[0])
                   yi[k,i]   = float(x[1])
                   mui[k,i]  = float(x[2])
                   ti[k,i]   = float(x[3])
                   # print "i = ",i," xi,yi,mui,ti = ",xi[k,i],yi[k,i],mui[k,i],ti[k,i]
                   i += 1

        file.close()

        # Check numbers:
        if (count != size):
            print "ERROR: found %d lenses in logfile and %d in result file" % count,size

        if vb: print "om10.DB: read in image parameters for",size,"lenses"

        # Now package into a pyfits table:

        self.lenses = []
        self.lenses.append(pyfits.Column(name='LENSID  ',format='J       ',array=id))
        self.lenses.append(pyfits.Column(name='FLAGTYPE',format='I       ',array=type))
        self.lenses.append(pyfits.Column(name='NIMG    ',format='I       ',array=nim))
        self.lenses.append(pyfits.Column(name='ZLENS   ',format='D       ',array=zd))
        self.lenses.append(pyfits.Column(name='VELDISP ',format='D       ',array=sigma))
        self.lenses.append(pyfits.Column(name='ELLIP   ',format='D       ',array=epsd))
        self.lenses.append(pyfits.Column(name='PHIE    ',format='D       ',array=phid))
        self.lenses.append(pyfits.Column(name='GAMMA   ',format='D       ',array=gammax))
        self.lenses.append(pyfits.Column(name='PHIG    ',format='D       ',array=phix))
        self.lenses.append(pyfits.Column(name='ZSRC    ',format='D       ',array=zs))
        self.lenses.append(pyfits.Column(name='XSRC    ',format='D       ',array=xs))
        self.lenses.append(pyfits.Column(name='YSRC    ',format='D       ',array=ys))
        self.lenses.append(pyfits.Column(name='MAGI_IN ',format='D       ',array=ms))
        self.lenses.append(pyfits.Column(name='MAGI    ',format='D       ',array=m3))
        self.lenses.append(pyfits.Column(name='IMSEP   ',format='D       ',array=dtheta))
        self.lenses.append(pyfits.Column(name='XIMG    ',format='4D      ',array=xi))
        self.lenses.append(pyfits.Column(name='YIMG    ',format='4D      ',array=yi))
        self.lenses.append(pyfits.Column(name='MAG     ',format='4D      ',array=mui))
        self.lenses.append(pyfits.Column(name='DELAY   ',format='4D      ',array=ti))
        self.lenses.append(pyfits.Column(name='KAPPA   ',format='4D      ',array=np.zeros([size,4])))
        self.lenses.append(pyfits.Column(name='FSTAR   ',format='4D      ',array=np.zeros([size,4])))
        self.lenses.append(pyfits.Column(name='DD      ',format='D       ',array=Dd))
        self.lenses.append(pyfits.Column(name='DDLUM   ',format='D       ',array=DLd))
        self.lenses.append(pyfits.Column(name='ABMAG_I ',format='D       ',array=absmd))
        self.lenses.append(pyfits.Column(name='APMAG_I ',format='D       ',array=md))
        self.lenses.append(pyfits.Column(name='KCORR   ',format='D       ',array=Kcorrd))
        self.lenses.append(pyfits.Column(name='DS      ',format='D       ',array=Ds))
        self.lenses.append(pyfits.Column(name='DDS     ',format='D       ',array=Dds))
        self.lenses.append(pyfits.Column(name='SIGCRIT ',format='D       ',array=Sigcrit))
        self.lenses.append(pyfits.Column(name='DSLUM   ',format='D       ',array=DLs))
        self.lenses.append(pyfits.Column(name='L_I     ',format='D       ',array=np.zeros([size])))
        self.lenses.append(pyfits.Column(name='REFF    ',format='D       ',array=np.zeros([size])))
        self.lenses.append(pyfits.Column(name='REFF_T  ',format='D       ',array=np.zeros([size])))

        return

    # ------------------------------------------------------------------

    def write_table(self,catalog):

        try: os.remove(catalog)
        except OSError: pass

        if self.sample is None:
            pyfits.writeto(catalog,self.lenses)
        else:
            pyfits.writeto(catalog,self.sample)

        if vb: print "om10.DB: wrote catalog of ",self.Nlenses," OM10 lenses to file at "+catalog

        return

    # ------------------------------------------------------------------

    def export_to_cpt(self,pars,cptfile):

        try: os.remove(cptfile)
        except OSError: pass

        if self.sample is None:
            data = self.lenses
        else:
            data = self.sample

        # Define labels and ranges for cpt files:
        labels = {}
        ranges = {}

        labels['ZLENS']           = '$z_{\\rm d}$,  '
        ranges['ZLENS']           = '0.0,2.5,    '

        labels['ZSRC']            = '$z_{\\rm s}$,  '
        ranges['ZSRC']            = '0.0,6.0,    '

        labels['IMSEP']           = '$\\Delta \\theta / "$,  '
        ranges['IMSEP']           = '0.0,5.0,    '

        labels['APMAG_I']         = '$i_{\\rm d}$ / mag,  '
        ranges['APMAG_I']         = '15.0,24.5,    '

        labels['MAGI']            = '$i_{\\rm 2/3}$ / mag,  '
        ranges['MAGI']            = '15.0,24.5,    '


        # Write header lines:

        hline1 = '# '
        hline2 = '# '
        for par in pars:
            hline1 += labels[par]
            hline2 += ranges[par]

        # Prepare data for writing:

        outbundle = np.zeros([self.Nlenses,len(pars)])

        for k,par in enumerate(pars):
            outbundle[:,k] = data[par]

        # The actual writing:

        output = open(cptfile,'w')
        output.write("%s\n" % hline1)
        output.write("%s\n" % hline2)
        output.close()
        np.savetxt('junk', outbundle)
        cat = subprocess.call("cat junk >> " + cptfile, shell=True)
        rm = subprocess.call("rm junk", shell=True)
        if cat != 0 or rm != 0:
          print "Error: write subprocesses failed in some way :-/"
          sys.exit()

        if vb: print "om10.DB: wrote a ",len(pars),"-column plain text file of ",len(data)," OM10 lenses to file at "+cptfile

        return

    # ------------------------------------------------------------------

    def get_lens(self,ID):

        try: rec = self.lenses[self.lenses.LENSID == ID]
        except: rec = None

        return rec

    # ------------------------------------------------------------------

    def select_random(self,Nlens=None,maglim=99.0,area=100000.0,IQ=0.0):

        # Select all lenses that meet the rough observing criteria:

        try:
            sample = self.lenses.copy()
            sample = sample[sample.MAGI < maglim]
            sample = sample[sample.IMSEP > 0.67*IQ]
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

    def get_sky_positions(self,dmag=0.2,dz=0.2):

        LRGfile = os.path.expandvars("$OM10_DIR/data/CFHTLS_LRGs.txt")
	try:
            d = np.loadtxt(LRGfile)
        except IOError:
            print "Cannot find LRG catalog!"        
        if vb: print "om10.DB: read in LRG sky position data from ",LRGfile

        # Put LRG parameters in LRG structure:

        self.LRGs = {}
        self.LRGs['RA']       = np.array(d[:, 0])
        self.LRGs['DEC']      = np.array(d[:, 1])
        self.LRGs['redshift'] = np.array(d[:, 2])
        self.LRGs['mag_i']    = np.array(d[:, 6])

        print "Mean LRG RA,DEC,z,i = ",np.average(self.LRGs['RA']),np.average(self.LRGs['DEC']),np.average(self.LRGs['redshift']),np.average(self.LRGs['mag_i']);

        # Bin LRGs in mag_i and redshift, and record bin numbers for each one:

        imin,imax = np.min(self.LRGs['mag_i']),np.max(self.LRGs['mag_i'])
        nibins = int((imax - imin)/dmag) + 1
        ibins = np.linspace(imin, imax, nibins)
        self.LRGs['ivals'] = np.digitize(self.LRGs['mag_i'],ibins)
        self.LRGs['ibins'] = ibins

        zmin,zmax = np.min(self.LRGs['redshift']),np.max(self.LRGs['redshift'])
        nzbins = int((zmax - zmin)/dz) + 1
        zbins = np.linspace(zmin, zmax, nzbins)
        self.LRGs['zvals'] = np.digitize(self.LRGs['redshift'],zbins)
        self.LRGs['zbins'] = zbins

        if vb: print "om10.DB: number of LRGs stored = ",len(self.LRGs['redshift'])

        return

    # ------------------------------------------------------------------

    def assign_sky_positions(self):

        reallyverbose = False

        # Prepare new columns for LRG properties:
        size = len(self.lenses.LENSID)
        dummy = np.zeros(size)
        # For some reason the following lines don't work - we need to
        # switch to a better dataframe system, like pandas, so that we can
        # add table columns in the way we need to...
        self.lenses['RA'] = dummy
        self.lenses['DEC'] = dummy

        # First digitize the lenses to the same bins as the LRGs:

        ii = np.digitize(self.lenses.APMAG_I,self.LRGs['ibins'])
        iz = np.digitize(self.lenses.ZLENS,self.LRGs['zbins'])

        # Loop over lenses, finding all LRGs in its bin:

        for k in range(len(self.lenses)):
            iindex = list(np.where(self.LRGs['ivals'] == ii[k])[0])
            zindex = list(np.where(self.LRGs['zvals'] == iz[k])[0])
            i = list(set(iindex).intersection(set(zindex)))

            if len(i) == 0:
                if reallyverbose: print "WARNING: Lens ",k," has no matching LRG..."
                self.lenses.RA[k]       = -99
                self.lenses.DEC[k]      = -99
                pass

            else:
                if reallyverbose: print "Lens ",k," has ",len(i)," neighbour LRGs... "

                # Compute distances to find nearest neighbour:
                i0 = self.lenses.APMAG_I[k]
                z0 = self.lenses.ZLENS[k]
                R = np.sqrt((self.LRGs['mag_i'][i]-i0)**2 + (self.LRGs['redshift'][i])**2)
                best = np.where(R == np.min(R))
                ibest = i[best[0][0]]
                if reallyverbose:
                    print "  LRG  i,z: ",self.LRGs['mag_i'][ibest],self.LRGs['redshift'][ibest]

                self.lenses.RA[k]       = self.LRGs['RA'][ibest]
                self.lenses.DEC[k]      = self.LRGs['DEC'][ibest]

            if reallyverbose:
                print "  Lens i,z: ",self.lenses.APMAG_I[k],self.lenses.ZLENS[k]
                print "  Lens RA,DEC: ",self.lenses.RA[k],self.lenses.DEC[k]

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
