# ======================================================================

import sys,os,subprocess
import numpy as np
import os
from numpy import *
import math
from astropy.table import Table, hstack

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
            #self.lenses = pyfits.getdat(self.catalog)
            self.lenses = Table.read(self.catalog,format='fits')
            #print self.lenses

        # No down-sampling has been done yet:
        self.sample = None

        # Count lenses:
        try: self.Nlenses = len(self.lenses['LENSID'])
        except: self.Nlenses = 0

        return

    # ------------------------------------------------------------------

    def make_table(self):
        """
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
        """
        return

    # ------------------------------------------------------------------

    def write_table(self,catalog):
        """
        try: os.remove(catalog)
        except OSError: pass

        if self.sample is None:
            pyfits.writeto(catalog,self.lenses)
        else:
            pyfits.writeto(catalog,self.sample)

        if vb: print "om10.DB: wrote catalog of ",self.Nlenses," OM10 lenses to file at "+catalog
        """
        return

    # ------------------------------------------------------------------

    def export_to_cpt(self,pars,cptfile):
        """
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
        """
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
            sample = self.lenses.copy()
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

    def assign_sky_positions(self,verbose=False):

        reallyverbose = verbose

        # Prepare new columns for LRG properties:
        size = len(self.lenses['LENSID'])
        dummy = np.zeros(size)
        print size
        #cols = []
        #cols.append(pyfits.Column(name='RA',format='D       ',array=dummy))
        #cols.append(pyfits.Column(name='DEC',format='D       ',array=dummy))

        self.lenses['RA'] = 0.0
        self.lenses['DEC'] = 0.0

        #self.lenses.columns.add_col(pyfits.Column(name='RA',format='D       ',array=dummy))
        #self.lenses.columns.add_col(pyfits.Column(name='DEC',format='D       ',array=dummy))
        #new_cols = pyfits.ColDefs(cols)
        #self.lenses = pyfits.FITS_rec.from_columns(orig_cols + new_cols)

        # First digitize the lenses to the same bins as the LRGs:

        ii = np.digitize(self.lenses['APMAG_I'],self.LRGs['ibins'])
        iz = np.digitize(self.lenses['ZLENS'],self.LRGs['zbins'])

        # Loop over lenses, finding all LRGs in its bin:

        for k in range(len(self.lenses)):
            iindex = list(np.where(self.LRGs['ivals'] == ii[k])[0])
            zindex = list(np.where(self.LRGs['zvals'] == iz[k])[0])
            i = list(set(iindex).intersection(set(zindex)))

            if len(i) == 0:
                if reallyverbose: print "WARNING: Lens ",k," has no matching LRG..."
                self.lenses['RA'][k]       = -99
                self.lenses['DEC'][k]      = -99
                pass

            else:
                if reallyverbose: print "Lens ",k," has ",len(i)," neighbour LRGs... "

                # Compute distances to find nearest neighbour:
                i0 = self.lenses['APMAG_I'][k]
                z0 = self.lenses['ZLENS'][k]
                R = np.sqrt((self.LRGs['mag_i'][i]-i0)**2 + (self.LRGs['redshift'][i])**2)
                best = np.where(R == np.min(R))
                ibest = i[best[0][0]]
                if reallyverbose:
                    print "  LRG  i,z: ",self.LRGs['mag_i'][ibest],self.LRGs['redshift'][ibest]

                self.lenses['RA'][k]       = self.LRGs['RA'][ibest]
                self.lenses['DEC'][k]      = self.LRGs['DEC'][ibest]

            if reallyverbose:
                print "  Lens i,z: ",self.lenses['APMAG_I'][k],self.lenses['ZLENS'][k]
                print "  Lens RA,DEC: ",self.lenses['RA'][k],self.lenses['DEC'][k]

        return

# ----------------------------------------------------------------------------

    def makeSimCatalog(self):
        raise NotImplementedError
        return

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

        ######
        zmin=0.
        zmaxq=5.
        zmaxg=1.
        imin=15.
        imax=22.
        sming=80.
        smaxg=400.
        binsz=50
        binsobs=20
        step1=float((zmaxg-zmin)/binsz)
        step2=float((smaxg-sming)/binsobs)
        step3=float((zmaxq-zmin)/binsz)
        step4=float((imax-imin)/binsobs)

        RSDSS=1.5 #aperture size for SDSS 1.5 arcseconds
        ng=zeros((500,200))
        nq=zeros((500,200))
        mlrg=zeros((50,20,14))
        slrg=zeros((50,20,14))
        mqso=zeros((50,20,9))
        sqso=zeros((50,20,9))

        ##########

        for i in range(lrg.shape[0]):

            sigma=lrg[i,1]*(2.45/2.0)**0.5*(lrg[i,3]/RSDSS)**-0.066
            kz=math.ceil((lrg[i,0]-zmin)/step1)-1
            ko=math.ceil((sigma-sming)/step2)-1

            if -1<kz<50 and -1<ko<20:
               ng[kz,ko]+=1
               mlrg[kz,ko,:]+=lrg[i,:]
               slrg[kz,ko,:]+=lrg[i,:]**2
               #slrg[kz,ko]+=


        #print ng[16,9]
        #meang=mlrg[16,9,:]/ng[16,9]
        #print meang
        #print sqrt(slrg[16,9,:]/ng[16,9]-meang**2)

        ###########

        for j in range(qso.shape[0]):

            kz=math.ceil((qso[j,0]-zmin)/step3)-1
            ko=math.ceil((qso[j,3]-imin)/step4)-1

            if -1<kz<50 and -1<ko<20:
               nq[kz,ko]+=1
               mqso[kz,ko,:]+=qso[j,:]
               sqso[kz,ko,:]+=qso[j,:]**2


        #print sum(nq)

        #meanq=mqso[21,14,:]/nq[21,14]
        #print meanq
        #print sqrt(sqso[21,14,:]/nq[21,14]-meanq**2)
        size = len(self.lenses)
        if verbose: print 'size is ', size

        lens_props = ['MAGG_LENS','MAGR_LENS','MAGI_LENS','MAGZ_LENS', \
        'MAGW1_LENS','MAGW2_LENS','MAGW3_LENS','MAGW4_LENS', 'SDSS_FLAG_LENS']

        src_props = ['MAGG_SRC','MAGR_SRC','MAGI_SRC','MAGZ_SRC', \
        'MAGW1_SRC','MAGW2_SRC','MAGW3_SRC','MAGW4_SRC', 'SDSS_FLAG_SRC']

        tmp_lens = Table(np.zeros((size,len(lens_props)),dtype='>f8'),names=lens_props)
        tmp_src = Table(np.zeros((size,len(src_props)),dtype='>f8'),names=src_props)

        if verbose: print 'setup done'

        if Nmax == None:
            Nmax = size
        for k in range(Nmax):

            if verbose: print 'analysing lens', k
            z_om10_g = self.lenses['ZLENS'][k]
            sigma_om10 = self.lenses['VELDISP'][k]
            z_om10_q = self.lenses['ZSRC'][k]
            i_q = self.lenses['MAGI'][k]

            kz_g=math.ceil((z_om10_g-zmin)/step1)-1
            ko_g=math.ceil((sigma_om10-sming)/step2)-1
            if kz_g<50 and -1<ko_g<20 and ng[kz_g,ko_g]>0:

               #print 'non-zero'
               tmp_lens['SDSS_FLAG_LENS'][k] = 0
               meang=mlrg[kz_g,ko_g,:]/ng[kz_g,ko_g]
               errorg=sqrt(slrg[kz_g,ko_g,:]/ng[kz_g,ko_g]-meang**2)
               om10g=meang+(-1.+2.*random.rand(14))*errorg
            else:
               if kz_g>49 or ko_g>19 or ko_g<0:
                  #print 'out of the region!'
                  tmp_lens['SDSS_FLAG_LENS'][k] = 2
               else:
                  tmp_lens['SDSS_FLAG_LENS'][k] = 1
                  #print 'zero'
               meangg=zeros((50,20,14))
               errorgg=zeros((50,20,14))
               disg=zeros((50,20))
               wmeang=zeros(14)
               werrorg=zeros(14)
               wg=0

               for i in range(50):
                   for j in range(20):
                       if ng[i,j]>0:
                          meangg[i,j,:]=mlrg[i,j,:]/ng[i,j]
                          errorgg[i,j,:]=sqrt(slrg[i,j,:]/ng[i,j]-meangg[i,j,:]**2)
                       else:
                          meangg[i,j,:]=mlrg[i,j,:]
                          errorgg[i,j,:]= meangg[i,j,:]
                       disg[i,j]=abs(z_om10_g-(i*step1+step1/2))+abs(sigma_om10-(j*step2+step2/2))
                       wmeang += meangg[i,j,:]*e**(-disg[i,j])
                       wg +=e**(-disg[i,j])
                       werrorg += errorgg[i,j,:]*e**(-disg[i,j])

               om10g=wmeang/wg+(-1.+2.*random.rand(14))*werrorg/wg

            [tmp_lens['MAGG_LENS'][k], tmp_lens['MAGR_LENS'][k], tmp_lens['MAGI_LENS'][k], tmp_lens['MAGZ_LENS'][k],
             tmp_lens['MAGW1_LENS'][k], tmp_lens['MAGW2_LENS'][k], tmp_lens['MAGW3_LENS'][k], tmp_lens['MAGW4_LENS'][k]] = om10g[6:]

            #print "Lens Properties"
            #print "redshift  sigma  g Reff  r   i    z   g mag   r mag   i mag   z   w1    w2    w3   w4"
            #print om10g


            ######################################3


            kz_q=math.ceil((z_om10_q-zmin)/step3)-1
            ko_q=math.ceil((i_q-imin)/step4)-1
            if kz_q<50 and -1<ko_q<20 and nq[kz_q,ko_q]>0:

               #print 'non-zero'
               tmp_src['SDSS_FLAG_SRC'][k] = 0
               meanq=mqso[kz_q,ko_q,:]/nq[kz_q,ko_q]
               errorq=sqrt(sqso[kz_q,ko_q,:]/nq[kz_q,ko_q]-meanq**2)
               om10q=meanq+(-1.+2.*random.rand(9))*errorq
            else:
               if kz_q>49 or ko_q>19 or ko_q<0:
                  tmp_src['SDSS_FLAG_SRC'][k] = 2
                  #print 'out of the region!'
               else:
                  tmp_src['SDSS_FLAG_SRC'][k] = 1
                  #print 'zero'
               meanqq=zeros((50,20,9))
               errorqq=zeros((50,20,9))
               disq=zeros((50,20))
               wmeanq=zeros(9)
               werrorq=zeros(9)
               wq=0

               for i in range(50):
                   for j in range(20):
                       if nq[i,j]>0:
                          meanqq[i,j,:]=mqso[i,j,:]/nq[i,j]
                          errorqq[i,j,:]=sqrt(sqso[i,j,:]/nq[i,j]-meanqq[i,j,:]**2)
                       else:
                          meanqq[i,j,:]=mqso[i,j,:]
                          errorqq[i,j,:]= meanqq[i,j,:]
                       disq[i,j]=abs(z_om10_q-(i*step3+step3/2))+abs(i_q-(j*step4+step4/2))
                       wmeanq += meanqq[i,j,:]*e**(-disq[i,j])
                       wq +=e**(-disq[i,j])
                       werrorq += errorqq[i,j,:]*e**(-disq[i,j])

               om10q=wmeanq/wq+(-1.+2.*random.rand(9))*werrorq/wq

            [tmp_src['MAGG_SRC'][k], tmp_src['MAGR_SRC'][k], tmp_src['MAGI_SRC'][k], tmp_src['MAGZ_SRC'][k],
            tmp_src['MAGW1_SRC'][k], tmp_src['MAGW2_SRC'][k], tmp_src['MAGW3_SRC'][k], tmp_src['MAGW4_SRC'][k]] = om10q[1:]

            #print "QSO Properties"
            #print "redshift   g     r     i     z     w1      w2       w3     w4"
            #print om10q

        self.lenses = hstack([self.lenses,tmp_lens,tmp_src])

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
