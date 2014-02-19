# ======================================================================

import numpy,pyfits,sys,os,subprocess,math

import om10

vb = True

# ======================================================================
# define useful functions for photometry

conv = 1./(numpy.log(2.)*(2**0.5))
from numpy import exp, sin, cos

def G(x,dx):
    gau = exp(-0.5*(x/dx)**2.)/(dx*(2.*numpy.pi)**0.5)
    return gau

def ks(n):
    kaser = 2.*n-1/3.
    return kaser

def sernorm(Re,n):
    sn = 2.*numpy.pi*(Re**2.)*(ks(n)**(-2*n))*math.gamma(2*n)
    return sn

def Sersic(R,n):
    SB = exp(-ks(n)*R**(1./n))
    return SB

def flaser(x,y,flat,pa,n):
    SB = Sersic(((cos(pa)*x+sin(pa)*y)**2. + (-sin(pa)*x+cos(pa)*y)**2./(flat**2.))**0.5,n)
    return SB

# ======================================================================

class Imager(object):

    """
    NAME
        Image

    PURPOSE
        Realize an OM10 lens system as a set of pixelated images, with 
        approximate properties for a given survey.

    COMMENTS

    INITIALISATION
        survey            Survey name (optional, default=LSST).
    
    METHODS
        reset canvas      in case the survey is not recognised, call after setting the new parameter values
        target            read in properties from the OM10 catalogue
        write             write image (sci, var) arrays to fits files       
        make              Generate the images and write them
        
    BUGS

    AUTHORS
      This file is part of the OM10 project, distributed under the
      GPL v2, by Phil Marshall (KIPAC) & Adriano Agnello (UCSB). 
      Please cite: Oguri & Marshall (2010), MNRAS, 405, 2579.

    HISTORY
      2013-11-15  Started Marshall & Agnello (UCSB)
    """
    # ------------------------------------------------------------------

    def __init__(self,survey=None):

        self.name = 'OM10 imager'
        # default initialisations
        self.pixscale = 0.25 # in arcseconds
        self.meanIQ = 1. # psf FWHM, in arcseconds
        self.meandepth = 20.0
        self.errdepth = 0.3

        if survey == None:
           self.survey = 'LSST'
        else:
           self.survey = survey
             
        # Set up observing parameters:   
        # they're accessible, in case one wants to bypass the survey choice and feed them:   
        if survey=='LSST':
            self.pixscale = 0.2 # in arcseconds
            self.meanIQ = 0.75 # psf FWHM, in arcseconds
            self.meandepth = 23.3
            self.errdepth = 0.3

        elif survey == 'PS1':
            self.pixscale = 0.25
            self.meanIQ = 1.0
            self.meandepth = 21.4
            self.errdepth = 0.3 

        else:
            raise "ERROR: unrecognised survey "+survey+": default config for parameters and canvas"       

        self.fov = 10.0 # arcsec
        self.imsize = int(self.fov/self.pixscale) +1
        self.midpsf = 6 #midpoint of the psf array, size set arbitrarily just for convenience!
        self.psfsize = 2*self.midpsf-1
        self.canvas = self.imsize +2*self.psfsize -2 #this is the 'canvas' size
        
        self.image = numpy.zeros([self.imsize,self.imsize])
        self.lens_galaxy_image = numpy.zeros([self.imsize,self.imsize])
        # magnitude fluctuation map, needed for noise later
        self.fluct = exp(-np.log(10.)*np.random.normal(0,self.errdepth,(self.imsize,self.imsize)))-1.

        self.x,self.y = numpy.mgrid[0:canvas,0:canvas] #coord.grid for the canvas
        self.center = int(self.canvas/2.)+1 # center of the raw image grid
        self.newcen = int(self.imsize/2)+1 # center of the blurred image grid
        self.xpsf,self.ypsf = numpy.mgrid[0:self.psfsize,0:self.psfsize] # coord.grid for the psf
        self.pixpsf=conv*self.meanIQ/self.pixscale # width of the Gaussian psf, given the FWHM
        
    # define interpolated psf grid, used in convolutions and point sources;
    # this is an array psfsize*psfsize!
        self.psf = (9./16.)*G(self.xpsf-self.midpsf,self.pixpsf)G(self.ypsf-self.midpsf,self.pixpsf)+
        +(3./32.)*(G(self.xpsf-self.midpsf-1,self.pixpsf)G(self.ypsf-self.midpsf,self.pixpsf)+G(self.xpsf-self.midpsf+1,self.pixpsf)G(self.ypsf-self.midpsf,self.pixpsf)
        +G(self.xpsf-self.midpsf,self.pixpsf)G(self.ypsf-self.midpsf-1,self.pixpsf)+G(self.xpsf-self.midpsf,self.pixpsf)G(self.ypsf-self.midpsf+1,self.pixpsf))
        +(1/64)*(G(self.xpsf-self.midpsf-1,self.pixpsf)G(self.ypsf-self.midpsf-1,self.pixpsf)+G(self.xpsf-self.midpsf-1,self.pixpsf)G(self.ypsf-self.midpsf+1,self.pixpsf)
        +G(self.xpsf-self.midpsf+1,self.pixpsf)G(self.ypsf-self.midpsf-1,self.pixpsf)+G(self.xpsf-self.midpsf+1,self.pixpsf)G(self.ypsf-self.midpsf+1,self.pixpsf))

        return
        

    def resetcanvas(self):     
        self.imsize = int(self.fov/self.pixscale) +1
        self.psfsize = 2*self.midpsf-1
        self.canvas = self.imsize +2*self.psfsize -2 
        
        self.image = numpy.zeros([self.imsize,self.imsize])
        self.lens_galaxy_image = numpy.zeros([self.imsize,self.imsize])
        self.fluct = exp(-np.log(10.)*np.random.normal(0,self.errdepth,(self.imsize,self.imsize)))-1.

        self.x,self.y = numpy.mgrid[0:canvas,0:canvas]
        self.center = int(self.canvas/2.)+1
        self.newcen = int(self.imsize/2)+1
        self.xpsf,self.ypsf = numpy.mgrid[0:self.psfsize,0:self.psfsize]
        self.pixpsf=conv*self.meanIQ/self.pixscale
        
        self.psf = (9./16.)G(self.xpsf-self.midpsf,self.pixpsf)G(self.ypsf-self.midpsf,self.pixpsf)+
        +(3./32.)(G(self.xpsf-self.midpsf-1,self.pixpsf)G(self.ypsf-self.midpsf,self.pixpsf)+G(self.xpsf-self.midpsf+1,self.pixpsf)G(self.ypsf-self.midpsf,self.pixpsf)
        +G(self.xpsf-self.midpsf,self.pixpsf)G(self.ypsf-self.midpsf-1,self.pixpsf)+G(self.xpsf-self.midpsf,self.pixpsf)G(self.ypsf-self.midpsf+1,self.pixpsf))
        +(1/64)(G(self.xpsf-self.midpsf-1,self.pixpsf)G(self.ypsf-self.midpsf-1,self.pixpsf)+G(self.xpsf-self.midpsf-1,self.pixpsf)G(self.ypsf-self.midpsf+1,self.pixpsf)
        +G(self.xpsf-self.midpsf+1,self.pixpsf)G(self.ypsf-self.midpsf-1,self.pixpsf)+G(self.xpsf-self.midpsf+1,self.pixpsf)G(self.ypsf-self.midpsf+1,self.pixpsf))

        return
        
    # ------------------------------------------------------------------

    def target(self,lens):

        # Read positions and set up WCS:
        self.RA = 0.0
        self.DEC = 36.0
        # BUG: these should be read from lens!
        # BUG: if preceded by self., they must go in the init! otherwise just make them local.
        
        # Set up WCS:
        # self.set_WCS() 

        # here we get Reff, band-magnitudes, p.a., flattening, positions...
#        logfile = os.path.expandvars("$OM10_DIR/data/qso_mock_log.dat") ...


        return


    # ------------------------------------------------------------------

    def make(self,filters=['i'],Nepochs=1):
        # sky flux in nano-maggies, needed in the simple recipe for noise-map;
        # it's just an internal variable!
        skyflux = exp(9-0.4*np.log(10.)*self.meandepth)

        # Loop over epochs, making images:
        # MIND THE GAP: initialisations missing at the moment
        for k in range(Nepochs):
            
            # Make raw image
            lgalflux = # in nano-maggies, set it from the bla-band magnitude read from target
            reff=self.Re/self.pixscale
            # normalisation: central flux in nanomaggies/pixscale^2
            S0 = lgalflux/sernorm(reff,n)
            self.sbraw = flaser((self.x-self.center)/reff,(self.y-self.center)/reff,self.flat,self.pa,4.)

            # Convolve lens galaxy image with PSF:
            # alternative without nested for loop highly desirable!
            # numpy.fft would avoid that and not invoke scipy

            for i in range(0,self.imsize-1)
                for j in range(0,self.imsize-1)
                    for i1 in range(0,psfsize-1)
                        for i2 in range(0,psfsize-1)
                            self.lens_galaxy_image += self.sbraw[i+i1,j+i2]*self.psf[i1,i2]
            
            # Paint in point-source PSFs:

            for kq in range(quasims-1)
                ipos = int(self.qim[kq][1]/self.pixscale + self.newcen)
                jpos = int(self.qim[kq][2]/self.pixscale + self.newcen)
                imin = max(0,ipos-self.midpsf)
                jmin = max(0,jpos-self.midpsf)
                imax = min(self.canvas-1,ipos+self.midpsf)
                jmax = min(self.canvas-1,jpos+self.midpsf)
                dx = self.qim[kq][1]/self.pixscale + self.newcen - ipos
                dy = self.qim[kq][2]/self.pixscale + self.newcen - jpos
                # unif.dither and drift the psf grid by (dx,dy), for each quasar image
                # can we avoid the nested loop?
                for i1 in range(imin,imax)
                    for i2 in range(jmin,jmax)
                    self.image += (self.psf[self.midpsf+i1-ipos,self.midpsf+i2-jpos]*(1.-abs(dx))*(1.-abs(dy))
                                  + self.psf[self.midpsf+i1-ipos-cmp(dx,0),self.midpsf+i2-jpos-cmp(dy,0)]*abs(dx*dy)
                                  + self.psf[self.midpsf+i1-ipos-cmp(dx,0),self.midpsf+i2-jpos]*abs(dy)*(1-abs(dx))
                                  + self.psf[self.midpsf+i1-ipos,self.midpsf+i2-jpos-cmp(dy,0)]*abs(dx)*(1-abs(dy)))*qflux[kq]


            self.image += self.lens_galaxy_image*lgalflux + skyflux

            # Image with noise:
            self.sci = self.image*self.fluct + self.image

            # Noise map:
            self.var = self.image*(exp(9.-0.4*np.log(10)*self.errdepth)-1.)

            # Write out image and weight map to file:
            self.write()


        return

    # ------------------------------------------------------------------

    def write(self):

        # Set filenames:
        scifile = 'test_sci.fits'
        varfile = 'test_var.fits'

        # Start an HDU object:
        hdus = pyfits.PrimaryHDU(self.sci)
        hduv = pyfits.PrimaryHDU(self.var)

        # Add WCS keywords to the header, e.g.:
        # hdus.header.set('CRVAL1',0.0,'Right ascension (J2000 degrees)')
        
        # Write them out
        hdus.writeto(scifile)   
        hduv.writeto(varfile)      

        return


# ======================================================================

if __name__ == '__main__':

# Some examples!
                
    db = om10.DB(catalog=os.path.expandvars("$OM10_DIR/data/qso_mock.fits"))

# Get one lens:
 
    id = 7176527
    lens = db.get_lens(id)
    
# Set up imager:

    imager = om10.Imager(survey='PS1')

    # imager.set('fov',6.0) # Needs to change imsize too!

# Make a single image of the lens:

    imager.target(lens)

    imager.make(filters=['i'],Nepochs=1)


# 10-sigma detection in a single epoch?
# surveys = PS1-3PI PS1-MDS DES-WL KIDS  HSC-WIDE HSC-DEEP LSST  SDSS-S82x100
# maglims = 21.4    23.3    23.6   22.9  24.9     25.3     23.3  21.3        
# areas   = 30000   84      5000   1500  1500     30       20000 30000        # survey area in sq deg
# psfs    = 1.0     1.0     0.9    0.7   0.75     0.75     0.75  1.4          # PSF FWHM in arcsec
# Note that these numbers give lower yields that OM10, by about a factor of 2:
# this is just due to the single epoch requirement, in the stacked images we 
# should be able to go much deeper.

# ======================================================================


