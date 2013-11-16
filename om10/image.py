# ======================================================================

import numpy,pyfits,sys,os,subprocess

import om10

vb = True

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
        make              Generate the images!
        
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
        if survey == None:
           self.survey = 'LSST'
        else:
           self.survey = survey
             
        # Set up observing parameters:   
        if survey=='LSST':
            self.pixscale = 0.2
            self.meanIQ = 0.75
            self.meandepth = 23.3
            self.errdepth = 0.3

        elif survey == 'PS1':
            self.pixscale = 0.25
            self.meanIQ = 1.0
            self.meandepth = 21.4
            self.errdepth = 0.3 

        else:
            raise "ERROR: unrecognised survey "+survey        
        
        self.fov = 10.0 # arcsec
        self.imsize = int(self.fov/self.pixscale)+1

        self.image = numpy.zeros([self.imsize,self.imsize])
        self.lens_galaxy_image = numpy.zeros([self.imsize,self.imsize])

        return
        
    # ------------------------------------------------------------------

    def target(self,lens):

        # Read positions and set up WCS:
        self.RA = 0.0
        self.DEC = 36.0
        # BUG: these should be read from lens!
        
        # Set up WCS:
        # self.set_WCS() 

        return


    # ------------------------------------------------------------------

    def make(self,filters=['i'],Nepochs=1):

        # Make arrays of noise rms and PSF FWHM:
        # self.sample_observing_conditions()        
        
        # Realize the lens galaxy image (as we'll need this many times):
        # self.make_lens_galaxy_image()

        # Loop over epochs, making images:
        for k in range(Nepochs):

            # Make PSF image:
            # self.make_psf_image(k)
        
            # Convolve lens galaxy image with PSF:
            # self.convolve_lens_galaxy_image()

            # Paint in point source PSFs:
            # self.add_quasar_images()

            # Add noise:
            # self.add_noise()

            # Write out image and weight map to file:
            self.write()


        return

    # ------------------------------------------------------------------

    def write(self):

        # Set filenames:
        scifile = 'test_sci.fits'
        varfile = 'test_var.fits'

        # Start an HDU object:
        hdu = pyfits.PrimaryHDU(self.image)

        # Add WCS keywords to the header:
        hdu.header.set('CRVAL1',0.0,'Right ascension (J2000 degrees)')
        
        # Write them out
        hdu.writeto(scifile)        

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


