from lenspop import population_functions, distances
from stellarpop import tools
from astropy.table import Table
import os

# seed
lens_redshift = 0.4
lens_veldisp = 220.0
source_redshift = 1.5
source_veldisp = 100.0

#veldisp = 200.0
#redshift = 1.5

# This function calculates magnitude for r, g, i, z filters
def CalculateMagnitude(dataPath, target):
	#dataFile = os.path.expandvars(dataPath)
	#data = Table.read(dataFile, format='fits')
	# call constructor of the distance class
	d = distances.Distance()
	#redshift = data['ZLENS']
	#veldisp = data['VELDISP']
	if target == 'source':
		# if target is lens, use appropriate SED
		sed = tools.getSED('QSO1_template_norm')
		veldisp = source_veldisp
                redshift = source_redshift
	elif target == 'lens':
		# if target is galaxy, use appropriate SED
		sed = tools.getSED('M82_template_norm')
                veldisp = lens_veldisp
                redshift = lens_redshift
	RF_Rmag_app, offset = CalculateRestFrameRMag(sed, veldisp, redshift, d)
	Gfilter = tools.filterfromfile('g_SDSS')
	Ifilter = tools.filterfromfile('i_SDSS')
	Zfilter = tools.filterfromfile('z_SDSS')
	# maybe explain it?	
	RF_Gmag_app = tools.ABFilterMagnitude(Gfilter, sed, redshift) + offset + d.distance_modulus(redshift)
	RF_Imag_app = tools.ABFilterMagnitude(Ifilter, sed, redshift) + offset + d.distance_modulus(redshift)
	RF_Zmag_app = tools.ABFilterMagnitude(Zfilter, sed, redshift) + offset + d.distance_modulus(redshift)
	return RF_Rmag_app, RF_Gmag_app, RF_Imag_app, RF_Zmag_app

# Decomposition to increase readability but the function is not readable at all I NEED TO COMMENT
def CalculateRestFrameRMag(sed, veldisp, redshift, d):
	# call constructor. Name should be changed
	lenspop_const = population_functions.LensPopulation_()
	# Reference Frame Absolute R magnitude
	RF_RMag_abs, _ = lenspop_const.EarlyTypeRelations(veldisp)
	Rfilter = tools.filterfromfile('r_SDSS')
	RMag_abs = tools.ABFilterMagnitude(Rfilter, sed, redshift)
	Rmag_app = RMag_abs + d.distance_modulus(redshift)
	offset_abs_app = RMag_abs - Rmag_app
	offset_RF_abs = RF_RMag_abs - RMag_abs
	RF_Rmag_app = RF_RMag_abs - offset_abs_app
	return RF_Rmag_app, offset_RF_abs

