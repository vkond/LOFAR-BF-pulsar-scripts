#
# Defines functions to calculate the effective area of
# the HBA station at a given freq, EL/ZA and time
# using the mscorpol package functions by Tobia Carozzi
# They were originally part of the lofar_gain.py file
# and were arranged in a separate file on 24.04.2015
#
# Vlad Kondratiev (c) - 24.04.2015
#
# 03.04.2015 - Vlad Kondratiev
#              added new model to calculate Aeff, called 'arisN'
#              It's using maximum theoretical value +
#              scaling with EL as sin(EL)^1.39 as in Noutsos et al. (2015)
#              The original model called 'arts' is default
# 06.04.2015 - Vlad Kondratiev
#              added new model to calculate Aeff, called 'hamaker_carozzi'
#	       It's using also maximum theoretical valueof Aeff in zenith
#	       and correcting it by beam correction coefficient calculated
#	       from Jones matrix based on Hamaker model and 'antennaJones.py'
#	       script written by Tobia Carozzi.
#	       This model only implemented for "..range" function to be called
#              from lofar_psrflux.py or lofar_fluxcal.py
# 16.04.2015 - Vlad Kondratiev
#              speed up hamaker_carozzi model by about 30 times by
#              calling proper functions from mscorpol directly instead of
#              running python script 'antennaJones.py'
# 24.04.2015 - Vlad Kondratiev
#              hamaker_carozzi functions are arranged in a separate file
# 05.08.2016 - Vlad Kondratiev
#              now casa_beamcorr is a dictionary for CS and INTL stations
#              and it gets imported from casa_beamcorr_pkg.py
"""
	Defines functions to calculate the effective area of
	the HBA station at a given freq, EL/ZA and time
	using the mscorpol package functions by Tobia Carozzi
	They were originally part of the lofar_gain.py file
	and were arranged in a separate file on 24.04.2015

	Vlad Kondratiev (c) - 24.04.2015

        03.04.2015 - Vlad Kondratiev
                     added new model to calculate Aeff, called 'arisN'
                     It's using maximum theoretical value +
                     scaling with EL as sin(EL)^1.39 as in Noutsos et al. (2015)
                     The original model called 'arts' is default
	06.04.2015 - Vlad Kondratiev
		     added new model to calculate Aeff, called 'hamaker_carozzi'
		     It's using also maximum theoretical valueof Aeff in zenith
		     and correcting it by beam correction coefficient calculated
		     from Jones matrix based on Hamaker model and 'antennaJones.py'
		     script written by Tobia Carozzi.
		     This model only implemented for "..range" function to be called
		     from lofar_psrflux.py or lofar_fluxcal.py
        16.04.2015 - Vlad Kondratiev
                     speed up hamaker_carozzi model by about 30 times by
                     calling proper functions from mscorpol directly instead of
                     running python script 'antennaJones.py'
	24.04.2015 - Vlad Kondratiev
		     hamaker_carozzi functions are arranged in a separate file
        05.08.2016 - Vlad Kondratiev
                     now casa_beamcorr is a dictionary for CS and INTL stations
                     and it gets imported from casa_beamcorr_pkg.py
"""

import os, sys, math
import numpy as np
from pyrap.measures import measures
from pyrap.quanta import quantity
from antennaJones import getJonesByAntFld

#
# Table of beam correction coefficients for CasA observation for normalization of pulsar obs
# This reference CasA observation was taken by Wijnholds (2011)
# at MJD=55159.77650462962963  that corresponds to 2009-11-24 18:38:10
# This table was generated with casa-print.py for center freq of a subband
# Jn=getJonesByAntFld("Hamaker", quantity([55159.77650462962963, 55159.77650462962963], 'd'), station, \
# 	measures().direction('J2000', str(6.123487681)+'rad', str(1.0265154)+'rad'), freq*1.e6)
# or the following command for the antennaJones.py script:
# antennaJones.py -- Hamaker CS001 "2009-11-24 18:38:10" 2.0 1.0 6.123487681 1.0265154 <freq in Hz>
#
from casa_beamcorr_pkg import *

# Calculate the Aeff theoretical maximum and scaling from Noutsos et al. (2015)
# Aeff here is in m^2 for a 48-tile HBA station or 48-dipole LBA station
# For LBA, the Aeff should be correct for LBA_OUTER, but for the LBA_INNER
# one must take into account the distance between nearest dipoles, as the 
# correct expression for LBA dipole is:
# aeff = min{ lambda^2/3 ; pi * d^2/4 }, where d - is the distance to the nearest
# dipole within the full array. Here I am using only first term, which is incorrect
# for LBA_INNER. Also, for LBA_INNER there are 46 dipoles in the station - not 48 
# freq is in MHz, el is in degrees
def get_lofar_aeff_max(freq, el):
	"""
	Calculate the Aeff using given frequency and EL
	"""
	wavelen = 300.0 / freq
	# HBA
	if freq >= 100.:
		aeff = 48. * 16. * np.minimum((wavelen * wavelen)/3., 1.5625)
	# LBA (LBA_OUTER)
	else:
		aeff = 48. * (wavelen * wavelen)/3.
	return aeff

# Function to return beam correction factor using Jones matrix
# returned by antennaJones.py script written by Tobia Carozzi
# freq is in MHz, direction is the value returned by pyrap.measures().direction
# based on pulsar RA and DEC
def get_beam_correction_factor_hamaker_carozzi(mjd, direction, freq, station="CS002"):
	# first we calculate jones matrix for our pulsar observation
	Jn=getJonesByAntFld("Hamaker", quantity([mjd, mjd], 'd'), station, direction, freq*1.e6)
	# calculate beam correction parameter
	bc_psr = 1./np.abs(0.5*(Jn[0,0,0]*np.conj(Jn[0,0,0]) + Jn[0,0,1]*np.conj(Jn[0,0,1]) + \
			Jn[0,1,0]*np.conj(Jn[0,1,0]) + Jn[0,1,1]*np.conj(Jn[0,1,1])))

	# get proper beam correction factor for CasA observation
	if freq < 10. or freq > 300.:
		print "Frequency %f MHz is out of range for hamaker_carozzi model!" % (freq)
		print "It should be between 10 and 300 MHz"
		sys.exit(1)

	# sorting by the difference between our freq and freq in the table
	facts=sorted([f for f in casa_beamcorr[station]], key=lambda arr: np.abs(freq-arr[0]))
	beamcorr = bc_psr / facts[0][1]
	return beamcorr

# Function to calculate the average Aeff for a frequency range f1-f2 and EL, given MJD, RA and DEC.
# Return value is either single average Aeff value.
# If frequency array 'freqs' is given, then average Aeff is calculated for each frequency range f0-f1, f1-f2, f2-f3,....
# in the array and returned value is the list of average Aeffs.
# The size of the returned array is smaller by 1 than the size of the input 'freqs' array
# f1 and f2 are in MHz, freqs in MHz, direction is the value returned by pyrap.measures().direction
# based on pulsar RA and DEC
def lofar_gain_range_hamaker_carozzi(f1, f2, el, mjd, direction, freqs=None):

	nparts=1 # 100
	# Calculate single average value of Aeff
	if freqs == None:
		tot=0
		for ii in xrange(nparts+1):
			freq = f1 + ii*(f2-f1)/nparts
			aeff = get_lofar_aeff_max(freq, el)  # get first the maximum theoretical Aeff value
			# Here we calculate the beam correction factor to correct Aeff
			beamcorr = get_beam_correction_factor_hamaker_carozzi(mjd, direction, freq)
			aeff /= beamcorr
			tot += aeff
		tot /= (nparts+1)
		return tot
	else:
		aeffs=[]
		for ff in xrange(1, len(freqs)):
			tot=0
			for ii in xrange(nparts+1):
				freq = freqs[ff-1] + ii*(freqs[ff]-freqs[ff-1])/nparts
				aeff = get_lofar_aeff_max(freq, el)  # get first the maximum theoretical Aeff value
				# Here we calculate the beam correction factor to correct Aeff
				beamcorr = get_beam_correction_factor_hamaker_carozzi(mjd, direction, freq)
				aeff /= beamcorr
				tot += aeff
			tot /= (nparts+1)
			aeffs.append(tot)
		return aeffs

# Function to calculate the average Aeff for a frequency range f1-f2 and EL, MJD, RA and DEC.
# Return value is either single average Aeff value.
# If frequency array 'freqs' is given, then average Aeff is calculated for each frequency range f0-f1, f1-f2, f2-f3,....
# in the array and returned value is the list of average Aeffs.
# The size of the returned array is smaller by 1 than the size of the input 'freqs' array
# f1 and f2 are in MHz, freqs in MHz, direction is the value returned by pyrap.measures().direction
# based on pulsar RA and DEC
def lofar_gain_range_hamaker_carozzi_center(f1, f2, el, mjd, direction, freqs=None):

	# Calculate single average value of Aeff
	if freqs == None:
		freq = (f1 + f2)/2.
		aeff = get_lofar_aeff_max(freq, el)  # get first the maximum theoretical Aeff value
		# Here we calculate the beam correction factor to correct Aeff
		beamcorr = get_beam_correction_factor_hamaker_carozzi(mjd, direction, freq)
		aeff /= beamcorr
		return aeff
	else:
		aeffs=[]
		for ff in xrange(1, len(freqs)):
			freq = (freqs[ff-1] + freqs[ff])/2.
			aeff = get_lofar_aeff_max(freq, el)  # get first the maximum theoretical Aeff value
			# Here we calculate the beam correction factor to correct Aeff
			beamcorr = get_beam_correction_factor_hamaker_carozzi(mjd, direction, freq)
			aeff /= beamcorr
			aeffs.append(aeff)
		return aeffs

# Function to calculate the average Aeff for a frequency range f1-f2 and EL, given MJD, RA and DEC.
# Return value is either single average Aeff value.
# If frequency array 'freqs' is given, then average Aeff is calculated for each frequency range f0-f1, f1-f2, f2-f3,....
# in the array and returned value is the list of average Aeffs.
# The size of the returned array is smaller by 1 than the size of the input 'freqs' array
# f1 and f2 are in MHz, freqs in MHz, direction is the value returned by pyrap.measures().direction
# based on pulsar RA and DEC
def lofar_gain_range_hamaker_carozzi_station(f1, f2, el, mjd, direction, station="CS002", freqs=None):

	nparts=1 # 100
	# Calculate single average value of Aeff
	if freqs == None:
		tot=0
		for ii in xrange(nparts+1):
			freq = f1 + ii*(f2-f1)/nparts
			aeff = get_lofar_aeff_max(freq, el)  # get first the maximum theoretical Aeff value
			# Here we calculate the beam correction factor to correct Aeff
			beamcorr = get_beam_correction_factor_hamaker_carozzi(mjd, direction, freq, station)
			aeff /= beamcorr
			tot += aeff
		tot /= (nparts+1)
		return tot
	else:
		aeffs=[]
		for ff in xrange(1, len(freqs)):
			tot=0
			for ii in xrange(nparts+1):
				freq = freqs[ff-1] + ii*(freqs[ff]-freqs[ff-1])/nparts
				aeff = get_lofar_aeff_max(freq, el)  # get first the maximum theoretical Aeff value
				# Here we calculate the beam correction factor to correct Aeff
				beamcorr = get_beam_correction_factor_hamaker_carozzi(mjd, direction, freq, station)
				aeff /= beamcorr
				tot += aeff
			tot /= (nparts+1)
			aeffs.append(tot)
		return aeffs

# Function to calculate the average Aeff for a frequency range f1-f2 and EL, given MJD, RA and DEC.
# Return value is either single average Aeff value.
# If frequency array 'freqs' is given, then average Aeff is calculated for each frequency range f0-f1, f1-f2, f2-f3,....
# in the array and returned value is the list of average Aeffs.
# The size of the returned array is smaller by 1 than the size of the input 'freqs' array
# f1 and f2 are in MHz, freqs in MHz, direction is the value returned by pyrap.measures().direction
# based on pulsar RA and DEC
def lofar_gain_range_hamaker_carozzi_center_station(f1, f2, el, mjd, direction, station="CS002", freqs=None):

	# Calculate single average value of Aeff
	if freqs == None:
		freq = (f1 + f2)/2.
		aeff = get_lofar_aeff_max(freq, el)  # get first the maximum theoretical Aeff value
		# Here we calculate the beam correction factor to correct Aeff
		beamcorr = get_beam_correction_factor_hamaker_carozzi(mjd, direction, freq, station)
		aeff /= beamcorr
		return aeff
	else:
		aeffs=[]
		for ff in xrange(1, len(freqs)):
			freq = (freqs[ff-1] + freqs[ff])/2.
			aeff = get_lofar_aeff_max(freq, el)  # get first the maximum theoretical Aeff value
			# Here we calculate the beam correction factor to correct Aeff
			beamcorr = get_beam_correction_factor_hamaker_carozzi(mjd, direction, freq, station)
			aeff /= beamcorr
			aeffs.append(aeff)
		return aeffs
