#!/usr/bin/env python
#
# Calculates the LOFAR HBA effective area for a single 48-tile station
# at the given frequency in the HBA frequency range 110-250 MHz. 
# The step in frequency is 5 MHz from 110 MHz up.
# Based on Stefan Wijnholds' model. Updated from original
# script from Joeri.
#
# Vlad Kondratiev (c) - 08.11.2014
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
#              moved hamaker_carozzi functions to another file
#              lofar_gain_hamaker_carozzi.py
"""
	Calculates the LOFAR HBA effective area for a single 48-tile station
	at the given frequency in the HBA frequency range 110-250 MHz. 
	The step in frequency is 5 MHz from 110 MHz up.
	Based on Stefan Wijnholds' model. Updated from original
	script from Joeri.

	Vlad Kondratiev (c) - 08.11.2014

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
		     moved hamaker_carozzi functions to another file
		     lofar_gain_hamaker_carozzi.py
"""

import os, sys, math
import optparse as opt
import numpy as np

# gain file
try:
        # Lofar software
        lofarsoft=os.environ['LOFARSOFT']
        # script that gets the status of processed data and also returns the size of processing directory
        gainfile="%s/release/share/pulsar/data/lofar_sensitivity_senstable-v02.txt" % (lofarsoft,)
except:
	gainfile="lofar_sensitivity_senstable-v02.txt" # current directory



# Calculate the Aeff using Gain's table for a given frequency and EL
def get_lofar_gain_arts(gtable, freq, el, flow=110, fhigh=250, fstep=5):
	"""
	Calculate the Aeff using Gain's table for a given frequency and EL
	"""
        # if frequency is equal exactly a value in the table
       	if freq in xrange(flow, fhigh+fstep, fstep):
		# if elevation is equal exactly a value in the table
		if el in xrange(0, gtable.shape[1], 1):
			aeff=gtable[int(((freq-flow)/fstep)+0.5)][el]
		else:
       			# if elevation is between table's values, then we need to interpolate (linearly)
		        for ii, val in enumerate(xrange(0, gtable.shape[1], 1)):
       			        if el < val:
					rel=val
                       			break
			aeff=np.interp(el, [rel-1, rel], [gtable[int(((freq-flow)/fstep)+0.5)][rel-1], gtable[int(((freq-flow)/fstep)+0.5)][rel]])
	else:
		# if elevation is equal exactly a value in the table
		if el in xrange(0, gtable.shape[1], 1):
       			# if frequency is between table's values, then we need to interpolate (linearly)
		        for ii, val in enumerate(xrange(flow, fhigh+fstep, fstep)):
       			        if freq < val:
               			        rii=ii
					rfreq=val
                       			break
			aeff=np.interp(freq, [rfreq-fstep, rfreq], [gtable[rii-1][el], gtable[rii][el]])
		else:
			import scipy.interpolate
		        for ii, val in enumerate(xrange(flow, fhigh+fstep, fstep)):
       			        if freq < val:
               			        rii=ii
					rfreq=val
                       			break
		        for ii, val in enumerate(xrange(0, gtable.shape[1], 1)):
       			        if el < val:
					rel=val
                       			break
			laeff=np.interp(freq, [rfreq-fstep, rfreq], [gtable[rii-1][rel-1], gtable[rii][rel-1]])
			raeff=np.interp(freq, [rfreq-fstep, rfreq], [gtable[rii-1][rel], gtable[rii][rel]])
			aeff=np.interp(el, [rel-1, rel], [laeff, raeff])
			# does not work - old version of Scipy on CEP2
			#aeff=scipy.interpolate.griddata([(rfreq-fstep, rel-1), (rfreq-fstep, rel), (rfreq, rel), (rfreq, rel-1)], \
			#		[gtable[rii-1][rel-1], gtable[rii-1][rel], gtable[rii][rel], gtable[rii][rel-1]], \
			#		[(freq, el)], method='linear')
	return aeff
		

# Main function to calculate the Aeff using Gain's table for a given frequency and EL
# Return value is either single Aeff value when both freq and elevation are given, 0 - when both EL and freq
# are not given, or list of aeffs for a range of freqs or elevations
# If frequency array 'freqs' is given, then Aeff is calculated for each frequency in the array
# and returned value is the list of Aeffs
def lofar_gain_arts(freq, el, is_quiet=False, freqs=None, is_plot=False, gfile=gainfile, flow=110, fhigh=250, fstep=5):

	nfreqs=len(xrange(flow, fhigh+fstep, fstep)) # 29
	# checking that given frequency is within the HBA range
	if freq != 0:
		if freq < flow:
			if is_quiet:
				freq = flow
			else:
				print "Given frequency of %g MHz is below %g MHz!" % (freq, flow)
				sys.exit(1)
		if freq > fhigh:
			if is_quiet:
				freq = fhigh
			else:
				print "Given frequency of %g MHz is above %g MHz!" % (freq, fhigh)
				sys.exit(1)

	# checking if given EL value is within the range
	if el != -100.:
		if el < 0 or el > 90.:
			print "Given elevation of %g deg is outside 0-90 deg range!" % (el)
			sys.exit(1)

	# reading gain's file

	# This is the data as provided by Stephan Wijnholds on 20120919 
	# These numbers are A/(T_sky+T_rec)

	# Dit tekstbestand bevat 91 (elevatiepunten) * 361 (azimuthpunten)
	# * 29 (frequentiepunten) = 952679 getallen. Deze zijn dus te
	# transformeren tot een 91 x 361 x 29 datakubus, waarbij de eerste
	# index de snelst lopende index is. De elevatie-as loopt van 0 to
	# 90 graden in stappen van 1 graad, de azimuth-as van 0 tot 360
	# graden in stappen van 1 graad en de frequentie-as van 110 tot
	# 250 MHz in stappen van 5 MHz.

	try:
		gtable=np.loadtxt(gfile) # first check if gainfile is in $LOFARSOFT/... or
					 # if LOFARSOFT was not set, then in the current dir
	except:
		# now we should check if gainfile is in the same directory as executable script
		try:
			gtable=np.loadtxt("%s/%s" % (os.path.dirname(os.path.realpath(sys.argv[0])), gfile.split("/")[-1]))
		except:
			# now we check again in the current dir
			# in case we were checking $LOFARSOFT/ in the first time...
			try:
				gtable=np.loadtxt("%s" % (gfile.split("/")[-1]))
			except:
				print "Can't open gain file '%s'. It's not present neither in \$LOFARSOFT/release/share/pulsar/data"
				print "nor in the same location as the calling script nor in the current directory."
				print "Use --file option to point to the correct location"
				sys.exit(1)

	naz = 361 # from 0 to 360 deg with 1-deg step
	nel = 91  # from 0 to 90 deg with 1-deg step
	gtable.shape=(nfreqs, naz, nel)  # 91 elevations * 361 azimuths * 29 frequencies
	gtable=np.mean(gtable, axis=1) # average over azimuths, as the stations are randomly rotated.

	# "Ik ben vandaag eens even weer in de code gedoken. Wat ik in de
	#  simulatie doe, is het uitrekenen van Aeff (apart dus). Vervolgens
	#  bepaal ik Aeff/Tsys door de gevonden Aeff te delen door de
	#  systeemtemperatuur. Die laatste heb ik als volgt bepaald (quote uit
	#  de code):"
	# Trec = 400;
	# Tsky = (lambda / 0.2008).^2.55 + (freq / 1e9).^1.8 + 2.7;
	# Tsys = Tsky + Trec;

	# correcting gain table to have proper values of Aeff
	# this is needed because in the input file, values Aeff are divided by T
	# using constant Trec and some Tsky
	for ii in xrange(gtable.shape[0]): # all rows
		f=(flow+fstep*ii)*1.0E6
		l=3.0E8/f
		# this is the T value used in the table on file
		t=400.0 + (l/0.2008)**2.55 + (f/1E9)**1.8 + 2.7
		# so, take out that division by t:
		gtable[ii,:]=gtable[ii,:] * t

	if is_plot:
		import matplotlib.pyplot as plt

	# Checking if both frequency and elevation values are given
	# If one of them or both are not given, then we will print out
	# the corresponding table
	if freq == 0 and el == -100.:
		print "#"
		print "# Freq (MHz)\tEL (deg)\tAeff (m^2)"
		print "#------------------------------------------"
		for ff in xrange(nfreqs):
			for ee in xrange(nel):
				print "%d\t\t%d\t\t%g" % (flow+fstep*ff, ee, gtable[ff][ee])

	        # making the plot
        	if is_plot:
	                fig = plt.figure()
        	        ax = fig.add_subplot(111)
                	plt.xlabel("Zenith angle (deg)")
	                plt.ylabel("Aeff (m^2)")
        	        ax.plot(xrange(0, nel, 1), gtable[-1][::-1], "b-", alpha=0.7, linewidth=2, label="%g MHz" % (fhigh))
        	        ax.plot(xrange(0, nel, 1), gtable[0][::-1], "-", color="orange", alpha=0.7, linewidth=2, label="%g MHz" % (flow))
			for ff in xrange(1, nfreqs-1, 1):
        	        	ax.plot(xrange(0, nel, 1), gtable[ff][::-1], "g-", alpha=0.3, linewidth=1)
	                plt.grid()
			plt.legend(loc=0)
        	        plt.show()
		return 0

	if freq == 0 and el != -100.:
		print "#"
		print "# Freq (MHz)\tAeff (m^2)"
		print "#--------------------------"
		aeffs=[]
	        # if elevation is equal exactly a value in the table
        	if el in xrange(0, nel, 1):
			for ff in xrange(nfreqs):
				aeff=gtable[ff][el]
				aeffs.append(aeff)
				print "%d\t\t%g" % (flow+ff*fstep, aeff)
		else:
        		# if elevation is between table's values, then we need to interpolate (linearly)
		        for ii, val in enumerate(xrange(0, nel, 1)):
        		        if el < val:
					rel=val
                        		break
			for ff in xrange(nfreqs):
				aeff=np.interp(el, [rel-1, rel], [gtable[ff][rel-1], gtable[ff][rel]])
				aeffs.append(aeff)
				print "%d\t\t%g" % (flow+ff*fstep, aeff)

	        # making the plot
        	if is_plot:
	                fig = plt.figure()
        	        ax = fig.add_subplot(111)
                	plt.xlabel("Frequency (MHz)")
	                plt.ylabel("Aeff (m^2)")
        	        ax.plot(xrange(flow, fhigh+fstep, fstep), aeffs, "g-", alpha=0.7, linewidth=2)
	                plt.grid()
        	        plt.show()
		return aeffs

	if freq != 0 and el == -100.:
		print "#"
		print "# EL (deg)\tAeff (m^2)"
		print "#--------------------------"
		aeffs=[]
	        # if frequency is equal exactly a value in the table
        	if freq in xrange(flow, fhigh+fstep, fstep):
			for ee in xrange(nel):
				aeff=gtable[int(((freq-flow)/fstep)+0.5)][ee]
				aeffs.append(aeff)
				print "%d\t\t%g" % (ee, aeff)
		else:
        		# if frequency is between table's values, then we need to interpolate (linearly)
		        for ii, val in enumerate(xrange(flow, fhigh+fstep, fstep)):
        		        if freq < val:
                		        rii=ii
					rfreq=val
                        		break
			for ee in xrange(nel):
				aeff=np.interp(freq, [rfreq-fstep, rfreq], [gtable[rii-1][ee], gtable[rii][ee]])
				aeffs.append(aeff)
				print "%d\t\t%g" % (ee, aeff)

	        # making the plot
        	if is_plot:
	                fig = plt.figure()
        	        ax = fig.add_subplot(111)
                	plt.xlabel("Zenith angle (deg)")
	                plt.ylabel("Aeff (m^2)")
        	        ax.plot(xrange(0, nel, 1), aeffs[::-1], "g-", alpha=0.7, linewidth=2)
	                plt.grid()
        	        plt.show()
		return aeffs

	# Calculate (interpolate) single value of Aeff and print it out
	if freq != 0 and el != -100.:
		if freqs == None:
			aeff = get_lofar_gain_arts(gtable, freq, el)
			if not is_quiet:
				print "#"
				print "# Freq (MHz)\tEL (deg)\tAeff (m^2)"
				print "#------------------------------------------"
				print "%g\t\t%g\t\t%g" % (freq, el, aeff)

			        # making the plot
        			if is_plot:
	        	        	fig = plt.figure()
	        	        	ax = fig.add_subplot(111)
		                	plt.xlabel("Zenith angle (deg)")
			                plt.ylabel("Aeff (m$^2$)")
        			        ax.plot(xrange(0, nel, 1), gtable[-1][::-1], "b-", alpha=0.7, linewidth=2, label="%g MHz" % (fhigh))
        	        		ax.plot(xrange(0, nel, 1), gtable[0][::-1], "-", color="orange", alpha=0.7, linewidth=2, label="%g MHz" % (flow))
					for ff in xrange(1, nfreqs-1, 1):
        			        	ax.plot(xrange(0, nel, 1), gtable[ff][::-1], "g-", alpha=0.3, linewidth=1)
					ax.plot([90-el], [aeff], "o", color="black", label="ZA=%g$^{\circ}$\nf=%g MHz\nA$_\mathrm{eff}$=%g m$^2$" % (90-el, freq, aeff))
	        		        plt.grid()
					plt.legend(loc=0)
		        	        plt.show()
			return aeff
		else:
			aeffs=[]
			for freq in freqs:
				if freq < flow: freq = flow
				if freq > fhigh: freq = fhigh
				aeff = get_lofar_gain_arts(gtable, freq, el)
				aeffs.append(aeff)
			return aeffs

# Function to calculate the average Aeff using Gain's table for a frequency range f1-f2 and EL
# Return value is either single average Aeff value.
# If frequency array 'freqs' is given, then average Aeff is calculated for each frequency range f0-f1, f1-f2, f2-f3,....
# in the array and returned value is the list of average Aeffs.
# The size of the returned array is smaller by 1 than the size of the input 'freqs' array
# f1 and f2 are in MHz, el in degrees
def lofar_gain_range_arts(f1, f2, el, freqs=None, gfile=gainfile, flow=110, fhigh=250, fstep=5):

	nfreqs=len(xrange(flow, fhigh+fstep, fstep)) # 29
	# checking that given frequency is within the HBA range
	if f1 < flow: f1 = flow
	if f2 < flow: f2 = flow
	if f1 > fhigh: f1 = fhigh
	if f2 > fhigh: f2 = fhigh

	# checking if given EL value is within the range
	if el < 0 or el > 90.:
		print "Given elevation of %g deg is outside 0-90 deg range!" % (el)
		sys.exit(1)

	# reading gain's file

	# This is the data as provided by Stephan Wijnholds on 20120919 
	# These numbers are A/(T_sky+T_rec)

	# Dit tekstbestand bevat 91 (elevatiepunten) * 361 (azimuthpunten)
	# * 29 (frequentiepunten) = 952679 getallen. Deze zijn dus te
	# transformeren tot een 91 x 361 x 29 datakubus, waarbij de eerste
	# index de snelst lopende index is. De elevatie-as loopt van 0 to
	# 90 graden in stappen van 1 graad, de azimuth-as van 0 tot 360
	# graden in stappen van 1 graad en de frequentie-as van 110 tot
	# 250 MHz in stappen van 5 MHz.

	try:
		gtable=np.loadtxt(gfile) # first check if gainfile is in $LOFARSOFT/... or
					 # if LOFARSOFT was not set, then in the current dir
	except:
		# now we should check if gainfile is in the same directory as executable script
		try:
			gtable=np.loadtxt("%s/%s" % (os.path.dirname(os.path.realpath(sys.argv[0])), gfile.split("/")[-1]))
		except:
			# now we check again in the current dir
			# in case we were checking $LOFARSOFT/ in the first time...
			try:
				gtable=np.loadtxt("%s" % (gfile.split("/")[-1]))
			except:
				print "Can't open gain file '%s'. It's not present neither in \$LOFARSOFT/release/share/pulsar/data"
				print "nor in the same location as the calling script nor in the current directory."
				print "Use --file option to point to the correct location"
				sys.exit(1)

	naz = 361 # from 0 to 360 deg with 1-deg step
	nel = 91  # from 0 to 90 deg with 1-deg step
	gtable.shape=(nfreqs, naz, nel)  # 91 elevations * 361 azimuths * 29 frequencies
	gtable=np.mean(gtable, axis=1) # average over azimuths, as the stations are randomly rotated.

	# "Ik ben vandaag eens even weer in de code gedoken. Wat ik in de
	#  simulatie doe, is het uitrekenen van Aeff (apart dus). Vervolgens
	#  bepaal ik Aeff/Tsys door de gevonden Aeff te delen door de
	#  systeemtemperatuur. Die laatste heb ik als volgt bepaald (quote uit
	#  de code):"
	# Trec = 400;
	# Tsky = (lambda / 0.2008).^2.55 + (freq / 1e9).^1.8 + 2.7;
	# Tsys = Tsky + Trec;

	# correcting gain table to have proper values of Aeff
	# this is needed because in the input file, values Aeff are divided by T
	# using constant Trec and some Tsky
	for ii in xrange(gtable.shape[0]): # all rows
		f=(flow+fstep*ii)*1.0E6
		l=3.0E8/f
		# this is the T value used in the table on file
		t=400.0 + (l/0.2008)**2.55 + (f/1E9)**1.8 + 2.7
		# so, take out that division by t:
		gtable[ii,:]=gtable[ii,:] * t


	nparts=100
	# Calculate (interpolate) single average value of Aeff
	if freqs == None:
		tot=0
		for ii in xrange(nparts+1):
			freq = f1 + ii*(f2-f1)/nparts
			aeff = get_lofar_gain_arts(gtable, freq, el)
			tot += aeff
		tot /= (nparts+1)
		return tot
	else:
		aeffs=[]
		for ff in xrange(1, len(freqs)):
			tot=0
			for ii in xrange(nparts+1):
				freq = freqs[ff-1] + ii*(freqs[ff]-freqs[ff-1])/nparts
				if freq < flow: freq = flow
				if freq > fhigh: freq = fhigh
				aeff = get_lofar_gain_arts(gtable, freq, el)
				tot += aeff
			tot /= (nparts+1)
			aeffs.append(tot)
		return aeffs

#
# Functions to get Aeff using maximum theoretical values for a station in HBA or LBA range
# scaled with EL as sin(EL)^1.39 as in Noutsos et al. (2015) - 'arisN' model

# Calculate the Aeff theoretical maximum and scaling from Noutsos et al. (2015)
# Aeff here is in m^2 for a 48-tile HBA station or 48-dipole LBA station
# For LBA, the Aeff should be correct for LBA_OUTER, but for the LBA_INNER
# one must take into account the distance between nearest dipoles, as the 
# correct expression for LBA dipole is:
# aeff = min{ lambda^2/3 ; pi * d^2/4 }, where d - is the distance to the nearest
# dipole within the full array. Here I am using only first term, which is incorrect
# for LBA_INNER. Also, for LBA_INNER there are 46 dipoles in the station - not 48 
# freq is in MHz, el is in degrees
def get_lofar_gain_arisN(freq, el):
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
	# scaling with elevation
	aeff *= ((math.sin((math.pi*el)/180.))**(1.39))
	return aeff

#
# Main function to calculate the Aeff using given frequency and EL
# Return value is either single Aeff value when both freq and elevation are given, 0 - when both EL and freq
# are not given, or list of aeffs for a range of freqs or elevations
# If frequency array 'freqs' is given, then Aeff is calculated for each frequency in the array
# and returned value is the list of Aeffs
# freq in MHz, elevation in degrees
def lofar_gain_arisN(freq, el, is_quiet=False, freqs=None, is_plot=False):

	# checking  the frequency range
	if freq != 0:
		if freq >= 100:
			flow = 110
			fhigh = 250
			fstep = 5			
		else:
			flow = 10
			fhigh = 90
			fstep = 5
	else:
		flow = 10
		fhigh = 250
		fstep = 5

	# getting this number for plotting purposes only
	nfreqs=len(xrange(flow, fhigh+fstep, fstep))
	# number of elevations (plotting purposes only)
	nel=len(xrange(0, 91, 1))

	# checking if given EL value is within the range
	if el != -100.:
		if el < 0 or el > 90.:
			print "Given elevation of %g deg is outside 0-90 deg range!" % (el)
			sys.exit(1)

	if is_plot:
		import matplotlib.pyplot as plt

	# Checking if both frequency and elevation values are given
	# If one of them or both are not given, then we will print out
	# the corresponding table
	if freq == 0 and el == -100.:
		print "#"
		print "# Freq (MHz)\tEL (deg)\tAeff (m^2)"
		print "#------------------------------------------"
		aeffs=np.zeros((nfreqs, nel), dtype=float)
		for ff in xrange(nfreqs):
			for ee in xrange(nel):
				aeffs[ff][ee] = get_lofar_gain_arisN(flow+fstep*ff, ee)
				print "%d\t\t%d\t\t%g" % (flow+fstep*ff, ee, aeffs[ff][ee])

	        # making the plot
        	if is_plot:
	                fig = plt.figure()
        	        ax = fig.add_subplot(111)
                	plt.xlabel("Zenith angle (deg)")
	                plt.ylabel("Aeff (m^2)")
        	        ax.plot(xrange(0, nel, 1), aeffs[-1][::-1], "b-", alpha=0.7, linewidth=2, label="%g MHz" % (fhigh))
        	        ax.plot(xrange(0, nel, 1), aeffs[0][::-1], "-", color="orange", alpha=0.7, linewidth=2, label="%g MHz" % (flow))
			for ff in xrange(1, nfreqs-1, 1):
        	        	ax.plot(xrange(0, nel, 1), aeffs[ff][::-1], "g-", alpha=0.3, linewidth=1)
	                plt.grid()
			plt.legend(loc=0)
        	        plt.show()
		return 0

	if freq == 0 and el != -100.:
		print "#"
		print "# Freq (MHz)\tAeff (m^2)"
		print "#--------------------------"
		aeffs=[]
		for ff in xrange(nfreqs):
			aeff=get_lofar_gain_arisN(flow+fstep*ff, el)
			aeffs.append(aeff)
			print "%d\t\t%g" % (flow+ff*fstep, aeff)

	        # making the plot
        	if is_plot:
	                fig = plt.figure()
        	        ax = fig.add_subplot(111)
                	plt.xlabel("Frequency (MHz)")
	                plt.ylabel("Aeff (m^2)")
        	        ax.plot(xrange(flow, fhigh+fstep, fstep), aeffs, "g-", alpha=0.7, linewidth=2)
	                plt.grid()
        	        plt.show()
		return aeffs

	if freq != 0 and el == -100.:
		print "#"
		print "# EL (deg)\tAeff (m^2)"
		print "#--------------------------"
		aeffs=[]
		for ee in xrange(nel):
			aeff=get_lofar_gain_arisN(freq, ee)
			aeffs.append(aeff)
			print "%d\t\t%g" % (ee, aeff)

	        # making the plot
        	if is_plot:
	                fig = plt.figure()
        	        ax = fig.add_subplot(111)
                	plt.xlabel("Zenith angle (deg)")
	                plt.ylabel("Aeff (m^2)")
        	        ax.plot(xrange(0, nel, 1), aeffs[::-1], "g-", alpha=0.7, linewidth=2)
	                plt.grid()
        	        plt.show()
		return aeffs

	# Calculate single value of Aeff and print it out
	if freq != 0 and el != -100.:
		if freqs == None:
			aeff = get_lofar_gain_arisN(freq, el)
			if not is_quiet:
				print "#"
				print "# Freq (MHz)\tEL (deg)\tAeff (m^2)"
				print "#------------------------------------------"
				print "%g\t\t%g\t\t%g" % (freq, el, aeff)

			        # making the plot
        			if is_plot:
					aeffs=np.zeros((nfreqs, nel), dtype=float)
					for ff in xrange(nfreqs):
						for ee in xrange(nel):
							aeffs[ff][ee] = get_lofar_gain_arisN(flow+fstep*ff, ee)

	        	        	fig = plt.figure()
	        	        	ax = fig.add_subplot(111)
		                	plt.xlabel("Zenith angle (deg)")
			                plt.ylabel("Aeff (m$^2$)")
        			        ax.plot(xrange(0, nel, 1), aeffs[-1][::-1], "b-", alpha=0.7, linewidth=2, label="%g MHz" % (fhigh))
        	        		ax.plot(xrange(0, nel, 1), aeffs[0][::-1], "-", color="orange", alpha=0.7, linewidth=2, label="%g MHz" % (flow))
					for ff in xrange(1, nfreqs-1, 1):
        			        	ax.plot(xrange(0, nel, 1), aeffs[ff][::-1], "g-", alpha=0.3, linewidth=1)
					ax.plot([90-el], [aeff], "o", color="black", label="ZA=%g$^{\circ}$\nf=%g MHz\nA$_\mathrm{eff}$=%g m$^2$" % (90-el, freq, aeff))
	        		        plt.grid()
					plt.legend(loc=0)
		        	        plt.show()
			return aeff
		else:
			aeffs=[]
			for freq in freqs:
				aeff = get_lofar_gain_arisN(freq, el)
				aeffs.append(aeff)
			return aeffs

# Function to calculate the average Aeff for a frequency range f1-f2 and EL
# Return value is either single average Aeff value.
# If frequency array 'freqs' is given, then average Aeff is calculated for each frequency range f0-f1, f1-f2, f2-f3,....
# in the array and returned value is the list of average Aeffs.
# The size of the returned array is smaller by 1 than the size of the input 'freqs' array
# f1 and f2 are in MHz, el in degrees
def lofar_gain_range_arisN(f1, f2, el, freqs=None):

	nparts=100
	# checking if given EL value is within the range
	if el < 0 or el > 90.:
		print "Given elevation of %g deg is outside 0-90 deg range!" % (el)
		sys.exit(1)

	# Calculate single average value of Aeff
	if freqs == None:
		tot=0
		for ii in xrange(nparts+1):
			freq = f1 + ii*(f2-f1)/nparts
			aeff = get_lofar_gain_arisN(freq, el)
			tot += aeff
		tot /= (nparts+1)
		return tot
	else:
		aeffs=[]
		for ff in xrange(1, len(freqs)):
			tot=0
			for ii in xrange(nparts+1):
				freq = freqs[ff-1] + ii*(freqs[ff]-freqs[ff-1])/nparts
				aeff = get_lofar_gain_arisN(freq, el)
				tot += aeff
			tot /= (nparts+1)
			aeffs.append(tot)
		return aeffs

#   m a i n
if __name__=="__main__":
        #
        # Parsing the command line options
        #
        usage = "Usage: %prog [options]"
        cmdline = opt.OptionParser(usage)
        cmdline.add_option('-f', '--freq', dest='freq', metavar='MHz', help="HBA frequency (in MHz) \
for which to calculate effective area. If frequency falls between table values, then effective area \
will be linearly interpolated. If frequency value is not given, then effective area values \
will be calculated for the whole HBA range and printed out.", default=0, type='float')
	cmdline.add_option('--el', '--elevation', dest='el', metavar='deg', help="Elevation of \
the source. If elevation falls between table values, then effective area will be linearly \
interpolated. If elevation is not given, the effective area values will be calculated for the whole \
elevation range.", default=-100, type='float')
	cmdline.add_option('--model', '--beam-model', dest='model', metavar='STRING', help="Beam model to use. Possible \
values: 'arts' for the beam model by Arts et al. (2013), and 'arisN' that uses theoretical maximum of Aeff in \
zenith and then scales it with elevation as sin(EL)^1.39 as in Noutsos et al. (2015). Default: %default", default="arts", type='str')
	cmdline.add_option('--file', dest='gtable', metavar='STRING', help="Default sensitivity \
file from Stefan Wijnholds having values of A/(T_sky+T_inst) for the HBA frequency range 110-250 MHz, \
wiht 5-MHz step, i.e. 29 frequencies, for 91 elevations and 361 azimuths. By default, file is found \
automatically in the default location: %default. This option is only relevant for the 'arts' model", default=gainfile, type='str')
	cmdline.add_option('--plot', dest='is_plot', action="store_true", help="To plot EL-vs-Aeff plot, or \
Freq-vs-Aeff, or EL-vs-Freq-vs-Aeff plot, depending on given values of frequency/elevation", default=False)
	cmdline.add_option('-q', '--quiet', dest='is_quiet', action="store_true", help="To print out only single Aeff value. Useful \
when run inside the script", default=False)

        # reading cmd options
        (opts,args) = cmdline.parse_args()

	# for the beam model by Arts et al. (2013)
	if opts.model == 'arts':
		lofar_gain_arts(opts.freq, opts.el, opts.is_quiet, None, opts.is_plot, opts.gtable)
	# for the EL-scaling ration used in Noutsos et al. (2015)
	elif opts.model == 'arisN':
		lofar_gain_arisN(opts.freq, opts.el, opts.is_quiet, None, opts.is_plot)
	else:
		print "Can't use model '%s' when run in a standalone mode!" % (opts.model)
		sys.exit(1)
