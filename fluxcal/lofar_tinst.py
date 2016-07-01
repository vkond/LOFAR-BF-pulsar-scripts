#!/usr/bin/env python
#
# Calculates the LOFAR HBA and LBA instrument temperature
# in the HBA frequency range 110-250 MHz and 
# in the LBA frequency range 10-90 MHz.
#
# Vlad Kondratiev (c) - 02.11.2014
#
"""
	Calculates the LOFAR HBA instrument temperature
	in the HBA frequency range 110-250 MHz.
	and LBA instrument temperature
	in the LBA frequency range 10-90 MHz.
	
	Vlad Kondratiev (c) - 02.11.2014 
"""

import os, sys
import numpy as np
import optparse as opt

# calculates the LOFAR HBA/LBA Tinst using polynomial expressions 
# for Tinst from fit to Wijnholds (2011).
# Argument freq is in MHz. Return value is Tinst in Kelvins.
# If given frequency is 0, then return value is a tuple with two lists:
# one is frequencies, and another is Tinsts
# If frequency array 'freqs' is given, then Tinst will be calculated for each
# frequency of the array and returned value is list of Tinst's
def lofar_tinst(freq, freqs=None):
	"""
	calculates the LOFAR HBA/LBA Tinst using polynomial expressions 
	for Tinst from fit to Wijnholds (2011).
	Argument freq is in MHz. Return value is Tinst in Kelvins.
	If given frequency is 0, then return value is a tuple with two lists:
	one is frequencies, and another is Tinsts
	If frequency array 'freqs' is given, then Tinst will be calculated for each
	frequency of the array and returned value is list of Tinst's
	"""

	if freq < 100:
		flow=10
		fhigh=90
		fstep=5
		# polynomial coefficients
		T_inst_poly = [6.2699888333e-05, -0.019932340239, 2.60625093843, -179.560314268, 6890.14953844, -140196.209123, 1189842.07708]
		# if freq is 0, then we calculate Tinst for all freqs in the Table and return tuple of lists
		if freq == 0:
			frequencies=range(flow, fhigh+fstep, fstep)
			temps=[]
			dpoly = len(T_inst_poly)
			for ff in frequencies:
				tinst = 0.0
				for ii in xrange(dpoly): tinst += T_inst_poly[ii]*(ff)**(dpoly-ii-1)
				temps.append(tinst)
	if freq >= 100 or freq == 0:
		flow=110
		fhigh=250
		fstep=5
		# polynomial coefficients
	        T_inst_poly = [6.64031379234e-08, -6.27815750717e-05, 0.0246844426766, -5.16281033712, 605.474082663, -37730.3913315, 975867.990312]
		# if freq is 0, then we calculate Tinst for all freqs in the Table and return tuple of lists
		if freq == 0:
			frequencies.extend(range(flow, fhigh+fstep, fstep))
			dpoly = len(T_inst_poly)
			for ff in xrange(flow, fhigh+fstep, fstep):
				tinst = 0.0
				for ii in xrange(dpoly): tinst += T_inst_poly[ii]*(ff)**(dpoly-ii-1)
				temps.append(tinst)

	tinst = 0.0
	dpoly = len(T_inst_poly)

	# if freq is 0, then we calculate Tinst for all freqs in the Table and
	# return tuple of lists
	if freq == 0:
		return (frequencies, temps)
	else:
		if freqs == None:
			for ii in xrange(dpoly): tinst += T_inst_poly[ii]*(freq)**(dpoly-ii-1)
			return tinst
		else:
			tinsts=[]
			for freq in freqs:
				tinst = 0.0
				for ii in xrange(dpoly): tinst += T_inst_poly[ii]*(freq)**(dpoly-ii-1)
				tinsts.append(tinst)
			return tinsts

# calculates the LOFAR HBA/LBA average Tinst using polynomial expressions 
# for Tinst from fit to Wijnholds (2011) between frequencies f1 and f2 (in MHz).
# Return value is Tinst in Kelvins.
# If frequency array 'freqs' is given, then average Tinst will be calculated for each
# frequency range f0-f1, f1-f2, f2-f2 of the array and returned value is list of average Tinst's.
# Size of the returned array is smaller by 1 than the size of the input freqs array
# Each pair of frequencies should be either above 100 MHz or below 100 MHz
def lofar_tinst_range(f1, f2, freqs=None):
	"""
	calculates the LOFAR HBA/LBA average Tinst using polynomial expressions 
	for Tinst from fit to Wijnholds (2011) between frequencies f1 and f2 (in MHz).
	Return value is Tinst in Kelvins.
	If frequency array 'freqs' is given, then average Tinst will be calculated for each
	frequency range f0-f1, f1-f2, f2-f2 of the array and returned value is list of average Tinst's.
	Size of the returned array is smaller by 1 than the size of the input freqs array
	Each pair of frequencies should be either above 100 MHz or below 100 MHz
	"""

	if (f1 < 100 and f2 >= 100) or (f1 >= 100 and f2 < 100):
		print "Wrong pair of frequencies: (%g, %g)! Both should be either above 100 MHz of below 100 MHz." % (f1, f2)
		sys.exit(1)

	if f1 < 100 and f2 < 100:
		flow=10
		fhigh=90
		fstep=5
		# polynomial coefficients
		T_inst_poly = [6.2699888333e-05, -0.019932340239, 2.60625093843, -179.560314268, 6890.14953844, -140196.209123, 1189842.07708]

	if f1  >= 100 and f2 >= 100:
		flow=110
		fhigh=250
		fstep=5
		# polynomial coefficients
	        T_inst_poly = [6.64031379234e-08, -6.27815750717e-05, 0.0246844426766, -5.16281033712, 605.474082663, -37730.3913315, 975867.990312]
	dpoly = len(T_inst_poly)

	if freqs == None:
		tot=0
		for ff in xrange(101):
			freq = f1 + ff*(f2-f1)/100.
			tinst = 0.0
			for ii in xrange(dpoly): tinst += T_inst_poly[ii]*(freq)**(dpoly-ii-1)
			tot += tinst
		tot /= 100.
		return tot
	else:
		tinsts=[]
		for ff in xrange(1, len(freqs)):
			if (freqs[ff-1] < 100 and freqs[ff] >= 100) or (freqs[ff-1] >= 100 and freqs[ff] < 100):
				print "Wrong pair of frequencies: (%g, %g)! Both should be either above 100 MHz of below 100 MHz." % (freqs[ff-1], freqs[ff])
				tinsts.append(-1) # adding -1 to the array
				continue
			tot = 0
			for ii in xrange(101):
				freq = freqs[ff-1] + ii*(freqs[ff]-freqs[ff-1])/100.
				tinst = 0.0
				for jj in xrange(dpoly): tinst += T_inst_poly[jj]*(freq)**(dpoly-jj-1)
				tot += tinst
			tot /= 100.
			tinsts.append(tot)
		return tinsts
	
#   m a i n
if __name__=="__main__":
        #
        # Parsing the command line options
        #
        usage = "Usage: %prog"
        cmdline = opt.OptionParser(usage)
	cmdline.add_option('-f', '--freq', dest='freq', metavar='MHz', help="HBA or LBA frequency (in MHz) \
for which to calculate instrument temperature (Tinst) in K. If frequency value is not given, then \
the temperatures will be calculated for the whole LBA range 10-90 MHz and HBA range 110-250 MHz with 5-MHz steps \
and printed out. With --plot option you still get the plot of Tinst(f) dependence.", default=0, type='float')
	cmdline.add_option('--plot', dest='is_plot', action="store_true", help="To make Tinst-vs-frequency \
plot. If --freq is used, the given frequency value will be also marked", default=False)

        # reading cmd options
        (opts,args) = cmdline.parse_args()

	# getting the list
	if opts.is_plot or opts.freq == 0:
		(freqs, temps) = lofar_tinst(0)

	# getting Tinst for specific frequency
	if opts.freq != 0:
		tinst = lofar_tinst(opts.freq)
		print "Tinst at %g MHz = %g K" % (opts.freq, tinst)

	if opts.freq == 0:
		print "#"
		print "# Freq (MHz)\tTinst (K)"
		print "#-------------------------"
		for (f,t) in zip(freqs, temps):
			print "%g\t\t%g" % (f, t)
		print

	# making the plot
	if opts.is_plot:
		import matplotlib.pyplot as plt	
		fig = plt.figure()
		ax = fig.add_subplot(111)
		plt.xlabel("Frequency (MHz)")
		plt.ylabel("Tinst (K)")
		ax.plot(freqs, temps, "go", alpha=0.7, linewidth=2)
		ax.plot(freqs, temps, "-", color="orange", alpha=0.5, linewidth=1)
#		(xmin, xmax, ymin, ymax) = ax.axis()
		if opts.freq != 0 and sign == "":
#			ax.axvline(x=opts.freq, ymax=(tinst-ymin)/(ymax-ymin), linestyle="--", color="black", alpha=0.5)
#			ax.axhline(y=tinst, xmax=(opts.freq-xmin)/(xmax-xmin), linestyle="--", color="black", alpha=0.5)
			ax.plot([opts.freq], [tinst], "yo")

		plt.grid()
		plt.show()
