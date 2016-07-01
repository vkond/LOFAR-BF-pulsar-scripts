#!/usr/bin/env python
#
import numpy as np
import os, os.path, stat, glob, sys, getopt, re
import optparse as opt
import psrchive as pc
import warnings
warnings.simplefilter('ignore', np.RankWarning)
warnings.simplefilter('ignore', RuntimeWarning)

# Main
if __name__=="__main__":
        #
        # Parsing the command line options
        #
        usage = "Usage: %prog [options] .ar"
        cmdline = opt.OptionParser(usage)
        cmdline.add_option('-b', '--bscrunch', dest='bscr', metavar='FACTOR', help="Bscrunch factor, \
default: %default", default=1, type='int')
        cmdline.add_option('-r', '--rotate', dest='rot_bins', metavar='#BIN|PHASE', help="Rotate profile by this number of bins (if absolute value >= 1). \
If absolute value is < 1, then the value is treated as pulse phase in turns. Negative values are to move right, default: %default", default=0, type='float')
        cmdline.add_option('--off-left', dest='off_left', metavar='BIN#', help="Left edge of the off-pulse window to calculate \
mean/rms in the 'Off' mode (inclusive). Default: %default", default=0, type='int')
        cmdline.add_option('--off-right', dest='off_right', metavar='BIN#', help="Right edge of the off-pulse window to calculate \
mean/rms in the 'Off' mode (exclusive). Default: 10%-bin of the profile", default=-1, type='int')
        cmdline.add_option('--fit-left', dest='fit_left', metavar='BIN#', help="Left edge of the window used for a fit. \
Default: %default", default=0, type='int')
        cmdline.add_option('--fit-right', dest='fit_right', metavar='BIN#', help="Right edge of the window used for a fit. \
Default: last bin of the profile", default=-1, type='int')
	cmdline.add_option('--saveonly', dest='is_saveonly', action="store_true", help="Save diagnostic plot \
to png-file instead of GUI", default=False)

        # reading cmd options
        (opts,args) = cmdline.parse_args()

        # check if input file is given
        if len(args) == 0:
                cmdline.print_usage()
                sys.exit(0)

	infile = args[0]
	raw = pc.Archive_load(infile)
       	if not(raw.get_dedispersed()):
		raw.dedisperse()
	# raw.remove_baseline() - we subtract it outselves
        raw.pscrunch()
        nchan = raw.get_nchan()
        nsubint = raw.get_nsubint()
	target = raw.get_source()
        if nchan > 1: raw.fscrunch()
        if nsubint > 1: raw.tscrunch()
        if opts.bscr > 1: raw.bscrunch(opts.bscr)
        nbins = raw.get_nbin()
	if opts.rot_bins != 0:
		if abs(opts.rot_bins) < 1:
			raw.rotate_phase(opts.rot_bins)
		else:
			raw.rotate_phase(opts.rot_bins/nbins)

       	r = raw.get_data()
       	#time stokes f phase
	data = r[0,0,0,:]
        weights = raw.get_weights()
	data[(weights[0]==0)] = 0.0

	if opts.off_right == -1:
		opts.off_right = int(nbins*0.1)
	if opts.off_right-opts.off_left<=1:
		opts.off_right = nbins

	if opts.fit_right == -1:
		opts.fit_right = nbins
	if opts.fit_right-opts.fit_left<=1:
		opts.fit_left = 0

	mean = np.mean(data[opts.off_left:opts.off_right])
	rms = np.std(data[opts.off_left:opts.off_right])
	prof = (data - mean)/rms

	# Fitting
	fit_range = range(opts.fit_left, opts.fit_right)
	ylog_fit=np.log(prof[opts.fit_left:opts.fit_right])
	crit=np.isfinite(ylog_fit)
	ylog_fit[-crit] = 0.0
	crit=np.isnan(ylog_fit)
	ylog_fit[crit] = 0.0
	polynom_coeffs = np.polyfit(fit_range, ylog_fit, 1)
        polynom_bline = np.polyval(polynom_coeffs, fit_range)
	tau = -1./polynom_coeffs[0]
	ampl = np.exp(polynom_coeffs[1])
	print "Coeffs: %s" % (", ".join(["%f" % ii for ii in polynom_coeffs]))
	print "A = %f  tau = %f" % (ampl, tau)
	scat_tail = [np.exp(polynom_coeffs[0]*ii + polynom_coeffs[1]) for ii in fit_range]

	# Plotting
	if opts.is_saveonly:
		pngname = ".".join(infile.split(".")[:-1]) + ".png"
	       	import matplotlib
		matplotlib.use("Agg")

	import matplotlib.pyplot as plt
	import matplotlib.cm as cm
        from matplotlib.ticker import *

	x=range(len(prof))

	fig = plt.figure(figsize=(8,6))
	ax = fig.add_subplot(211)
	plt.xlabel("Bin")
	plt.ylabel("S/N")
	ax.plot(x, prof, "g-", alpha=0.5)
	ax.plot(fit_range, scat_tail, "b--", alpha=0.7, lw=2, label=r'$\tau = %.2f$ bins' % tau)
	ax.axvspan(opts.off_left, opts.off_right, color="yellow", alpha=0.2)
	ax.axvline(x=opts.off_left, linestyle=":", color="black")
	ax.axvline(x=opts.off_right, linestyle=":", color="black")
	ax.axvline(x=np.argmax(prof), linestyle="--", color="black")
	ax.axvline(x=opts.fit_left, linestyle=":", color="black")
	ax.axvline(x=opts.fit_right, linestyle=":", color="black")
	ax.set_xlim(xmin=0, xmax=len(prof)-1)
	ax.set_ylim(ymin=np.min(prof), ymax=np.max(prof))
	ax.axhline(y=0.0, linestyle="--", color="black")
	ax.legend(loc=0)

	ax2 = fig.add_subplot(212)
	plt.xlabel("Bin")
	plt.ylabel("Log(S/N)")
	y=np.log(prof)
	crit=np.isfinite(y)
	y[-crit] = 0.0
	crit=np.isnan(y)
	y[crit] = 0.0
	ax2.plot(x, y, "g-", alpha=0.5)
	ax2.plot(fit_range, polynom_bline, "b--", alpha=0.7, lw=2, label=r'$\tau = %.2f$ bins' % tau)
	ax2.axvline(x=np.argmax(y), linestyle="--", color="black")
	ax2.axvline(x=opts.fit_left, linestyle=":", color="black")
	ax2.axvline(x=opts.fit_right, linestyle=":", color="black")
	ax2.set_xlim(xmin=np.min(x), xmax=np.max(x))
	ax2.set_ylim(ymin=np.min(y), ymax=np.max(y))
	ax2.axhline(y=0.0, linestyle="--", color="black")
	ax2.legend(loc=0)

	if opts.is_saveonly:
		plt.savefig(pngname)
	else:
		plt.show()
