#!/usr/bin/env python
#
# Determine the S/N of the profile using diffent methods
# and print out noise/signal statistics for comparison
# Can be used with both Presto .bestprof files and PSRCHIVE
# archive files.
#
# Vlad Kondratiev (c) - 26.11.2014 
#
# 13.07.2016 - Vlad Kondratiev
#	       added --auto-off option to automatically 
#              determine the OFF-pulse window 
#
import numpy as np
import os, os.path, stat, glob, sys, getopt, re
import scipy.stats as sc
import optparse as opt
import warnings
warnings.simplefilter('ignore', np.RankWarning)
warnings.simplefilter('ignore', RuntimeWarning)

# rotate profile by rot_bins
# negative value of rot_bins is to move right
def bestprof_rotate (x, rot_bins):
	if rot_bins == 0: return x
	if abs(rot_bins) < 1.0: # it means that rot_bins is in pulse phase
		rot_bins = (rot_bins/abs(rot_bins))*int(abs(rot_bins)*len(x)+0.5)
	if abs(rot_bins) >= len(x):
		rot_bins = (rot_bins/abs(rot_bins))*(int(abs(rot_bins)+0.5)%len(x))
	if rot_bins > 0:
		out=np.append(x[rot_bins:],x[:rot_bins])
	else:
		out=np.append(x[len(x)-abs(rot_bins):],x[:len(x)-abs(rot_bins)])
	return out

# bin-scrunch the profile
def bestprof_bscrunch (x, bscr):
	out=np.array(x)
	out.shape=(len(x)/bscr, bscr)
	out=np.mean(out, axis=1)
	return out

# normalizing ind channel
def get_mean_rms (x, osm_min, osm_max):
        osm, osr = sc.probplot(x, sparams=(), dist='norm', fit=0)
	if np.size(np.where(osm > osm_max)) == 0:
		q_max = np.size(osm)
	else:
	        q_max = np.min(np.where(osm > osm_max))
	if np.size(np.where(osm < osm_min)) == 0:
		q_min = 0
	else:
	        q_min = np.max(np.where(osm < osm_min))
        rms, mean = np.polyfit(osm[q_min:q_max], osr[q_min:q_max], 1)
        return (mean, rms, osm, osr, q_min, q_max)

# exclude single bins representating 1-bin outliers
def trim_bins(x):
	x_diffs=[x[ii]-x[ii-1] for ii in xrange(1, len(x))]
	# trim left side
	cut_bin = 0
	for ii in xrange(0, len(x_diffs)/2):
		if x_diffs[ii] == 1: 
			if cut_bin != 0 : x=x[cut_bin:]
			break
		else: cut_bin += 1
	# trim right side
	cut_bin = 0
	for ii in xrange(len(x_diffs)-1, len(x_diffs)/2, -1):
		if x_diffs[ii] == 1: 
			if cut_bin != 0: x=x[:-cut_bin]
			break
		else: cut_bin += 1
	# trim in the middle
	x_diffs=[x[ii]-x[ii-1] for ii in xrange(1, len(x))]
	ii_to_trim=[]
	prev = 1
	for ii in xrange(0, len(x_diffs)):
		if x_diffs[ii] != 1 and prev == 1:
			prev = x_diffs[ii]
		elif x_diffs[ii] != 1 and prev != 1:
			ii_to_trim.append(ii)
			prev = x_diffs[ii]
		else: prev = 1
	x=np.delete(x, ii_to_trim, axis=0)
	x_diffs=[x[ii]-x[ii-1] for ii in xrange(1, len(x))]
	return x

# automatic search for the off-pulse window
# input profile will be rotated as necessary
# return tuple (data, rotphase, off-left, off-right)
def auto_find_off_window(data, rot_bins, nbins):
	# find first the bin with maximum value
	maxbin = np.argmax(data)
	# exclude the area of 60% of all bins around the maxbin
	# make the 60%-area the even number
	exclsize=int(nbins*0.6)+int(nbins*0.6)%2
	le=maxbin-exclsize/2
	re=maxbin+exclsize/2
	# extra rotation by "le" bins, so left edge will be at 0
	data = bestprof_rotate(data, le)
	# total rotation in phase
	if abs(rot_bins) < 1:
		rot_bins += float(le)/nbins
	else:
		rot_bins = float(rot_bins + le)/nbins
	amean = np.mean(data[re-le:nbins])
	arms = np.std(data[re-le:nbins])
	aprof = (data - amean)/arms
       	abins=np.arange(0,nbins)[(aprof>2.5)]
	abins=trim_bins(abins) # trimming bins
	# updating pulse window
	exclsize=abins[-1]-abins[0]
	# to be extra-cautious, increase it by 10% of the pulse window on both sides
	le=abins[0]-int(0.1*exclsize)
	re=abins[-1]+1+int(0.1*exclsize)
	# extra rotation by "le" bins again, so left edge will be at 0
	data = bestprof_rotate(data, le)
	# total rotation in phase
	rot_bins += float(le)/nbins
	return (data, rot_bins, re-le, nbins)


# Main
if __name__=="__main__":
        #
        # Parsing the command line options
        #
        usage = "Usage: %prog [options] .bestprof|.ar [.bestprof|.ar] [...]"
        cmdline = opt.OptionParser(usage)
	cmdline.add_option('--presto', dest='is_presto', action="store_true", help="if input file is Presto's .bestprof file,\
otherwise scripts assumes Prschive's archive file", default=False)
        cmdline.add_option('--snrmethod', dest='method', metavar='STRING', help="Method to calculate the mean/rms. Possible values \
are: 'Off', 'QQ', 'Polynom', 'Psrstat'. 'Off' is when user specifies off-pulse window to use, 'QQ' is to use Q-Q probability plot \
to determine range of quantiles that satisfies Gaussian distribution, 'Polynom' is like 'Off' but all profile bins are used to calculate \
mean/rms after subtracting polynomial fit to the pulse profile, and 'Psrstat' is when values of off:avg, off:rms, etc. reported by \
psrstat command are used. 'Psrstat' method is valid only when --presto is \
not used. Default: %default", default="QQ", type='str')
        cmdline.add_option('-b', '--bscrunch', dest='bscr', metavar='FACTOR', help="Bscrunch factor, \
default: %default", default=1, type='int')
        cmdline.add_option('-t', '--threshold', dest='thres', metavar='FACTOR', help="Threshold sigma \
value of the profile bins, default: %default", default=5, type='float')
        cmdline.add_option('-r', '--rotate', dest='rot_bins', metavar='#BIN|PHASE', help="Rotate profile by this number of bins (if absolute value >= 1). \
If absolute value is < 1, then the value is treated as pulse phase in turns. Negative values are to move right, default: %default", default=0, type='float')
        cmdline.add_option('--off-left', dest='off_left', metavar='BIN#', help="Left edge of the off-pulse window to calculate \
mean/rms in the 'Off' mode (inclusive). Default: %default", default=0, type='int')
        cmdline.add_option('--off-right', dest='off_right', metavar='BIN#', help="Right edge of the off-pulse window to calculate \
mean/rms in the 'Off' mode (exclusive). Default: 10%-bin of the profile", default=-1, type='int')
	cmdline.add_option('--auto-off', dest='is_auto_off', action="store_true", help="Automatic determination of OFF region", default=False)
        cmdline.add_option('--osm-min', dest='osm_min', metavar='MIN PROB', help="Minimum probability value to be used \
to calculate the mean and rms, default: min possible", type='float')
        cmdline.add_option('--osm-max', dest='osm_max', metavar='MAX PROB', help="Maximum probability value to be used \
to calculate the mean and rms, default: %default", default=+0.95, type='float')
        cmdline.add_option('--polynom-deg', dest='polynom_ndeg', metavar='DEGREE', help="Degree of a polynomial for 'Polynom' method. \
Default is maximum fittable starting from number of bins minus 1", default=-1, type='int')
        cmdline.add_option('--histbins', dest='histbins', metavar='#HISTBINS', help="Number of histogram bins to see samples distribution, \
default: %default", default=200, type='int')
	cmdline.add_option('--plot', dest='is_plot', action="store_true", help="Make diagnostic plots", default=False)
	cmdline.add_option('--saveonly', dest='is_saveonly', action="store_true", help="Save diagnostic plot \
to png-file instead of GUI", default=False)
        cmdline.add_option('--oneliner', dest='oneliner', metavar='PSR', help="To print info in one line comprising all interesting \
information. The argument is the pulsar name or any other label. If argument is 'meta', then source name will be taken from the file header", default="", type='str')

        # reading cmd options
        (opts,args) = cmdline.parse_args()

        # check if input file is given
        if len(args) == 0:
                cmdline.print_usage()
                sys.exit(0)

	if not opts.is_presto:
		import psrchive as pc
	else:
		import bestprof as bp

	if opts.is_presto and opts.method == "Psrstat":
		print "You can't use 'Psrstat' method when option --presto is used!"
		sys.exit(1)

	for infile in args:
		if opts.is_presto: # Presto's .bestprof
			bpobj = bp.bestprof(infile)
			data = bpobj.profile
			target = bpobj.psr
			if opts.bscr > 1:
				data = bestprof_bscrunch(data, opts.bscr)
			nbins = len(data)

			# auto-find of the OFF-pulse window
			if opts.is_auto_off:
				(data, opts.rot_bins, opts.off_left, opts.off_right) = auto_find_off_window(data, opts.rot_bins, nbins)
			else:
				if opts.rot_bins != 0:
					data = bestprof_rotate(data, opts.rot_bins)

		else: # Psrchive's file
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

	        	r = raw.get_data()
	        	#time stokes f phase
        		data = r[0,0,0,:]
		        weights = raw.get_weights()
			data[(weights[0]==0)] = 0.0

			# auto-find of the OFF-pulse window
			if opts.is_auto_off:
				(data, opts.rot_bins, opts.off_left, opts.off_right) = auto_find_off_window(data, opts.rot_bins, nbins)
				# and do total rotation for the input file as well (for psrstat)
				raw.rotate_phase(opts.rot_bins)
			else:
				if opts.rot_bins != 0:
					if abs(opts.rot_bins) < 1:
						raw.rotate_phase(opts.rot_bins)
					else:
						raw.rotate_phase(opts.rot_bins/nbins)
					data = bestprof_rotate(data, opts.rot_bins)

			# Psrstat
			cmd="psrstat -Q -c all:sum %s" % (infile) + " | awk '{printf \"%.20f\", $2}'"
                        allsum=float(os.popen(cmd).readlines()[0][:-1])
			cmd="psrstat -Q -c on:sum %s | awk '{print $2}' -" % (infile)
                        onsum=float(os.popen(cmd).readlines()[0][:-1])
			cmd="psrstat -Q -c on:max %s | awk '{print $2}' -" % (infile)
                        onmax=float(os.popen(cmd).readlines()[0][:-1])
			cmd="psrstat -Q -c snr %s" % (infile) + " | awk '{printf \"%.5f\", $2}'"
                        snr=float(os.popen(cmd).readlines()[0][:-1])
			cmd="psrstat -Q -c on:count %s | awk '{print $2}' -" % (infile)
                        oncount=int(os.popen(cmd).readlines()[0][:-1])
			cmd="psrstat -Q -c nbin %s | awk '{print $2}' -" % (infile)
                        nbin=int(os.popen(cmd).readlines()[0][:-1])
			cmd="psrstat -Q -c off:avg %s | awk '{print $2}' -" % (infile)
                        offavg=float(os.popen(cmd).readlines()[0][:-1])
			cmd="psrstat -Q -c off:rms %s | awk '{print $2}' -" % (infile)
                        offrms=float(os.popen(cmd).readlines()[0][:-1])
			if offavg == 0 or offrms == 0:
				cmd="psrstat -Q -c off=minimum -c off:avg %s | awk '{print $2}' -" % (infile)
	                        offavg=float(os.popen(cmd).readlines()[0][:-1])
				cmd="psrstat -Q -c off=minimum -c off:rms %s | awk '{print $2}' -" % (infile)
                	        offrms=float(os.popen(cmd).readlines()[0][:-1])
			if onsum == 0 or onmax == 0 or oncount == 0:
				cmd="psrstat -Q -c on=above -c on:sum %s | awk '{print $2}' -" % (infile)
	                        onsum=float(os.popen(cmd).readlines()[0][:-1])
				cmd="psrstat -Q -c on=above -c on:max %s | awk '{print $2}' -" % (infile)
                	        onmax=float(os.popen(cmd).readlines()[0][:-1])
				cmd="psrstat -Q -c on=above -c on:count %s | awk '{print $2}' -" % (infile)
	                        oncount=int(os.popen(cmd).readlines()[0][:-1])
			weff_bin=(allsum-nbin*offavg)/(onmax-offavg)
			psrstat_snrpeak = (onmax-offavg)/offrms
			psrstat_snrmean = (allsum-nbin*offavg)/(offrms*nbin)
			psrstat_profsign = (allsum-nbin*offavg)/(offrms*np.sqrt(weff_bin))
			psrstat_prof = (data - offavg)/offrms
			psrstat_binpeak = np.argmax(psrstat_prof)

		if opts.off_right == -1:
			opts.off_right = int(nbins*0.1)
		if opts.off_right-opts.off_left<=1:
			opts.off_right = nbins

		# Q-Q
		(qq_mean, qq_rms, osm, osr, qmin, qmax) = get_mean_rms(data, opts.osm_min, opts.osm_max)
		qq_prof = (data - qq_mean)/qq_rms
	        qq_snrpeak = np.max(qq_prof)
		qq_binpeak = np.argmax(qq_prof)
        	qq_snrmean = np.mean(qq_prof)
		qq_weq = np.sum(qq_prof)/qq_snrpeak
		qq_profsign = np.sum(qq_prof)/np.sqrt(qq_weq)
		qq_chi2 = np.sum(qq_prof*qq_prof)/(nbins-1)

	        qq_crit=(qq_prof>opts.thres)
	        qq_bins=np.arange(0,nbins)[qq_crit]
		qq_bins=trim_bins(qq_bins) # trimming bins
	        qq_profcrit=qq_prof[qq_bins]
        	qq_oncount=np.size(qq_bins)
			

		# Range
		range_mean = np.mean(data[opts.off_left:opts.off_right])
		range_rms = np.std(data[opts.off_left:opts.off_right])
		range_prof = (data - range_mean)/range_rms
	        range_snrpeak = np.max(range_prof)
		range_binpeak = np.argmax(range_prof)
        	range_snrmean = np.mean(range_prof)
		range_weq = np.sum(range_prof)/range_snrpeak
		range_profsign = np.sum(range_prof)/np.sqrt(range_weq)
		range_chi2 = np.sum(range_prof*range_prof)/(nbins-1)

	        range_crit=(range_prof>opts.thres)
	        range_bins=np.arange(0,nbins)[range_crit]
		range_bins=trim_bins(range_bins) # trimming bins
	        range_profcrit=range_prof[range_bins]
        	range_oncount=np.size(range_bins)

		# Polynom
		if opts.polynom_ndeg == -1:
			for polynom_ndeg in xrange(int(float(len(data)-1)),1,-1):
				try:
					tmp=float(len(data)-1)**polynom_ndeg
					break
				except: pass
		else:
			polynom_ndeg = opts.polynom_ndeg
		polynom_coeffs = np.polyfit(range(len(data)), data, polynom_ndeg)
		polynom_bline = np.polyval(polynom_coeffs, range(len(data)))
		if opts.method == "Polynom":
			polynom_mean = np.mean(data-polynom_bline)+np.mean(sorted(polynom_bline)[:int(0.2*len(polynom_bline))])
		elif opts.method == "Off":
			polynom_mean = range_mean
		elif opts.method == "QQ":
			polynom_mean = qq_mean
		elif opts.method == "Psrstat":
			polynom_mean = offavg
		polynom_rms = np.std(data-polynom_bline)
		polynom_prof = (data - polynom_mean)/polynom_rms
		polynom_fit = (polynom_bline - polynom_mean)/polynom_rms
	        polynom_snrpeak = np.max(polynom_prof)
		polynom_binpeak = np.argmax(polynom_prof)
        	polynom_snrmean = np.mean(polynom_prof)
		polynom_weq = np.sum(polynom_prof)/polynom_snrpeak
		polynom_profsign = np.sum(polynom_prof)/np.sqrt(polynom_weq)
		polynom_chi2 = np.sum(polynom_prof*polynom_prof)/(nbins-1)

	        polynom_crit=(polynom_prof>opts.thres)
	        polynom_bins=np.arange(0,nbins)[polynom_crit]
		polynom_bins=trim_bins(polynom_bins) # trimming bins
	        polynom_profcrit=polynom_prof[polynom_bins]
        	polynom_oncount=np.size(polynom_bins)

		if not opts.is_presto:
			if opts.oneliner == "":
				print infile
				print "-----------------------------------------------------------------------------------------"
				print "\t\t\t|\tQQ\t|    Off\t| Polynom (%d)\t| Psrstat" % (polynom_ndeg)
				print "-----------------------------------------------------------------------------------------"
				print "Mean:\t\t\t| %-13g | %-13g | %-13g | %-13g" % (qq_mean, range_mean, polynom_mean, offavg)
				print "RMS:\t\t\t| %-13g | %-13g | %-13g | %-13g" % (qq_rms, range_rms, polynom_rms, offrms)
				print "Peak S/N:\t\t| %-13g | %-13g | %-13g | %-13g" % (qq_snrpeak, range_snrpeak, polynom_snrpeak, psrstat_snrpeak)
				print "Peak bin:\t\t| %-13d | %-13d | %-13d | %-13d" % (qq_binpeak, range_binpeak, polynom_binpeak, psrstat_binpeak)
				print "Peak phase:\t\t| %-13g | %-13g | %-13g | %-13g" % (float(qq_binpeak)/float(nbins), float(range_binpeak)/float(nbins), \
										float(polynom_binpeak)/float(nbins), float(psrstat_binpeak)/float(nbins))
				print "Mean S/N:\t\t| %-13g | %-13g | %-13g | %-13g" % (qq_snrmean, range_snrmean, polynom_snrmean, psrstat_snrmean)
				print "Eff width (bins):\t| %-13d | %-13d | %-13d | %-13d" % (qq_weq, range_weq, polynom_weq, weff_bin)
				print "Weff/P ratio:\t\t| %-13g | %-13g | %-13g | %-13g" % (qq_weq/nbins, range_weq/nbins, polynom_weq/nbins, weff_bin/nbins)
				print "-----------------------------------------------------------------------------------------"
				print "S/N:\t\t\t| %-13.2f | %-13.2f | %-13.2f | %-13.2f" % (qq_profsign, range_profsign, polynom_profsign, psrstat_profsign)
				print "Chi^2/dof:\t\t| %-13.2f | %-13.2f | %-13.2f | %.2f (-c snr)" % (qq_chi2, range_chi2, polynom_chi2, snr)
				print "#bins with S/N>%.2f:\t| %-13d | %-13d | %-13d | %-13d" % (opts.thres, qq_oncount, range_oncount, polynom_oncount, oncount)
				print "-----------------------------------------------------------------------------------------"
				if opts.is_auto_off:
					print "AUTO-OFF: --off-left %d --off-right %d -r %f -b %d" % (opts.off_left, opts.off_right, opts.rot_bins, opts.bscr)	
			else:
				if opts.oneliner == "meta": psr = target
				else: psr = opts.oneliner
				print "# PSR\t\tS/Npeak (off, polynom, psrstat)\tS/Nmean (off, polynom, psrstat)\tS/N (off, polynom, psrstat)\tEff widths in bins \
(off, polynom, psrstat)\tW/P (off, polynom, psrstat)"
				print "#--------------------------------------------------------------------------------------------------------------------------\
-----------------------------------------------------------------"
				print "%-10s\t%g\t%g\t%g\t\t%g\t%g\t%g\t\t%.2f\t%.2f\t%.2f\t\t%d\t%d\t%d\t\t%g\t%g\t%g" % (psr, range_snrpeak, polynom_snrpeak, psrstat_snrpeak, \
					range_snrmean, polynom_snrmean, psrstat_snrmean, range_profsign, polynom_profsign, psrstat_profsign, range_weq, polynom_weq, weff_bin, \
					range_weq/nbins, polynom_weq/nbins, weff_bin/nbins)

		else: # for .bestprof files (the same but without printing values from psrstat)
			if opts.oneliner == "":
				print infile
				print "-----------------------------------------------------------------------"
				print "\t\t\t|\tQQ\t|    Off\t| Polynom (%d)" % (polynom_ndeg)
				print "-----------------------------------------------------------------------"
				print "Mean:\t\t\t| %-13g | %-13g | %-13g" % (qq_mean, range_mean, polynom_mean)
				print "RMS:\t\t\t| %-13g | %-13g | %-13g" % (qq_rms, range_rms, polynom_rms)
				print "Peak S/N:\t\t| %-13g | %-13g | %-13g" % (qq_snrpeak, range_snrpeak, polynom_snrpeak)
				print "Peak bin:\t\t| %-13d | %-13d | %-13d" % (qq_binpeak, range_binpeak, polynom_binpeak)
				print "Peak phase:\t\t| %-13g | %-13g | %-13g" % (float(qq_binpeak)/float(nbins), float(range_binpeak)/float(nbins), \
										float(polynom_binpeak)/float(nbins))
				print "Mean S/N:\t\t| %-13g | %-13g | %-13g" % (qq_snrmean, range_snrmean, polynom_snrmean)
				print "Eff width (bins):\t| %-13d | %-13d | %-13d" % (qq_weq, range_weq, polynom_weq)
				print "Weff/P ratio:\t\t| %-13g | %-13g | %-13g" % (qq_weq/nbins, range_weq/nbins, polynom_weq/nbins)
				print "-----------------------------------------------------------------------"
				print "S/N:\t\t\t| %-13.2f | %-13.2f | %-13.2f" % (qq_profsign, range_profsign, polynom_profsign)
				print "Chi^2/dof:\t\t| %-13.2f | %-13.2f | %-13.2f" % (qq_chi2, range_chi2, polynom_chi2)
				print "#bins with S/N>%.2f:\t| %-13d | %-13d | %-13d" % (opts.thres, qq_oncount, range_oncount, polynom_oncount)
				print "-----------------------------------------------------------------------"
				if opts.is_auto_off:
					print "AUTO-OFF: --off-left %d --off-right %d -r %f -b %d" % (opts.off_left, opts.off_right, opts.rot_bins, opts.bscr)	
			else:
				if opts.oneliner == "meta": psr = target
				else: psr = opts.oneliner
				print "# PSR\t\tS/Npeak (off, polynom)\tS/Nmean (off, polynom)\tS/N (off, polynom)\tEff widths in bins (off, polynom)\tW/P (off, polynom)"
				print "#--------------------------------------------------------------------------------------------------------------------------------------------------"
				print "%-10s\t%g\t%g\t\t%g\t%g\t\t%.2f\t%.2f\t\t%d\t%d\t\t\t%g\t%g" % (psr, range_snrpeak, polynom_snrpeak, \
					range_snrmean, polynom_snrmean, range_profsign, polynom_profsign, range_weq, polynom_weq, range_weq/nbins, polynom_weq/nbins)

		# Plotting
		if opts.is_plot:
			if opts.is_saveonly:
				pngname = ".".join(infile.split(".")[:-1]) + ".png"
               			import matplotlib
		                matplotlib.use("Agg")

		        import matplotlib.pyplot as plt
		        import matplotlib.cm as cm
		        from matplotlib.ticker import *

			if opts.method == "QQ":
				fig = plt.figure(figsize=(8,9))
				ax1 = fig.add_subplot(311)
			else:
				fig = plt.figure(figsize=(8,6))
				ax1 = fig.add_subplot(211)
			plt.xlabel("Bin")
			plt.ylabel("S/N")
			if opts.method == "QQ":
				ax1.plot(range(len(qq_prof)), qq_prof, "g-", alpha=0.7)
				ax1.axvline(x=qq_binpeak, linestyle="--", color="black")
				if len(qq_bins) != 0:
					on_start = qq_bins[0]	
					for ii in xrange(1, len(qq_bins)):
						if qq_bins[ii]-qq_bins[ii-1] == 1: continue
						else:
							ax1.axvspan(on_start, qq_bins[ii-1], color="green", alpha=0.1)
							on_start=qq_bins[ii]
					ax1.axvspan(on_start, qq_bins[-1], color="green", alpha=0.1)
				ax1.set_xlim(xmin=0, xmax=len(qq_prof)-1)
			if opts.method == "Off":
				ax1.plot(range(len(range_prof)), range_prof, "g-", alpha=0.7)
				ax1.axvline(x=range_binpeak, linestyle="--", color="black")
				ax1.axvspan(opts.off_left, opts.off_right, color="yellow", alpha=0.2)
				if len(range_bins) != 0:
					on_start = range_bins[0]	
					for ii in xrange(1, len(range_bins)):
						if range_bins[ii]-range_bins[ii-1] == 1: continue
						else:
							ax1.axvspan(on_start, range_bins[ii-1], color="green", alpha=0.1)
							on_start=range_bins[ii]
					ax1.axvspan(on_start, range_bins[-1], color="green", alpha=0.1)
				ax1.axvline(x=opts.off_left, linestyle=":", color="black")
				ax1.axvline(x=opts.off_right, linestyle=":", color="black")
				ax1.set_xlim(xmin=0, xmax=len(range_prof)-1)
			if opts.method == "Polynom":
				ax1.plot(range(len(polynom_prof)), polynom_prof, "g-", alpha=0.7)
				ax1.plot(range(len(polynom_prof)), polynom_fit, "b-", label="fit, ndeg=%d" % (polynom_ndeg), alpha=0.5)
				ax1.axvline(x=polynom_binpeak, linestyle="--", color="black")
				if len(polynom_bins) != 0:
					on_start = polynom_bins[0]	
					for ii in xrange(1, len(polynom_bins)):
						if polynom_bins[ii]-polynom_bins[ii-1] == 1: continue
						else:
							ax1.axvspan(on_start, polynom_bins[ii-1], color="green", alpha=0.1)
							on_start=polynom_bins[ii]
					ax1.axvspan(on_start, polynom_bins[-1], color="green", alpha=0.1)
				ax1.set_xlim(xmin=0, xmax=len(polynom_prof)-1)
			if opts.method == "Psrstat":
				ax1.plot(range(len(psrstat_prof)), psrstat_prof, "g-", alpha=0.7)
				ax1.axvline(x=psrstat_binpeak, linestyle="--", color="black")
				if len(psrstat_bins) != 0:
					on_start = psrstat_bins[0]	
					for ii in xrange(1, len(psrstat_bins)):
						if psrstat_bins[ii]-psrstat_bins[ii-1] == 1: continue
						else:
							ax1.axvspan(on_start, psrstat_bins[ii-1], color="green", alpha=0.1)
							on_start=psrstat_bins[ii]
					ax1.axvspan(on_start, psrstat_bins[-1], color="green", alpha=0.1)
				ax1.set_xlim(xmin=0, xmax=len(psrstat_prof)-1)
			ax1.axhline(y=0.0, linestyle="--", color="black")
			ax1.axhspan(-1, +1, color="grey", alpha=0.1)
			#plt.grid()
			plt.gca().minorticks_on()

			if opts.method == "QQ":
				ax2 = fig.add_subplot(312)
			else:
				ax2 = fig.add_subplot(212)
			plt.xlabel("Data sample values")
			plt.ylabel("Number")
			if opts.method == "QQ":
				histmean = qq_mean
				histrms = qq_rms
			if opts.method == "Off":
				histmean = range_mean
				histrms = range_rms
			if opts.method == "Polynom":
				histmean = polynom_mean
				histrms = polynom_rms
			if opts.method == "Psrstat":
				histmean = offavg
				histrms = offrms
			from scipy.optimize import leastsq
			fitfunc  = lambda p, x: p[0]*np.exp(-0.5*((x-histmean)/histrms)**2)+p[1]
			errfunc  = lambda p, x, y: (y - fitfunc(p, x))			
			init  = [1.0, 0.5]
			n, bins, patches = plt.hist(data, bins=opts.histbins, color='skyblue', alpha=0.4)
			out   = leastsq(errfunc, init, args=(bins[:opts.histbins], n[:opts.histbins]))
			plt.plot(bins, fitfunc(out[0], bins), color='green', linewidth=2, label='$\mu$ = %g\n$\sigma$ = %g' % (histmean, histrms))
			plt.legend()

			if opts.method == "QQ":
				ax3 = fig.add_subplot(313)
				plt.xlabel("OSM (quantiles)")
				plt.ylabel("OSR (ordered values)")
				osr=(osr-qq_mean)/qq_rms
				ax3.plot(osm, osr, "b-", alpha=0.7)
				ax3.axhline(y=osr[qmin], linestyle="--", color="orange")
				ax3.axhline(y=osr[qmax], linestyle="--", color="orange")
				ax3.axvline(x=osm[qmin], linestyle="--", color="orange")
				ax3.axvline(x=osm[qmax], linestyle="--", color="orange")
				plt.gca().minorticks_on()

	                if opts.is_saveonly:
        	                plt.savefig(pngname)
               		else:
                       		plt.show()
