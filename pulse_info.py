#!/usr/bin/env python
#
# This script does 1D search of pulse components in phase
# and then search of separate islands of emission in frequency
# domain with friends-of-friends algorithm
#
# Vlad Kondratiev, Jan 21, 2012
#
import numpy as np
import psrchive as pc
import os, os.path, stat, glob, sys, getopt, re
import scipy.stats as sc
import scipy.signal as ss
import optparse as opt
import gc

# normalizing each channel
def normalize(x, nch):
        for i in range(nch):
                osm, osr = sc.probplot(x[i], sparams=(), dist='norm', fit=0)
                q_max = np.min(np.where(osm > 1.0))
                q_min = np.max(np.where(osm < -1.0))
                rms, mean = np.polyfit(osm[q_min:q_max], osr[q_min:q_max], 1)
                x[i] -= mean
                x[i] /= rms
                crit=np.isfinite(x[i])
                x[i][-crit] = 0.0
	return x

# normalize profile
def prof_norm(x):
        osm, osr = sc.probplot(x, sparams=(), dist='norm', fit=0)
        q_max = np.min(np.where(osm > 1.0))
        q_min = np.max(np.where(osm < -1.0))
        rms, mean = np.polyfit(osm[q_min:q_max], osr[q_min:q_max], 1)
        x -= mean
        x /= rms
        crit=np.isfinite(x)
        x[-crit] = 0.0
        return x


# Main
if __name__=="__main__":

        #
        # Parsing the command line options
        #
        usage = "Usage: %prog [options] <PSRCHIVE fits-files>"
        cmdline = opt.OptionParser(usage)
        cmdline.add_option('-s', '--start', dest='start_bin', metavar='BIN',
                           help="start bin of ON-pulse window to search (default: %default)", default=0, type='int')
        cmdline.add_option('-w', '--win', dest='win_bin', metavar='#BINS',
                           help="size of ON-pulse window (default: number of bins in the files)", default=-1, type='int')
        cmdline.add_option('--off-start', dest='noise_start_bin', metavar='BIN',
                           help="start bin of OFF-pulse window to search (default: %default)", default=0, type='int')
        cmdline.add_option('--out', dest='outfile', metavar='FILE',
                           help="output file with search info for each input file (default: %default)", default="pulse_info.txt", type='str')
        cmdline.add_option('-n', '--numcomp', dest='ncomp', metavar='#COMP',
                           help="number of components to search (default: %default)", default=1, type='int')
        cmdline.add_option('-f', '--fscrunch', dest='fscr', metavar='FACTOR',
                           help="scrunch in frequency by a factor (default: %default)", default=1, type='int')
        cmdline.add_option('-b', '--bscrunch', dest='binscr', metavar='FACTOR',
                           help="scrunch in phase by a factor (default: %default)", default=1, type='int')
        cmdline.add_option('--threshold', dest='threshold', metavar='SNR',
                           help="cut-off threshold (in spectrum's sigmas) in frequency-domain search (default: %default)", default=7, type='float')
        cmdline.add_option('--td', dest='tdstr', metavar='LIST_OF_DOWNFACTS',
                           help="comma-separated list of time downfacts in bins for a search (default: %default)", default="2,3,4,6,9,14,20", type='str')
        cmdline.add_option('-i', '--indexing', action="store_true", dest="is_indexing",
                           help="use local indexing of pulses (files) instead of pulse number from the filename", default=False)
        cmdline.add_option('--saveonly', action="store_true", dest="is_saveonly", help="only saves png-file and exits", default=False)

        group = opt.OptionGroup(cmdline, "Graphics Options")
        group.add_option('--fontsize', dest='fs', metavar='SIZE',
                           help="font size for labels (default: %default)", default=10, type='int')
	cmdline.add_option_group(group)

        # reading cmd options
        (opts,args) = cmdline.parse_args()

        # check if input file is given
        if len(args) == 0:
                cmdline.print_usage()
                sys.exit(0)

	# create downfacts lists
	time_downfacts = [int(s) for s in opts.tdstr.split(',')]

	# checking the number of given components
	if opts.ncomp < 1:
		print "Number of components to search should be at least 1!"
		sys.exit(1)

	# we have to check that the last (maximum) time downfact is not larger than the size of the window
	if opts.win_bin != -1:
		if time_downfacts[-1] > opts.win_bin:
			print "Error: the largest time downfactor %d is larger than the size of ON-pulse window %d!" % (time_downfacts[-1], opts.win_bin)
			sys.exit(1)

	outf = open(opts.outfile, "wt")
	outf.write("# Pulse#\tPeak Prof SNR\tAverage Prof SNR\tToffset (bins)\tTwidth (bins)\tPeak SNR\tAverage SNR\tFoffset (bins)\tFwidth (bins)\n")
	count = 0

	if opts.is_saveonly:
		import matplotlib
		matplotlib.use("Agg")

	import matplotlib.pyplot as plt
	from matplotlib.patches import Rectangle
	import matplotlib.cm as cm
	from matplotlib.ticker import *

	# list of unit arrays
	psfs = []
	for td in range(len(time_downfacts)):
		psfs.append(np.ones((time_downfacts[td]), dtype=np.float32))

	# loop on input files
	for fitsfile in args:
		print fitsfile
		if opts.is_indexing:
			pulseN = count
		else:
			pulseN = int(fitsfile.split("/")[-1].split(".")[0].split("_")[1])
		pngname = ".".join(fitsfile.split("/")[-1].split(".")[:-1]) + ".png"
		# read input file using PSRCHIVE module
        	raw = pc.Archive_load(fitsfile)
	        raw.pscrunch()  # add all polarizations together
	        if opts.fscr > 1: raw.fscrunch(opts.fscr)   # fscrunch by fscr factor
	        if opts.binscr > 1: raw.bscrunch(opts.binscr) # bscrunch by binscr factor

		nphbin = raw.get_nbin()
		if opts.win_bin == -1:
			opts.win_bin = nphbin
			if time_downfacts[-1] > opts.win_bin:
				print "Error: the largest time downfactor %d is larger than the size of ON-pulse window %d!" % (time_downfacts[-1], opts.win_bin)
				sys.exit(1)
		nchan = raw.get_nchan()
		cfreq = raw.get_centre_frequency()
		bw = raw.get_bandwidth()
		chan_bw = bw / nchan
		if not(raw.get_dedispersed()):
			raw.dedisperse()

		r = raw.get_data()
		#time stokes f phase
        	data = r[0,0,:,:]

		# normalizing...
		data = normalize(data, nchan)

		weights = raw.get_weights()
		data[(weights[0]==0)] = 0.0
	        globmax=np.max(np.abs(data))
		globmin=-globmax

		# getting profile
		profile = np.sum(data, axis=0)
		profile = prof_norm(profile)
		# backup profile for plotting
		backup = np.copy(profile)
		
		# array of phases
		phase = np.arange(nphbin)*1.0/nphbin
		# patches of emission to plot
		rects = nrects = loffs = roffs = []

		# searching opts.ncomp times
		for tt in range(opts.ncomp):
			# array with max values for each box-function
			maxes = np.zeros(len(time_downfacts), dtype=np.float32)
			# array with time indexes for maximum values 
			maxtime = np.zeros(len(time_downfacts), dtype=int)

			# searching....
			for td in np.arange(len(time_downfacts)):
				out = ss.convolve(profile[opts.start_bin:opts.start_bin+opts.win_bin], psfs[td], mode='same')
				maxes[td] = np.max(out) / np.sqrt(time_downfacts[td])
				maxtime[td] = np.argmax(out) + opts.start_bin

			td = np.argmax(maxes)
			tw = time_downfacts[td]
			aver_prof_snr = np.max(maxes) / np.sqrt(tw) # this is really average snr, to calculate energy - just multiply by N
			toffset = int(maxtime[td]-tw/2)
			peak_prof_snr = np.max(profile[toffset:toffset+tw])
			loffs = np.append(loffs, float(toffset)/nphbin)
			roffs = np.append(roffs, float(toffset+tw)/nphbin)

			# getting pulse spectrum
			spectrum = np.sum(data[:,toffset:toffset+tw], axis=1)
			peak_snr = np.max(spectrum) / np.sqrt(tw)
			# if max value is less than threshold, then no search and we take the whole spectrum as single large patch
			if peak_snr < opts.threshold:
				fw = nchan
				foffset = 0
				aver_snr = np.sum(data[:,toffset:toffset+tw]) / (tw * fw)
				rects.append(Rectangle((float(toffset)/nphbin, cfreq-bw/2.+foffset*chan_bw), float(tw)/nphbin, fw*chan_bw, fill=False, linewidth=1, edgecolor='grey'))
				print count, pulseN, peak_prof_snr, aver_prof_snr, toffset, tw, peak_snr, aver_snr, foffset, fw
				str = "%d\t%f\t%f\t\t%d\t\t%d\t\t%f\t%f\t%d\t\t%d" % (pulseN, peak_prof_snr, aver_prof_snr, toffset, tw, peak_snr, aver_snr, foffset, fw)
				outf.write(str + "\n")
			else: # friends-of-friends searching...
				while peak_snr >= opts.threshold:	
					ind_max = np.argmax(spectrum)
					ledge = ind_max - 1 < 0 and 0 or ind_max - 1
					while ledge >= 0 and spectrum[ledge] >= opts.threshold: ledge -= 1
					ledge += 1
					redge = ind_max + 1 >= nchan and nchan - 1 or ind_max + 1
					while redge < nchan and spectrum[redge] >= opts.threshold: redge += 1
					foffset = ledge
					fw = redge - ledge
					aver_snr = np.sum(data[foffset:foffset+fw:,toffset:toffset+tw]) / (tw * fw)
					rects.append(Rectangle((float(toffset)/nphbin, cfreq-bw/2.+foffset*chan_bw), float(tw)/nphbin, fw*chan_bw, fill=False, linewidth=1, edgecolor='grey'))
					print count, pulseN, peak_prof_snr, aver_prof_snr, toffset, tw, peak_snr, aver_snr, foffset, fw
					str = "%d\t%f\t%f\t\t%d\t\t%d\t\t%f\t%f\t%d\t\t%d" % (pulseN, peak_prof_snr, aver_prof_snr, toffset, tw, peak_snr, aver_snr, foffset, fw)
					outf.write(str + "\n")
					spectrum[ledge:redge] = 0.0
					peak_snr = np.max(spectrum) / np.sqrt(tw)
		
			# noise searching...
			nmaxes = np.zeros(len(time_downfacts), dtype=np.float32)
			nmaxtime = np.zeros(len(time_downfacts), dtype=int)

			for ntd in np.arange(len(time_downfacts)):
				nout = ss.convolve(profile[opts.noise_start_bin:opts.noise_start_bin+opts.win_bin], psfs[ntd], mode='same')
				nmaxes[ntd] = np.max(nout) / np.sqrt(time_downfacts[ntd])
				nmaxtime[ntd] = np.argmax(nout) + opts.noise_start_bin

			ntd = np.argmax(nmaxes)
			ntw = time_downfacts[ntd]
			aver_prof_nsnr = np.max(nmaxes)/ np.sqrt(ntw) # this is really average snr, to calculate energy - just multiply by N
			ntoffset = int(nmaxtime[ntd]-ntw/2)
			peak_prof_nsnr = np.max(profile[ntoffset:ntoffset+ntw])

			# getting noise spectrum
			nspectrum = np.sum(data[:,ntoffset:ntoffset+ntw], axis=1)
			peak_nsnr = np.max(nspectrum) / np.sqrt(ntw)
			# if max value is less than threshold, then no search and we take the whole spectrum as single large patch
			if peak_nsnr < opts.threshold:
				nfw = nchan
				nfoffset = 0
				aver_nsnr = np.sum(data[:,ntoffset:ntoffset+ntw]) / (ntw * nfw)
				nrects.append(Rectangle((float(ntoffset)/nphbin, cfreq-bw/2.+nfoffset*chan_bw), float(ntw)/nphbin, nfw*chan_bw, fill=False, linewidth=1, edgecolor='white'))
				print count, pulseN, peak_prof_nsnr, aver_prof_nsnr, ntoffset, ntw, peak_nsnr, aver_nsnr, nfoffset, nfw
				str = "%d\t%f\t%f\t\t%d\t\t%d\t\t%f\t%f\t%d\t\t%d" % (pulseN, peak_prof_nsnr, aver_prof_nsnr, ntoffset, ntw, peak_nsnr, aver_nsnr, nfoffset, nfw)
				outf.write(str + "\n")
			else: # friends-of-friends searching...
				while peak_nsnr >= opts.threshold:	
					ind_max = np.argmax(nspectrum)
					ledge = ind_max - 1 < 0 and 0 or ind_max - 1
					while ledge >= 0 and nspectrum[ledge] >= opts.threshold: ledge -= 1
					ledge += 1
					redge = ind_max + 1 >= nchan and nchan - 1 or ind_max + 1
					while redge < nchan and nspectrum[redge] >= opts.threshold: redge += 1
					nfoffset = ledge
					nfw = redge - ledge
					aver_nsnr = np.sum(data[nfoffset:nfoffset+nfw:,ntoffset:ntoffset+ntw]) / (ntw * nfw)
					nrects.append(Rectangle((float(ntoffset)/nphbin, cfreq-bw/2.+nfoffset*chan_bw), float(ntw)/nphbin, nfw*chan_bw, fill=False, linewidth=1, edgecolor='white'))
					print count, pulseN, peak_prof_nsnr, aver_prof_nsnr, ntoffset, ntw, peak_nsnr, aver_nsnr, nfoffset, nfw
					str = "%d\t%f\t%f\t\t%d\t\t%d\t\t%f\t%f\t%d\t\t%d" % (pulseN, peak_prof_nsnr, aver_prof_nsnr, ntoffset, ntw, peak_nsnr, aver_nsnr, nfoffset, nfw)
					outf.write(str + "\n")
					nspectrum[ledge:redge] = 0.0
					peak_nsnr = np.max(nspectrum) / np.sqrt(ntw)

			# zero'ing the part of the period already found
			# zero'ing a bit more in phase than the size of rectangle
			profile[toffset-1:toffset+tw+1] = 0.0
			# zero'ing on noise as well
			profile[ntoffset-1:ntoffset+ntw+1] = 0.0

			outf.flush()


		# Plotting

	        # dimensions of the central plot
        	left = 0.07
        	bottom = 0.08
        	width = 0.97
        	height = 0.65

		fig = plt.figure()
		ax = plt.axes([left, bottom, width, height])
		ax.xaxis.tick_bottom()
		ax.xaxis.set_label_position("bottom")
                for label in ax.get_xticklabels(): label.set_fontsize(opts.fs)
		ax.yaxis.tick_left()
		ax.yaxis.set_label_position("left")
		for label in ax.get_yticklabels(): label.set_fontsize(opts.fs)
                plt.ylabel('Frequency (MHz)', fontsize=opts.fs)
                plt.xlabel('Pulse phase', fontsize=opts.fs)
	
		colormap = cm.get_cmap("BrBG")
		cax = plt.imshow(data, aspect='auto', interpolation='nearest', origin='lower', cmap=colormap, extent=(0.0, 1.0, cfreq-bw/2., cfreq+bw/2.), vmin=globmin, vmax=globmax)
		cbar = fig.colorbar(cax, orientation='vertical', spacing='proportional')
		cbar.ax.set_ylabel("Flux density ($\sigma$)", fontsize=opts.fs)
		for label in cbar.ax.get_yticklabels(): label.set_fontsize(opts.fs)
		for r in rects: ax.add_artist(r)
		for r in nrects: ax.add_artist(r)

	        # profile
                ax3 = plt.axes([left, bottom+height+0.005, width-0.194, height-0.41])
                ax3.xaxis.tick_bottom()
                ax3.xaxis.set_label_position("bottom")
                ax3.xaxis.set_major_formatter(NullFormatter())
                ax3.yaxis.tick_left()
                ax3.yaxis.set_major_locator(MaxNLocator(5))
                for label in ax3.get_yticklabels(): label.set_fontsize(opts.fs)
                ax3.yaxis.set_label_position("left")
                for label in ax3.get_yticklabels():
                        label.set_rotation(90)
                        label.set_fontsize(opts.fs)
                plt.ylabel('Flux density ($\sigma$)', fontsize=opts.fs)
                ax3.plot(phase, backup, color="darkgreen", linewidth=1, label="pulse %d" % (pulseN))
                ax3.set_xlim(xmin=0, xmax=1)
                ax3.set_ylim(ymin=0.9*np.min(backup), ymax=1.1*np.max(backup))
		for ll, rr in zip(loffs, roffs):
			ax3.axvspan(ll, rr, facecolor='yellow', alpha=0.1)

                plt.legend(loc=1) # upper right

		if opts.is_saveonly: plt.savefig(pngname)
		else: plt.show()

		count += 1
		gc.collect() # calling garbage collection trying to help to free the memory

	outf.close()
