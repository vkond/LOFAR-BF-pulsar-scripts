#!/usr/bin/env python
#
import numpy as np
import psrchive as pc
import os, os.path, stat, glob, sys, getopt
import scipy.stats as sc
import optparse as opt

# normalizing each channel
def normalize(x, first_chan, nch, first_sub, nsubs):
	for s in range(first_sub, first_sub + nsubs):
        	for i in range(first_chan, first_chan + nch):
                	osm, osr = sc.probplot(x[s,i], sparams=(), dist='norm', fit=0)
                	q_max = np.min(np.where(osm > 1.0))
                	q_min = np.max(np.where(osm < -1.0))
                	rms, mean = np.polyfit(osm[q_min:q_max], osr[q_min:q_max], 1)
                	x[s,i] -= mean
			if rms == 0.0: x[s,i] = 0.0
                	else: x[s,i] /= rms
			crit=np.isfinite(x[s,i])
			x[s,i][-crit] = 0.0
        return x

# normalize profile
def prof_norm(x):
	osm, osr = sc.probplot(x, sparams=(), dist='norm', fit=0)
      	q_max = np.min(np.where(osm > 1.0))
      	q_min = np.max(np.where(osm < -1.0))
       	rms, mean = np.polyfit(osm[q_min:q_max], osr[q_min:q_max], 1)
       	x -= mean
	if rms == 0.0: x = 0.0
      	else: x /= rms
	crit=np.isfinite(x)
	x[-crit] = 0.0
        return x

# normalize profiles
def profs_norm(x, nsubs):
	for s in range(nsubs):
		osm, osr = sc.probplot(x[s], sparams=(), dist='norm', fit=0)
      		q_max = np.min(np.where(osm > 1.0))
      		q_min = np.max(np.where(osm < -1.0))
       		rms, mean = np.polyfit(osm[q_min:q_max], osr[q_min:q_max], 1)
       		x[s] -= mean
		if rms == 0.0: x[s] = 0.0
      		else: x[s] /= rms
		crit=np.isfinite(x[s])
		x[s][-crit] = 0.0
        return x


# Main
if __name__=="__main__":

	# parsing command line options
	cmd = opt.OptionParser("Usage: %prog <arfile> [-h|--help] [OPTIONS]")
	cmd.add_option('-f', '--fscrunch', dest='fscr', metavar='FACTOR', help="Fscrunch factor, default: %default", default=4, type='int')
	cmd.add_option('-b', '--bscrunch', dest='bscr', metavar='FACTOR', help="Bscrunch factor, default: %default", default=2, type='int')
	cmd.add_option('-p', '--pscrunch', action="store_true", dest='is_pscr', help="Polarization scrunch to total intensity ", default=False)
	cmd.add_option('-t', '--tscrunch', dest='tscr', metavar='FACTOR', help="Tscrunch factor, default: %default", default=1, type='int')
	cmd.add_option('--first-subint', dest='first_subint', metavar='NUMBER', help="First subint to use, default: %default", default=0, type='int')
	cmd.add_option('--total-nsubints', dest='nsubints', metavar='NUMBER', help="Number of subints to use, default: all", default=-1, type='int')
	cmd.add_option('-n', '--group-nsubints', dest='ngroups', metavar='NUMBER', help="Number of subints in a group to show at once, default: %default", default=100, type='int')
	cmd.add_option('--left-sp', dest='phleft', metavar='PHASE', help="Left phase edge for ON-pulse spectra, default: %default", default=0.0, type='float')
	cmd.add_option('--right-sp', dest='phright', metavar='PHASE', help="Right phase edge for ON-pulse spectra, default: %default", default=1.0, type='float')
	cmd.add_option('--off-pulse', action="store_true", dest='is_off_pulse', help="spectra will be calculated in the OFF-pulse window, i.e. phase ranges complementary to window given by --left-sp and --right-sp options", default=False)
	cmd.add_option('--phase-min', dest='lview', metavar='PHASE', help="Left phase edge for T-vs.-Phase subplot, default: %default", default=0.0, type='float')
	cmd.add_option('--phase-max', dest='rview', metavar='PHASE', help="Right phase edge for T-vs.-Phase subplot, default: %default", default=1.0, type='float')
	cmd.add_option('--fmin', dest='fmin', metavar='FREQ', help="Minimum frequency in MHz for T-vs.-Freq subplot", default=-1.0, type='float')
	cmd.add_option('--fmax', dest='fmax', metavar='FREQ', help="Maximum frequency in MHz for T-vs.-Freq subplot", default=-1.0, type='float')
	cmd.add_option('--no-global-norm', action="store_true", dest='is_noglobalnorm', help="turn off global normalization using all subints, otherwise normalization is done separately for every group of subints", default=False)
	cmd.add_option('-s', '--saveonly', action="store_true", dest='is_saveonly', help="do not show GUI window", default=False)
	cmd.add_option('--save', action="store_true", dest='is_save', help="turn on saving png files in GUI mode", default=False)
	cmd.add_option('--legend', dest='legend', metavar='STRING', help="additional info to put to the plot header, usually ObsID, etc.", default="", type='str')
	cmd.add_option('-o', '--output', dest='outstem', metavar='STRING', help="output pngfile prefix, default: %default", default="sp", type='str')
	cmd.add_option('--fontsize', dest='fs', metavar='SIZE', help="Font size, default: %default", default=10, type='int')

	# reading cmd options
	(opts, args) = cmd.parse_args()

        # check if input file is given
        if len(args) != 0:
                infile = args[0]
        else:
                cmd.print_usage()
                sys.exit(0)

        if opts.is_saveonly:
                import matplotlib
                matplotlib.use("Agg")

	import matplotlib.pyplot as plt
	import matplotlib.cm as cm
	from matplotlib.ticker import *
	import matplotlib.font_manager as fm

	print "loading..."
       	raw = pc.Archive_load(infile)
	if not(raw.get_dedispersed()):
		raw.dedisperse()
	if opts.is_pscr: raw.pscrunch()
	if opts.fscr > 1: raw.fscrunch(opts.fscr)
	if opts.bscr > 1: raw.bscrunch(opts.bscr)
	if opts.tscr > 1: raw.tscrunch(opts.tscr)
	nphbin = raw.get_nbin()
	nchan = raw.get_nchan()
	cfreq = raw.get_centre_frequency()
	bw = raw.get_bandwidth()
	nsubint = raw.get_nsubint()
	chan_bw = bw / nchan
	target = raw.get_type()
	source = raw.get_source()
	telescope = raw.get_telescope()
	# update number of subints
	if opts.nsubints != -1:
		if opts.first_subint + opts.nsubints > nsubint:
			opts.nsubints -= (opts.first_subint + opts.nsubints - nsubint) 
	else: opts.nsubints = nsubint - opts.first_subint
	# update number of subints in a group
	if opts.ngroups > opts.nsubints:
		opts.ngroups = opts.nsubints

	# values of bin edges for profile
	profleft = int(nphbin * opts.lview)
	profright = int(nphbin * opts.rview) + 1
	phase = np.arange(profleft, profright, 1)*1.0/nphbin
	# and for spectra
	binleft = int(nphbin * opts.phleft)
	binright = int(nphbin * opts.phright) + 1
	freqs = np.arange(nchan)*chan_bw + cfreq - bw/2.
	if opts.fmin == -1.:
		fminchan = 0
		opts.fmin = cfreq - bw/2.
	else:
		if opts.fmin < freqs[0] or opts.fmin >= freqs[-1]: opts.fmin = freqs[0]
		fminchan = int((opts.fmin - cfreq + bw/2.)/chan_bw)
	if opts.fmax == -1.:
		fmaxchan = nchan
		opts.fmax = cfreq + bw/2.
	else:
		if opts.fmax > freqs[-1]+chan_bw or opts.fmax <= opts.fmin: opts.fmax = freqs[-1]+chan_bw
		fmaxchan = int((opts.fmax - cfreq + bw/2.)/chan_bw)

	# number of channels actually to show
	nshowchan = fmaxchan - fminchan
	# recalculate freqs array
	freqs = np.arange(fminchan,fmaxchan)*chan_bw + cfreq - bw/2.

	r = raw.get_data()
	#time stokes f phase
       	data = r[opts.first_subint:opts.first_subint+opts.nsubints,0,fminchan:fmaxchan,:]
	weights = raw.get_weights()
	for s in np.arange(opts.nsubints): data[s,(weights[s]==0)] = 0.0

        # normalizing...
	if not opts.is_noglobalnorm:
		print "normalizing..."
        	data = normalize(data, 0, nshowchan, 0, opts.nsubints)

		# getting average profile
		averprof = prof_norm(np.sum(np.sum(data, axis=1), axis=0))
		averprof = averprof[profleft:profright]

		# table of profiles
		proflist=profs_norm(np.sum(data, axis=1), opts.nsubints)
		proflist=proflist[:,profleft:profright]
		profmin=np.min(proflist)
		profmax=np.max(proflist)

		# table of spectra
		if opts.is_off_pulse: # calculate OFF-pulse spectrum
			splist=np.mean(data[:,:,range(0,binleft)+range(binright,nphbin)], axis=2)
		else: # calculate ON-pulse spectrum
			splist=np.mean(data[:,:,binleft:binright], axis=2)
		spmin=np.min(splist)
		spmax=np.max(splist)
		# average spectrum
		aversp = np.mean(splist, axis=0)

	# Plotting
	print "plotting..."

        # dimensions of the central plot
        left = 0.08
        bottom = 0.08
        width = 0.8
        height = 0.65

        # loop over many groups of pulses/spectra
        for ii in np.arange(0, opts.nsubints, opts.ngroups):

		if ii + opts.ngroups > opts.nsubints: jj = opts.nsubints
		else: jj = ii + opts.ngroups

		print "%d - %d" % (ii, jj-1)

		fig = plt.figure()
		fig.suptitle("%s %s  %s  Cfreq=%.3f MHz  BW=%.3f MHz  %s\nNchan=%d  Nbin=%d  Nsubint=%d  Current subints=%d-%d" % \
			(target, source, telescope, cfreq, bw, opts.legend, nchan, nphbin, opts.nsubints, ii, jj-1), fontsize=opts.fs)
		ax = plt.axes([left, bottom, width/2, height])
		ax.xaxis.tick_bottom()
		ax.xaxis.set_label_position("bottom")
		for label in ax.get_xticklabels(): label.set_fontsize(opts.fs)
		ax.yaxis.set_label_position("left")
		for label in ax.get_yticklabels(): label.set_fontsize(opts.fs)
		plt.ylabel('Pulse number', fontsize=opts.fs)
		plt.xlabel('Frequency (MHz)', fontsize=opts.fs)

		# pulse spectra
		if opts.is_noglobalnorm:
			normdata = normalize(data, 0, nshowchan, ii, jj-ii)
			if opts.is_off_pulse: # calculate OFF-pulse spectrum
				currsplist=np.mean(normdata[ii:jj,:,range(0,binleft)+range(binright,nphbin)], axis=2)
			else: # calculate ON-pulse spectrum
				currsplist=np.mean(normdata[ii:jj,:,binleft:binright], axis=2)
		else: 
			if opts.is_off_pulse: # calculate OFF-pulse spectrum
				currsplist=np.mean(data[ii:jj,:,range(0,binleft)+range(binright,nphbin)], axis=2)
			else: # calculate ON-pulse spectrum
				currsplist=np.mean(data[ii:jj,:,binleft:binright], axis=2)
		if opts.is_noglobalnorm:
			spmin=np.min(currsplist)
			spmax=np.max(currsplist)
		colormap = cm.get_cmap("PRGn_r") # other good ones: BrBG, PuOr, gist_stern, hot, jet, pink, gist_heat, gist_gray, copper, gist_earth
		cax = plt.imshow(currsplist, aspect='auto', interpolation='nearest', origin='lower', cmap=colormap, extent=(opts.fmin, opts.fmax, ii, jj), vmin=spmin, vmax=spmax)
		cbar = fig.colorbar(cax, orientation='horizontal', spacing='proportional')
		cbar.ax.set_xlabel("Spectral density ($\sigma$)", fontsize=opts.fs)
		for label in cbar.ax.get_xticklabels(): label.set_fontsize(opts.fs)

		# average spectrum
		ax3 = plt.axes([left, bottom+height+0.005, width/2, height-0.46])
		ax3.xaxis.set_major_formatter(NullFormatter())
		for label in ax3.get_yticklabels(): label.set_fontsize(opts.fs)
		for label in ax3.get_yticklabels(): label.set_fontsize(opts.fs)

		currsp = np.mean(currsplist, axis=0)
		if opts.is_noglobalnorm:
			ymin=np.min(currsp)
			ymax=np.max(currsp)
		else:
			ymin=np.min([np.min(aversp), np.min(currsp)])
			ymax=np.max([np.max(aversp), np.max(currsp)])

		ymin=0.9*ymin
		ymax=1.1*ymax

		if not opts.is_noglobalnorm:
			ax3.plot(freqs, aversp, color="black", alpha=0.3, linewidth=3, label="%s" % (opts.is_off_pulse and "OFF (total)" or "total"))
		ax3.plot(freqs, currsp, color="green", alpha=1, label="%s" % (opts.is_off_pulse and "OFF" or "current"))
		ax3.set_ylim(ymin=ymin, ymax=ymax)
		ax3.set_xlim(xmin=freqs[0], xmax=freqs[-1])
		plt.ylabel('Average spectrum ($\sigma$)', fontsize=opts.fs)
		prop=fm.FontProperties(size=opts.fs)
		if opts.is_off_pulse: plt.legend(prop=prop, loc=1)

		# pulse profiles	
		ax4 = plt.axes([left+width/2+0.05, bottom, width/2, height])
		ax4.xaxis.tick_bottom()
		ax4.xaxis.set_label_position("bottom")
		for label in ax4.get_xticklabels(): label.set_fontsize(opts.fs)
		ax4.yaxis.set_major_formatter(NullFormatter())
		plt.xlabel('Pulse phase', fontsize=opts.fs)

		if opts.is_noglobalnorm:
			currproflist=profs_norm(np.sum(normdata[ii:jj,:,:], axis=1), jj-ii)
		else: currproflist=profs_norm(np.sum(data[ii:jj,:,:], axis=1), jj-ii)
		currprof = prof_norm(np.sum(currproflist, axis=0))
		currproflist=currproflist[:,profleft:profright]
		currprof = currprof[profleft:profright]
		if opts.is_noglobalnorm:
			profmin=np.min(currproflist)
			profmax=np.max(currproflist)
		colormap = cm.get_cmap("YlOrBr") # other good ones: BrBG, PuOr, gist_stern, hot, jet, pink, gist_heat, gist_gray, copper, gist_earth
		cax = plt.imshow(currproflist, aspect='auto', interpolation='nearest', origin='lower', cmap=colormap, extent=(opts.lview, opts.rview, ii, jj), vmin=profmin, vmax=profmax)
		cbar = fig.colorbar(cax, orientation='horizontal', spacing='proportional')
		cbar.ax.set_xlabel("Flux density ($\sigma$)", fontsize=opts.fs)
		for label in cbar.ax.get_xticklabels(): label.set_fontsize(opts.fs)
	
		# average profile
		ax2 = plt.axes([left+width/2+0.05, bottom+height+0.005, width/2, height-0.461])
		ax2.xaxis.set_major_formatter(NullFormatter())
		ax2.yaxis.set_label_position("right")
		for label in ax2.get_yticklabels(): label.set_fontsize(opts.fs)
		for label in ax2.get_yticklabels(): label.set_fontsize(opts.fs)

		if opts.is_noglobalnorm:
			ymin=np.min(currprof)
			ymax=np.max(currprof)
		else:
			ymin=np.min([np.min(averprof), np.min(currprof)])
			ymax=np.max([np.max(averprof), np.max(currprof)])

		ymin=0.9*ymin
		ymax=1.1*ymax

		if not opts.is_noglobalnorm:
			ax2.plot(phase, averprof, color="black", alpha=0.3, linewidth=3, label="total")
		ax2.plot(phase, currprof, color="green", alpha=1, label="current")
		ax2.set_ylim(ymin=ymin, ymax=ymax)
		ax2.set_xlim(xmin=phase[0], xmax=phase[-1])
		plt.ylabel('Average profile ($\sigma$)', fontsize=opts.fs)
		prop=fm.FontProperties(size=opts.fs)
		plt.legend(prop=prop, loc=0)

		pngname="%s-%05d-%05d.png" % (opts.outstem, ii, jj-1)
		if opts.is_saveonly:
			plt.savefig(pngname)
		else:
			if opts.is_save: plt.savefig(pngname)
			plt.show()
