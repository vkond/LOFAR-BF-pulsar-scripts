#!/usr/bin/env python
#
import numpy as np
import psrchive as pc
import os, os.path, stat, glob, sys, getopt, re
import matplotlib.pyplot as plt
from matplotlib.ticker import *
import matplotlib.cm as cm
import scipy.signal as ss
import scipy.stats as sc
import cPickle
import math

fs=10 # fontsize
fscr_factor=1
tscr_factor=1
bscr_factor=1
is_acf=False
is_quick=False
is_load=False
dynspfile="dynsp"
tmin=0
tmax=0
fmin=0
fmax=0
is_secsp=False
is_off_pulse=False

# Help
def usage (prg):
        """ prints the usage info about the current program
        """
        print "Program %s makes dynamic spectrum or 2D ACF of it\n" % (prg,)
        print "Usage: %s [options] <PSRCHIVE-style fits-file>\n\
         --fscrunch N           - scrunching in frequency by a factor of N, default = 1\n\
         --tscrunch N           - scrunching in time by a factor of N, default = 1\n\
         --bscrunch N           - scrunching in phase by a factor of N, default = 1\n\
         --acf                  - to calculate 2D ACF of dynamic spectrum\n\
         --secsp                - to calculate secondary spectrum (not implemented)\n\
         --tmin tmin            - minimum subintegration to use\n\
         --tmax tmax            - maximum subintegration to use\n\
         --fmin fmin            - minimum frequency channel to use\n\
         --fmax fmax            - maximum frequency channel to use\n\
         -q, --quick            - quick normalizing, mean and rms for each channel\n\
                                  are calculated on tscrunched data\n\
         --off-pulse            - using only OFF-pulse region for mean/rms calculation\n\
         -l, --load             - load file with dynamic spectrum\n\
         --file DYNSPECTRUM     - name of the file with dynamic spectrum to load or save\n\
                                  default: 'dynsp'\n\
         -h, --help             - print this message\n" % (prg,)

# Parsing the command line
def parsecmdline (prg, argv):
        """ parse the command line arguments
        """
        if not argv:
                usage (prg)
                sys.exit()
        else:
                try:
                    opts, args = getopt.getopt (argv, "hql", ["help", "tscrunch=", "fscrunch=", "bscrunch=", "acf", "off-pulse", "quick", "load", "secsp", "tmin=", "tmax=", "fmin=", "fmax=", "file="])
                    for opt, arg in opts:

                                if opt in ("-h", "--help"):
                                        usage (prg)
                                        sys.exit()

                                if opt in ("--tscrunch"):
					global tscr_factor
					tscr_factor = int(arg)

                                if opt in ("--fscrunch"):
					global fscr_factor
					fscr_factor = int(arg)

                                if opt in ("--bscrunch"):
					global bscr_factor
					bscr_factor = int(arg)

                                if opt in ("--acf"):
					global is_acf
					is_acf = True

                                if opt in ("--off-pulse"):
					global is_off_pulse
					is_off_pulse = True

                                if opt in ("--secsp"):
					global is_secsp
					is_secsp = True

                                if opt in ("-q", "--quick"):
					global is_quick
					is_quick = True

                                if opt in ("-l", "--load"):
					global is_load
					is_load = True

                                if opt in ("--file"):
					global dynspfile
					dynspfile = arg

                                if opt in ("--tmin"):
					global tmin
					tmin = int(arg)

                                if opt in ("--tmax"):
					global tmax
					tmax = int(arg)

                                if opt in ("--fmin"):
					global fmin
					fmin = int(arg)

                                if opt in ("--fmax"):
					global fmax
					fmax = int(arg)

                    if not args:
                                print ".ar file is not given!\n"
                                sys.exit(2)
                    else:
                                global fitsfile
                                fitsfile = args[0]

                except getopt.GetoptError:
                        print "Wrong option!"
                        usage (prg)
                        sys.exit(2)


# Main
if __name__=="__main__":
        parsecmdline (sys.argv[0].split("/")[-1], sys.argv[1:])

	if is_load:
        	dp = open(dynspfile, "rb")
        	(sp, tscr_factor, fscr_factor, bscr_factor)=cPickle.load(dp)
        	dp.close()

	print fitsfile
	print "Tscrunch factor:", tscr_factor
	print "Fscrunch factor:", fscr_factor
	print "Bscrunch factor:", bscr_factor
	if is_acf:
		print "ACF"
	if is_load: print "Loaded ", dynspfile

       	raw = pc.Archive_load(fitsfile)
	if not(raw.get_dedispersed()):
		raw.dedisperse()
	raw.pscrunch()
	if fscr_factor > 1: raw.fscrunch(fscr_factor)
	if tscr_factor > 1: raw.tscrunch(tscr_factor)
	if bscr_factor > 1: raw.bscrunch(bscr_factor)
	nphbin = raw.get_nbin()
	noff = int(0.1 * nphbin)  # size of 10% window (for half of the OFF-pulse window)
	nchan = raw.get_nchan()
	cfreq = raw.get_centre_frequency()
	bw = raw.get_bandwidth()
	nsubint = raw.get_nsubint()
	chan_bw = bw / nchan
	tobs = raw.integration_length()
	tsubint = tobs / nsubint
	psrname = raw.get_source()

	if tmin < 0 or tmin >= nsubint: tmin = 0
	if tmax <= 0 or tmax>nsubint or tmax <= tmin: tmax = nsubint
	if fmin < 0 or fmin >= nchan: fmin = 0
	if fmax <= 0 or fmax>nchan or fmax <= fmin: fmax = nchan
	twin = tmax - tmin
	fwin = fmax - fmin

	if not is_load:
		r = raw.get_data()
       		data = r[:,0,:,:]
	        weights = raw.get_weights()
		for s in np.arange(nsubint): data[s,(weights[s]==0)] = 0.0

		#time stokes f phase
		sp=np.zeros((nsubint, nchan), dtype=float)

		print "Normalizing..."
		if not is_quick:
			for s in range(nsubint):
				for i in range(nchan):
					if is_off_pulse:
						binmin1=np.argmin(data[s,i,:nphbin/2])
						binmin2=np.argmin(data[s,i,nphbin/2:])
						mean = np.median(data[s,i,range(binmin1-noff/2,binmin1+noff/2)+range(binmin2-noff/2,binmin2+noff/2)])
						rms = np.std(data[s,i,range(binmin1-noff/2,binmin1+noff/2)+range(binmin2-noff/2,binmin2+noff/2)])
					else: 
						osm, osr = sc.probplot(data[s,i], sparams=(), dist='norm', fit=0)
						q_max = np.min(np.where(osm > 1.0))
						q_min = np.max(np.where(osm < -1.0))
						rms, mean = np.polyfit(osm[q_min:q_max], osr[q_min:q_max], 1)
					data[s,i] -= mean
					if rms == 0.0: data[s,i] = 0.0
					else: data[s,i] /= rms
				        crit=np.isfinite(data[s,i])
				        data[s,i][-crit] = 0.0
					sp[s,i] = np.max(data[s,i])
		else:
			scr = np.sum(data, axis=0)	
			scr /= nsubint
			for i in range(nchan):
				if is_off_pulse:
					binmin1=np.argmin(scr[i,:nphbin/2])
					binmin2=np.argmin(scr[i,nphbin/2:])	
					mean = np.mean(scr[i,range(binmin1-noff/2,binmin1+noff/2)+range(binmin2-noff/2,binmin2+noff/2)])
					rms = np.std(scr[i,range(binmin1-noff/2,binmin1+noff/2)+range(binmin2-noff/2,binmin2+noff/2)])
				else:
					osm, osr = sc.probplot(scr[i], sparams=(), dist='norm', fit=0)
					q_max = np.min(np.where(osm > 1.0))
					q_min = np.max(np.where(osm < -1.0))
					rms, mean = np.polyfit(osm[q_min:q_max], osr[q_min:q_max], 1)
				for s in range(nsubint): 
					data[s,i] -= mean
					if rms == 0.0: data[s,i] = 0.0
					else: data[s,i] /= rms
				        crit=np.isfinite(data[s,i])
				        data[s,i][-crit] = 0.0
					sp[s,i] = np.max(data[s,i])

		# saving dynamic spectrum to binary file
		dp = open (dynspfile, "wb")
		record = (sp, tscr_factor, fscr_factor, bscr_factor)
		cPickle.dump(record, dp, True)
		dp.close()

        if is_acf:
		if twin%2 != 0:
			tmax -= 1
			twin -= 1
		if fwin%2 != 0:
			fmax -= 1
			fwin -= 1
                def energy (x): return np.mean(x.ravel()**2)
                acf = ss.correlate2d(sp[tmin:tmax,fmin:fmax], sp[tmin:tmax,fmin:fmax], mode='same', boundary='fill', fillvalue = 0)
                acf /= np.sqrt (energy(sp[tmin:tmax,fmin:fmax]) * energy(sp[tmin:tmax,fmin:fmax]))
                for s in np.arange(twin/2):
                        for i in np.arange(fwin/2):
                                acf[s,i] /= (4*(s+1+twin/2)*(i+1+fwin/2))
                                acf[s,i+fwin/2] /= (4*(s+1+twin/2)*(fwin-1-i))
                                acf[s+twin/2,i] /= (4*(twin-1-s)*(i+1+fwin/2))
                                acf[s+twin/2,i+fwin/2] /= (4*(twin-1-s)*(fwin-1-i))
                (tm, fm) = np.unravel_index(np.argmax(acf), np.shape(acf))
                tcut=acf[tm,:]
                fcut=acf[:,fm]


	# Plotting
	fig = plt.figure()
	if not is_acf:
		str = "Pulsar: %s    File: %s    Fc: %.3f MHz    BW: %.3f MHz\nNchan: %d    chan_bw: %.3f MHz    Nsubint: %d    T_subint: %.3f s" % (psrname, fitsfile.split("/")[-1], cfreq, bw, nchan, chan_bw, nsubint, tsubint)
#		ax = plt.axes([0.08, 0.2, 0.8, 0.8])
		ax = plt.subplot(111)
		ax.set_title(str, fontsize=fs, ha='left', x=0)
		ax.xaxis.tick_bottom()
		ax.xaxis.set_label_position("bottom")
		for label in ax.get_xticklabels(): label.set_fontsize(fs)
		ax.yaxis.tick_left()
		ax.yaxis.set_label_position("left")
		for label in ax.get_yticklabels(): label.set_fontsize(fs)
		plt.xlabel('Frequency (MHz)', fontsize=fs)
		plt.ylabel('Time (s)', fontsize=fs)
	
                colormap = cm.get_cmap("YlGnBu")
		cax = plt.imshow(sp[tmin:tmax,fmin:fmax], aspect='auto', interpolation='nearest', origin='lower', cmap=colormap, extent=(cfreq-bw/2.+fmin*chan_bw, cfreq-bw/2.+fmax*chan_bw, tmin*tsubint, tmax*tsubint))
		cbar = fig.colorbar(cax, orientation='vertical', spacing='proportional')
		cbar.ax.set_ylabel("Flux density ($\sigma$)", fontsize=fs)
		for label in cbar.ax.get_xticklabels(): label.set_fontsize(fs)
		
	else: # acf

		# dimensions of the central plot
		left = 0.30
		bottom = 0.08
		width = 0.39
		height = 0.55
		
		# acf
		str = "Pulsar: %s    File: %s    Fc: %.3f MHz    BW: %.3f MHz\nNchan: %d    chan_bw: %.3f MHz    Nsubint: %d    T_subint: %.3f s" % (psrname, fitsfile.split("/")[-1], cfreq, bw, nchan, chan_bw, nsubint, tsubint)
		ax = plt.axes([left, bottom, width, height])
		for label in ax.get_yticklabels(): label.set_fontsize(fs)
		plt.xlabel('Frequency lag (MHz)', fontsize=fs)
		plt.ylabel('Time lag (s)', fontsize=fs)
		for label in ax.get_xticklabels(): label.set_fontsize(fs)
	
		colormap = cm.get_cmap("YlGnBu")
		cax = plt.imshow(acf, aspect='auto', interpolation='nearest', origin='lower', cmap=colormap, extent=(-fwin*chan_bw/2.+chan_bw, fwin*chan_bw/2., -0.5*twin*tsubint, 0.5*twin*tsubint))
		cbar = fig.colorbar(cax, orientation='horizontal', spacing='proportional')
		cbar.ax.set_xlabel("Correlation coefficient", fontsize=fs)
		for label in cbar.ax.get_xticklabels(): label.set_fontsize(fs)
                axl = plt.twinx()
                axl.yaxis.tick_right()
                axl.yaxis.set_label_position("right")
                for label in axl.get_yticklabels(): 
			label.set_rotation(90)
			label.set_fontsize(fs)
                axl.set_ylim(ymin=-0.5*twin*tsubint,ymax=0.5*twin*tsubint)
                plt.ylabel('Time lag (s)', fontsize=fs)

		# tcut
		ax3 = plt.axes([left, bottom+height, width, height-0.3])
		ax3.set_title(str, fontsize=fs, ha='left', x=0, y=1.2)
		ax3.xaxis.set_major_formatter(NullFormatter())
		for label in ax3.get_yticklabels(): label.set_fontsize(fs)
		plt.ylabel('ACF', fontsize=fs)
		freqlag = np.arange(fwin)
		freqlag = freqlag - fwin/2. + 1.
		freqlag *= chan_bw
		ax3.plot(freqlag, tcut)
		ax3.grid()
		halfpower = np.median(tcut)
		halfpower = halfpower + (np.max(tcut) - halfpower) / 2.
		ax3.plot(freqlag, halfpower * np.ones(fwin, dtype=float),'g--')
		ax3.set_xlim(xmin=-fwin*chan_bw/2.+chan_bw,xmax=fwin*chan_bw/2.)
		ax3.set_ylim(ymin=np.min(tcut), ymax=np.max(tcut))
		axb3 = plt.twiny()
		axb3.xaxis.tick_top()
		axb3.xaxis.set_label_position("top")
		axb3.set_xlim(xmin=-fwin*chan_bw/2.+chan_bw,xmax=fwin*chan_bw/2.)
		for label in axb3.get_xticklabels(): label.set_fontsize(fs)

		# fcut
		ax4 = plt.axes([left-0.14, bottom+0.165, width-0.25, width-0.005])
		for label in ax4.get_xticklabels(): label.set_fontsize(fs)
		plt.xlabel('ACF', fontsize=fs)
		timelag = np.arange(twin)
		timelag = timelag - twin/2. + 1.
		timelag *= tsubint
		ax4.plot(fcut, timelag)
		ax4.grid()
		halfpower = np.min(fcut) + (np.max(fcut) - np.min(fcut)) / math.e 
		halfline=np.ones(twin, dtype=float)
		ax4.plot(halfpower * np.ones(twin, dtype=float), timelag, 'g--')
		ax4.set_ylim(ymin=-0.5*twin*tsubint,ymax=0.5*twin*tsubint)
		ax4.set_xlim(xmin=np.min(fcut), xmax=np.max(fcut))
		plt.ylabel('Time lag (s)', fontsize=fs)
		for label in ax4.get_yticklabels():
        		label.set_rotation(90)
        		label.set_fontsize(fs)
		ax4.invert_xaxis()
		axb4 = plt.twinx()
		axb4.yaxis.set_label_position("right")
		axb4.yaxis.set_major_locator(NullLocator())

	plt.show()
