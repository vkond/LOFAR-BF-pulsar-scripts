#!/usr/bin/env python
#
# plots the dynamic spectrum for a given *.sub files. In interactive mode, 
# allows to manually select channels corrupted by RFI, view time series and 
# samples' histogram of selected channels.
#
# Vlad Kondratiev (c)
#
import numpy as np
import array as ar
import os, os.path, stat, glob, sys, getopt
import infodata as inf
import math
import scipy.stats as sc

is_saveonly = False      # if True, script will save the dynamic spectrum in png file
is_excludeonly = False   # if True, only completely bad subbands will be excluded
is_psrfits = False       # if True, the input file is the single fits-file rather than *.sub
threshold = 6 # threhold in sigmas to clip RFIs
rfilimit = 10  # (in percents) if more, whole subband will be excluded
samples2show = 0 #  size of window to show in seconds (if 0, show the whole file(s))
samples_offset = 0 # offset from the beginning to skip (if 0 , show the whole file(s))
Nbins = 7630 # number of bins to average (usually corresponds to ~10s)
histbins = 100 # number of bins in the histogram
xlow=-32768 # lowest sample value (for signed short)
xhigh=32767 # hoghest sample value
samplesize = 2 # 2 bytes, because we have 16-bit data
fs=10 # fontsize
is_statistics = False # if True prints statistics about large positive/negative samples
# number of occasions when number of positive large peaks is larger than number of negative ones
numposlarger=0
# number of occasions when number of positive large peaks is smaller than number of negative ones
numpossmaller=0
user_title=""  # Extra user info to plot on the title


def usage (prg):
        """ prints the usage info about the current program
        """
        print "Program %s plots the dynamic spectrum based on raw subbands data\n" % (prg,)
        print "Usage: %s [-n, --nbins <value>] [-t, --threshold <value>] [-w, --win <value>]\n\
			[-l, --offset <value>] [--rfilimit <value>] [--excludeonly] [--statistics]\n\
                        [--saveonly] [--psrfits] [--title <string>] [-h, --help] <*sub|*.fits>\n\
         -n, --nbins <value>     - number of samples to average (default: 7630)\n\
	 -w, --win   <value>     - size of the window (in seconds) to show\n\
	 -l, --offset  <value>   - offset from the beginning (in seconds)\n\
         -t, --threshold <value> - threshold (in sigma) to clip RFI (default: 6)\n\
	 --excludeonly           - only exclude completely junk subbands\n\
	 --rfilimit <value>      - percent of RFI per subband allowed not to exclude whole subband (default: 10)\n\
         --statistics            - print info about fraction of large positive and negative samples\n\
	 --saveonly              - only saves png file and exits\n\
         --psrfits               - input data is the single fits-file instead of *.sub\n\
         --title <string>        - set user info to show on the title of the plot (e.g. ObsID, station, etc.)\n\
         -h, --help              - print this message\n" % (prg,)


def parsecmdline (prg, argv):
        """ parse the command line arguments
        """
        if not argv:
                usage (prg)
                sys.exit()
        else:
                try:
                        opts, args = getopt.getopt (argv, "hn:t:w:l:", ["help", "nbins=", "threshold=", "rfilimit=", "excludeonly", "saveonly", "win=", "offset=", "statistics", "psrfits", "title="])
                        for opt, arg in opts:
                                if opt in ("-h", "--help"):
                                        usage (prg)
                                        sys.exit()

                                if opt in ("-n", "--nbins"):
                                        global Nbins
                                        Nbins = long(arg)

                                if opt in ("-t", "--threshold"):
                                        global threshold
                                        threshold = float(arg)

                                if opt in ("--rfilimit"):
                                        global rfilimit
                                        rfilimit = float(arg)

                                if opt in ("--excludeonly"):
                                        global is_excludeonly
					is_excludeonly = True

                                if opt in ("--saveonly"):
					global is_saveonly
					is_saveonly = True

				if opt in ("-w", "--win"):
					global samples2show
					samples2show = float(arg)

				if opt in ("-l", "--offset"):
					global samples_offset
					samples_offset = float(arg)

                                if opt in ("--statistics"):
					global is_statistics
					is_statistics = True

                                if opt in ("--psrfits"):
					global is_psrfits
					is_psrfits = True
					# we change xlow and xhigh because PSRFITS is 8-bit data
					global xlow, xhigh
					xlow = -256
					xhigh = 255

                                if opt in ("--title"):
                                        global user_title
                                        user_title = arg

                        if not args:
                                print "No subband files!\n"
                                usage (prg)
                                sys.exit(2)
                        else:
                                global subfiles
				subfiles = args
				subfiles.sort()

                except getopt.GetoptError:
                        print "Wrong option!"
                        usage (prg)
                        sys.exit(2)
	

def setup_plot(arr, title, colormap):
	""" initializing dynamic spectrum plot
	"""
	global cax, cbar
	plt.clf()
	ax = fig.add_subplot (111)
	ax.set_zorder(0.1)   # to make moise y-cursor to get values from this 'Subband' axis, rather then from axr 'Freq' axis
	cax = ax.imshow(arr, norm=colors.Normalize(arr.min(), arr.max()), interpolation='nearest', aspect='auto', origin='lower', cmap=colormap)
	cbar = fig.colorbar (cax, orientation='horizontal', spacing='uniform', pad=0.1)
	fig.suptitle ("%s %s" % (user_title, title), fontsize=fs, y=0.94)
        def printsub (val, pos=None): return '%d' % (subband_offset + val)
	ax.yaxis.set_major_formatter(ticker.FuncFormatter(printsub))
	for label in ax.get_yticklabels(): label.set_fontsize(fs)
	plt.ylabel ("Channels", fontsize=fs)
	plt.xlabel ("Time (s)", fontsize=fs)

	axr = plt.twinx()
	axr.yaxis.tick_right()
	axr.yaxis.set_label_position("right")
	axr.set_ylim(ymin=lofreq-chanbw/2., ymax=lofreq+totalbw-chanbw/2.)
	for label in axr.get_yticklabels(): label.set_fontsize(fs)
	plt.ylabel("Frequency (MHz)", fontsize=fs, rotation=-90)

	def printtime (val, pos=None): return '%1.0f'%(float(samples_offset)*tsamp + float(val) * tsamp * Nbins)
	ax.xaxis.set_major_formatter(ticker.FuncFormatter(printtime))
	for label in ax.get_xticklabels(): label.set_fontsize(fs)

	plt.draw()

def plot_update(arr):
	""" updating dynamic spectrum
	"""
	cax.set_data(arr)
	cbar.set_array(arr)
	cbar.autoscale()
	cax.changed()
	plt.draw()

def press(event):
        """ event handler for pressed keys in the main window
        """
        if event.key == 'm' or event.key == 'M': setup_plot (mask, "Clipping mask", cm.Reds)
        if event.key == 'c' or event.key == 'C': setup_plot (clipped, "Cleaned dynamic spectrum", cm.jet)
        if event.key == 'd' or event.key == 'D': setup_plot (spectrum, "Original dynamic spectrum", cm.jet)
	if event.key == ' ':
		subband_to_exclude = int(event.ydata + 0.5)
		badbands.append(subband_to_exclude)
		print "subband %d will be excluded" % (subband_to_exclude + subband_offset)
		clipped[subband_to_exclude] = np.zeros(int(size/Nbins))
	if event.key == 'enter': 
		print "Clipping RFI from subbands files ...", 
		sys.stdout.flush()
		clipsubband()
		print "done. Original files moved to 'backup' directory"
	if event.key == 'q' or event.key == 'Q': plt.close()
	if event.key == '?': plot_help()
	if event.key == 'b' or event.key == 'B':
		global selected_subband
		selected_subband = int(event.ydata + 0.5)
		plot_help_subband()
		setup_plottime(spectrum, selected_subband, "Time series with mask ovelapped", overlap=True)
	plt.show()

def plot_help():
        """ prints the info about keys available for usage
        """
        print
        print "Press following keys"
        print " m or M - to plot the clipping mask"
        print " d or D - to plot the dynamic spectrum"
	print " c or C - to plot the clipped dynamic spectrum"
	print " b or B - to plot time series of the selected subband in a separate window"
	print " enter  - clip RFI from subbands files"
	print " space  - select subband to be excluded"
	print "   ?    - print this help"
        print " q or Q - quit"
        print

def clipsubband():
	""" move subband files to backup directory and create new clipped subbands
	"""
	path="/".join(subfiles[0].split("/")[0:-1])
	if path == "": path = "./.orig"
	else: path += "/" + ".orig"
	if os.path.exists(path):  # change the name of backup directory to .origN (if .orig exists)
		origpath = path
		origcounter = 0
		while os.path.exists(path):
			origcounter += 1
			path = origpath + str(origcounter)
	cmd = "mkdir -p " + path
	os.system (cmd)
	# excluding first very bad subbands
	for s in np.unique(badbands):
		# move original subband file to backup (.orig) directory
		cmd = "cp " + subfiles[s] + " " + path
		os.system (cmd)
		# create a zero-sized subband file
		cmd = "touch " + subfiles[s]
		os.system (cmd)

	if is_excludeonly == False: # correcting other subbands for RFI
		indarr = np.arange (0, nfiles, 1)
		for i in np.unique(badbands): indarr.remove(i)
		for i in indarr:
			f = open(subfiles[i], "rb")
			data = ar.array('h')  # 'h' - for signed short
			data.read(f, sizes[i])
			f.close()
			cmd = "cp " + subfiles[i] + " " + path
			os.system (cmd)
			ndata = np.array(data)
			# clipping
			for elem in clipindices[i]: ndata[elem*Nbins:(elem+1)*Nbins] = 0
			# writing the new subband file
			f = open(subfiles[i], "wb")
			out = ar.array('h')
			out.fromlist(ndata.tolist())
			out.tofile(f)
			f.close()

def press_in_band(event):
        """ event handler for pressed keys in the subband window
        """
        if event.key == 'q' or event.key == 'Q': plt.close()
	if event.key == '?': plot_help_subband()
	if event.key == 'm' or event.key == 'M': plottime(mask, selected_subband, "Clipping mask", clr="red", overlap=False)
	if event.key == 'c' or event.key == 'C': plottime(clipped, selected_subband, "Flagged time series", overlap=False)
	if event.key == 'b' or event.key == 'B': plottime(spectrum, selected_subband, "Time series with mask ovelapped", overlap=True)
	if event.key == 't' or event.key == 'T': plottime(spectrum, selected_subband, "Original time series", overlap=False)
	if event.key == 'h' or event.key == 'H': plothist(spectrum, selected_subband, "Samples histogram")

def plot_help_subband():
        """ prints the info about keys available for usage in a subband window
        """
        print
        print "Following keys are available in the separate subband window:"
        print " m or M - to plot the clipping mask"
	print " c or C - to plot flagged time series"
	print " b or B - to plot original time series with the mask overlapped"
	print " t or T - to plot original time series"
	print " h or H - to plot samples histogram"
	print "   ?    - print this help"
        print " q or Q - quit"
        print

def setup_plottime(series, selband, title, clr="blue", overlap=False):
        """ setup the separate window with the time series for selected subband
        """
	global figt
        figt = plt.figure()
        figt.canvas.mpl_connect('key_press_event', press_in_band)
	plottime (series, selband, title, clr, overlap)

def plottime(series, selband, title, clr="blue", overlap=False):
	""" plots the separate window with the time series for selected subband
	"""
	plt.clf()
	axt = figt.add_subplot(111)
	def printtime (x, pos=None): return '%1.0f'%(float(samples_offset)*tsamp + float(x) * tsamp * Nbins)
	plt.xlabel("Time (s)", fontsize=fs)
	plt.ylabel("Flux density (arb. units)", fontsize=fs)

	axt.plot (range(0,int(size/Nbins)), series[selband], color=clr)
	if overlap: 
		axt.plot (range(0,int(size/Nbins)), mask[selband], color='red')
	axt.set_xlim(xmin=0, xmax=int(size/Nbins))
	axt.xaxis.set_major_formatter(ticker.FuncFormatter(printtime))
	axt.yaxis.set_major_locator(ticker.MaxNLocator(4))
	axt.annotate(subfiles[selband].split("/")[-1].split(".")[1] + ",   " + title, xy=(0,0), xytext=(0.0, 1.02), xycoords='axes fraction', fontsize=fs)
	for label in axt.get_xticklabels(): label.set_fontsize(fs)
	for label in axt.get_yticklabels():
        	label.set_fontsize(fs)

	plt.show()

def plothist (series, selband, title):
	""" plots the histogram of the chosen subband
	"""
	plt.clf()
	axt = figt.add_subplot(111)
	def percent (x, pos=None): return '%1.2f'%(float(x) * 100./int(size/Nbins))
	plt.xlabel("Sample Value", fontsize=fs)
	plt.ylabel("Fraction (%)", fontsize=fs)

	axt.hist (series[selband], histbins, range=(xlow, xhigh))
	axt.set_xlim(xmin=xlow, xmax=xhigh)
	axt.yaxis.set_major_formatter(ticker.FuncFormatter(percent))
	axt.xaxis.set_major_locator(ticker.AutoLocator())
	axt.annotate(subfiles[selband].split("/")[-1].split(".")[1] + ",   " + title, xy=(0,0), xytext=(0.0, 1.02), xycoords='axes fraction', fontsize=fs)
	for label in axt.get_xticklabels(): label.set_fontsize(fs)
	for label in axt.get_yticklabels():
	        label.set_fontsize(fs)

	plt.show()

# ============================== M A I N =============================================
if __name__=="__main__":
        parsecmdline (sys.argv[0].split("/")[-1], sys.argv[1:])
	if is_saveonly:
		import matplotlib
		matplotlib.use("Agg")
	else:
		import matplotlib

	import matplotlib.pyplot as plt
	import matplotlib.ticker as ticker
	import matplotlib.colors as colors
	import matplotlib.cm as cm

	if is_psrfits:
		import pyfits as py


	# if input files are *.sub
	if not is_psrfits:
	
		nfiles=len(subfiles)
		# array of file sizes
		sizes = [os.stat(file)[stat.ST_SIZE] / samplesize for file in subfiles]
		# maximum size
		size=max(sizes)

        	# reading inf-file to get corresponding info
        	inffile = subfiles[0].split(".sub")[0] + ".sub.inf"
        	id = inf.infodata(inffile)
        	lofreq = id.lofreq
        	chanbw = id.chan_width
		totalbw = id.BW
		tsamp = id.dt

		# handle offset from the beginning
		if samples_offset > 0:
			samples_offset = int (samples_offset / tsamp)
			if samples_offset > size - 2:
				samples_offset = 0
		# will show only samples2show number of bins (if chosen)
		if samples2show > 0:
			samples2show = int(samples2show / tsamp)	
			if size-samples_offset > samples2show: 
				size = samples2show
				sizes = [value-samples_offset > samples2show and samples2show or value for value in sizes]
			else:
				if samples_offset > 0:
					sizes = [value - samples_offset for value in sizes]

		# first subband number
		subband_offset = int(subfiles[0].split(".sub")[-1])

		if is_saveonly:
			pngname = subfiles[0].split(".sub")[0] + ".sub" + str(subband_offset) + "-" + str(subband_offset+nfiles-1) + ".png"

		rfirepname = subfiles[0].split(".sub")[0] + ".sub" + str(subband_offset) + "-" + str(subband_offset+nfiles-1) + ".rfirep"

		isize=int(size/Nbins)

		# forming the array for having the dynamic spectrum
		if isize == 0:
			print "Number of bins %d is larger than given size %d!" % (Nbins, size) 
			sys.exit(1)	
		spectrum=np.zeros((nfiles, isize))

		# forming a mask file with samples to reject
		mask=np.zeros((nfiles, isize))
		clipped=np.zeros((nfiles, isize))  # clipped data (only for plotting)
		mean=np.zeros((nfiles, isize))     # 2D array of means
		rms=np.zeros((nfiles, isize))      # 2D array of rms's
		levels=np.zeros((nfiles, isize))   # 2D array of levels (samples in sigma)

		# dictionary with indices of RFI'ish samples to be clipped
		clipindices = {}  # for not _very_ bad subbands
		badbands = []     # list of completely bad subbands

		if is_saveonly == False:
			plt.ion() # interactive plotting
			fig=plt.figure()
			fig.canvas.mpl_connect('key_press_event', press)
			setup_plot(spectrum, "Original dynamic spectrum", cm.jet)
			plot_help()

		for i in np.arange(0, nfiles, 1):
			f = open(subfiles[i], "rb")
			data = ar.array('h')  # 'h' - for signed short
			f.seek (samples_offset * samplesize)  # position to the first sample to read
			data.read(f, sizes[i])
			f.close()
			ndata = np.array(data)
			if sizes[i] == 0:
				rfi_fraction = 100.
				clipped[i] = np.zeros(isize)
				print "subband %d will be excluded (blanked)" % (i+subband_offset,)
				continue

			spectrum[i] = [np.average(ndata[k*Nbins:(k+1)*Nbins]) for k in np.arange(0, isize, 1)]
			# calculate mean and rms in the windows of Nbins size
			# to exclude outliers we sort the values in Nbins interval first and then use only first half of it
			mean[i] = [np.mean(np.sort(ndata[k*Nbins:(k+1)*Nbins])[0:Nbins/2]) for k in np.arange(0, isize, 1)]
			rms[i] = [np.std(np.sort(ndata[k*Nbins:(k+1)*Nbins])[0:Nbins/2]) for k in np.arange(0, isize, 1)]
			# getting statistics about fraction of positive/negative samples
			if is_statistics == True:
				condition=np.zeros(int(size/Nbins), dtype=bool)
				lev = np.array([(spectrum[i][k]-mean[i][k])/rms[i][k] for k in np.arange(0, isize, 1)])
				levsize=np.size(lev)
				condition = condition | (lev > threshold)
				pospeak = np.size(lev.compress(condition))
				condition=np.zeros(isize, dtype=bool)
				condition = condition | (lev < -threshold)
				negpeak = np.size(lev.compress(condition))
				if pospeak > negpeak:
					numposlarger += 1
				if pospeak < negpeak:
					numpossmaller += 1

				pospeak = float((pospeak * 100.)/levsize)
				negpeak = float((negpeak * 100.)/levsize)

			# check if all values in rms are zeros. If so, then exclude this subband
			if np.size(np.trim_zeros(np.sort(rms[i]), 'bf')) != 0:
#				sig = np.trim_zeros(np.sort(rms[i]), 'bf')[0]
				# I should remove np.abs at some point, because it's not really correct. I am using it here
				# to avoid really huge negative spikes
#				av = np.trim_zeros(np.sort(np.abs(mean[i])), 'bf')[0]

#				levels[i] = [np.abs((float(spectrum[i][k] - av))/sig) for k in np.arange(0, isize, 1)]
				levels[i] = [np.abs((spectrum[i][k] - mean[i][k])/rms[i][k]) for k in np.arange(0, isize, 1)]
				mask[i] = [(levels[i][k] > threshold and spectrum[i][k] or 0) for k in np.arange(0, isize, 1)]
				clipped[i] = [(spectrum[i][k] - mask[i][k]) for k in np.arange(0, isize, 1)]
				# are there many zeros?
				condition=np.zeros(isize, dtype=bool)
				condition=condition | (levels[i] > threshold)
				rfi_fraction = (float(np.size(clipped[i].compress(condition)))/np.size(clipped[i]))*100.
			else:
				rfi_fraction = 100.

			if rfi_fraction >= rfilimit:  # bad subband  
				clipped[i] = np.zeros(isize)
				badbands.append(i)
				print "subband %d will be excluded (rfi fraction = %.2f%%)" % (i+subband_offset, rfi_fraction)
			else:
				clipindices[i] = np.where(levels[i] > threshold)[0]

			if is_saveonly == False: plot_update(spectrum)
		
	# input file is single fits-file
	else:
		subband_offset = 0

		fitsfile = subfiles[0]
		hdu=py.open(fitsfile, 'readonly', memmap=1)
		cfreq=hdu[0].header['obsfreq']
		nchan=hdu[1].header['nchan']
		chanbw=hdu[1].header['chan_bw']
		totalbw=nchan*chanbw
		lofreq=cfreq-0.5*totalbw+0.5*chanbw   # central freq of lowest channel
		tsamp=hdu[1].header['tbin']
		nsize=np.size(hdu[1].data[0]['data'])
		nrows=hdu[1].header['naxis2']
		size=(nsize/nchan) * nrows

                # handle offset from the beginning
                if samples_offset > 0:
                        samples_offset = int (samples_offset / tsamp)
                        if samples_offset > size - 2:
                                samples_offset = 0
		row_offset=int(samples_offset/(nsize/nchan) + 0.5)
		samples_offset=row_offset*(nsize/nchan)
                # will show only samples2show number of bins (if chosen)
                if samples2show > 0:
                        samples2show = int(samples2show / tsamp)
			print samples2show
                        if size-samples_offset > samples2show:
                                nrows = int (samples2show/(nsize/nchan))
				size=(nsize/nchan) * nrows

		# updating Nbins 
		nrows_step=int(Nbins/(nsize/nchan) + 0.5)
		Nbins = nrows_step * (nsize/nchan)
		if Nbins == 0:
			Nbins = 1
			nrows_step = 1
		isize=int(nrows/nrows_step)

		print "CFreq: ", cfreq, " MHz   BW: ", totalbw, " MHz   Flow: ", lofreq, " MHz"
		print "Nchan: ", nchan, "   Channel bw: ", chanbw, " MHz   Sampling time: ", tsamp, " s"
		print "Nrows: ", nrows, "   Number of samples: ", size, "   Averaged: ", isize, " (by %d bins = %.2f s" % (Nbins, Nbins * tsamp), ")"
		print "Row offset: ", row_offset, "  Rows to use: ", nrows


                if is_saveonly:
			pngname = fitsfile.split(".fits")[0] + ".sp.png"

		rfirepname = fitsfile.split(".fits")[0] + ".rfirep"

		# forming the array for having the dynamic spectrum
		if isize == 0:
			print "Number of bins %d is larger than given size %d!" % (Nbins, size) 
			sys.exit(1)	
		spectrum=np.zeros((nchan, isize))
		spectrum8=np.zeros((nchan, isize))

		# forming a mask file with samples to reject
		mask=np.zeros((nchan, isize))
		clipped=np.zeros((nchan, isize))  # clipped data (only for plotting)
		mean=np.zeros((nchan, isize))     # 2D array of means
		rms=np.zeros((nchan, isize))      # 2D array of rms's
		levels=np.zeros((nchan, isize))   # 2D array of levels (samples in sigma)
		lev=np.zeros(nchan, dtype=float)

		# dictionary with indices of RFI'ish samples to be clipped
		clipindices = {}  # for not _very_ bad subbands
		badbands = []     # list of completely bad subbands

		scales = hdu[1].data.field('dat_scl')[row_offset:row_offset+nrows]
		offsets = hdu[1].data.field('dat_offs')[row_offset:row_offset+nrows]
		chandata = hdu[1].data.field('data')[row_offset:row_offset+nrows]

		for ii in xrange(nrows / nrows_step):
			print "row block: %d (%d)" % (ii, nrows / nrows_step)
			rb=ii*nrows_step
			re=(ii+1)*nrows_step
			for ch in np.arange(nchan):
				spectrum8[ch][ii] = np.average(chandata[rb:re,ch::nchan])
#				mean[ch][ii] = np.mean(np.sort(chandata[rb:re,ch::nchan], axis=None, kind='quicksort')[0:Nbins/2])
#				rms[ch][ii] = np.std(np.sort(chandata[rb:re,ch::nchan], axis=None, kind='quicksort')[0:Nbins/2])
				spectrum[ch][ii] = np.mean(chandata[rb:re,ch::nchan] * scales[rb:re,ch][:,np.newaxis] + offsets[rb:re,ch][:,np.newaxis])
				mean[ch][ii] = np.mean(np.sort(chandata[rb:re,ch::nchan] * scales[rb:re,ch][:,np.newaxis] + offsets[rb:re,ch][:,np.newaxis], axis=None, kind='quicksort')[0:Nbins/2])
				rms[ch][ii] = np.std(np.sort(chandata[rb:re,ch::nchan] * scales[rb:re,ch][:,np.newaxis] + offsets[rb:re,ch][:,np.newaxis], axis=None, kind='quicksort')[0:Nbins/2])

		# close input fits-file
		hdu.close()

		# getting statistics about fraction of positive/negative samples
		if is_statistics == True:
			for i in np.arange(nchan):
				lev[i] = [(spectrum[i][k] - mean[i][k])/rms[i][k] for k in np.arange(0, isize, 1)]	
				levsize=np.size(lev[i])
				condition = (lev[i] > threshold)
				pospeak = np.size(lev[i].compress(condition))
				condition = (lev[i] < -threshold)
				negpeak = np.size(lev[i].compress(condition))
				if pospeak > negpeak: numposlarger += 1
				if pospeak < negpeak: numpossmaller += 1

		# check if all values in rms are zeros. If so, then exclude this subband
		for i in np.arange(nchan):
			if np.size(np.trim_zeros(np.sort(rms[i]), 'bf')) != 0:
				levels[i] = [(spectrum[i][k] - mean[i][k])/rms[i][k] for k in np.arange(0, isize, 1)]
                        	mask[i] = [(np.abs(levels[i][k]) > threshold and spectrum[i][k] or 0) for k in np.arange(0, isize, 1)]
                        	clipped[i] = [(spectrum[i][k] - mask[i][k]) for k in np.arange(0, isize, 1)]
				# are there many zeros?
				condition=(np.abs(levels[i]) > threshold)
				rfi_fraction = (float(np.size(clipped[i].compress(condition)))/np.size(clipped[i]))*100.
			else:
				rfi_fraction = 100.

			if rfi_fraction >= rfilimit:  # bad subband  
				clipped[i] = np.zeros(isize)
				badbands.append(i)
				print "channel %d will be excluded (rfi fraction = %.2f%%)" % (i, rfi_fraction)
			else:
				clipindices[i] = np.where(np.abs(levels[i]) > threshold)[0]

		if is_saveonly == False:
			plt.ion() # interactive plotting
			fig=plt.figure()
			fig.canvas.mpl_connect('key_press_event', press)
			# normalizing
			setup_plot(spectrum8, "Original dynamic spectrum", cm.jet)
			plot_help()
			plot_update(spectrum8)

	if is_statistics == True:
		numposlarger = "%.1f" % (float(numposlarger)/nfiles, )
		numpossmaller = "%.1f" % (float(numpossmaller)/nfiles, )
		print "Number of pos>neg AND neg>pos occasions:  %s   %s" % (numposlarger, numpossmaller)
	badchanfreqs = [lofreq + float(sb) * chanbw for sb in badbands]
	np.savetxt("." + rfirepname, np.transpose((badbands, badchanfreqs)), fmt="%d\t\t%.6g")
	rfirep = open (rfirepname, 'w')	
	rfirep.write("# Subband	Freq (MHz)\n")
	rfirep.close()
	os.system("cat " + "." + rfirepname + " >> " + rfirepname)
	os.system("rm -f " + "." + rfirepname)

	if is_saveonly:
		fig=plt.figure()
		if not is_psrfits: 
			setup_plot(spectrum, "Original dynamic spectrum", cm.jet)
		else:
			setup_plot(spectrum8, "Original dynamic spectrum", cm.jet)
		plt.savefig(pngname)

	plt.show()
