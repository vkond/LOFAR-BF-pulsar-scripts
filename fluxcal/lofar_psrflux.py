#!/usr/bin/env python
#
# Measures pulsar flux density in Jansky's for LOFAR HBA data
# for the input PSRFITS (Psrchive) file.
# Fluxes are calculated for the total band or the wanted
# number of output channels. Diagnostic plots are made
# if requested
#
# Vlad Kondratiev (c) - 26.11.2014
#
# 26.03.2015 - Vlad Kondratiev 
#              minor/cosmetic changes, improved log output
# 03.04.2015 - Vlad Kondratiev
#              added new model to calculate Aeff, called 'arisN'
#              It's using maximum theoretical value +
#              scaling with EL as sin(EL)^1.39 as in Noutsos et al. (2015)
#              The original model called 'arts' is default
# 06.04.2015 - Vlad Kondratiev
#              added new model to calculate Aeff, called 'hamaker_carozzi'
#              It's using also maximum theoretical valueof Aeff in zenith
#              and correcting it by beam correction coefficient calculated
#              from Jones matrix based on Hamaker model and 'antennaJones.py'
#              script written by Tobia Carozzi.
# 16.04.2015 - Vlad Kondratiev
#              speed up hamaker_carozzi model by about 30 times by
#	       calling proper functions from mscorpol directly instead of
#	       running python script 'antennaJones.py'
# 24.04.2015 - Vlad Kondratiev
#              re-arranged hamaker_carozzi functions in a separate file,
#              so changed the module import calls
# 19.04.2016 - Vlad Kondratiev
#              added calculating "real" SEFD without taking into account
#              npol, integration time and bandwidth. In the output, this
#              is now called "SEFD", and previous "SEFD" is renamed to
#              "Sensit."
# 05.08.2016 - Vlad Kondratiev
#              added --station option to be used for calibrating observations
#              where only single station was used (e.g. FE). casa_beamcorr is 
#              a dictionary with CasA correction factors for all the CS and INTL 
#              stations. For RS stations mscorpol is currently failing. 
# 19.04.2017 - Vlad Kondratiev
#              changed default model to 'hamaker_carozzi';
#              checking the bw value now, if it's negative or not, indicating the
#              reversed frequency order. If it's negative, then calibration is
#              aborted and user is warned to reverse the frequency order with
#              the pam command; Also printing a warning now that Freq-averaged
#              values should be used with caution.
# 27.11.2017 - Vlad Kondratiev
#              improved calibration for international stations; correct number of
#              tiles/dipoles is used now (96). This is valid only for hamaker_carozzi
#              and arisN models. The arts model in its essence was implemented
#              only for core stations
"""
	Measures pulsar flux density in Jansky's for LOFAR HBA data
	for the input PSRFITS (Psrchive) file.
	Fluxes are calculated for the total band or the wanted
	number of output channels. Diagnostic plots are made
	if requested

	Vlad Kondratiev (c) - 26.11.2014

        26.03.2015 - Vlad Kondratiev 
                     minor/cosmetic changes, improved log output
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
	16.04.2015 - Vlad Kondratiev
		     speed up hamaker_carozzi model by about 30 times by
		     calling proper functions from mscorpol directly instead of
		     running python script 'antennaJones.py'
	24.04.2015 - Vlad Kondratiev
		     re-arranged hamaker_carozzi functions in a separate file,
		     so changed the module import calls
        19.04.2016 - Vlad Kondratiev
                     added calculating "real" SEFD without taking into account
                     npol, integration time and bandwidth. In the output, this
                     is now called "SEFD", and previous "SEFD" is renamed to "Sensit."
        05.08.2016 - Vlad Kondratiev
                     added --station option to be used for calibrating observations
                     where only single station was used (e.g. FE). casa_beamcorr is 
                     a dictionary with CasA correction factors for all the CS and INTL 
                     stations. For RS stations mscorpol is currently failing. 
        19.04.2017 - Vlad Kondratiev
                     changed default model to 'hamaker_carozzi';
                     checking the bw value now, if it's negative or not, indicating the
                     reversed frequency order. If it's negative, then calibration is
                     aborted and user is warned to reverse the frequency order with
                     the pam command; Also printing a warning now that Freq-averaged
                     values should be used with caution.
        27.11.2017 - Vlad Kondratiev
                     improved calibration for international stations; correct number of
                     tiles/dipoles is used now (96). This is valid only for hamaker_carozzi
                     and arisN models. The arts model in its essence was implemented
                     only for core stations
"""

import numpy as np
import os, os.path, stat, glob, sys, getopt, re
import scipy.stats as sc
import optparse as opt
import psrchive as pc
import ephem
import h5py
from tsky import *
from lofar_tinst import *
from lofar_gain import lofar_gain_range_arts, lofar_gain_range_arisN
import warnings
warnings.simplefilter('ignore', np.RankWarning)
warnings.simplefilter('ignore', RuntimeWarning)


# function to read HDF5 .h5 file to get info about number of stations
# and flagged tiles/dipoles (in the future)
# return a tuple with number of Core stations and 
# fraction of flagged tiles/dipoles (currently we are returning None)
def read_meta(h5file):
	try:
		f5 = h5py.File(h5file, 'r')
		bandFilter = f5.attrs['FILTER_SELECTION']
		antenna = bandFilter.split("_")[0]
		stations=f5.attrs['OBSERVATION_STATIONS_LIST']
		ncorestations = len([s for s in stations if s[0:2] == "CS"])
		# because in the list for HBA there are sub-stations
		if antenna == "HBA": ncorestations /= 2
                # in case only one station (remote or international) was used
                # if several non-core stations were used, then flux calibration won't be correct as
                # here we assume coherent stations addition
                if ncorestations == 0:
                        ncorestations = len(stations)
		# in the future here we should also add reading info
		# of the flagged tiles/dipoles (broken antenna info)
		flagged_fraction = None
	except:
		ncorestations = 0
		flagged_fraction = None
		stations = []
	return (ncorestations, flagged_fraction, stations)

# Tom's function to get number of bad tiles/dipoles using imaging
# data close in time to given pulsar obs
# "obsid" is assumed to have leading "L"
def badDipoles(obsid, band):
	if band == "HBA":
		f = open('/home/astron/kondratiev/pulsar/etc/flagged-antennas-tom-hassall/flaggedAntennasHBA','r')
	elif band == "LBA":
		f = open('/home/astron/kondratiev/pulsar/etc/flagged-antennas-tom-hassall/flaggedAntennasLBA','r')
	else:
		print "Band not recognised, assuming all antennas are active"
		return 0.0
	# checking the ObsIDs if they are higher than the last entries
	# in Tom's file and return 0, if they are more recent
	if float(obsid[1:]) > 181858 and band == "HBA":
		return 0
	if float(obsid[1:]) > 180236 and band == "LBA":
		return 0
	x = np.inf
	for ll in f.readlines():
		line = ll.split()
		if (float(obsid[1:]) - float(line[0]))**2 < x:
			x = (float(obsid[1:]) - float(line[0]))**2
			nFlagged = float(line[1])
	return nFlagged

# returning mean/rms using Q-Q probability
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

# calculating mean/rms of Calibrated Profile
# method is one of the 3 predetermined strings: "QQ", "Off", or "Polynom"
def mean_rms (method, prof, osm_min=0.0, osm_max=1.0, off_left=0, off_right=1, polynom_ndeg=-1):
	if method == "QQ": # Q-Q
               	(mean, rms, osm, osr, qmin, qmax) = get_mean_rms(prof, osm_min, osm_max)
	elif method == "Off": # Off window
               	mean = np.mean(prof[off_left:off_right])
               	rms = np.std(prof[off_left:off_right])
	elif method == "Polynom": # Polynom
		if polynom_ndeg == -1:
	               	for polynom_ndeg in xrange(int(float(len(prof)-1)),1,-1):
        	               	try:
                	               	tmp=float(len(prof)-1)**polynom_ndeg
                        	        break
	       	                except: pass
		polynom_coeffs = np.polyfit(range(len(prof)), prof, polynom_ndeg)
                polynom_bline = np.polyval(polynom_coeffs, range(len(prof)))
               	mean = np.mean(prof-polynom_bline)+np.mean(sorted(polynom_bline)[:int(0.2*len(polynom_bline))])
               	rms = np.std(prof-polynom_bline)
	else:
		print "Wrong S/N method is used!"
		sys.exit(1)
	return (mean, rms)

# convert "dd:mm:ss" to degrees
def ddmmss2deg(pres):
	vals=pres.split(":")
	val = abs(float(vals[0])) + float(vals[1])/60. + float(vals[2])/3600.
	if vals[0].strip()[0] == '-':
		val *= -1.
	return val

# Main
if __name__=="__main__":
        #
        # Parsing the command line options
        #
        usage = "Usage: %prog [options] .ar"
        cmdline = opt.OptionParser(usage)
	cmdline.add_option('--meta', dest='metafile', metavar='.h5', help="HDF5 .h5 file with the meta info about the observation. \
Info about number of stations and bad tiles/dipoles (in the future) will be extracted. When used, options --nstat and --flagged are \
ignored", default="", type='str')
	cmdline.add_option('--flagged', '--badtiles', '--baddipoles', dest='badtiles', metavar='FRACTION | ObsID', help="Fraction of \
bad tiles/dipoles. For older observations one can use ObsId, in this case the fraction of bad tiles/dipoles will be looked up in the \
Tom Hassal's list based on metadata from imaging observations close in time. An ObsID should have leading letter before the number. \
Default = %default.", default="0.05", type='str')
        cmdline.add_option('--nstat', dest='nstations', metavar='#STATIONS', help="Number of 48-tile stations used in the observation. \
Default = %default", default=24, type='int')
        cmdline.add_option('--cohfactor', '--coherence-factor', dest='cohfactor', metavar='FACTOR', help="Coherence scaling factor 'gamma' in \
S/N~Nstat^gamma, where Nstat is number of 48-tile stations. Default = %default", default=0.85, type='float')
        cmdline.add_option('--snrmethod', dest='method', metavar='STRING', help="Method to calculate the mean/rms. Possible values \
are: 'Off', 'QQ', and 'Polynom'. 'Off' is when user specifies off-pulse window to use, 'QQ' is to use Q-Q probability plot \
to determine range of quantiles that satisfies Gaussian distribution, and 'Polynom' is like 'Off' but all profile bins are used to calculate \
mean/rms after subtracting polynomial fit to the pulse profile. Default: %default", default="QQ", type='str')
	cmdline.add_option('--model', '--beam-model', dest='model', metavar='STRING', help="Beam model to use. Possible \
values: 'arts' for the beam model by Arts et al. (2013), 'arisN' that uses theoretical maximum of Aeff in \
zenith and then scales it with elevation as sin(EL)^1.39 as in Noutsos et al. (2015), and 'hamaker_carozzi' that uses \
also theoretical maximum of Aeff in zenith and then calculates beam correcting coefficient based on Jones matrix \
calculated using Hamaker model and antennaJones.py script by Tobia Carozzi. Default: %default", default="hamaker_carozzi", type='str')
        cmdline.add_option('-b', '--bscrunch', dest='bscr', metavar='FACTOR', help="Bscrunch factor, \
default: %default", default=1, type='int')
	cmdline.add_option('-F', dest='to_Fscrunch', action="store_true", help="To Fscrunch _all_ frequency channels together", default=False)
        cmdline.add_option('-f', '--fscrunch', dest='fscr', metavar='FACTOR', help="Fscrunch factor, \
default: %default", default=1, type='int')
	cmdline.add_option('-r', '--rotate', dest='rot_bins', metavar='#BIN|PHASE', help="Rotate profile by this number of bins (if absolute value >= 1). \
If absolute value is < 1, then the value is treated as pulse phase in turns. Negative values are to move right, default: %default", default=0, type='float')
        cmdline.add_option('--spectrum', dest='spchan', metavar='#OUTCHAN', help="Calculate flux density for #OUTCHAN of output channels", default=-1, type='int')
        cmdline.add_option('--spectrum-skip-first-channels', dest='spskip_first', metavar='#INCHAN', help="For spectrum calculation do not take into account \
_first_ #INCHAN input channels (possibly fscrunched if --fscrunch was used). Default: %default", default=0, type='int')
        cmdline.add_option('--spectrum-skip-last-channels', dest='spskip_last', metavar='#INCHAN', help="For spectrum calculation do not take into account \
_last_ #INCHAN input channels (possibly fscrunched if --fscrunch was used). Default: %default", default=0, type='int')
	cmdline.add_option('--max-weight', dest='str_weight', metavar='WEIGHT,NSUBS,NCHANS', help="Provide a string with 3 values separated by commas. First value \
is the maximum weight in the original PSRFITS file, second and third values are number of subints and channels in the this original file, respectively. \
For the expert use only. Without this option, all weights will be normalized by the maximum weight in the input file", default="", type='str')
        cmdline.add_option('--off-left', dest='off_left', metavar='BIN#', help="Left edge of the off-pulse window to calculate \
mean/rms in the 'Off' mode (inclusive). Default: %default", default=0, type='int')
        cmdline.add_option('--off-right', dest='off_right', metavar='BIN#', help="Right edge of the off-pulse window to calculate \
mean/rms in the 'Off' mode (exclusive). Default: 10%-bin of the profile", default=-1, type='int')
        cmdline.add_option('--osm-min', dest='osm_min', metavar='MIN PROB', help="Minimum probability value to be used \
to calculate the mean and rms, default: min possible", default=-9999., type='float')
        cmdline.add_option('--osm-max', dest='osm_max', metavar='MAX PROB', help="Maximum probability value to be used \
to calculate the mean and rms, default: %default", default=+0.95, type='float')
	cmdline.add_option('--polynom-deg', dest='polynom_ndeg', metavar='DEGREE', help="Degree of a polynomial for 'Polynom' method. \
Default is maximum fittable starting from number of bins minus 1", default=-1, type='int')
        cmdline.add_option('--latitude', dest='latitude', metavar='DD:MM:SS.SS', help="specify the latitude of a site. \
Default is for LOFAR (%default)", default="52:52:59.88", type='str')
        cmdline.add_option('--longitude', dest='longitude', metavar='DD:MM:SS.SS', help="specify the longitude of a site. \
Default is for LOFAR (%default)", default="06:52:00.12", type='str')
	cmdline.add_option('--station', dest='onlystation', metavar='NAME', help="Name of the station for this observation to be used for calibration. \
Without this option the default station is CS002 when many core stations were used in observation. If only one station was used, this station will be the \
reference one unless this option is used", default="", type='str')
	cmdline.add_option('-v', '--verbose', dest='is_verbose', action="store_true", help="Verbose output, print extra info", default=False)
	cmdline.add_option('-V', '--vv', '--very-verbose', dest='is_very_verbose', action="store_true", help="Print related info for every calibrated channel", default=False)
	cmdline.add_option('--no-print-total', dest='is_no_print_total', action="store_true", help="Do no print flux values in the total band when --spectrum is used", default=False)
	cmdline.add_option('--plot', dest='is_plot', action="store_true", help="Make diagnostic plots", default=False)
	cmdline.add_option('--plot-saveonly', dest='is_saveonly', action="store_true", help="Save diagnostic plot \
to png-file instead of GUI", default=False)

        # reading cmd options
        (opts,args) = cmdline.parse_args()

        # check if input file is given
        if len(args) == 0:
                cmdline.print_usage()
                sys.exit(0)

	# loading pyrap.measures as it required by mscorpol
	if opts.model == 'hamaker_carozzi':
		from pyrap.measures import measures
		from lofar_gain_hamaker_carozzi import lofar_gain_range_hamaker_carozzi_center_station, lofar_gain_range_hamaker_carozzi_station

	if opts.is_very_verbose:
		opts.is_verbose = True	

	# reading input ar-file
	infile = args[0]
	try:
		raw = pc.Archive_load(infile)
	except:
		print "Can't open input file: '%s'" % (infile)
		sys.exit(1)
	if not raw.get_dedispersed():
		raw.dedisperse()
	# raw.remove_baseline() - we subtract it outselves
        raw.pscrunch()
	npol = 2 # as we scrunched all pols together
        nsubint = raw.get_nsubint()
	target = raw.get_source()
        orig_nchan = raw.get_nchan()

	# checking the weights and re-normalizing them if needed
	weights=raw.get_weights()
	subint_durs = [raw.get_Integration(ii).get_duration() for ii in xrange(nsubint)]
	if opts.str_weight == "":
		max_weight = np.max(weights)
	else:
		try:
			str_warr = opts.str_weight.split(",")
			max_weight = float((float(str_warr[0]) * int(str_warr[1]) * int(str_warr[2])) / (nsubint * orig_nchan))
		except:
			max_weight = np.max(weights)
			print "Error in reading --max-weight option!"
			print "Warning!!! Weights will be calibrated using the weights' maximum value from the input file!"
	max_subint_dur = np.max(subint_durs)
	for ii in xrange(nsubint):
		isub=raw.get_Integration(ii)
		dur=isub.get_duration()
		for ch in xrange(orig_nchan):
			w=isub.get_weight(ch)
			if w<max_weight:
				w *= (max_subint_dur/dur)
				if w>max_weight: w=max_weight
			w /= max_weight
			isub.set_weight(ch, float(w))

	# re-reading weights again and calculate the RFI fraction
	weights=raw.get_weights()
	rfi_fraction=1.-np.sum(weights)/(nsubint*orig_nchan)

	# tscrunching
        if nsubint > 1: raw.tscrunch()

	# bscrunching
        if opts.bscr > 1: raw.bscrunch(opts.bscr)
        nbins = raw.get_nbin()

	# fscrunching
	if opts.to_Fscrunch: # then fscrunching _all_ freq channels
		raw.fscrunch()
	else:
        	if opts.fscr > 1: raw.fscrunch(opts.fscr)
        nchan = raw.get_nchan()

	# rotate profile if necessary
	if opts.rot_bins != 0:
		if abs(opts.rot_bins) < 1:
			raw.rotate_phase(opts.rot_bins)
		else:
			raw.rotate_phase(opts.rot_bins/nbins)

	# getting weights of our possibly-fscrunched data
	weights = raw.get_weights()

	cfreq = raw.get_centre_frequency()  # center freq in MHz
	bw = raw.get_bandwidth()            # bandwidth in MHz
        # checking if bw is negative. If it < 0, then it means that frequency order is reversed
        # in this case we need to abort, as lowest frequency is expected to have lower index
        # so, we abort and let user that (s)he needs to reverse the frequency order with the pam
        # pam --reverse_freqs
        if bw < 0:
                print "Frequency order is reversed (highest frequency goes first, bandwidth value is negative)!"
                print "You must reverse the frequency order before continuing with the calibration."
                print "You can do that with the pam command:"
                print "pam --reverse_freqs [-e|-m|-u] <infile>"
                sys.exit(1)
	lowfreq = cfreq - bw/2.0            # lowest freq in MHz
	chan_bw = bw/nchan                  # channel width in MHz
	tobs = raw.integration_length()     # obs duration (in seconds)

	# reading h5 file to get info about Nstat and in the future flagged tiles/dipoles
	badtiles = None
	if opts.metafile != "":
		(nstat, badtiles, stations) = read_meta(opts.metafile)
		if nstat != 0:
			opts.nstations = nstat
		if len(stations) == 1:
			if opts.onlystation == "": opts.onlystation = stations[0][0:5]
			opts.nstations = 1
		else:
			if opts.onlystation == "": opts.onlystation = "CS002"
	else:
		if opts.nstations > 1:
			if opts.onlystation == "": opts.onlystation = "CS002"
		elif opts.nstations <= 0:
			print "Number of stations should be more than 0!"
			sys.exit(1)
		else:
			if opts.onlystation == "":
				print "Only one station is used, but it's not known what station."
				print "Either is .h5 metafile or provide reference station using --station option!"
				sys.exit(1)
		
	# checking if user gave ObsId for the badtiles fraction
	if badtiles is None:
		if opts.badtiles[0].isalpha():
			if lowfreq > 100:
				badtiles=badDipoles(opts.badtiles, "HBA")
			else:
				badtiles=badDipoles(opts.badtiles, "LBA")
		else:
			badtiles=float(opts.badtiles)

	# getting pulsar RA/DEC
	coord=raw.get_coordinates().getHMSDMS()
	raj=coord.split()[0]
	decj=coord.split()[-1]

	# getting Galactic coordinates
	(gl, gb) = eq2gal(raj, decj)

	# getting AZ/Elevation
	star = ephem.FixedBody()
	star._ra = raj
	star._dec = decj
	if opts.model == 'hamaker_carozzi':
		direction=measures().direction('J2000', str(float(repr(star._ra)))+'rad', str(float(repr(star._dec)))+'rad')

	# getting start date/time and mid-point of obs
	start_mjd = float(raw.start_time().strtempo())
	midpoint_mjd = start_mjd + 0.5*(tobs/86400.)
	midpoint_dublin_day = midpoint_mjd - 15019.5 # in Dublin Days, Dublin day = JD - 2415020
	midpoint = ephem.Date(midpoint_dublin_day)

	observer = ephem.Observer()
	observer.date = midpoint
	try:
		observer.long = opts.longitude
	except:
		# for now old versions of ephem has longitude attribute "long"
		# and newer version have both "long" and "lon". I am checking this here
		# in case, "long" will become deprecated
		observer.lon = opts.longitude
	observer.lat = opts.latitude

	star.compute(observer)
	azimuth = ddmmss2deg(str(star.az))
	elevation = ddmmss2deg(str(star.alt))



	# Print extra info in verbose mode
	if opts.is_verbose:
		print "#"
		print "# Source: %s" % (target)
		print "# RA: %s   DEC: %s   GL(deg): %g   GB(deg): %g" % (raj, decj, gl, gb)
		print "# Start MJD: %.15g   Mid-obs MJD: %.15g" % (start_mjd, midpoint_mjd)
		print "# AZ(deg): %g   EL(deg): %g" % (azimuth, elevation)
		print "# Tobs(s): %g   Cfreq(MHz): %g   BW(MHz): %g   OrigNchan: %d   ChanWidth(MHz): %g" % (tobs, cfreq, bw, orig_nchan, bw/orig_nchan)
		if opts.to_Fscrunch:
			print "# Nsub: %d (tscrunched)   Nchan: %d (Fscrunched)   ChanWidth(MHz): %g   Nbins: %d%s" % (nsubint, nchan, \
				chan_bw, nbins, opts.bscr>1 and " (bscrunched: %d)" % (opts.bscr) or "")
		else:
			print "# Nsub: %d (tscrunched)   Nchan: %d%s   ChanWidth(MHz): %g   Nbins: %d%s" % (nsubint, nchan, \
				opts.fscr>1 and " (fscrunched: %d)" % (opts.fscr) or "", chan_bw, nbins, opts.bscr>1 and " (bscrunched: %d)" % (opts.bscr) or "")
		print "# Npol: %d   Data state: %s" % (raw.get_npol(), raw.get_state())
		print "# Nstations: %d   Coherence factor: %g" % (opts.nstations, opts.cohfactor)
		print "# Bad dipoles/tiles(%%): %g" % (badtiles * 100.)
		print "# RFI fraction(%%): %g" % (rfi_fraction * 100.)
		sn_params = ""
		if opts.method == "QQ":
			if opts.osm_min != -9999.: sn_params += " --osm_min %g" % (opts.osm_min)
			sn_params += " --osm_max %g" % (opts.osm_max)
		if opts.method == "Off":
			sn_params += " --off-left %d" % (opts.off_left)
			if opts.off_right != -1: sn_params += " --off-right %d" % (opts.off_right)
		if opts.method == "Polynom":
			if opts.polynom_ndeg != -1: sn_params += " --polynom-deg %d" % (opts.polynom_ndeg)
		print "# S/N method: %s%s" % (opts.method, sn_params != "" and " [%s]" % sn_params or "")
		print "# Beam model: %s" % (opts.model)
		print "# File: %s" % (infile)
	if opts.is_very_verbose:
		print "#"
		print "# Cmdline: %s %s" % (sys.argv[0].split("/")[-1], " ".join(sys.argv[1:]))

	# Boltzmann constant
	kb = 1.3806488*1e-16 # erg/K
	beta = 1 # digitization factor
	skytemp_index = -2.55

       	r = raw.get_data()
       	#time stokes f phase
	data = r[0,0,:,:]

	if opts.off_right == -1:
		opts.off_right = int(nbins*0.1)
	if opts.off_right-opts.off_left<=1:
		opts.off_right = nbins
	if opts.off_right > nbins:
		opts.off_right = nbins

	prof=np.zeros((nchan,nbins),dtype=float) # 2-D array with calibrated profiles for each channel
	prof_err=np.zeros(nchan, dtype=float)    # flux errors for each channel
	true_sefd=np.zeros(nchan, dtype=float)   # True SEFD for each channel

	# RMS->Jy scaling factor 
	# here it's without T/aeff, only part that does not depend on frequency
	# 1e4 -> this is conversion from m^2 to cm^2
	# 1e23 -> this is conversion from  erg/s/cm^2/Hz to Jy
	# 1e3 -> this is conversion from Jy to mJy
	delta_s = 1e3*1e23*(2.*beta*kb)/((1.-badtiles)*(1e4*opts.nstations**opts.cohfactor)*np.sqrt(npol*(tobs/nbins)*chan_bw*1e6)) # in mJy
	# same as delta_s but without radiometer "advantage", i.e. without sqrt(npol*T*BW) - for true SEFD calculation
	true_delta_s = 1e23*(2.*beta*kb)/((1.-badtiles)*(1e4*opts.nstations**opts.cohfactor)) # in Jy

	# array of freqs
	freqs=[lowfreq + ch * chan_bw for ch in xrange(nchan + 1)]

	# array of sky temperatures for each freq channel
	skytemps=tsky_range(gl, gb, lowfreq, lowfreq+nchan*chan_bw, skytemp_index, freqs)

	# array of Tinst for each frequency channel
	tinsts=lofar_tinst_range(lowfreq, lowfreq+nchan*chan_bw, freqs)

	# array of Aeffs for each freq channel
	if opts.model == 'arts':  # beam model by Arts et al. (2013)
		aeffs=lofar_gain_range_arts(lowfreq, lowfreq+nchan*chan_bw, elevation, freqs)
	elif opts.model == 'arisN':  # here we are using scaling with EL from Noutsos et al. (2015)
		aeffs=lofar_gain_range_arisN(lowfreq, lowfreq+nchan*chan_bw, elevation, opts.onlystation, freqs)
	elif opts.model == 'hamaker_carozzi':  # here we are using Hamaker/Carozzi model
		aeffs=lofar_gain_range_hamaker_carozzi_center_station(lowfreq, lowfreq+nchan*chan_bw, elevation, midpoint_mjd, direction, opts.onlystation, freqs)
	else: # by default the beam model is Hamaker/Carozzi model
		aeffs=lofar_gain_range_hamaker_carozzi_center_station(lowfreq, lowfreq+nchan*chan_bw, elevation, midpoint_mjd, direction, opts.onlystation, freqs)

	if opts.is_very_verbose:
		print "#"
		print "# Ch# Freq%6smean%7srms%8sS/N%5sS/N%6sSEFD%5sSensit.%2sTsky%5sTinst%4sTsys%5sAeff%5sGain%4sRFI Fr.%3sFlux%4sFlux Err" % (" ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ")
		print "#%5s(MHz)%27smean%4speak%5s(Jy)%6s(mJy)%5s(K)%6s(K)%6s(K)%5s(m^2)%3s(K/Jy)%4s(%%)%6s(mJy)%4s(mJy)" % (" ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ")
		print "#--------------------------------------------------------------------------------------------------------------------------------------------"

	# Loop over frequency channels
	for ch in xrange(nchan):


		# getting RFI fraction, or actually Good fraction
		good_fraction=weights[0,ch]/(nsubint*(orig_nchan/nchan))
		if good_fraction == 0: 
			if opts.is_very_verbose:
				print "#%-4d %-8g   zapped" % (ch, freqs[ch]+chan_bw/2.)
			continue

		# Full system temperature
		Tsys = skytemps[ch]+tinsts[ch]
		# SEFD
		sefd = delta_s*Tsys/(aeffs[ch]*np.sqrt(good_fraction)) # 1-sigma sensitivity in mJy
		true_sefd[ch] = true_delta_s*Tsys/(aeffs[ch])          # real SEFD in Jy
		# Full effective area (m^2)
		Aeff = aeffs[ch] * (opts.nstations**opts.cohfactor) * (1.-badtiles)
		# Gain (K/Jy)
		Gain = (Aeff * 1e4 * 1e-23) / (2.* kb)

		# getting mean/rms
		(mean, rms) = mean_rms(opts.method, data[ch], opts.osm_min, opts.osm_max, opts.off_left, opts.off_right, opts.polynom_ndeg)
		prof[ch] = (data[ch] - mean)/rms
		snrmean = np.mean(prof[ch])
		snrpeak = np.max(prof[ch])

		prof[ch] *= sefd
		prof_err[ch] = (sefd * sefd)

		if opts.is_very_verbose:
			print "%-4d %-8g %-8g %-12g %-8g %-8g %-8g %-8g %-8g %-8g %-8g %-8g %-8g %-8g %-8g %-8g" % (ch, freqs[ch]+chan_bw/2., mean, rms, snrmean, snrpeak, \
				true_sefd[ch], sefd, skytemps[ch], tinsts[ch], Tsys, Aeff, Gain, (1.-good_fraction)*100., np.mean(prof[ch]), sefd/np.sqrt(nbins))
		
	if opts.is_very_verbose:
		print "#--------------------------------------------------------------------------------------------------------------------------------------------"

	# checking if more than 1 output channels are specified
	if opts.spskip_first < 0: opts.spskip_first = 0
	if opts.spskip_first >= nchan: opts.spskip_first = 0
	if opts.spskip_last < 0: opts.spskip_last = 0
	if opts.spskip_last >= nchan: opts.spskip_last = 0
	if opts.spskip_first+opts.spskip_last >= nchan:
		opts.spskip_last -= (1+opts.spskip_first+opts.spskip_last-nchan)
	if opts.spchan != -1:
		if opts.spchan < 1: opts.spchan = nchan
		if opts.spskip_first + opts.spchan > nchan - opts.spskip_last: opts.spchan = nchan - opts.spskip_last - opts.spskip_first
		inch_per_outch = int(float(nchan-opts.spskip_first-opts.spskip_last)/opts.spchan) # number of input channels per 1 output spectrum channel

		# 2-D array of spectrum mean flux densities and their errors, [0] - value, [1] - error, [2] - frequency
		spectrum=np.zeros((3, opts.spchan), dtype=float)

		print "#"
		print "# Freq%4sSEFD%5sS/N%6sS/N%6sProf%4sChi^2/%4sWeff%6sDC%6sSpeak%3sSensit.%3sSmean  Smean Err" % (" ", " ", " ", " ", " ", " ", " ", " ", " ", " ")
		print "# (MHz)%3s(Jy)%5smean%5speak%5sSign.%4sdof%5s(bins)%5s(%%)%5s(mJy)%4s(mJy)%4s(mJy)%4s(mJy)" % (" ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ")
		print "#-----------------------------------------------------------------------------------------------------------"
		
		for spch in xrange(0,opts.spchan):

			nzapped=0
			# normalize by a number of used channels

			spleft = opts.spskip_first + spch*inch_per_outch
			spright = opts.spskip_first + (spch+1)*inch_per_outch

			totprof = np.sum(prof[spleft:spright], axis=0)

			# width of output spectrum channel
			spbw = inch_per_outch * chan_bw
			# lowest freq of the output channel
			splow = freqs[spleft]

			for ii in xrange(spleft, spright):
				if not np.any(prof[ii]): nzapped += 1
			if inch_per_outch-nzapped == 0:
				print "#%-8g   zapped" % (splow+spbw/2.)
				continue

			totprof /= (inch_per_outch-nzapped)
			flux_err = np.sum(prof_err[spleft:spright])
			# 1-sigma sensitivity
			SEFD = np.sqrt(flux_err) / (inch_per_outch-nzapped) # 1-sigma sensitivity or peak flux err in mJy
			# true SEFD in Jy
			real_sefd = np.sqrt(np.sum([s*s for s in true_sefd[spleft:spright]]) / (inch_per_outch-nzapped))

			# getting mean/rms
			(mean, rms) = mean_rms(opts.method, totprof, opts.osm_min, opts.osm_max, opts.off_left, opts.off_right, opts.polynom_ndeg)
			flux_peak = np.max(totprof)                # peak flux in mJy
			binpeak = np.argmax(totprof)               # bin with the max flux
			snrpeak = (flux_peak - mean)/rms           # S/N peak

			flux_mean = np.mean(totprof)               # mean flux in mJy
			flux_mean_err = SEFD / np.sqrt(nbins)      # mean flux err in mJy
			snrs = (totprof-mean)/rms
			snrmean = np.mean(snrs)                    # S/N mean
			# filling in spectrum array
			spectrum[0][spch] = flux_mean
			spectrum[1][spch] = flux_mean_err
			spectrum[2][spch] = splow+spbw/2. # center freq of the output channel

			# pulse effective width
			weq = np.sum(totprof)/flux_peak            # in bins
			dc = 100.*weq/nbins                        # pulse duty cycle in %
			profsign = np.sum(snrs)/np.sqrt(weq)       # profile significance
			chi2 = np.sum(snrs*snrs)/(nbins-1)         # chi^2 of the profile

			# output info
			print "%-8g %-8g %-8g %-8g %-8g %-8g %-8g %-8g %-8g %-8g %-8g %-8g" % (splow+spbw/2., real_sefd, snrmean, snrpeak, profsign, chi2, weq, dc, flux_peak, SEFD, flux_mean, flux_mean_err)



	# Flux measurement in the total band

	# normalize by a number of used channels
	totprof = np.sum(prof, axis=0)
	nzapped=0
	for ii in xrange(nchan):
		if not np.any(prof[ii]): nzapped += 1
	totprof /= (nchan-nzapped)

	# getting mean/rms
	(mean, rms) = mean_rms(opts.method, totprof, opts.osm_min, opts.osm_max, opts.off_left, opts.off_right, opts.polynom_ndeg)
	flux_peak = np.max(totprof)                # peak flux in mJy
	binpeak = np.argmax(totprof)               # bin with the max flux
	flux_err = np.sum(prof_err)
	# 1-sigma sensitivity
	SEFD = np.sqrt(flux_err) / (nchan-nzapped) # 1-sigma sensitivity or peak flux err in mJy
	# true SEFD in Jy
	real_sefd = np.sqrt(np.sum([s*s for s in true_sefd]) / (nchan-nzapped))
	snrpeak = (flux_peak - mean)/rms           # S/N peak

	flux_mean = np.mean(totprof)               # mean flux in mJy
	flux_mean_err = SEFD / np.sqrt(nbins)      # mean flux err in mJy
	snrs = (totprof-mean)/rms
	snrmean = np.mean(snrs)                    # S/N mean

	# pulse effective width
	weq = np.sum(totprof)/flux_peak            # in bins
	dc = 100.*weq/nbins                        # pulse duty cycle in %
	profsign = np.sum(snrs)/np.sqrt(weq)       # profile significance
	chi2 = np.sum(snrs*snrs)/(nbins-1)         # chi^2 of the profile

	# Calculating the average system parameters over the whole freq range
	Tsky = tsky_range(gl, gb, lowfreq, lowfreq+nchan*chan_bw, skytemp_index)
	Tinst = lofar_tinst_range(lowfreq, lowfreq+nchan*chan_bw)
	Tsys = Tsky + Tinst  # full system temperature in K
	# calculating Aeff
	if opts.model == 'arts':  # beam model by Arts et al. (2013)
		Aeff = lofar_gain_range_arts(lowfreq, lowfreq+nchan*chan_bw, elevation)
	elif opts.model == 'arisN':  # here we are using scaling with EL from Noutsos et al. (2015)
		Aeff = lofar_gain_range_arisN(lowfreq, lowfreq+nchan*chan_bw, elevation, opts.onlystation)
	elif opts.model == 'hamaker_carozzi':  # here we are using Hamaker/Carozzi model
		Aeff = lofar_gain_range_hamaker_carozzi_center_station(lowfreq, lowfreq+nchan*chan_bw, elevation, midpoint_mjd, direction, opts.onlystation)
	else: # by default the beam model is Hamaker/Carozzi model
		Aeff = lofar_gain_range_hamaker_carozzi_center_station(lowfreq, lowfreq+nchan*chan_bw, elevation, midpoint_mjd, direction, opts.onlystation)
	Aeff *= ((opts.nstations**opts.cohfactor) * (1.-badtiles))
	# Gain (K/Jy)
	Gain = (Aeff * 1e4 * 1e-23) / (2.* kb)
	# average SEFD (in Jy)
	average_sefd = 1e23 * (2. * beta * kb * Tsys)/(Aeff * (1.-badtiles) * (1e4 * opts.nstations**opts.cohfactor)) 

	if opts.spchan == -1 or not opts.is_no_print_total:
		# Print extra info in verbose mode
		if opts.is_verbose:
			print "#"
                        print "# ===== WARNING! Use this Freq-averaged values with the great caution! ====="
                        print "# =====          They are only meant for the rough estimation          ====="
                        print "#"
			print "# Center Frequency(MHz): %g   BW(MHz): %g" % (cfreq, bw)
			print "# Freq-averaged Tsky(K): %g" % (Tsky)
			print "# Freq-averaged Tinst(K): %g" % (Tinst)
			print "# Freq-averaged Tsys (K): %g" % (Tsys)
			print "# Freq-averaged Aeff(m^2): %g" % (Aeff)
			print "# Freq-averaged Gain(K/Jy): %g" % (Gain)
			print "# Freq-averaged SEFD(Jy): %g" % (average_sefd)
			print "#"

		print "#"
		print "#%3sPSR%6sSEFD%5sS/N%6sS/N%6sProf%4sChi^2/%4sWeff%6sDC%6sSpeak%4sSensit.%2sSmean  Smean Err" % (" ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ")
		print "#%12s(Jy)%5smean%5speak%5sSign.%4sdof%5s(bins)%5s(%%)%5s(mJy)%4s(mJy)%4s(mJy)%4s(mJy)" % (" ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ")
		print "#--------------------------------------------------------------------------------------------------------------"
		print "%-10s  %-8g %-8g %-8g %-8g %-8g %-8g %-8g %-8g %-8g %-8g %-8g" % (target, real_sefd, snrmean, snrpeak, profsign, chi2, weq, dc, flux_peak, SEFD, flux_mean, flux_mean_err)

	# Plotting
	if opts.is_plot:
		if opts.is_saveonly:
			pngname = ".".join(infile.split(".")[:-1]) + "-flux.png"
			import matplotlib
			matplotlib.use("Agg")

		import matplotlib.pyplot as plt
		import matplotlib.cm as cm
		from matplotlib.ticker import *

		fig = plt.figure(figsize=(8,6))
		if opts.spchan != -1:
			ax1 = fig.add_subplot(211)
		else:
			ax1 = fig.add_subplot(111)
		plt.xlabel("Bin", fontsize=12)
		plt.ylabel("Flux density (mJy)", fontsize=12)

		ax1.plot(range(len(totprof)), totprof, "g-", alpha=0.7)
		ax1.axvline(x=binpeak, linestyle="--", color="black")
		ax1.set_xlim(xmin=0, xmax=len(totprof)-1)

		if opts.method == "Off":
			ax1.axvspan(opts.off_left, opts.off_right, color="yellow", alpha=0.2)
			ax1.axvline(x=opts.off_left, linestyle=":", color="black")
			ax1.axvline(x=opts.off_right, linestyle=":", color="black")

		if opts.method == "Polynom":
			if opts.polynom_ndeg == -1:
				for polynom_ndeg in xrange(int(float(nbins-1)),1,-1):
					try:
						tmp=float(nbins-1)**polynom_ndeg
						break
					except: pass
			else:
				polynom_ndeg = opts.polynom_ndeg
			polynom_coeffs = np.polyfit(range(nbins), totprof, polynom_ndeg)
			polynom_bline = np.polyval(polynom_coeffs, range(nbins))
			ax1.plot(range(len(totprof)), polynom_bline, "b-", alpha=0.5)

		ax1.axhline(y=0.0, linestyle="--", color="black")
		ax1.axhspan(-SEFD, SEFD, color="grey", alpha=0.1)
		plt.gca().minorticks_on()

		# plotting the spectrum
		if opts.spchan != -1:
			ax2 = fig.add_subplot(212)
			plt.xlabel("Frequency (MHz)", fontsize=12)
			plt.ylabel("Flux density (mJy)", fontsize=12)
			ax2.set_xscale('log')
			ax2.set_yscale('log')
			ax2.set_xlim(xmin=100*int((cfreq+bw/2.)/100.), xmax=100*(1+int((cfreq+bw/2.)/100.)))
			ax2.set_ylim(ymin=10**(-0.5+np.log10(np.min(spectrum[0]))), ymax=10**(0.5+np.log10(np.max(spectrum[0]))))
			ax2.xaxis.set_major_locator(MaxNLocator(10))
			ax2.xaxis.set_minor_locator(MaxNLocator(20))
			ax2.xaxis.set_major_formatter(FormatStrFormatter("%g"))
			ax2.xaxis.set_minor_formatter(FormatStrFormatter(""))
			# fit straight line in log-log plot
			lf=[np.log10(ff) for ff in spectrum[2]]
			coeffs = np.polyfit(lf, [np.log10(ff) for ff in spectrum[0]], 1)
	                spfit = np.polyval(coeffs, lf)
			ax2.plot(spectrum[2], [pow(10., ff) for ff in spfit], "--", color="blue", label=r"$\alpha=%g$" % (coeffs[0]), alpha=0.8)
			ax2.errorbar(spectrum[2], spectrum[0], yerr=spectrum[1], color='darkgreen', fmt="o", alpha=0.9)
			plt.legend()
			plt.grid()

		plt.tight_layout()
		if opts.is_saveonly:
			plt.savefig(pngname)
		else:
			plt.show()
