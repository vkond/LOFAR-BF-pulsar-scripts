#!/usr/bin/env python
#
# Reads the ascii file with info about flagged tiles per station per date
# that was e.g. created by getState.py script that reads hardware_states_latest.txt
# file. Calculates the overall flagged tiles fraction to be used for flux-calibration
# software.
# 
# Vlad Kondratiev (c) - 23.06.2016
#
import os, sys
import re
import numpy as np
import optparse as opt
import h5py
import time
from datetime import *

# Main
if __name__=="__main__":
        #
        # Parsing the command line options
        #
        usage = "Usage: %prog [options] [.h5 | -t]"
        cmdline = opt.OptionParser(usage)
        cmdline.add_option('-f', '--file', dest='infofile', metavar='FILE', help="Ascii file with info about flagged tiles per station per date", default="", type='str')
        cmdline.add_option('-t', '--time', dest='startTime', metavar='YYYY-MM-DD hh:mm:ss', help="specify the UT start date/time of an observation", default="", type='str')
        cmdline.add_option('-a', '--antenna', dest='antenna', metavar='HBA | LBA', help="specify either HBA or LBA antenna", default="", type='str')
        cmdline.add_option('-s', '--stations', dest='stations', metavar='CS001HBA0,CS017HBA1,...', help="specify the comma-separated list of stations used", default="", type='str')
	cmdline.add_option('-v', '--verbose', dest='is_verbose', action="store_true", help="Verbose output, print extra info", default=False)

        # reading cmd options
        (opts,args) = cmdline.parse_args()

        # check if input file is given
        if len(args) == 0 and (opts.startTime == "" or opts.antenna == "" or opts.stations == ""):
                cmdline.print_usage()
                sys.exit(0)

	if opts.infofile == "":
		print "No ascii file is given with --file option to read info about flagged tiles!"
		sys.exit(1)

	# reading info file
	f = open(opts.infofile, 'r')
	# ignoring also comments and empty lines
	comments_and_empty=re.compile(r"(^\s*#+.*$)|(^\s*$)")
	lines = [ff for ff in f.read().splitlines() if comments_and_empty.search(ff) is None]
	f.close()

	if len(args) != 0:
	        # reading input .h5 file
        	h5file = args[0]
		try:
			f5 = h5py.File(h5file, 'r')

			startdate = f5.attrs['OBSERVATION_START_UTC']
			startdate = startdate.split("T")[0]
			bandFilter = f5.attrs['FILTER_SELECTION']
			antenna = bandFilter.split("_")[0]
			stations=f5.attrs['OBSERVATION_STATIONS_LIST']
			f5.close()
		except:
			print "Can't open input .h5 file: %s" % (h5file)
			sys.exit(1)

	# then we will be getting necessary info from the command line	
	else:
		try:
			startdate = opts.startTime.split()[0]
			antenna = opts.antenna
			stations=opts.stations.split(",")
		except:
			print "Some of the input parameters are wrong..."
			sys.exit(1)
		
	try:

		# getting the list of stations/flagged-tiles for this date
		res=[ii for ii in lines if "%s" % (startdate) in ii]
		if len(res) == 0:
			print "No info found for the observation date: %s" % (startdate)
			sys.exit(1)


		# loop over used stations to collect info about flagged tiles
		total_tiles = 0
		nflagged = 0
		st_worst = []
		worst = -1
		for st in stations:
			out=[ii for ii in res if "%s" % (st) in ii]
			if len(out) == 0:
				print "Missing info about flagged tiles for the station: %s" % (st)
				sys.exit(1)
			if len(out) > 1:
				if out[0] != out[1]:
					print "There are more than 1 line with flagged tiles for the station: %s" % (st)
					sys.exit(1)
			(dd, tt, ss, ntiles) = out[0].split()
			if int(ntiles) == worst:
				st_worst.append(ss)
			if int(ntiles) > worst:
				worst = int(ntiles)
				st_worst = [ss]
			nflagged += int(ntiles)
			
			
		ncorestations = len([s for s in stations if s[0:2] == "CS"])
		# because in the list for HBA there are sub-stations
		if antenna == "HBA": 
			if opts.is_verbose:
				print "Number of core sub-stations: %d" % (ncorestations)
			ncorestations /= 2
			total_tiles=24*len([s for s in stations if s[0:2] == "CS"]) + 48*len([s for s in stations if s[0:2] == "RS"]) + \
				96*len([s for s in stations if s[0:2] != "CS" and s[0:2] != "RS"])
			fraction = float(nflagged) / total_tiles
			if any([s for s in st_worst if s[0:2] == "CS"]):
				worst_fraction = float(worst) / 24.
			elif any([s for s in st_worst if s[0:2] == "RS"]):
				worst_fraction = float(worst) / 48.
			else:
				worst_fraction = float(worst) / 96.
		else:
			total_tiles=48*len([s for s in stations if s[0:2] == "RS" or s[0:2] == "CS"]) + \
				96*len([s for s in stations if s[0:2] != "CS" and s[0:2] != "RS"])
			fraction = float(nflagged) / total_tiles
			if any([s for s in st_worst if s[0:2] == "CS" or s[0:2] == "RS"]):
				worst_fraction = float(worst) / 48.
			else:
				worst_fraction = float(worst) / 96.
		if opts.is_verbose:
			print "Number of core stations: %d" % (ncorestations)
			print "Number of total flagged tiles: %d" % (nflagged)
			print "Fraction of bad tiles: %g [%g%%]" % (fraction, fraction * 100.)
			if worst != 0:
				print "Worst (sub-)station(s) [#tiles=%d, fraction=%g%%]: %s" % (worst, worst_fraction * 100., ",".join(st_worst))
		else:
			print "%g %d %d" % (fraction, nflagged, total_tiles)
	except:
		print "Crashed..."
		sys.exit(1)

