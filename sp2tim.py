#!/usr/bin/env python
#
# converts the .singlepulse file (output from single_pulse_search.py script 
# from Presto) to the tim-file in either Tempo or Tempo2 format.
#
# Vlad Kondratiev (c)
#
import numpy as np
import os, os.path, stat, glob, sys, getopt
import infodata as inf

is_phase = False   # if True then calculate the phase of the pulses
inffile=""   # inf-file
spfile=""    # .singlepulse file
polycofile="polyco.dat"  # polyco-file
obscode="1"  # GBT
psrname=""
is_tempo2 = False  # if True then output tim-file is in Tempo2 format
extra=""
is_col6 = False    # if True then also add 6th column from the .singlepulse file (if it exists there)
                   # before extra-field

def usage (prg):
        """ prints the usage info about the current program
        """
        print "Program %s converts the .singlepulse file to tim-file\n" % (prg,)
        print "Usage: %s [-p, --polyco <polyco-file>] [-o, --obs <obscode>]\n\
                         [-s, --source <psrname>] [--tempo2] [--extra <field>]\n\
                         [--col6] [-h, --help] <inf-file> <singlepulse-file>\n\
         -o, --obs <obscode>        - set the observatory code (default = 1 for GBT)\n\
         -s, --source <psrname>     - pulsar name to be used for polyco and output tim-file\n\
         -p, --polyco <polyco-file> - use polyco-file to calculate the phase of the pulses\n\
         --tempo2                   - convert to Tempo2 format, default - Princeton\n\
                                      in Tempo2 format the first column is sigma, downfact, phase\n\
                                      (when -p option used), and extra field (when --extra option used)\n\
                                      separated by comma with no spaces\n\
         --extra <field>            - extra field to pass to tim-file. It will be added to the first column\n\
                                      in Tempo2 format, or be the last column in Princeton format\n\
         --col6                     - add 6th column from .singlepulse if it exists there\n\
         -h, --help                 - print this help\n" % (prg,)

def parsecmdline (prg, argv):
        """ parse the command line arguments
        """
        if not argv:
                usage (prg)
                sys.exit()
        else:
                try:
                        opts, args = getopt.getopt (argv, "hs:p:o:", ["help", "source=", "polyco=", "obs=", "tempo2", "extra=", "col6"])
                        for opt, arg in opts:
                                if opt in ("-h", "--help"):
                                        usage (prg)
                                        sys.exit()

                                if opt in ("-p", "--polyco"):
                                        global is_phase
					is_phase = True
					global polycofile
					polycofile = arg

                                if opt in ("-s", "--source"):
                                        global psrname
                                        psrname = arg

                                if opt in ("-o", "--obs"):
                                        global obscode
                                        obscode = arg

                                if opt in ("--tempo2"):
                                        global is_tempo2
					is_tempo2 = True

                                if opt in ("--extra"):
                                        global extra
                                        extra = arg

                                if opt in ("--col6"):
                                        global is_col6
					is_col6 = True

                        if not args:
                                print "inf-file and .singlepulse files are not given!\n"
                                sys.exit(2)
                        else:
				if len(args) < 2:
					print ".singlepulse file is not given!\n"
					sys.exit(2)
				else:
					global inffile
					inffile = args[0]
                                	global spfile
                                	spfile = args[1]

                except getopt.GetoptError:
                        print "Wrong option!"
                        usage (prg)
                        sys.exit(2)

if __name__=="__main__":
        parsecmdline (sys.argv[0].split("/")[-1], sys.argv[1:])
	# reading inf-file
	id = inf.infodata(inffile)

	if is_phase == True:
		import polycos as poly

	startmjd=id.epoch
	tres = id.dt
	unc = "%9f" % (tres * 1000000., )
	if psrname != "":
		source = psrname
	else:
		source = id.object
	cfreq = id.lofreq
	totbw = id.BW
	chanbw = id.chan_width
	freq = cfreq + totbw - chanbw     # central freq of the highest channel
                                          # prepdata is dedispersing to higher freq (the highest freq channel gets zero delay)

	# reading .singlepulse file
	dm, sigma, secs = np.loadtxt(spfile, usecols=(0,1,2), comments='#', dtype=float, unpack=True)
	offset, downfact = np.loadtxt(spfile, usecols=(3,4), comments='#', dtype=int, unpack=True)
	if is_col6:
		col6 = np.loadtxt(spfile, usecols=(5,5), comments='#', dtype=str, unpack=True)[0]
	toa = ["%.13f" % (startmjd + (offset[i] * tres)/86400.,) for i in np.arange(np.size(offset))]

	# calculating the phases of pulses 
	if is_phase == True:
		pid=poly.polycos(source, polycofile)
		phase=[pid.get_phs_and_freq(float(t.split(".")[0]), float("0." + t.split(".")[1]))[0] for t in toa]

	# writing the tim-file
	# Princeton format (+ additional extra field is for sigma)
	timfile=inffile.split(".inf")[0] + ".tim"
	if is_tempo2:  # output tim-file is in Tempo2 format
		if is_phase == True:
			if extra != "" or is_col6:
				if not is_col6:
					lines=["%s,%d,%f,%s   %8.3f   %s   %s   %s" % (str(sigma[i]), downfact[i], phase[i], extra, freq, str(toa[i]), str(unc), obscode) for i in np.arange(np.size(offset))]
				else:
					if extra == "":
						lines=["%s,%d,%f,%s   %8.3f   %s   %s   %s" % (str(sigma[i]), downfact[i], phase[i], col6[i], freq, str(toa[i]), str(unc), obscode) for i in np.arange(np.size(offset))]
					else:
						lines=["%s,%d,%f,%s,%s   %8.3f   %s   %s   %s" % (str(sigma[i]), downfact[i], phase[i], col6[i], extra, freq, str(toa[i]), str(unc), obscode) for i in np.arange(np.size(offset))]
			else:
				lines=["%s,%d,%f   %8.3f   %s   %s   %s" % (str(sigma[i]), downfact[i], phase[i], freq, str(toa[i]), str(unc), obscode) for i in np.arange(np.size(offset))]
		else:
			if extra != "" or is_col6:
				if not is_col6:
					lines=["%s,%d,,%s   %8.3f   %s   %s   %s" % (str(sigma[i]), downfact[i], extra, freq, str(toa[i]), str(unc), obscode) for i in np.arange(np.size(offset))]
				else:
					if extra == "":
						lines=["%s,%d,,%s   %8.3f   %s   %s   %s" % (str(sigma[i]), downfact[i], col6[i], freq, str(toa[i]), str(unc), obscode) for i in np.arange(np.size(offset))]
					else:
						lines=["%s,%d,,%s,%s   %8.3f   %s   %s   %s" % (str(sigma[i]), downfact[i], col6[i], extra, freq, str(toa[i]), str(unc), obscode) for i in np.arange(np.size(offset))]
			else:
				lines=["%s,%d   %8.3f   %s   %s   %s" % (str(sigma[i]), downfact[i], freq, str(toa[i]), str(unc), obscode) for i in np.arange(np.size(offset))]
	else: # Princeton format
		if is_phase == True:
			if extra != "" or is_col6:
				if not is_col6:
					lines=["%1s %-12s %8.3f %-20s%9s%10s   %s   %d   %f   %s" % (obscode, source, freq, str(toa[i]), str(unc), "0", str(sigma[i]), downfact[i], phase[i], extra) for i in np.arange(np.size(offset))]
				else:
					if extra == "":
						lines=["%1s %-12s %8.3f %-20s%9s%10s   %s   %d   %f   %s" % (obscode, source, freq, str(toa[i]), str(unc), "0", str(sigma[i]), downfact[i], phase[i], col6[i]) for i in np.arange(np.size(offset))]
					else:
						lines=["%1s %-12s %8.3f %-20s%9s%10s   %s   %d   %f   %s   %s" % (obscode, source, freq, str(toa[i]), str(unc), "0", str(sigma[i]), downfact[i], phase[i], col6[i], extra) for i in np.arange(np.size(offset))]
			else:
				lines=["%1s %-12s %8.3f %-20s%9s%10s   %s   %d   %f" % (obscode, source, freq, str(toa[i]), str(unc), "0", str(sigma[i]), downfact[i], phase[i]) for i in np.arange(np.size(offset))]
		else:
			if extra != "" or is_col6:
				if not is_col6:
					lines=["%1s %-12s %8.3f %-20s%9s%10s   %s   %d   0.0   %s" % (obscode, source, freq, str(toa[i]), str(unc), "0", str(sigma[i]), downfact[i], extra) for i in np.arange(np.size(offset))]
				else:
					if extra == "":
						lines=["%1s %-12s %8.3f %-20s%9s%10s   %s   %d   0.0   %s" % (obscode, source, freq, str(toa[i]), str(unc), "0", str(sigma[i]), downfact[i], col6[i]) for i in np.arange(np.size(offset))]
					else:
						lines=["%1s %-12s %8.3f %-20s%9s%10s   %s   %d   0.0   %s   %s" % (obscode, source, freq, str(toa[i]), str(unc), "0", str(sigma[i]), downfact[i], col6[i], extra) for i in np.arange(np.size(offset))]
			else:
				lines=["%1s %-12s %8.3f %-20s%9s%10s   %s   %d" % (obscode, source, freq, str(toa[i]), str(unc), "0", str(sigma[i]), downfact[i]) for i in np.arange(np.size(offset))]

	timp = open(timfile, 'w')	
	if is_tempo2:
		timp.write("FORMAT 1\n")
	np.savetxt(timp, np.transpose((lines)), fmt="%s")
	timp.close()
