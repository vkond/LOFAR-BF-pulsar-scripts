#!/usr/bin/env python
#
import numpy as np
import psrchive as pc
import os, os.path, stat, glob, sys, getopt
import scipy.stats as sc
import optparse as opt

# Main
if __name__=="__main__":

	# parsing command line options
	cmd = opt.OptionParser("Usage: %prog <arfile> [-h|--help] [OPTIONS]")
	cmd.add_option('-f', '--fscrunch', dest='fscr', metavar='FACTOR', help="Fscrunch factor, default: %default", default=1, type='int')
	cmd.add_option('-b', '--bscrunch', dest='bscr', metavar='FACTOR', help="Bscrunch factor, default: %default", default=1, type='int')
	cmd.add_option('-t', '--tscrunch', dest='tscr', metavar='FACTOR', help="Tscrunch factor, default: %default", default=1, type='int')
        cmd.add_option('-o', '--out', dest='outname', metavar='STRING', help="Name of the output file. Default is \
input name without .ar + \".flat.ar\". This option takes precedence over -m and -e", default="", type='str')
        cmd.add_option('-u', dest='outpath', metavar='STRING', help="Write output file to this location", default="", type='str')
        cmd.add_option('-e', dest='extension', metavar='STRING', help="Write output file with this extension. Take precedence over -m", default="", type='str')
        cmd.add_option('-m', dest='is_modify', action="store_true", help="Modify the original file on disk", default=False)

	# reading cmd options
	(opts, args) = cmd.parse_args()

        # check if input file is given
        if len(args) != 0:
                infile = args[0]
        else:
                cmd.print_usage()
                sys.exit(0)

        # forming the name of the output file
        if opts.outname == "":
                if opts.extension != "":
                        opts.outname = "%s%s" % (".".join(infile.split(".")[:-1]), opts.extension)
                        if opts.outpath != "":
                                opts.outname = "%s/%s" % (opts.outpath, opts.outname.split("/")[-1])
                elif opts.is_modify:
                        opts.outname = infile
                else:
                        opts.outname = "%s.flat.ar" % (infile.split(".ar")[0])
                        if opts.outpath != "":
                                opts.outname = "%s/%s" % (opts.outpath, opts.outname.split("/")[-1])
        else:
                if opts.outpath != "":
                        opts.outname = "%s/%s" % (opts.outpath, opts.outname.split("/")[-1])
                else:
                        # getting location of the input file
                        loc="/".join(infile.split("/")[:-1])
                        if loc == "": loc = "."
                        opts.outname = "%s/%s" % (loc, opts.outname.split("/")[-1])

	print "loading..."
       	raw = pc.Archive_load(infile)
	raw.pscrunch() # scrunching to total intensity
	if opts.fscr > 1: raw.fscrunch(opts.fscr)
	if opts.bscr > 1: raw.bscrunch(opts.bscr)
	if opts.tscr > 1: raw.tscrunch(opts.tscr)
	nchan = raw.get_nchan()
	nsubint = raw.get_nsubint()

        # normalizing...
	print "normalizing..."
	for ss in xrange(nsubint):
		for ff in xrange(nchan):  
			print "(%d, %d) out of (%d, %d)" % (ss, ff, nsubint, nchan)
			weight = raw.get_Integration(ss).get_weight(ff)
			prof = raw.get_Integration(ss).get_Profile(0,ff)
			if weight == 0.0:
				raw.get_Integration(ss).get_Profile(0,ff).scale(0)
				continue
			osm, osr = sc.probplot(prof.get_amps(), sparams=(), dist='norm', fit=0)
      			q_max = np.min(np.where(osm > 1.0))
      			q_min = np.max(np.where(osm < -1.0))
       			rms, mean = np.polyfit(osm[q_min:q_max], osr[q_min:q_max], 1)
       			prof.offset(-mean)
			if rms == 0.0: prof.scale(0)
      			else: prof.scale(1./rms)

	# saving (unload) data to *.flat.ar file
	raw.unload(opts.outname)
