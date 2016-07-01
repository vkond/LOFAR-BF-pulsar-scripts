#!/usr/bin/env python
#
# Script that determines the raw data on CEP2/4 clusters
# for a given ObsID, maps these data for the DRAGNET cluster
# and then copies there there as fast as possible
#
# Vlad Kondratiev, Aug 25, 2015 (c)
#
# Aug 30, 2015 - added new copying scheme: copying only one raw
#                file per dragnet node at once; and of course
#                copying the raw files for different dragnet
#                nodes in parallel. We have to see which scheme
#                is faster, and also possibly use either other
#                scp options, and maybe use rsync/bscp/.. instead
# Sep  1, 2015 - added the check if there are no files to copy
# Sep  2, 2015 - added removing empty directories in the case, when
#                there are no files to copy; also added catching
#                the SIGTERM signal and KeyboardInterruption, so
#                all open scp processed could be killed properly.
# Oct  6, 2015 - added option to run bbcp copy instead of scp
#                (much faster) 
# Nov 17, 2015 - added specific files' mapping for Ziggy's Pleunis
#                project (--ziggy). For this mapping, all polarisation
#                products for the same part go to the same dragnet
#                node. Also, it is the same for all TABs.
# Nov 18, 2015 - added the option to specify the first dragnet node to use;
#                fixed bug related to mkdir/rmdir only for those nodes
#                that were specified with --dragnet_nodes, --first_dragnet_node
# Nov 20, 2015 - added option to specify maximum files to copy at a time;
#                useful to limit the number of files when copying to only
#                one dragnet node and yet not too slow when using the option
#                of copying just one file at a time per dragnet node
# Dec  3, 2015 - added new bbcp command (still testing); added new options
#                --h5-only and --raw-only; changed code for copying .h5 files - 
#                now only copy one .h5 file at a time, it is slower, but at
#                least does not crush (do not know the reason for crashing);
#                also now use fully-qualified hostnames for dragnet nodes, otherwise
#                target is an NFS mount
# Dec 11, 2015 - rewrote the code to copy .h5 files the same time with .raw files;
#                added loop for many given ObsIDs
# Jun 30, 2016 - added option --ziggy2 to map files when there are twice as many
#                parts then for --ziggy with part's number exceeding potentially
#                the drgnode number. Mapping puts parts 0,20 on the drg01, parts 1,21
#                on the node drg02, etc.
#
import os, sys, glob, copy
import numpy as np
import optparse as opt
import re, time
import signal
import subprocess, shlex
from subprocess import PIPE, STDOUT, Popen

# cexec command to run. Using this mapfile makes keep mapping of the locus to be always the same
cexeccmd="cexec -f /etc/c3.conf.full"

cep2nodes=100 # we have 100 locus nodes
cep2node_stem="locus"

# Dragnet cluster setup
#drgnodes=23  # we have 23 nodes on Dragnet # now it is set as an option (default is 23)
drgnode_stem="drg"
network_suffix="-10g.online.lofar"
#network_suffix=""


# printing info about files on CEP2 cluster
def print_locus(rawfiles):
	nodes=rawfiles.keys()
	nfiles = 0
	print
	for node in sorted(nodes):
		print "%s [%d] -> " % (node, len(rawfiles[node]))
		nfiles += len(rawfiles[node])
		for ff in rawfiles[node]:
			print "   %s" % (ff)
	print "=" * 52
	print "Total number of files: %d (on %d nodes)" % (nfiles, len(nodes))
	return nfiles



# start subprocess
def start_and_go(cmd, shell=False):
	"""
	Execute the command 'cmd' after logging the command
	This function start the cmd and leaves the function
	returning the Popen object, it does not wait for process to finish
	"""
	status=1
	try:
		process = Popen(shlex.split(cmd), stdout=PIPE, stderr=STDOUT, shell=shell)
		status=process.returncode
		return (process, cmd)
	except Exception:
		print "Error! Job has crashed!\n%s\nStatus=%s" % (re.sub("\n", "\\\\n", cmd), status)
		sys.exit(1)

# waiting for one process to finish
# popen_tuple is the tuple of 2 elements: first - popen object, 2nd - cmd string
def waiting(popen_tuple):
	"""
	Waiting for process to finish
	"""
	try:
		(popen, cmd) = popen_tuple
		job_start = time.time()
		#print "Waiting for process to finish, pid=%d" % (popen.pid)
		(sout, serr) = popen.communicate()
		job_end = time.time()
		job_total_time = job_end - job_start
		#print "Process pid=%d has finished at UTC %s, status=%d, Waiting time: %.1f s (%.2f min)" % \
		#	(popen.pid, time.asctime(time.gmtime()), popen.returncode, job_total_time, job_total_time/60.)
		# if job is not successful
		if popen.poll() != 0 and popen.poll() is not None:
			print cmd
			print sout, serr
			raise Exception
	except Exception:
		print cmd
		print("Error! Job has crashed!\npid=%d, Status=%s" % (popen.pid, popen.returncode))
		sys.exit(1)

# waiting for processes to finish
def waiting_list(popen_list):
	"""
	Waiting for list of processes to finish
	"""
	try:
		job_start = time.time()
		for unit in popen_list:
			waiting(unit)
			finished_units = [u for u in popen_list if u[0].poll() is not None]
			for fu in finished_units:
				if fu[0].returncode != 0:
					print fu[1]
					print "Error! Job has crashed!\npid=%d, Status=%s" % (fu[0].pid, fu[0].returncode)
					raise Exception
		job_end = time.time()
		job_total_time = job_end - job_start
		#print "Processes have finished at UTC %s, Waiting time: %.1f s (%.2f min)" % \
		#	(time.asctime(time.gmtime()), job_total_time, job_total_time/60.)
	except Exception:
		print "Error! Job has crashed!\npids = %s" % (",".join(["%d" % (fu[0].pid) for fu in popen_list if fu[0].poll() is not None]))
		sys.exit(1)


# function that checks all processes in the list and kill them if they are still running
def kill(popen_list):
	print "Killing all open processes..."
	for popen in popen_list:
		if popen[0].poll() is None: # process is still running
			popen[0].kill()
			#if popen[0] != None: popen[0].communicate()
			#if popen[0] != None: popen[0].poll()


# main
if __name__=="__main__":

        cmdline=opt.OptionParser("Usage: %prog [options] <ObsID 1> <ObsID 2> ... <ObsID n>")
        cmdline.add_option('-d', '--data', dest='datadir', metavar='DATADIR',
                           help="Datadir in the DRAGNET destination (/data1 or /data2), default = %default", default="/data1", type='str')
        cmdline.add_option('--bbcp', action="store_true", dest='is_bbcp',
                           help="use bbcp to copy the data", default=False)
	cmdline.add_option('--one-file-per-dragnet-node-at-once', action="store_true", dest='is_one_file_per_drg', help="turn on the scheme of copying \
only one file for each dragnet node at once", default=False)
	cmdline.add_option('--h5-only', action="store_true", dest='is_h5_only', help="copy only .h5 files", default=False)
	cmdline.add_option('--raw-only', action="store_true", dest='is_raw_only', help="copy only .raw files", default=False)
	cmdline.add_option('-m', '--max-number-of-files-at-a-time', dest='max_number_of_files', metavar='N', help="set the maximum number of files \
to copy at a time (or less), default - all", default=-1, type='int')
        cmdline.add_option('--ziggy', action="store_true", dest='is_ziggy',
                           help="to copy complex-voltage data putting all polarisation products of the same part to the same dragnet node, \
same for all the TABs", default=False)
        cmdline.add_option('--ziggy2', action="store_true", dest='is_ziggy2',
                           help="to copy complex-voltage data putting all polarisation products of the same part to the same dragnet node, \
same for all the TABs. Same as --ziggy, but used when there are twice as many parts, i.e. 2 parts per node are stored", default=False)
        cmdline.add_option('--dragnet-nodes', dest='drgnodes', metavar='N',
                           help="Number of dragnet nodes to use, default = %default", default=23, type='int')
        cmdline.add_option('--first-dragnet-node', dest='first_dragnet_node', metavar='N',
                           help="First dragnet node to use, default = %default", default=1, type='int')
        (opts, args) = cmdline.parse_args()

	if len(args) == 0:
		print "Error: no ObsID is given!"
		sys.exit(1)

	# setting the number of dragnet nodes
	drgnodes = opts.drgnodes
	first_dragnet_node = opts.first_dragnet_node

	# copy command to use
	scp_copy_command = "scp -B -c arcfour"
	if opts.is_bbcp:
		#copy_command="/home/kondratiev/bbcp/bin/amd64_linux31/bbcp -P 2 -s 20 -w 2M"
		#copy_command="/home/kondratiev/bbcp/bin/amd64_linux31/bbcp -s 1 -w 2M"
		# command from Alexander
		#copy_command="/home/kondratiev/bbcp/bin/amd64_linux31/bbcp -y dd -e -s 8 -B 4M -r -g -@ follow -v -- "
		copy_command="/home/kondratiev/bbcp/bin/amd64_linux31/bbcp -a -y dd -e -s 8 -B 4M -r -g -@ follow -v -- "
	else:
		copy_command=scp_copy_command

	locus_nodes=["%s%03d" % (cep2node_stem, num+1) for num in range(cep2nodes)]
	cexec_nodes={}
	for num in range(cep2nodes):
		key="%s%03d" % (cep2node_stem, num+1)
		cexec_nodes[key] = "%s:%d" % (cep2node_stem, num)
	# we are using 10g names here
	dragnet_nodes=["%s%02d%s" % (drgnode_stem, num+1, network_suffix) for num in range(drgnodes)]

	# first define signal handler for SIGTERM (if system or user sends kill signal), "kill -9" can _not_ be caught 
	def sigterm_handler(signum, frame):
		raise Exception("SIGTERM signal got caught!")
	signal.signal(signal.SIGTERM, sigterm_handler)

	try:
		popens=[]
		all_obsid_start = time.time()

		# loop over Obsids
		for obsid in args:

			print "*" * 20
			print "* ObsID = %s" % (obsid)
			print "*" * 20

			rawfiles={}  # dictionary that keeps info about all raw files
        	        	     # key - locus node, value - list of rawfiles on this node (with full path)

			# first, we have to create /data?/obsid directories on dragnet nodes
			total_start = time.time()
			job1_start = time.time()
			print "Creating %s/%s directories on Dragnet nodes..." % (opts.datadir, obsid)
			popens=[]
			for node_num in xrange(first_dragnet_node, first_dragnet_node + drgnodes, 1):
				cmd="ssh %s%02d%s 'mkdir -m 775 -p %s/%s'" % (drgnode_stem, node_num, network_suffix, opts.datadir, obsid)
				popen = start_and_go(cmd)
				popens.append(popen)

			# Determining where the raw data are....
        		# forming string with all locus nodes needed to check in one cexec command
			# here we are using only nodes that are alive
			job_start = time.time()
			print
			print "Collecting info about files for the ObsID=%s..." % (obsid)
			print
			cexeclocus=cexec_nodes[locus_nodes[0]]
			for s in locus_nodes[1:]:
				cexeclocus += ",%s" % (cexec_nodes[s].split(":")[1])

			cmd="%s %s 'ls -1L /data/%s/*_bf.raw' | grep -v such | grep -v match | grep -v xauth | grep -v connect | grep -v closed | egrep -v \'\\*\\*\\*\\*\\*\'" % (cexeccmd, cexeclocus, obsid)
       			cexec_output=[line.strip() for line in os.popen(cmd).readlines()]
			for l in xrange(len(cexec_output)):
				if re.match("^-----", cexec_output[l]) is not None:
					loc=cexec_output[l].split(" ")[1]
				else:
					# first we are checking that this file (exactly the same with full path) DO NOT already present
					# in the current dictionary 'rawfiles'.
					rawfile_exist=False
					for cloc in rawfiles:
						if cexec_output[l] in rawfiles[cloc]:
							rawfile_exist=True
							break
					# adding the file if it's not already added
					if not rawfile_exist:
						if loc in rawfiles: rawfiles[loc].append(cexec_output[l])
						else: rawfiles[loc]=[cexec_output[l]]	

			# printing info about files on CEP2 for a given ObsID
			# and getting total number of files
			nfiles = print_locus(rawfiles)

			job_end = time.time()
			job_time = job_end - job_start
			print("Running time: %.1f s (%.2f min)" % (job_time, job_time/60.))

			# checking the given value of max number of files to copy
			if opts.max_number_of_files <= 0 or opts.max_number_of_files >= nfiles:
				opts.max_number_of_files = -1

			# check if no files are there
			if nfiles == 0:
				print "\nThere are no files to copy! Skipping this ObsID=%s..." % (obsid)
				job_start = time.time()
				print "  Removing %s/%s empty directories on Dragnet nodes..." % (opts.datadir, obsid)
				popens=[]
				for node_num in xrange(first_dragnet_node, first_dragnet_node + drgnodes, 1):
					cmd="ssh %s%02d%s 'rmdir %s/%s'" % (drgnode_stem, node_num, network_suffix, opts.datadir, obsid)
					popen = start_and_go(cmd)
					popens.append(popen)
				waiting_list(popens)
				job_end = time.time()
				job_time = job_end - job_start
				print("  Running time (rmdir %s/%s): %.1f s (%.2f min)" % (opts.datadir, obsid, job_time, job_time/60.))
				total_end = time.time()
				total_time = total_end - total_start
				print("  Total running time for ObsID=%s: %.1f s (%.2f min)" % (obsid, total_time, total_time/60.))
				continue

			nmin=int(nfiles/drgnodes) # min number of files per Dragnet nodes
			nmaxnodes=nfiles - nmin*drgnodes # first 'nmaxnodes' will have nmin+1 files

			#
			# Files mapping
			#
			print "Files mapping between the clusters..."
			cep2files={} # key - the file name, value - tuple of locus node and dragnet node
			drgmap={}  # key - dragnet node, value - list of tuples of (locus node, cep2file)

			# default mapping
			if not opts.is_ziggy and not opts.is_ziggy2:
				dnode=1 # index for dragnet node
				curfile=0 # index for a file number
				for node in sorted(rawfiles.keys()):
					for ff in rawfiles[node]:
						drgnode = "%s%02d%s" % (drgnode_stem, first_dragnet_node - 1 + dnode, network_suffix)
						cep2files[ff]=(node, drgnode)
						# forming another map (to be used when option --one-file-per-dragnet-node-at-once is used)
						if drgnode in drgmap:
							drgmap[drgnode].append((node, ff))
						else:
							drgmap[drgnode] = [(node, ff)]
						#print "%s    %s  ->  %ss" % (ff, node, drgnode)
						curfile += 1
						if nmaxnodes == 0:
							if curfile == nmin:
								curfile=0
								dnode += 1
						else:
							if dnode <= nmaxnodes:
								if curfile == nmin+1:
									curfile=0
									dnode += 1
							else:
								if curfile == nmin:
									curfile=0
									dnode += 1
			# Mapping for the specific Ziggy's case
			else:
				for node in sorted(rawfiles.keys()):
					for ff in rawfiles[node]:
						# based on the name of the file (specifically the part number) we select the dragnet node
						# For Ziggy's project there should be 20 parts (from 0 to 19)
						# Part 0 goes to drg01, etc
						part=int(ff.split("/")[-1].split("_P")[-1].split("_")[0])
						if opts.is_ziggy:
							drgnode = "%s%02d%s" % (drgnode_stem, first_dragnet_node - 1 + part + 1, network_suffix)
						if opts.is_ziggy2:
							if part > 19:
								drgnode = "%s%02d%s" % (drgnode_stem, first_dragnet_node - 1 + part%20 + 1, network_suffix)
							else:
								drgnode = "%s%02d%s" % (drgnode_stem, first_dragnet_node - 1 + part + 1, network_suffix)
						cep2files[ff]=(node, drgnode)
						if drgnode in drgmap:
							drgmap[drgnode].append((node, ff))
						else:
							drgmap[drgnode] = [(node, ff)]
	
			# Now, copying...
	
			# now waiting for the end of processes that should create directories on Dragnet nodes
			waiting_list(popens)
			job_end = time.time()
			job_time = job_end - job1_start
			print("Running time (mkdir -p %s/%s): %.1f s (%.2f min)" % (opts.datadir, obsid, job_time, job_time/60.))

			# Copying all .raw and .h5 files
			job_start = time.time()
			if opts.is_h5_only:
				print "Copying .h5 files..."
			elif opts.is_raw_only:
				print "Copying .raw files..."
			else:
				print "Copying .raw and .h5 files..."

			# copying all raw|h5 files at once
			if not opts.is_one_file_per_drg:
				if opts.is_h5_only:
					print "   Copying all h5 files at once..."
				elif opts.is_raw_only:
					print "   Copying all raw files at once..."
				else:
					print "   Copying all raw and h5 files at once..."
				
				popens=[]
				# copying everything at once
				if opts.max_number_of_files == -1:
					for ff in cep2files.keys(): # loop over the files
						if opts.is_raw_only or (not opts.is_h5_only and not opts.is_raw_only):
							cmd="%s %s:%s %s:%s/%s" % (copy_command, cep2files[ff][0], ff, cep2files[ff][1], opts.datadir, obsid)
							popen = start_and_go(cmd)
							popens.append(popen)
						if opts.is_h5_only or (not opts.is_h5_only and not opts.is_raw_only):
							cmd="%s %s:%s.h5 %s:%s/%s" % (copy_command, cep2files[ff][0], ff.split(".raw")[0], cep2files[ff][1], opts.datadir, obsid)
							popen = start_and_go(cmd)
							popens.append(popen)
					waiting_list(popens)
				else:
					cepfiles=cep2files.keys() # list of cep2 files to copy
					while len(popens) <= opts.max_number_of_files:
						if opts.is_raw_only or (not opts.is_h5_only and not opts.is_raw_only):
							cmd="%s %s:%s %s:%s/%s" % (copy_command, cep2files[cepfiles[-1]][0], cepfiles[-1], cep2files[cepfiles[-1]][1], opts.datadir, obsid)
							popen = start_and_go(cmd)
							popens.append(popen)
						if opts.is_h5_only or (not opts.is_h5_only and not opts.is_raw_only):
							cmd="%s %s:%s.h5 %s:%s/%s" % (copy_command, cep2files[cepfiles[-1]][0], cepfiles[-1].split(".raw")[0], cep2files[cepfiles[-1]][1], opts.datadir, obsid)
							popen = start_and_go(cmd)
							popens.append(popen)
						cepfiles.pop()
						if len(cepfiles) == 0: break
					while True:
						if len(popens) == 0: break
						ii_to_delete = -1
						for ii in xrange(len(popens)):
							popen = popens[ii][0]
							cmd = popens[ii][1]
							(sout, serr) = popen.communicate()
							if popen.poll() is not None:
								if popen.returncode != 0:
									print cmd
									print sout, serr
									print "Error! Job has crashed!\npid=%d, Status=%s" % (popen.pid, popen.returncode)	
									raise Exception
								# starting new scp to the same drg node
								if len(cepfiles) != 0:
									if opts.is_raw_only or (not opts.is_h5_only and not opts.is_raw_only):
										cmd="%s %s:%s %s:%s/%s" % (copy_command, cep2files[cepfiles[-1]][0], cepfiles[-1], cep2files[cepfiles[-1]][1], opts.datadir, obsid)
										newpopen = start_and_go(cmd)
										popens.append(newpopen)
									if opts.is_h5_only or (not opts.is_h5_only and not opts.is_raw_only):
										cmd="%s %s:%s.h5 %s:%s/%s" % (copy_command, cep2files[cepfiles[-1]][0], cepfiles[-1].split(".raw")[0], cep2files[cepfiles[-1]][1], opts.datadir, obsid)
										newpopen = start_and_go(cmd)
										popens.append(newpopen)
									cepfiles.pop()
								# removing this finished process from the popens list
								ii_to_delete = ii
								break
						if ii_to_delete != -1: del popens[ii_to_delete]

			else: # copying only one raw file per dragnet node at once
				if opts.is_h5_only:
					print "   Copying only one h5 file per dragnet node at once..."
				elif opts.is_raw_only:
					print "   Copying only one raw file per dragnet node at once..."
				else:
					print "   Copying only one raw and h5 file per dragnet node at once..."
				popens=[]
				for ff in drgmap.keys(): # loop over dragnet nodes
					(node, cep2file) = drgmap[ff][-1]
					if opts.is_raw_only or (not opts.is_h5_only and not opts.is_raw_only):
						cmd="%s %s:%s %s:%s/%s" % (copy_command, node, cep2file, ff, opts.datadir, obsid)
						popen = start_and_go(cmd)
						popen = popen + (ff,)
						popens.append(popen)
					if opts.is_h5_only or (not opts.is_h5_only and not opts.is_raw_only):
						cmd="%s %s:%s.h5 %s:%s/%s" % (copy_command, node, cep2file.split(".raw")[0], ff, opts.datadir, obsid)
						popen = start_and_go(cmd)
						popen = popen + (ff,)
						popens.append(popen)
				while True:
					if len(popens) == 0: break
					for (popen, cmd, drg) in popens:
						(sout, serr) = popen.communicate()
						if popen.poll() is not None:
							if popen.returncode != 0:
								print cmd
								print sout, serr
								print "Error! Job has crashed!\npid=%d, Status=%s" % (popen.pid, popen.returncode)	
								raise Exception
							# removing last tuple from the drgmap as it is finished
							drgmap[drg].pop()
							# removing this finished process from the popens list
							for ii in xrange(len(popens)):
								if popens[ii][2] == drg:
									del popens[ii]
									break
							# starting new scp to the same drg node
							if len(drgmap[drg]) != 0:
								(node, cep2file) = drgmap[drg][-1]
								if opts.is_raw_only or (not opts.is_h5_only and not opts.is_raw_only):
									cmd="%s %s:%s %s:%s/%s" % (copy_command, node, cep2file, drg, opts.datadir, obsid)
									newpopen = start_and_go(cmd)
									newpopen = newpopen + (drg,)
									popens.append(newpopen)
								if opts.is_h5_only or (not opts.is_h5_only and not opts.is_raw_only):
									cmd="%s %s:%s.h5 %s:%s/%s" % (copy_command, node, cep2file.split(".raw")[0], drg, opts.datadir, obsid)
									newpopen = start_and_go(cmd)
									newpopen = newpopen + (drg,)
									popens.append(newpopen)
							break

			job_end = time.time()
			job_time = job_end - job_start
			print("Running time: %.1f s (%.2f min)" % (job_time, job_time/60.))

			total_end = time.time()
			total_time = total_end - total_start
			print("Total running time for ObsID=%s: %.1f s (%.2f min)" % (obsid, total_time, total_time/60.))

		all_obsid_end = time.time()
		all_obsid_time = all_obsid_end - all_obsid_start
		print("Total running time for all ObsIDs: %.1f s (%.2f min)" % (all_obsid_time, all_obsid_time/60.))

	except Exception as ex:
		print "Exception has been caught: %s" % (", ".join(ex.args))
		kill(popens)
		sys.exit(1)
	except KeyboardInterrupt as ex:
		print "User interruption..."
		kill(popens)
		sys.exit(1)
