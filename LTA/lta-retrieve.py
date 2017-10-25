#!/usr/bin/env python
#
# Script to retrieve the processed data of LOFAR PWG
# from the LTA
#
# Vlad Kondratiev, Feb 12, 2014 (c)
#
# Jul  4, 2015 - checking srm link now to see if it points
#                to Julich, and then select proper html-prefix
#                for download link with wget, either for SARA
#                or Julich
# Jul 21, 2015 - added new option --csvfile to provide
#                full csv-file from lta-query.py
# Jul 22, 2015 - added three more options:
#                --sap, --tab, and --part to retrieve only
#                tarballs for a specific SAP, TAB, or PART
# Jul 23, 2015 - added --query option and --project option
#                to query directly LTA for the specific data
#                for a given project. 
# Aug 18, 2015 - retrieving data first to tempdir and then
#                moving one level up. This is done to avoid
#                clash with the same files (like parfile)
#                in different tarballs. 
# Mar 17, 2016 - added downloadng from Poznan (PL) LTA site
# Jul 27, 2016 - added --stage-only option
# Nov 26, 2016 - added query to retrieve the data from projects
#                that became public and you were not the member
# Nov 26, 2016 - added logging handling when communicating with the
#                database; there is no logfile now when there are
#                errors; output on the screen can be turned off
#                with the cmdline option -q
# Mar 20, 2017 - added new options: --skip-summary and --skip-staging
# Oct 25, 2017 - domain name sara.nl has changed to surfsara.nl
#
import numpy as np
import time
import xmlrpclib
import os, os.path, sys, re
import optparse as opt
import subprocess, shlex
from subprocess import PIPE, STDOUT, Popen

# html prefix to wrap up srm links
sara_html_prefix="https://lofar-download.grid.surfsara.nl/lofigrid/SRMFifoGet.py?surl="
juelich_html_prefix="https://lofar-download.fz-juelich.de/webserver-lofar/SRMFifoGet.py?surl="
poznan_html_prefix="https://lta-download.lofar.psnc.pl/lofigrid/SRMFifoGet.py?surl="

# grid directory
griddir="/home/kondratiev/Lofar/grid"

#
# SQL "super-query" that allows to collect data from different datatypes, namely:
# beamformed data products, pipeline (pulp) data products, pipeline summary data products,
# and unspecified data products (plenty for pulsar data)
# It uses UNION to merge queries together
# The last column is ObservationID, however due to LTA metadata  is completely
# messed up it is often the same as Pipeline ID.
#
superquery="SELECT fo.FILENAME, fo.FILESIZE, fo.CREATION_DATE, fo.URI, obs.ObservationID \
FROM AWOPER.BeamformedDataProduct dp, \
     AWOPER.FileObject fo, \
     AWOPER.Observation obs, \
     AWOPER.FORMEDDATAPRODUCT$OBSERVATIONS dp_obs \
WHERE fo.data_object = dp.object_id \
  AND fo.isValid > 0 \
  AND dp.isValid > 0 \
  AND dp.\"+PROJECT\" = SYS_CONTEXT('AWCONTEXT','PROJECTID') \
  AND dp_obs.object_id = dp.object_id \
  AND dp_obs.column_value = obs.object_id \
UNION \
SELECT fo.FILENAME, fo.FILESIZE, fo.CREATION_DATE, fo.URI, obs.ObservationID \
FROM AWOPER.PulpDataProduct dp, \
     AWOPER.FileObject fo, \
     AWOPER.Observation obs, \
     AWOPER.PulpDataProduct$Observations dp_obs \
WHERE fo.data_object = dp.object_id \
  AND fo.isValid > 0 \
  AND dp.isValid > 0 \
  AND dp.\"+PROJECT\" = SYS_CONTEXT('AWCONTEXT','PROJECTID') \
  AND dp_obs.object_id = dp.object_id \
  AND dp_obs.column_value = obs.object_id \
UNION \
SELECT fo.FILENAME, fo.FILESIZE, fo.CREATION_DATE, fo.URI, obs.ObservationID \
FROM AWOPER.PulpSummaryDataProduct dp, \
     AWOPER.FileObject fo, \
     AWOPER.Observation obs, \
     AWOPER.UMMARYDATAPRODUCT$OBSERVATIONS dp_obs \
WHERE fo.data_object = dp.object_id \
  AND fo.isValid > 0 \
  AND dp.isValid > 0 \
  AND dp.\"+PROJECT\" = SYS_CONTEXT('AWCONTEXT','PROJECTID') \
  AND dp_obs.object_id = dp.object_id \
  AND dp_obs.column_value = obs.object_id \
UNION \
SELECT fo.FILENAME, fo.FILESIZE, fo.CREATION_DATE, fo.URI, pr.ObservationID \
FROM AWOPER.UnspecifiedDataProduct dp, \
     AWOPER.FileObject fo, \
     AWOPER.UnspecifiedProcess pr \
WHERE fo.data_object = dp.object_id \
  AND fo.isValid > 0 \
  AND dp.isValid > 0 \
  AND dp.\"+PROJECT\" = SYS_CONTEXT('AWCONTEXT','PROJECTID') \
  AND dp.UnspecifiedProcess = pr.object_id"

#
# This superquery is for the case when data already became public and you are not the original
# member of this project. Then, these data can be retrieved from the LTA by specifying "All public data" as project
# and then choosing appropriate project within these All data
# Luciano's email:
# Hi Vlad,
# if you select as project "LC3_024" (for example), the system will always check if you are a member or not and will not let you get data.
# You need to select "All public data" as your project, then in the advanced window for Observation search, you can select the project within "All public data".
# 
# To do it from Python:
# Hi Vlad,
# Willem Jan clarified that when you use with the python scripts, you need to set project "ALL" first.
#
# Then you need to construct a query and, if you know the project the data is in, limit 
# the query to that project data. An example:
# 
# python
# context.set_project('ALL')
# q = CorrelatedDataProduct.select_all()
# q &= q.project_only('LC0_017')
# print(len(q))
# -> 1800
# 
# Before running this superquery, you have set context as "ALL" and then 
# substitute %s with the corresponding project code
#
public_superquery="SELECT fo.FILENAME, fo.FILESIZE, fo.CREATION_DATE, fo.URI, obs.ObservationID \
FROM AWOPER.BeamformedDataProduct dp, \
     AWOPER.FileObject fo, \
     AWOPER.Observation obs, \
     AWOPER.FORMEDDATAPRODUCT$OBSERVATIONS dp_obs \
WHERE fo.data_object = dp.object_id \
  AND fo.isValid > 0 \
  AND dp.isValid > 0 \
  AND dp.\"+PROJECT\" = (SELECT ID FROM AWOPER.AWEPROJECTS WHERE NAME='%s') \
  AND dp_obs.object_id = dp.object_id \
  AND dp_obs.column_value = obs.object_id \
UNION \
SELECT fo.FILENAME, fo.FILESIZE, fo.CREATION_DATE, fo.URI, obs.ObservationID \
FROM AWOPER.PulpDataProduct dp, \
     AWOPER.FileObject fo, \
     AWOPER.Observation obs, \
     AWOPER.PulpDataProduct$Observations dp_obs \
WHERE fo.data_object = dp.object_id \
  AND fo.isValid > 0 \
  AND dp.isValid > 0 \
  AND dp.\"+PROJECT\" = (SELECT ID FROM AWOPER.AWEPROJECTS WHERE NAME='%s') \
  AND dp_obs.object_id = dp.object_id \
  AND dp_obs.column_value = obs.object_id \
UNION \
SELECT fo.FILENAME, fo.FILESIZE, fo.CREATION_DATE, fo.URI, obs.ObservationID \
FROM AWOPER.PulpSummaryDataProduct dp, \
     AWOPER.FileObject fo, \
     AWOPER.Observation obs, \
     AWOPER.UMMARYDATAPRODUCT$OBSERVATIONS dp_obs \
WHERE fo.data_object = dp.object_id \
  AND fo.isValid > 0 \
  AND dp.isValid > 0 \
  AND dp.\"+PROJECT\" = (SELECT ID FROM AWOPER.AWEPROJECTS WHERE NAME='%s') \
  AND dp_obs.object_id = dp.object_id \
  AND dp_obs.column_value = obs.object_id \
UNION \
SELECT fo.FILENAME, fo.FILESIZE, fo.CREATION_DATE, fo.URI, pr.ObservationID \
FROM AWOPER.UnspecifiedDataProduct dp, \
     AWOPER.FileObject fo, \
     AWOPER.UnspecifiedProcess pr \
WHERE fo.data_object = dp.object_id \
  AND fo.isValid > 0 \
  AND dp.isValid > 0 \
  AND dp.\"+PROJECT\" = (SELECT ID FROM AWOPER.AWEPROJECTS WHERE NAME='%s') \
  AND dp.UnspecifiedProcess = pr.object_id"

# funtion to download files from LTA
def retrieve(data, rf=""):
	# downloading
	curdir=os.getcwd()
	job_start = time.time()
	print "Retrieving..."
	for ii, ff in enumerate(data[0]):
		print " %s..." % (ff),
		sys.stdout.flush()
		srmlink=data[4][ii]
		if "_locus" in ff:
			obsid=ff.split("_locus")[-1].split("_")[1]
		else:
			obsid=ff.split("_")[0]
		cmd="mkdir -p %s" % (obsid)
		os.system(cmd)
		os.chdir(obsid)
		try:
			retrieve_start = time.time()
			retrieve_end = 0
			extract_end = -1

			# first, we create temporary directory, where we do download, and extraction
			# and only then we move the downloaded stuff above
			workdir = os.getcwd()
			tp = Popen(shlex.split("mktemp -d -p %s .XXXXXX" % (workdir)), stdout=PIPE, stderr=STDOUT, cwd=workdir)
			(tpout, tperr) = tp.communicate()
			tempdir = tpout.strip()
			if tperr is not None: raise Exception
			# go to this temp dir
			os.chdir(tempdir)

			if "juelich.de" in srmlink:
				cmd="wget -c -O %s -t 0 %s%s" % (ff, juelich_html_prefix, srmlink)
			elif "psnc.pl" in srmlink:
				cmd="wget -c -O %s -t 0 %s%s" % (ff, poznan_html_prefix, srmlink)
			else:
				cmd="wget -c -O %s -t 0 %s%s" % (ff, sara_html_prefix, srmlink)
			p = Popen(shlex.split(cmd), stdout=PIPE, stderr=STDOUT)
			if opts.is_log:
				while True:
					line = p.stdout.readline()
					if not line: break
					print line.rstrip()
			p.communicate()

			if p.poll() !=0: raise Exception
			retrieve_end = time.time()
			retrieve_time = retrieve_end - retrieve_start
			print "done"
			# get file size in GB
			fsize=os.path.getsize(ff)/1000./1000./1000.
			print "\tTarball size (GB): %s (LTA), %.1f (Retrieved)" % (data[3][ii], fsize)
			print "\tRetrieve time: %.1f s (%.1f min)" % (retrieve_time, retrieve_time/60.)
			extract_start = time.time()
			extract_end = 0
			print " %s extracting..." % (ff),
			sys.stdout.flush()
			# checking if file is gzipped or not
			if re.search(".tgz$", ff) or re.search(".gz$", ff):
				cmd="tar xvfz %s" % (ff)
			else:
				cmd="tar xvf %s" % (ff)
			p = Popen(shlex.split(cmd), stdout=PIPE, stderr=STDOUT)
			p.communicate()
			if p.poll() !=0: raise Exception
			extract_end = time.time()
			extract_time = extract_end - extract_start
			print "done"
			print "\tExtract time: %.1f s (%.1f min)" % (extract_time, extract_time/60.)
			print " %s removing..." % (ff),
			sys.stdout.flush()
			cmd="rm -f %s" % (ff)
			os.system(cmd)
			print "done"
			if re.search("%s_red" % (obsid), ff):	
				locus=ff.split("ARCHIVE_")[-1].split("_")[0]
				is_is=""
				if re.search("%s_redIS" % (obsid), ff): 
					is_is="IS"
				cmd="mv %s_red%s %s_red%s_%s" % (obsid, is_is, obsid, is_is, locus)
				print " renaming %s_red%s  to  %s_red%s_%s" % (obsid, is_is, obsid, is_is, locus)
				os.system(cmd)

			# now we have to the extracted files one dir above
			# and then remove the temporary directory
			os.system("cp -r * ../")
			os.chdir(workdir)
			os.system("rm -rf %s" % (tempdir))
		except:
			print "failed!"
			if retrieve_end == 0:
				retrieve_end = time.time()
				retrieve_time = retrieve_end - retrieve_start
				print "\tRetrieve time till the crash: %.1f s (%.1f min)" % (retrieve_time, retrieve_time/60.)
			if extract_end == 0:
				extract_end = time.time()
				extract_time = extract_end - extract_start
				print "\tExtract time till the crash: %.1f s (%.1f min)" % (extract_time, extract_time/60.)
			if rf!="":
				retryout=open(rf, "at")
				retryout.write("%s   %s   %s   %s   %s\n" % (data[0][ii], data[1][ii], data[2][ii], data[3][ii], data[4][ii]))
				retryout.close()

		os.chdir(curdir)

	job_end = time.time()	
	job_time = job_end - job_start
	print "Total retrieve/extract time: %.1f s (%.1f min)" % (job_time, job_time/60.)

# main
if __name__=="__main__":

	cmdline=opt.OptionParser("Usage: %prog <ObsID.txt1> <ObsID.txt2>...")
	cmdline.add_option('--sap', dest='sap', metavar='SAP#', help="retrieve data only for the given SAP", default=-1, type='int')
	cmdline.add_option('--tab', dest='tab', metavar='TAB#', help="retrieve data only for the given TAB", default=-1, type='int')
	cmdline.add_option('--part', dest='part', metavar='PART#', help="retrieve data only for the given PART", default=-1, type='int')
	cmdline.add_option('--summary-only', action="store_true", dest='is_summary_only',
                           help="retrieve only summary directories (CSplots, or CVplots, or redIS). This option has priority over \
options --sap, --tab, --part, and --skip-summary", default=False)
	cmdline.add_option('--skip-summary', action="store_true", dest='is_skip_summary',
                           help="Do not retrieve summaries", default=False)
	cmdline.add_option('--stage-only', action="store_true", dest='is_stage_only',
                           help="Only stage the data without downloading", default=False)
	cmdline.add_option('--skip-staging', action="store_true", dest='is_skip_staging',
                           help="Skip staging and start downloading right away", default=False)
	cmdline.add_option('--obsids', action="store_true", dest='is_obsids',
                           help="input arguments are ObsIDs instead of ascii files. Based on given ObsIDs \
corresponding files will be looked at designated area on CEP2", default=False)
	cmdline.add_option('-f', '--format', dest='format', metavar='FORMAT',
                           help="column format of input ascii files. By default (%default), it is the same as from \
web-summary pages. Other format is 'manual', it's csv format from manual LTA query (expert mode)", default="websummary", type='str')
	cmdline.add_option('--csvfile', dest='csvfile', metavar='CSV-FILE',
                           help="specify single csv-file (comma-separated-values) with srm-links for all given ObsIDs. \
With this option, it is assumed that you give the list of ObsIDs instead of ascii files, therefore this option automatically \
sets --obsids and --format='manual'. Only lines for given ObsIDs will be used from this csv-file", default="", type='str')
	cmdline.add_option('--query', action="store_true", dest='is_query', help="as --csvfile but runs SQL query instead of \
using given csv file. One must specify project as well with --project option. If both --csvfile and --query are given, then \
--csvfile option has the priority", default=False)
        cmdline.add_option('-p', '--project', dest='project', metavar='PROJECT',
                           help="specify the project to query. Only to be used with --query option", default="", type='str')
	cmdline.add_option('-u', '--username', dest='user', metavar='USERNAME',
                           help="specify the LTA username. By default, it's the same as your current login name", default="", type='str')
	cmdline.add_option('-l', '--log', action="store_true", dest='is_log',
                           help="optional parameter to turn on wget output", default=False)
        cmdline.add_option('-q', '--quiet', action="store_true", dest='is_quiet',
                           help="turn off logging from the communication with the LTA database", default=False)
        (opts, args) = cmdline.parse_args()

	# canceling --sap, --tab, and --part options when --summary-only is given
	if opts.is_summary_only:
		opts.sap = -1
		opts.tab = -1
		opts.part = -1
		opts.is_skip_summary = False

	if opts.csvfile != "":
		opts.is_obsids = True
		opts.format = "manual"
		opts.is_query = False

	if opts.is_query:
		opts.is_obsids = True
		opts.format = "manual"
		import awlofar
		from common.database.Database import database
		from common.database.Context import context
		from common.config.Environment import Env

		# re-directing log-messages from log-files in the current dir to /dev/null
		Env['log'].filename = "/dev/null"
		# to turn off logging on the screen
		if opts.is_quiet:
			Env['logprint'] = 0

	        # setting specific project
	        if opts.project != "":
        	        context.set_project(opts.project)
	        else:
        	        print "You should specify project to query when using --query option!"
                	sys.exit(1)

		# running the super-query
		query_result=database.execute_select(superquery)
		if len(query_result) == 0: # this potentially means that you should run query for public data
			context.set_project("ALL")
			query_result=database.execute_select(public_superquery % (opts.project, opts.project, opts.project, opts.project))

	if not opts.is_obsids:
		inputfiles = args
	else:
		if opts.csvfile == "":
			inputfiles = ["%s/%s.txt" % (griddir, ff) for ff in args]

	info=np.asarray([[]]*5)
	if opts.format == "websummary":
		for infile in inputfiles:
			info=np.row_stack((info.T,np.loadtxt(infile, comments='#', usecols=(0,1,2,3,4), dtype=str, unpack=True).T)).T
	elif opts.format == "manual":
		if opts.csvfile != "": # if csv-file is given
			tmp=np.asarray([[]]*5)
			tmp=np.loadtxt(opts.csvfile, comments='#', delimiter=",", usecols=(0,2,2,1,3), dtype=str, unpack=True)
			indices=[ii for ii in xrange(len(tmp[0])) if re.search("|".join(args), tmp[0][ii])]
			info=np.vstack((np.array([ii.strip('"') for ii in tmp[0][indices]]), \
				np.array([ii.split(" ")[0] for ii in tmp[1][indices]]), \
				np.array([ii.split(" ")[-1] for ii in tmp[2][indices]]), \
				np.array([float(ii)/1e9 for ii in tmp[3][indices]]), \
				np.array([ii.strip('"') for ii in tmp[4][indices]]) ))
		elif opts.is_query:
			cond=re.compile(r"%s" % ("|".join(args)))
			indices=[ii for ii in xrange(len(query_result)) if cond.search(query_result[ii][0]) or cond.search("L%s" % (query_result[ii][4]))]
			info=np.vstack((np.array([ii for ii in map(lambda j: query_result[j][0], indices)]), \
				np.array([str(ii).split(" ")[0] for ii in map(lambda j: query_result[j][2], indices)]), \
				np.array([str(ii).split(" ")[-1] for ii in map(lambda j: query_result[j][2], indices)]), \
				np.array([float(ii)/1e9 for ii in map(lambda j: query_result[j][1], indices)]), \
				np.array([ii for ii in map(lambda j: query_result[j][3], indices)]) ))
		else: # if neither --csvfile not --query are given
			tmp=np.asarray([[]]*5)
			for infile in inputfiles:
				tmp=np.loadtxt(infile, comments='#', delimiter=",", usecols=(0,2,2,1,3), dtype=str, unpack=True)
				tmp[0] = np.array([ii.strip('"') for ii in tmp[0]])
				tmp[1] = np.array([ii.split(" ")[0] for ii in tmp[1]])
				tmp[2] = np.array([ii.split(" ")[-1] for ii in tmp[2]])
				tmp[3] = np.array([float(ii)/1e9 for ii in tmp[3]])
				tmp[4] = np.array([ii.strip('"') for ii in tmp[4]])
				info=np.row_stack((info.T, tmp.T)).T
	else:
		print "Wrong format of the input ascii files. Available formats are: websummary (default), and manual"
		sys.exit(1)

	if np.size(info[0]) == 0:
		print "Nothing to retrieve. Exiting..."
		sys.exit(0)

	# file with all srm links that failed to be retrieved/staged
	if opts.csvfile != "":
		retryfile="%s/%s.retry" % (os.getcwd(), opts.csvfile.split("/")[-1])
	elif opts.is_query:
		retryfile="%s/%s.retry" % (os.getcwd(), opts.project.lower())
	else:
		retryfile="%s/%s.retry" % (os.getcwd(), inputfiles[-1].split("/")[-1])

	# excluding files that we don't need
	indices=[]
	print "Skipping..."
	for ii,ff in enumerate(info[0]): 
		srmlink=info[4][ii]
		if re.search("^srm://", srmlink) is None:
			print " %s\t - use 'grid-pulsar-retrieve-file'" % (ff)
			continue # this is not LTA link, but Joeri's grid link
		if opts.is_summary_only:
			if re.search("CSplots", ff) is None and re.search("CVplots", ff) is None and \
					(re.search("locus094", ff) and re.search("redIS", ff)) is None and \
					re.search("summaryCV", ff) is None and re.search("summaryCS", ff) is None and \
					re.search("summaryIS", ff) is None:
				print " %s\t - not summaries" % (ff)
				continue	
		if opts.is_skip_summary:
			if re.search("CSplots", ff) is not None or re.search("CVplots", ff) is not None or \
					(re.search("locus094", ff) and re.search("redIS", ff)) is not None or \
					re.search("summaryCV", ff) is not None or re.search("summaryCS", ff) is not None or \
					re.search("summaryIS", ff) is not None:
				print " %s\t - summaries are skipped" % (ff)
				continue	
		if opts.sap != -1:
			if re.search("_SAP%03d_" % (opts.sap), ff) is None:
				print " %s\t - data are not for a given SAP=%d" % (ff, opts.sap)
				continue	
		if opts.tab != -1:
			if re.search("_B%03d_" % (opts.tab), ff) is None:
				print " %s\t - data are not for a given TAB=%d" % (ff, opts.tab)
				continue	
		if opts.part != -1:
			if re.search("_P%03d_" % (opts.part), ff) is None:
				print " %s\t - data are not for a given PART=%d" % (ff, opts.part)
				continue	
		# list of good indices
		indices.append(ii)

	# extra check in case eventually there is nothing to retrieve
	if np.size(info[0]) == 0 or len(indices) == 0:
		print "Nothing to retrieve. Exiting..."
		sys.exit(0)

	# taking only selected files
	info=info[:,indices]
	# then, we remove duplicating files (e.g. if same ObsID was given more than once)
	tmp, uniq_indices = np.unique(info[0], return_index=True)
	info=info[:,uniq_indices]
	totsize=sum([float(ii) for ii in info[3]])
	print
	print "Files to retrieve (Total = %s):" % (totsize < 1000 and "%.1f GB" % (totsize) or "%.1f TB" % (totsize/1000.))
	for ii, ff in enumerate(info[0]):
		print " %s\t\t(%s GB)" % (ff, info[3][ii])

	# staging
	if opts.user == "":
		opts.user=os.environ['USER']

	total_start = time.time()
	print

	if not opts.is_skip_staging:  # if we are not skipping staging...
		print "Staging..."
		staging_start = time.time()
		proxy=xmlrpclib.ServerProxy("https://webportal.astron.nl/service/xmlrpc")
		staging_id=proxy.LtaStager.add_getid(opts.user, ["%s" % ss for ss in info[4]])
		staging_status=""

		download_urls=set()
		while True:
			time.sleep(5)
			status=proxy.LtaStager.getstatus(opts.user, staging_id).strip().lower()
			if status == "failed":
				staging_end = time.time()
				staging_time = staging_end - staging_start
				print "failed!"
				print "Staging time till the crash: %.1f s (%.1f min)" % (staging_time, staging_time/60.)
				sys.exit(1)
			elif status in ["new", "scheduled"]:
				if status != staging_status:
					staging_status=status
					print "Staging status: %s" % (staging_status)
				continue
			elif status == "in progress":
				if status != staging_status:
					staging_status=status
					print "Staging status: %s" % (staging_status)
				if not opts.is_stage_only:
					ready_urls=set(proxy.LtaStager.getstagedurls(opts.user, staging_id))
					if set(ready_urls-download_urls): # if we have new files to retrieve
						staging_end = time.time()
						staging_time = staging_end - staging_start
						print "Staging time since the staging start or after the previous retrieval: %.1f s (%.1f min)" % (staging_time, staging_time/60.)
						staging_start = time.time()
						indices=[ii for ii in np.arange(np.size(info[4])) if info[4,ii] in set(ready_urls - download_urls)]
						# downloading
						retrieve(info[:, indices], retryfile)
						download_urls = ready_urls.copy()
						staging_status=""
			elif status == "success":
				print "Staging status: %s" % (status)
				staging_end = time.time()
				staging_time = staging_end - staging_start
				print "Staging time (includes possible overhead due to retrieval): %.1f s (%.1f min)" % (staging_time, staging_time/60.)
				if not opts.is_stage_only:
					print "Retrieving the rest of the files:"
					ready_urls=set(proxy.LtaStager.getstagedurls(opts.user, staging_id))
					if set(ready_urls-download_urls): # if we have new files to retrieve
						indices=[ii for ii in np.arange(np.size(info[4])) if info[4,ii] in set(ready_urls - download_urls)]
						# downloading
						retrieve(info[:, indices], retryfile)
					else:
						print "All files have been retrieved!"
				break
			elif status == "partial success":
				print "Staging status: %s" % (status)
				staging_end = time.time()
				staging_time = staging_end - staging_start
				print "Staging time (includes possible overhead due to retrieval): %.1f s (%.1f min)" % (staging_time, staging_time/60.)
				if not opts.is_stage_only:
					print "Retrieving the rest of the staged files:"
					ready_urls=set(proxy.LtaStager.getstagedurls(opts.user, staging_id))
					if set(ready_urls-download_urls): # if we have new files to retrieve
						indices=[ii for ii in np.arange(np.size(info[4])) if info[4,ii] in set(ready_urls - download_urls)]
						# downloading
						retrieve(info[:, indices], retryfile)
					else:
						print "All staged files have been retrieved!"
				failed_indices=[ii for ii in np.arange(np.size(info[4])) if info[4,ii] not in ready_urls]
				retryout=open(retryfile, "at")
				print "\nFiles failed to be staged:"
				for ff in info[0,failed_indices]:
					print " %s" % ff
					retryout.write("%s   %s   %s   %s   %s\n" % (info[0][ff], info[1][ff], info[2][ff], info[3][ff], info[4][ff]))
				retryout.close()
				break
			else:
				print "unknown staging status: %s. Exiting..." % (status)
				break

	else: # start retrieving immediately
		print "Staging is skipped"
		print
		retrieve(info[:,:], retryfile)

	total_end = time.time()
	total_time = total_end - total_start
	print "Total retrieve time: %.1f s (%.1f min)" % (total_time, total_time/60.)
