#!/usr/bin/env python
#
# Script that connects to LTA db, runs SQL query
# for a given project and saves query results to csv file
#
# Vlad Kondratiev, Jul 21, 2015 (c)
#
# Nov 25, 2016 - added query to retrieve the data from projects
#                that became public and you were not the member
# Nov 25, 2016 - added logging handling when communicating with the
#                database; there is no logfile now when there are
#                errors; output on the screen can be turned off
#                with the cmdline option -q
#
import os, os.path, sys, re
import optparse as opt
import numpy as np
# LTA related
import awlofar
from common.database.Database import database
from common.database.Context import context
from common.config.Environment import Env

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




# main
if __name__=="__main__":

	cmdline=opt.OptionParser("Usage: %prog -p <project>")
	cmdline.add_option('-p', '--project', dest='project', metavar='PROJECT',
                           help="specify the project to query (mandatory)", default="", type='str')
        cmdline.add_option('-q', '--quiet', action="store_true", dest='is_quiet',
                           help="turn off logging from the communication with the LTA database", default=False)
        (opts, args) = cmdline.parse_args()

	# check if any input parameters are given
	if len(sys.argv[1:]) == 0:
		cmdline.print_usage()
		sys.exit(0)

	# re-directing log-messages from log-files in the current dir to /dev/null
	Env['log'].filename = "/dev/null"
	# to turn off logging on the screen
	if opts.is_quiet:
		Env['logprint'] = 0

	# setting specific project
	if opts.project != "":
		context.set_project(opts.project)
	else:
		print "You should specify project to query!"
		sys.exit(1)

	# running the super-query
	result=database.execute_select(superquery)
	if len(result) == 0: # this potentially means that you should run query for public data
		print "Now will try to query public data..."
		context.set_project("ALL")
		result=database.execute_select(public_superquery % (opts.project, opts.project, opts.project, opts.project))

	if len(result) != 0:
		print "Success!"
		# writing results to the output csv file
		out = open("%s.csv" % (opts.project.lower()), "w")
		out.write("FILENAME,FILESIZE,CREATION_DATE,URI,OBSERVATIONID\n")
		for row in result:
			line = "\"%s\",%s,%s,\"%s\",\"%s\"\n" % (row[0], row[1], row[2], row[3], row[4])
			out.write(line)
		out.close()
	else:
		print "No data available!"
