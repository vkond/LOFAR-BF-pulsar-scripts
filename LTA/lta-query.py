#!/usr/bin/env python
#
# Script that connects to LTA db, runs SQL query
# for a given project and saves query results to csv file
#
# Vlad Kondratiev, Jul 21, 2015 (c)
#
import os, os.path, sys, re
import optparse as opt
import numpy as np
# LTA related
import awlofar
from common.database.Database import database
from common.database.Context import context

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


# main
if __name__=="__main__":

	cmdline=opt.OptionParser("Usage: %prog -p <project>")
	cmdline.add_option('-p', '--project', dest='project', metavar='PROJECT',
                           help="specify the project to query (mandatory)", default="", type='str')
        (opts, args) = cmdline.parse_args()

	# check if any input parameters are given
	if len(sys.argv[1:]) == 0:
		cmdline.print_usage()
		sys.exit(0)

	# setting specific project
	if opts.project != "":
		context.set_project(opts.project)
	else:
		print "You should specify project to query!"
		sys.exit(1)

	# running the super-query
	result=database.execute_select(superquery)

	# writing results to the output csv file
	out = open("%s.csv" % (opts.project.lower()), "w")
	out.write("FILENAME,FILESIZE,CREATION_DATE,URI,OBSERVATIONID\n")
	for row in result:
		line = "\"%s\",%s,%s,\"%s\",\"%s\"\n" % (row[0], row[1], row[2], row[3], row[4])
		out.write(line)
	out.close()
