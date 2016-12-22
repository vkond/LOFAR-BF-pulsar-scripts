#!/usr/bin/env python
#
#
# Version 1.0: uses old TAB addressing
# Version 1.1: updated TABs but doesn't support multiple TABs
# Version 1.2: TABs use offsets from the templates
# Version 1.3 (current): changed default template directory to LOFARSOFT/...
import optparse as opt
import numpy as np
import random, string
import warnings
import time, datetime
from datetime import *
import ephem
import os, sys, math, re
from xml.dom import minidom

# This suprresses some DeprecationWarnings from pyephem
warnings.filterwarnings("ignore", category=DeprecationWarning)

##############################################################################################

# define class for xml generation of the schedule
class xmlSched:
	# dur in minutes
	# ra as HH:MM:SS.SS
	# dec as [+/-]DD:MM:SS.SS
	# tempdir - dir for output xml files
	def __init__(self, template, project, index, psr, dur, start, end, ra, dec, bad_stations_str, tempdir=None):

		# list of all core stations (only numbers), real names should be preceded with CS%03d
		self.corestations=[1, 2, 3, 4, 5, 6, 7, 11, 13, 17, 21, 24, 26, 28, 30, 31, 32, 101, 103, 201, 301, 302, 401, 501]
		# template dir
		try:
			# Lofar software
			lofarsoft=os.environ['LOFARSOFT']
			self.template_dir = "%s/release/share/pulsar/data/templates" % (lofarsoft,)
		except:
			self.template_dir = "/home/kondratiev/scheduling/templates"
		self.template = template
		self.project = project
		self.obsindex = index
		self.psr = psr
		self.length = str(dur)
		self.duration = int(self.length) * 60 # duration in secs
		self.start = start
		self.end = end
		self.ra = ra
		self.dec = dec
		self.update_stations = True

		self.radeg = self.ra2deg(self.ra)
		self.decdeg = self.dec2deg(self.dec)

		if bad_stations_str == "-1": # don't update stations in xml at all
			self.update_stations = False
		elif bad_stations_str == "0" or bad_stations_str == "": # use all core stations (don't exclude anything)
			self.good_stations = ["CS%03d" % cs for cs in self.corestations]
		else: # considering this column as a list of stations to exclude
			self.good_stations=[]
			bads = [int(s) for s in bad_stations_str.split(",")]
			for cs in self.corestations:
				if cs in bads: continue
				else:
					self.good_stations.append("CS%03d" % cs)

		# initializing the xml-doc
		if not os.path.exists(self.template):
			# now checking the templates/ sub-directory in the current directory
			self.template = "templates/%s" % (template)
			if not os.path.exists(self.template):
				self.template = "%s/%s" % (self.template_dir, template)
				if not os.path.exists(self.template):
					print "Can't find template '%s' neither in the current directory \
nor in the templates/ sub-directory nor in the '%s' directory!" % (template, self.template_dir)
					sys.exit(1)
		try:
			self.xmldoc=minidom.parse(self.template)
		except:
			print "Can't open/read XML template for the pulsars %s!" % (self.psr)
			sys.exit(1)

	# update the xml
	def update(self):
		self.updateProjectName(self.xmldoc.firstChild, self.project, "PROJECT_NAME")
		self.updatePointingIndex(self.xmldoc.firstChild, self.obsindex, "OBSINDEX")
		self.updatePipelineIndex(self.xmldoc.firstChild, self.obsindex, "PIPEINDEX")
		self.updatePointingName(self.xmldoc.firstChild, self.psr, "PSRNAME")
		self.updateDescription(self.xmldoc.firstChild, self.length, "LENGTHMIN")
		self.updateStartEndTime(self.xmldoc.firstChild, self.start, self.end)
		self.updateSAPcoordinates(self.xmldoc.firstChild, self.radeg, self.decdeg)
		self.updateTABcoordinates(self.xmldoc.firstChild, (self.radeg/180.)*math.pi, (self.decdeg/180.)*math.pi)
		self.updateObsDuration(self.xmldoc.firstChild, self.duration, "LENGTHSEC")
		self.updatePipelineDuration(self.xmldoc.firstChild, 10 * self.duration, "PIPELENGTH") # assume that pipeline runtime is 10x longer than obslength
		if self.update_stations:
			self.updateStationsList(self.xmldoc.firstChild, self.good_stations)

	# write xml to a file
	def write(self, filename):
		out=open(filename, 'wb')
		self.xmldoc.writexml(out, encoding='UTF-8')
		out.close()

	# convert RA string to DEG
	def ra2deg(self, rastr):
		hh, mm, ss = rastr.split(":")
		return 15.*(float(hh)+(float(mm)+float(ss)/60.)/60.)

	# convert DEC string to DEG
	def dec2deg(self, decstr):
		sign = 1.
		dd, mm, ss = decstr.split(":")
		if dd[0] == '-': sign = -1.
		return sign*(float(abs(int(dd)))+(float(mm)+float(ss)/60.)/60.)

	# update the pointing index from "OBSINDEX" in the template
	def updatePointingIndex(self, root, obsindex, stem_index):
        	if root.nodeName == "item" and "index" in root.attributes.keys() and root.attributes['index'].value == stem_index:
                	root.attributes['index'].value = "%d" % (obsindex)
        	if root.nodeName == "name"  or root.nodeName == "topology":
			root.childNodes[0].data=re.sub(stem_index, "%d" % (obsindex), root.childNodes[0].data)
		if root.nodeName == "predecessor_topology" and len(root.childNodes) != 0:
			root.childNodes[0].data=re.sub(stem_index, "%d" % (obsindex), root.childNodes[0].data)
		if root.nodeName == "lofar:bfDataProduct" and "topology" in root.attributes.keys():
			root.attributes['topology'].value=re.sub(stem_index, "%d" % (obsindex), root.attributes['topology'].value)
	        for child in root.childNodes:
        	        self.updatePointingIndex(child, obsindex, stem_index)

	# update the pipeline index from "PIPEINDEX" in the template
	def updatePipelineIndex(self, root, obsindex, stem_index):
        	if root.nodeName == "item" and "index" in root.attributes.keys() and root.attributes['index'].value == stem_index:
                	root.attributes['index'].value = "%d" % (10000+obsindex)
	        for child in root.childNodes:
        	        self.updatePipelineIndex(child, obsindex, stem_index)

	# update the name of the project
	def updateProjectName(self, root, projectname, projectname_stem):
        	if root.nodeName == "name" or root.nodeName == "projectName":
                	if root.childNodes[0].data == projectname_stem:
                        	root.childNodes[0].data=projectname
	        for child in root.childNodes:
        	        self.updateProjectName(child, projectname, projectname_stem)

	# update the name of the pointing from "PSRNAME" to proper pulsar name
	def updatePointingName(self, root, psrname, psrname_stem):
        	if root.nodeName == "name"  or root.nodeName == "targetName" or root.nodeName == "description" or root.nodeName == "topology":
			root.childNodes[0].data=re.sub(psrname_stem, psrname, root.childNodes[0].data)
		if root.nodeName == "predecessor_topology" and len(root.childNodes) != 0:
			root.childNodes[0].data=re.sub(psrname_stem, psrname, root.childNodes[0].data)
		if root.nodeName == "lofar:bfDataProduct" and "topology" in root.attributes.keys():
			root.attributes['topology'].value=re.sub(psrname_stem,psrname,root.attributes['topology'].value)
	        for child in root.childNodes:
        	        self.updatePointingName(child, psrname, psrname_stem)

	# update the pointing name in the description from descr_stem to descr
	def updateDescription(self, root, descr, descr_stem):
        	if root.nodeName == "description":
			root.childNodes[0].data=re.sub(descr_stem, descr, root.childNodes[0].data)
	        for child in root.childNodes:
        	        self.updateDescription(child, descr, descr_stem)

	# set start and stop times of observation
	def updateStartEndTime(self, root, st, et):
        	if root.nodeName == "startTime": root.childNodes[0].data=st
	        if root.nodeName == "endTime": root.childNodes[0].data=et
        	for child in root.childNodes:
                	self.updateStartEndTime(child, st, et)

	# update RA and DEC for SAP
	def updateSAPcoordinates(self, root, ra, dec):
        	if root.nodeName == "ra": root.childNodes[0].data=ra
	        if root.nodeName == "dec": root.childNodes[0].data=dec
        	for child in root.childNodes:
                	self.updateSAPcoordinates(child, ra, dec)

	# update RA and DEC for TABs (values in template are *offsets* from the center)
	def updateTABcoordinates(self, root, ra, dec):
        	if root.nodeName == "angle1": 
			d = root.childNodes[0].data
			root.childNodes[0].data = str(float(d)+float(ra))
	        if root.nodeName == "angle2": 
			d = root.childNodes[0].data
			root.childNodes[0].data = str(float(d)+float(dec))
        	for child in root.childNodes:
                	self.updateTABcoordinates(child, ra, dec)

	# update duration
	def updateObsDuration(self, root, dur, durstem):
        	if root.nodeName == "duration":
	               	if root.childNodes[0].data == durstem:
				root.childNodes[0].data=dur
	        for child in root.childNodes:
        	        self.updateObsDuration(child, dur, durstem)

	# update pipeline duration
	def updatePipelineDuration(self, root, dur, pipelength_stem):
        	if root.nodeName == "duration":
	               	if root.childNodes[0].data == pipelength_stem:
				root.childNodes[0].data=dur
	        for child in root.childNodes:
        	        self.updatePipelineDuration(child, dur, pipelength_stem)

	# update stations list
	def updateStationsList(self, root, good_stations):
        	if root.nodeName == "stations":
			for st in good_stations:
				newnode=self.xmldoc.createElement("station")
				newnode.setAttribute('name', st)
				root.appendChild(newnode)
	        for child in root.childNodes:
        	        self.updateStationsList(child, good_stations)

##############################################################################################

# define class to schedule a pulsars for a given slot
class psrSched:
	# targets is list of list of fields from input source list
	def __init__(self, targets, start, end, latitude, longitude, horizon, transit_tolerance, gap, project, is_plot, blocked=""):
		self.targets = targets
		self.start = self.round_time(start)
		self.end = self.round_time(end)
		self.mstart = 0 # start of the slot in mins relative to the starting point
		self.mend = int(1440.0*(self.round_time(self.end - self.start))) # end of slot in mins after the start
		self.transit_tolerance = transit_tolerance
		self.horizon = horizon
		self.gap = gap
		self.project = project
		self.is_plot = is_plot
		self.blocked = blocked

		# The coordinates and horizon of LOFAR
		self.observer = ephem.Observer()
		self.observer.lat = latitude
		self.observer.long = longitude
		self.observer.horizon = (self.horizon * math.pi)/180.0

		self.toobserve = [] # For holding a subset of pulsars to observe
		self.transits  = [] # For holding the transit time of each pulsar
		self.scores    = [] # For holding the scores for each pulsar

		for ii in np.arange(len(self.targets)):
			# Set up this source in pyephem
			star      = ephem.FixedBody()
		        star._ra  = self.targets[ii][1]
	        	star._dec = self.targets[ii][2]
	        	star.compute(self.observer)

		        visible = False
        		# Determine if the source is visible at some point during the observation
	        	if star.circumpolar:
        			visible = True
		        elif star.neverup:
				visible = False
		        else:
        			for mm in np.arange(self.mstart, self.mend+1, 1):
					self.observer.date = ephem.Date(self.start + mm*ephem.minute)
	                		star.compute(self.observer)
        	        		if star.alt > self.observer.horizon: 
						visible = True
						break

	        	if visible:
				self.observer.date = self.start
	        	        star.compute(self.observer)
				next_transit = self.observer.next_transit(star)
				previous_transit = self.observer.previous_transit(star)
				score = 0.0

				# if transit time was within 1h from the start or will be 
				# within 1h after the end of observing slot, then these sources
				# get a much higher priority to schedule as close to transit
				# as possible
				approach = 1440*(self.start-previous_transit)
				if approach > 0:
					if 1440*(next_transit-self.end) > 0:
						approach = np.min([approach, 1440*(next_transit-self.end)])
				else:
					approach = 1440*(next_transit-self.end)		
				if approach > 0 and approach <= self.transit_tolerance:
					score += 100.0
					# we also change the score depending on how far from transit the time is (within this 1h)
					# further from transit, score is higher
					score += approach/10. # divide by 10 here in order to have low EL similar (or higher) impact

				# depending on source's maximum EL we change the score, the smaller DEC get larger score
				elev=(math.pi/2.-np.abs(self.observer.lat-star.dec))*(180./math.pi)
				score += 90/elev

				# incriment the score by 1 if the source transits during the observation
			        if next_transit >= self.start and next_transit <= self.end:
					score += 10.0
                			# Append the transit time
	                		self.transits.append(next_transit)
				# If the source doesn't transit during the obervation, store
				# the time of the most closest transit (may be in past or future)
				elif (self.start - previous_transit) < (next_transit - self.end):
					self.transits.append(previous_transit)
				else:
					self.transits.append(next_transit)

				self.toobserve.append(targets[ii])
				# check if there is a Priority field in the source list and it's not 0 (use default)
				if len(self.targets[ii]) > 4 and int(self.targets[ii][4]) != 0:
					self.scores.append(float(self.targets[ii][4]))
				else:
					self.scores.append(score)

		# Now sort everything by the scores
		args = list(reversed(np.argsort(self.scores)))
		self.toobserve = np.asarray(self.toobserve)[args]
		self.transits = np.asarray(self.transits)[args]
		self.scores = np.asarray(self.scores)[args]

		# The following portion of the script schedules sources in order of
		# their score as close to their transits as possible.  It makes use of
		# python sets to find overlapping observing slots.  The sets store
		# obsering times incremented every minute in string format to avoid
		# problems with floating point comparison.  Some sets have equivalent
		# times as floats for compuation.

		# we need to put 1-min reservations at the start and end so we can calculate time differences till ends of the slot as well
		reserved = set([t for t in (self.mstart-1, self.mend+1)])

		# if we have some time intervals pre-reserved within the main slot, then we parse
		# the reserved-blocks string here and add these reserved blocks to our set
		if self.blocked != "":
			blocked_slots = self.blocked.split(",")
			for block in blocked_slots:
				if not re.search("/", block):
					print "Wrong format of reserved-blocks option! Missing /!"
					sys.exit(1)
				blockstart = ephem.Date(block.split("/")[0].replace("-", "/").replace("T", " "))
				blockend = ephem.Date(block.split("/")[1].replace("-", "/").replace("T", " "))
				blockstart_min = int(1440.0*(self.round_time(blockstart - self.start)))
				blockend_min = int(1440.0*(self.round_time(blockend - self.start)))
				reserved = reserved.union(set(list(np.arange(blockstart_min, blockend_min + 1, 1))))
				reserved_sorted = np.sort(list(reserved))

		self.schedpsrs = [] # list of scheduled pulsars
		self.start_times = []
		self.end_times = []
		self.alts = [] # list of elevations
		self.advanced = [] # list of diffs between observation midpoint and transit

		if self.is_plot:
			import matplotlib
			matplotlib.use("Agg")
			import matplotlib.pyplot as plt
			import matplotlib.gridspec as gridspec
			import matplotlib.dates as mdates
			fig = plt.figure(figsize=(10, 7.5), dpi=300)
			gs  = gridspec.GridSpec(2, 1, height_ratios=[5,1])
			ax1 = plt.subplot(gs[0])
			ax2 = plt.subplot(gs[1], sharex=ax1)

		for ii in np.arange(len(self.toobserve)):
			source=self.toobserve[ii][0]
			transit=self.transits[ii]
			obstime=int(self.toobserve[ii][3])
			# Start off assuming this source can be scheduled
			schedule = True

			star = ephem.FixedBody()
			star._ra = self.toobserve[ii][1]
			star._dec = self.toobserve[ii][2]

			# Try setting the mid-point of a scan at the transit, rounded to the closest minute
			scanstart = int(1440.0*(self.round_time(transit - (0.5*obstime)*ephem.minute) - self.start))
			# The extra minute accounts for setup time in between scans
			scanend = scanstart + obstime

			# Adjust the observing times as needed if we scheduled a scan outside our observing window
			if scanstart < self.mstart:
				shift = self.mstart - scanstart
				if self.mstart - (scanstart + 0.5*obstime) > self.transit_tolerance:
					schedule = False 
				scanstart += shift
				scanend += shift
			if scanend > self.mend:
				shift = scanend - self.mend
				if scanend - 0.5*obstime - self.mend > self.transit_tolerance:
					schedule = False 
				scanstart -= shift
				scanend -= shift
	
			# This set holds the times that our scan takes up
			scan = set(list(np.arange(scanstart, scanend + 1, 1)))

			# This if statement will trigger if a scan overlaps with a previously scheduled one
			if not reserved.isdisjoint(scan) and schedule:
				# This will check to see if there any open slots greater than the scan time
				if any(np.diff(reserved_sorted) >= obstime + 2*self.gap):
					# If we have room to schedule this scan, find the open time slot closest to transit
					smallestdiff = 24.*ephem.hour
					for arg,slot in enumerate(np.diff(reserved_sorted)):
						# This finds the slots that are big enough for our observation
						if slot >= obstime + 2*self.gap:
							slotstart = np.array(reserved_sorted)[arg] + self.gap
							slotend = np.array(reserved_sorted)[arg+1] - self.gap
							if ephem.Date(self.start+ephem.minute*slotend) <= self.round_time(transit):
								trialstart = slotend - obstime
								trialend = slotend
							elif ephem.Date(self.start + ephem.minute*slotstart) >= self.round_time(transit):
								trialstart = slotstart
								trialend = slotstart + obstime
							else:
								trialstart = int(1440.0*(self.round_time(transit - 0.5*obstime*ephem.minute) - self.start))
								if trialstart < slotstart:
									trialstart = slotstart
								trialend = trialstart + obstime
								if trialend > slotend:
									trialend = slotend
									trialstart = trialend - obstime

							difference = abs(ephem.Date(self.start + (trialstart + 0.5*obstime)*ephem.minute) - transit)
							# Check to see if this is the best slot so far
							if difference < smallestdiff:
								smallestdiff = difference
								scanstart = trialstart
								scanend = trialend

					# check if scheduled slot is far from the transit
					if 1440*abs(ephem.Date(self.start + (scanstart + 0.5*obstime)*ephem.minute) - transit) > self.transit_tolerance:
						schedule = False

					# here we check that for chosen times, the elevation at any given minute is above the limit
					elevations = []
					for mm in np.arange(scanstart, scanend + 1, 1):
						self.observer.date = ephem.Date(self.start + mm*ephem.minute)
						star.compute(self.observer)
						elevations.append(star.alt)
					if np.min(elevations) < self.observer.horizon:
						schedule = False
				else:
					schedule = False


			if schedule:
				star = ephem.FixedBody()
				star._ra = self.toobserve[ii][1]
				star._dec = self.toobserve[ii][2]
				reserved = reserved.union(set(list(np.arange(scanstart, scanend + 1, 1))))
				reserved_sorted = np.sort(list(reserved))
				self.observer.date = ephem.Date(self.start + (scanstart + 0.5*obstime)*ephem.minute)
				star.compute(self.observer)
				alt = (star.alt/math.pi)*180.
				transit_advanced=1440.0*(ephem.Date(self.start + (scanstart + 0.5*obstime)*ephem.minute) - transit)

				# adding source to lists
				self.schedpsrs.append(self.toobserve[ii])
				self.start_times.append(self.round_time(self.start + scanstart*ephem.minute))
				self.end_times.append(self.round_time(self.start + scanend*ephem.minute))
				self.alts.append(alt)
				self.advanced.append(transit_advanced)

				# Plotting part
				if self.is_plot:
					ts = np.arange(self.mstart, self.mend + 1, 1)
					elevs = []
					for t in ts:
						self.observer.date = ephem.Date(self.start + t*ephem.minute)
						star.compute(self.observer)
						elevs.append(star.alt*180/math.pi)
					ax1.plot([datetime.strptime(str(ephem.Date(self.start + t*ephem.minute)), "%Y/%m/%d %H:%M:%S") for t in ts], elevs, ":", color='black', alpha=0.7)
					ts = np.arange(scanstart, scanstart + obstime + 1, 1)
					elevs = []
					for t in ts:
						self.observer.date = ephem.Date(self.start + t*ephem.minute)
						star.compute(self.observer)
						elevs.append(star.alt*180/math.pi)
					ax1.plot([datetime.strptime(str(ephem.Date(self.start +    t*ephem.minute)), "%Y/%m/%d %H:%M:%S") for t in ts], elevs, "r-", lw=2)
					ax2.axvspan(datetime.strptime(str(ephem.Date(self.start + scanstart*ephem.minute)), "%Y/%m/%d %H:%M:%S"), \
						datetime.strptime(str(ephem.Date(self.start + (scanstart + obstime + 1)*ephem.minute)), "%Y/%m/%d %H:%M:%S"), facecolor="darkgreen", alpha=0.3)
					ax2.text(datetime.strptime(str(ephem.Date(self.start + (scanstart + 0.5*obstime)*ephem.minute)), "%Y/%m/%d %H:%M:%S"), \
						0.5, source, fontsize=8, ha="center", va="center", rotation=90)
					# drawing a hatched region when time is blocked by other observations
					if self.blocked != "":
						blocked_slots = self.blocked.split(",")
						for block in blocked_slots:
							if not re.search("/", block):
								print "Wrong format of reserved-blocks option! Missing /!"
								sys.exit(1)
							blockstart = ephem.Date(block.split("/")[0].replace("-", "/").replace("T", " "))
							blockend = ephem.Date(block.split("/")[1].replace("-", "/").replace("T", " "))
							ax2.add_patch(plt.Rectangle([datetime.strptime(str(blockstart), "%Y/%m/%d %H:%M:%S"), 0], \
								ephem.Date(blockend-blockstart), 10, fill=True, edgecolor="None", facecolor="None", hatch="x", alpha=0.1))

		# Now sort everything by the start time
		args = np.argsort(self.start_times)
		self.schedpsrs = np.asarray(self.schedpsrs)[args]
		self.start_times = np.asarray(self.start_times)[args]
		self.end_times = np.asarray(self.end_times)[args]
		self.alts = np.asarray(self.alts)[args]
		self.advanced = np.asarray(self.advanced)[args]

		# Plotting part
		if self.is_plot:
			# '7' - width, i.e. 7 days, could be any large number
			rect = plt.Rectangle([datetime.strptime(str(self.start), "%Y/%m/%d %H:%M:%S"), 0], 7, \
				180.*self.observer.horizon/math.pi, fill=True, facecolor="gray", hatch="x", alpha=0.1)
			ax1.add_patch(rect)
			# finding local time for exact LST hours 
			self.observer.date = self.start
			ss = ephem.Date(self.start + ephem.second*(60.-float(str(self.observer.sidereal_time()).split(":")[-1])))
			ax1.text(datetime.strptime(str(ephem.Date(self.end+10*ephem.minute)), "%Y/%m/%d %H:%M:%S"), 2, "LST", color="orange", fontsize=9, ha="left")
			nn = 0
			for tt in np.arange(ss, self.end, ephem.minute):
				corrected = ephem.Date(tt-ephem.second*0.163889*nn) # every second LST gets larger by 0.163889 s
				self.observer.date = ephem.Date(corrected)
				nn += 1
				lst = str(self.observer.sidereal_time())
				lsth = int(lst.split(":")[0])
				lstm = int(lst.split(":")[1])
				if lstm == 0:
					ax1.axvline(datetime.strptime(str(corrected), "%Y/%m/%d %H:%M:%S"), color='orange', ls='--', alpha=0.5)
					ax1.text(datetime.strptime(str(ephem.Date(corrected+5*ephem.minute)), "%Y/%m/%d %H:%M:%S"), 2, \
						"%dh" % (lsth), color="orange", fontsize=9, ha="left")

			star = ephem.Sun()
			self.observer.date = self.start
			self.observer.horizon = 0
			star.compute(self.observer)
			ax1.axvline(datetime.strptime(str(star.rise_time), "%Y/%m/%d %H:%M:%S"), color='black', ls='-.')
			if star.rise_time >= self.start and star.rise_time <= self.end:
				ax1.text(datetime.strptime(str(ephem.Date(star.rise_time+3*ephem.minute)), "%Y/%m/%d %H:%M:%S"), 6, "sunrise", \
					color="red", fontsize=8, va='bottom', rotation=90)
			ax1.axvline(datetime.strptime(str(star.set_time), "%Y/%m/%d %H:%M:%S"), color='black', ls='-.')
			if star.set_time >= self.start and star.set_time <= self.end:
				ax1.text(datetime.strptime(str(ephem.Date(star.set_time+3*ephem.minute)), "%Y/%m/%d %H:%M:%S"), 6, "sunset", \
					color="red", fontsize=8, va='bottom', rotation=90)

			ax1.set_xlim(datetime.strptime(str(self.start), "%Y/%m/%d %H:%M:%S"), datetime.strptime(str(ephem.Date(self.end+ephem.second)), "%Y/%m/%d %H:%M:%S"))
			ax1.set_ylim(0, 90)
			ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
			ax1.set_ylabel("Elevation (deg)")

			ax2.set_xlim(datetime.strptime(str(self.start), "%Y/%m/%d %H:%M:%S"), datetime.strptime(str(ephem.Date(self.end+ephem.second)), "%Y/%m/%d %H:%M:%S"))
			ax2.set_yticklabels([], visible=False)
			ax2.set_yticks([])
			ax2.set_xlabel("Time (UTC)")

			fig.suptitle("\n\n%sObserving Schedule for %s" % (self.project != "" and "%s " % self.project or "", \
				self.start.datetime().strftime("%Y-%m-%d")))
			fig.savefig("%ssched-%s.png" % (self.project != "" and ("%s-" % self.project) or "", \
				self.start.datetime().strftime("%Y-%m-%d")), bbox_inches="tight", pad_inches=0.5)

	# convert DEC string to RAD
	def dec2rad(self, decstr):
        	sign = 1.
	        dd, mm, ss = decstr.split(":")
        	if dd[0] == '-': sign = -1.
	        return (math.pi*sign*(float(abs(int(dd)))+(float(mm)+float(ss)/60.)/60.))/180.

	# return the list of scheduled pulsar, times, etc as a tuple of lists
	def get_schedule(self):
		return (self.schedpsrs, self.start_times, self.end_times, self.alts, self.advanced)

	# print the schedule
	def print_schedule(self):
		print
		print "Pulsar\t\tDEC\t\tSTART\t\t\tLENGTH (min)\tEND\t\t\tTRANSIT (MIN)\tEL/MaxEL (deg)"
		print "-------------------------------------------------------------------------------------------------------------------------------"
		for ii in np.arange(len(self.schedpsrs)):
	        	start_scan_str = self.round_time(self.start_times[ii]).datetime().strftime("%Y-%m-%dT%H:%M:%S")
	        	end_scan_str = self.round_time(self.end_times[ii]).datetime().strftime("%Y-%m-%dT%H:%M:%S")
			maxelev=(math.pi/2.-np.abs(self.observer.lat-self.dec2rad(self.schedpsrs[ii][2])))*(180./math.pi)
		        print "%s\t%s\t%s\t%d\t\t%s\t%.1f\t\t%.1f (%.1f)" % (self.schedpsrs[ii][0], self.schedpsrs[ii][2], \
				start_scan_str, int(self.schedpsrs[ii][3]), end_scan_str, self.advanced[ii], self.alts[ii], maxelev)


	# Round a pyephem date to the nearest minute and return it in string form
	def round_time(self, date, round_to=60):
		date     = ephem.Date(date).datetime()
		seconds  = (date - date.min).seconds
		# The // operator does floor diviison
		rounding = (seconds + 0.5*round_to) // round_to * round_to
		date     = date + timedelta(0, rounding - seconds, -date.microsecond)

                if int(ephem.Date(date).tuple()[-1]) == 59:
		#if str(ephem.Date(date)).split(":")[-1] == "59":
			return ephem.Date(ephem.Date(date) + 0.000001*ephem.second)
		else:
			return ephem.Date(date)



###  M A I N ###
if __name__ == "__main__":

        version = "v1.3"
        usage = "Usage: %prog [-h|--help] [OPTIONS]"
        cmd = opt.OptionParser(usage, version="%prog " + version)
        cmd.add_option('-f', '--infile', dest='targets_file', metavar='FILE',
                           help="input ascii list of sources to observe. Columns are: name, \
RA, DEC, Duration in minutes. Optionally one can add also Priority, Template, Excluded (core) stations, \
and Extra script in the next columns for specific pulsars. If priority is 0 (or not given), then auto-assigned \
priorities will be used. To significantly exceed any auto-priorities use values above 1000. Template column can \
have name of an xml file located in ~kondratiev/scheduling/templates. This template has precedence over the template \
given by -t option. If template is 0 (or not given) template from the command line will be used. The column for excluded \
core stations was meant to provide a way to override the excluded stations from the command line for the case when \
different pulsars need different stations setup (e.g. LBA vs. HBA, for LBA exclyding CS013 is not necessary, etc.) \
Or, if all stations are set up already in the template and are not needed to be updated. If value in this column is -1 \
then stations won't be updated in xml. If value is 0, then use _all_ core stations, any other number or comma-separated \
field of numbers will be considered as particular core stations to exclude. Extra script can be used to update already pre-made \
xml to do other extra stuff, e.g. add manual TABs, etc. This script should take only one argument which is the name of the \
xml-file. Body of the script will likely be some other program to do the necessary tweaks.", default="", type='str')
        cmd.add_option('-s', '--start', dest='startTime', metavar='YYYY-MM-DDThh:mm:ss',
                           help="specify the UT start date/time of a given time slot", default="", type='str')
        cmd.add_option('-e', '--end', dest='endTime', metavar='YYYY-MM-DDThh:mm:ss',
                           help="specify the UT end date/time of a given time slot", default="", type='str')
        cmd.add_option('-p', '--project', dest='project', metavar='PROJECT',
                           help="specify project code. It's obligatory when --sched is used", default="", type='str')
        cmd.add_option('--folder', dest='folder', metavar='MoM-FOLDER',
                           help="specify MoM folder for scheduled observations. By default, \
all observations will be put in the root tree", default="", type='str')
        cmd.add_option('-t', '--template', dest='template', metavar='XML',
                           help="specify the XML template. It's obligatory when --sched is used, \
unless template is given in the input ascii source list", default="", type='str')
        cmd.add_option('--excluded-stations', dest='bad_stations', metavar='LIST',
                           help="specify comma-separated list (no spaces) of excluded Core stations. \
Give only numbers without CS prefix, e.g. 13,21,30. By default excluded are: %default", default="13", type='str')
        cmd.add_option('--reserved-blocks', dest='reserved_blocks', metavar='LIST',
                           help="specify a list of reserved blocks within the main slot (that is specified by -s and -e) to be \
excluded from the schedule. This might be useful if other urgent observation from other project to be scheduled, etc. The format is: \
START1/END1,START2/END..., where START# and END# are start and end time of reserved block given in YYYY-MM-DDThh:mm:ss format", default="", type='str')
        cmd.add_option('--plot', action="store_true", dest='is_plot',
                help="make Time-vs-Elevation plot for scheduled pulsars", default=False)
        cmd.add_option('--sched', action="store_true", dest='is_sched',
                help="create XML files for scheduled pulsars", default=False)
        cmd.add_option('-g', '--gap', dest='gap', metavar='MIN',
                           help="specify gap (in min) between scans. Default = %default", default="1", type='int')
        cmd.add_option('--horizon', dest='horizon', metavar='DEG',
                           help="specify minimum elevation (in deg) for sources. Default = %default", default="20", type='float')
        cmd.add_option('--transit-tolerance', dest='transit_tolerance', metavar='MIN',
                           help="specify the maximum possible delay/advance (in mins) of scan midtime \
from the transit. Default = %default", default="120", type='int')
        cmd.add_option('--latitude', dest='latitude', metavar='DD:MM:SS.SS',
                           help="specify the latitude of a site. Default is for LOFAR (%default)", default="52:52:59.88", type='str')
        cmd.add_option('--longitude', dest='longitude', metavar='DD:MM:SS.SS',
                           help="specify the longitude of a site. Default is for LOFAR (%default)", default="06:52:00.12", type='str')

        # reading cmd options
        (opts, args) = cmd.parse_args()

        # check if any input parameters are given
        if len(sys.argv[1:]) == 0:
                cmd.print_usage()
                sys.exit(0)

        # check if start and end date/times are given
        if (opts.startTime == "" or opts.endTime == ""):
                print "Either start or end date/time of a time slot is not specified!"
                sys.exit(1)

	# check if input fule is given
	if opts.targets_file == "":
		print "Need to provide input source list to make a schedule. Use option -f for this."
		sys.exit(1)

	if opts.project == "" and opts.is_sched:
		print "Need to give the project code before making an XML schedule!"
		sys.exit(1)

	if opts.gap < 1: opts.gap = 1
	if opts.transit_tolerance < 1: opts.transit_tolerance = 1

	# Get the start and end time of the observation from the command line
	start = ephem.Date(opts.startTime.replace("-", "/").replace("T", " "))
	end = ephem.Date(opts.endTime.replace("-", "/").replace("T", " "))

#	targets = np.loadtxt(opts.targets_file, dtype=str, unpack=True, comments='#')
#	targets = targets.T # lines are for different pulsars now
	targets = []
	f = open(opts.targets_file, 'r')
	plines = f.read().splitlines()	
	f.close()
	for line in plines:
		if line[0] == "#": continue
		targets.append(line.split())

	# Get the Schedule
	sched = psrSched(targets, start, end, opts.latitude, opts.longitude, opts.horizon, \
			opts.transit_tolerance, opts.gap, opts.project, opts.is_plot, opts.reserved_blocks)
	sched.print_schedule()
	(schedpsrs, start_times, end_times, alts, advanced) = sched.get_schedule()

	# XML Part
	if opts.is_sched:
		# reading targets file to make an updated one with scheduled pulsars commented out
		if re.search("-after-", opts.targets_file):
			upd_targets_file = "%s-after-%s" % (opts.targets_file.split("-after-")[0], start.datetime().strftime("%Y-%m-%d"))
		else:
			upd_targets_file = "%s-after-%s" % (opts.targets_file, start.datetime().strftime("%Y-%m-%d"))
                f = open(upd_targets_file, 'w')
		for line in plines:
			if line.split()[0] in [s[0] for s in schedpsrs]:
				f.write("#%s\n" % line)
			else: f.write("%s\n" % line)
		f.close()

		# will be writing first the individual xml-files to temporary directory
		tempdir=".tmp-%s" % (''.join(random.choice(string.lowercase) for iii in range(7)))
		cmd="mkdir -p %s" % (tempdir)
		os.system(cmd)
		
		# making first a separate xml-file for each scheduled pulsars
		for ii in np.arange(len(schedpsrs)):
		        start_scan_str = sched.round_time(start_times[ii]).datetime().strftime("%Y-%m-%dT%H:%M:%S")
        		end_scan_str = sched.round_time(end_times[ii]).datetime().strftime("%Y-%m-%dT%H:%M:%S")

			if len(schedpsrs[ii]) > 5:
				if schedpsrs[ii][5] != "0": template="%s" % (schedpsrs[ii][5])
				else: template = opts.template
			else: template = opts.template
			
			if template == "":
				print "No XML template is given!\nEither use -t option or specify template in the input source list"
				sys.exit(1)

			if len(schedpsrs[ii]) > 6: 
				if schedpsrs[ii][6] == "-2": bad_stations_str = opts.bad_stations
				else: bad_stations_str = schedpsrs[ii][6]
			else: bad_stations_str = opts.bad_stations

			xmltje = xmlSched(template, opts.project, ii, schedpsrs[ii][0], schedpsrs[ii][3], start_scan_str, end_scan_str, \
					schedpsrs[ii][1], schedpsrs[ii][2], bad_stations_str)
			# update xml
			xmltje.update()

			# writing out XML file for the pulsar "psr"
			xmltje.write("%s/%s-%d.xml" % (tempdir, schedpsrs[ii][0], ii))

			# run extra command when needed to tweak the xml-file
			if len(schedpsrs[ii]) > 7:
				cmd="%s/%s %s/%s.xml" % (os.getcwd(), schedpsrs[ii][7], tempdir, schedpsrs[ii][0])
				os.system(cmd)

		# making full single XML file that includes all scheduled pulsars
		xmlsched="%s-sched-%s.xml" % (opts.project, start.datetime().strftime("%Y-%m-%d"))		
		if len(schedpsrs) > 0:
			out = open(xmlsched, 'w')
			for ii in np.arange(len(schedpsrs)):
				psr = schedpsrs[ii][0]
				px = open("%s/%s-%d.xml" % (tempdir, psr, ii), 'r')
				xmllines = px.read().splitlines()
				px.close()
				if ii == 0: # in the first file exclude only 2 last lines (if no folder is used)
					if opts.folder == "":
						if len(schedpsrs) == 1: outline="\n".join(xmllines)
						else: outline="\n".join(xmllines[:-2])
					else:
						outline="\n".join(xmllines[:5])
						# for stupid reason MoM does not support descriptions longer than 256 symbols, so we first merge all pulsars together
						# then we just take first 253 symbols, then exclude last truncated pulsar and then add ",..." in the end
						descr=", ".join([p[0] for p in schedpsrs])
						if len(descr) > 256: descr=", ".join(descr[:253].split(", ")[:-1]) + ",..."
						outline+="\n\t<item index=\"0\">\n\t<lofar:folder topology_parent=\"true\">\n\t\t<topology>0</topology>\n\t\t<name>%s</name>\n\t\t<description>%s</description>\n\t\t<children>\n" % \
							(opts.folder, descr)
						outline+="\n".join(xmllines[5:-2])
						if len(schedpsrs) == 1:
							outline+="\n</children>\n</lofar:folder>\n</item>\n"
							outline+="\n".join(xmllines[-2:])
				elif ii == len(schedpsrs) - 1: # in the last file exclude only the first 5 lines (if no folder is used)
					if opts.folder == "":
						outline="\n".join(xmllines[5:])
					else:
						outline="\n".join(xmllines[5:-2])
						outline+="\n</children>\n</lofar:folder>\n</item>\n"
						outline+="\n".join(xmllines[-2:])
				else: # in all other files exclude both first 5 and last 2 lines
					outline="\n".join(xmllines[5:-2])
				out.write("%s\n" % outline)	
			out.close()		

		# removing temporary directory
		cmd="rm -rf %s" % (tempdir)
		os.system(cmd)
