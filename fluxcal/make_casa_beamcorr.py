#!/usr/bin/env python
#
# Accessory script to create the dictionary
# casa_beamcorr with CasA correction factors
# for different stations
#
# Vlad Kondratiev (c) - 05.08.2016
#
# 17.04.2017 - wrote info about the bug in antennaJones.py
#              and updated casa_beamcorr_pkg.py for RS stations;
#              also updated the stations' list
#
import os, sys, math
import numpy as np
from pyrap.measures import measures
from pyrap.quanta import quantity
import antennaJones as aj

# In order to calculate correction factors for Remote stations (RS)
# one needs to fix the bug in antennaJones.py script of mscorpol
# package in line 210 (see below).
stations=["CS001", "CS002", "CS003", "CS004", "CS005", "CS006", "CS007",
          "CS011", "CS013", "CS017", "CS021", "CS024", "CS026", "CS028",
          "CS030", "CS031", "CS032", "CS101", "CS103", "CS201", "CS301", 
          "CS302", "CS401", "CS501", "RS106", "RS205", "RS208", "RS210",
          "RS305", "RS306", "RS307", "RS310", "RS406", "RS407", "RS409", 
          "RS503", "RS508", "RS509", "RS511", "DE601", "DE602", "DE603",
          "DE604", "DE605", "FR606", "SE607", "UK608", "DE609", "FI609", 
          "PL610", "PL611", "PL612"]

# currently tables are done only for CS and INTL stations. Mscorpol package
# for whatever reasons fails to give Jones matrices for RS stations
# There is a bug in mscorpol package script antennaJones.py the line 210. 
# It has to be changed from "if stnLoc=='CS' or stnLoc=='RS':"  to  
# "if stnLoc=='CS': # or stnLoc=='RS':", i.e. to comment out 'RS' part. 
# It makes sense as all RS stations have just only one ear of 48 tiles similar 
# to International stations with 96 tiles, thus antenna must be "HBA" instead of "HBA0"
#stations=["CS001", "CS002", "CS003", "CS004", "CS005", "CS006", "CS007",
#          "CS011", "CS013", "CS017", "CS021", "CS024", "CS026", "CS028",
#          "CS030", "CS031", "CS032", "CS101", "CS103", "CS201", "CS301", 
#          "CS302", "CS401", "CS501", "DE601", "DE602", "DE603", "DE604", 
#          "DE605", "FR606", "SE607", "UK608", "DE609", "PL610", "PL611", "PL612"]

#   m a i n
if __name__=="__main__":

	subwidth=100./512.
	obstimes=quantity([55159.77650462962963, 55159.77650462962963], 'd')
	direction=measures().direction('J2000', str(6.123487681)+'rad', str(1.0265154)+'rad')

	out = open("casa_beamcorr_pkg.py", "wt")
	out.write("casa_beamcorr={")
	for ss in xrange(len(stations)):
		print "%s, #%d of %d" % (stations[ss], ss+1, len(stations))
		out.write("\"%s\" : [" % (stations[ss]))
		for ii in xrange(51, 1536, 6): # freq between 10 MHz and 300 MHz	
			outline=""
			for jj in xrange(ii, ii+6 > 1536 and 1536 or ii+6):
				if jj != ii: outline += ", "
				freq=jj*subwidth+subwidth/2. # in MHz
				Jn=aj.getJonesByAntFld("Hamaker", obstimes, stations[ss], direction, freq*1.e6)
				bc_casa = 1./np.abs(0.5*(Jn[0,0,0]*np.conj(Jn[0,0,0]) + Jn[0,0,1]*np.conj(Jn[0,0,1]) + \
						Jn[0,1,0]*np.conj(Jn[0,1,0]) + Jn[0,1,1]*np.conj(Jn[0,1,1])))
				outline+="[%lf, %lf]" % (freq, bc_casa)
			if ii+6<1536: outline += ",\n"
			else: 
				if ss + 1 == len(stations): outline += "]\n"
				else: outline += "],\n"
			out.write(outline)
	out.write("}\n")
	out.close()
