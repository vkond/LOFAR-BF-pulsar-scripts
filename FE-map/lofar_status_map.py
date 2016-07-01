#!/usr/bin/python
#
# Script to make a "status" maps of LOFAR stations
# based on FE observations. Scripts will create a few maps:
# for international, remote and core stations, with overlaid
# pulsar profiles for each
#
# Vlad Kondratiev, May 13, 2011 (c)
#
# Oct 16, 2015 - added settings for three new Polish stations;
#		 added printing of S/N for each used station to stdout
############################################################################

import matplotlib
matplotlib.use("Agg")
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
import re
import os, sys

# Getting LOFARSOFT location
lofarsoft=os.environ['LOFARSOFT']

chisq_file="chi-squared.txt"  # file from the pipeline
#stations_file="%s/release/share/pulsar/data/stations.txt" % (lofarsoft,) # location of George's stations.txt file
stations_file="/home/kondratiev/bin/stations.txt" # location of George's stations.txt file

# at the core location, how many meters per degree latitude,longitude?
# Assuming that http://www.csgnetwork.com/degreelenllavcalc.html is right,
# then at latitude = 52.91663 we have:
metersToLat = 1.0/111284.65883820018
metersToLong = 1.0/67266.45194985923

imw_orig=200 # original width of the png thumbnail image
imh_orig=140 # original height of the png thumbnail image
imw=143 # width of png image of the station
imh=100 # height of png image of the station
figx=12.8 # h-size of the fig in inches
figy=7.7  # v-size of the fig in inches
dpi=144
continents_color='lightgreen'

# dictionary of X-, Y- offsets of the stations thumnail plots left-bottom corners
# relative to the position of the station
# offsets are in fraction of the width and height of the plot in pixel
# (in fraction of imw and imh)
# if 3rd element = 1, then X-, Y- are actual coordinates, rather then offsets
offsets={'FR606': (-1, 0.1, 0), 'UK608': (-1, 0.1, 0), 'SE607': (0, 0.1, 0), 'DE601': (-1, -1.2, 0), 
         'DE602': (0.1, -1.1, 0), 'DE603': (-1, 0.1, 0), 'DE604': (0, 0.1, 0), 'DE605': (-1, 0.1, 0), 'DE609': (-1, 0.1, 0),
         'PL610': (-1, -1.2, 0), 'PL611': (0, 0.1, 0), 'PL612': (0, 0.1, 0),
	 'RS305': (10, 530, 1), 'RS106': (0, 0.1, 0), 'RS205': (0.3, -1, 0), 'RS208': (0, -1.2, 0), 'RS210': (-1, 0.1, 0), 
	 'RS306': (170, 530, 1), 'RS307': (-1, -1.3, 0), 'RS310': (-1, -1.3, 0), 
	 'RS404': (10, 290, 1), 'RS406': (10, 410, 1), 'RS407': (170, 410, 1), 'RS409': (-1, 0.1, 0), 'RS410': (-1, -1.3, 0), 
	 'RS503': (10, 650, 1), 'RS508': (0, 0.1, 0), 'RS509': (170, 650, 1),
	 'CS002': (1100, 650, 1), 'CS003': (1100, 530, 1), 'CS004': (1100, 410, 1),  # superterp
	 'CS005': (1100, 290, 1), 'CS006': (1100, 170, 1), 'CS007': (1100, 50, 1),   # superterp
	 'CS001': (900, 410, 1), 'CS011': (900, 530, 1), 'CS013': (0, 0.1, 0), 'CS017': (900, 650, 1),
	 'CS021': (10, 650, 1), 'CS024': (-0.4, -1.2, 0), 'CS026': (0, 0.1, 0), 'CS028': (-0.5, 0.1, 0),
	 'CS030': (-1, 0.1, 0), 'CS031': (10, 530, 1), 'CS032': (-1, -1.2, 0),
	 'CS101': (0, 0.1, 0), 'CS103': (0, -1.2, 0), 'CS201': (0, -1.2, 0), 'CS301': (-1, -1.2, 0), 'CS302': (-1.1, 0.1, 0),
	 'CS401': (-1, 0.1, 0), 'CS501': (-1, 0.4, 0)
	 }
# same for the text labels together with justification
text_offsets={'FR606': (-40000, 0, 'right'), 'UK608': (-40000, 0, 'right'), 'SE607': (40000, 0, 'left'),
              'DE601': (-40000, 0, 'right'), 'DE602': (-40000, 0, 'right'), 'DE603': (-40000, 0, 'right'), 
	      'DE604': (40000, 0, 'left'), 'DE605': (-40000, 0, 'right'), 'DE609': (-40000, 0, 'right'),
              'PL610': (-40000, 0, 'right'), 'PL611': (40000, 0, 'left'), 'PL612': (40000, 0, 'left'),
	      'RS305': (2500, 2500, 'right'), 'RS106': (-2000, 0, 'right'), 
	      'RS205': (-2000, 0, 'right'), 'RS208': (-2000, 0, 'right'), 'RS210': (-2000, 0, 'right'), 
	      'RS306': (2000, 200, 'left'), 'RS307': (-2000, 0, 'right'), 'RS310': (-2000, 0, 'right'), 
	      'RS404': (200, 2500, 'right'), 'RS406': (-2000, 0, 'right'), 'RS407': (2000, 0, 'left'),
	      'RS409': (-2000, 0, 'right'), 'RS410': (-2000, 0, 'right'),
	      'RS503': (0, -2700, 'left'), 'RS508': (2000, 0, 'left'), 'RS509': (-2000, 0, 'right'),
	      'CS002': (-50, 0, 'right'), 'CS003': (-50, 0, 'right'), 'CS004': (-50, 0, 'right'), # superterp
	      'CS005': (-30, -50, 'right'), 'CS006': (50, 0, 'left'), 'CS007': (50, 0, 'left'), # superterp
	      'CS001': (-70, 0, 'right'), 'CS011': (70, 0, 'left'), 'CS013': (70, 0, 'left'), 'CS017': (70, 0, 'left'),
	      'CS021': (0, -70, 'left'), 'CS024': (70, 0, 'left'), 'CS026': (70, 0, 'left'), 'CS028': (70, 0, 'left'),
	      'CS030': (-70, 0, 'right'), 'CS031': (0, 70, 'right'), 'CS032': (-70, 0, 'right'),
	      'CS101': (70, 0, 'left'), 'CS103': (70, 0, 'left'), 'CS201': (70, 0, 'left'), 
	      'CS301': (-70, 0, 'right'), 'CS302': (-70, 0, 'right'),
	      'CS401': (-70, 0, 'right'), 'CS501': (-70, 0, 'right')
	      }
# map parameters for different zooms
zooms={ 'intl': (-10.90, 43.4, 31.90, 60.4, 'lcc', 51, 12, 'i', 100.),
        'remote': (4.1, 52.2, 7.7, 53.5, 'lcc', 52.85, 5.9, 'h', 10.),
        'core': (6.83, 52.90, 6.93, 52.9362, 'lcc', 52.9181, 6.88, 'c', 10000.)
       }




# function that return pixel coordinates from the station coordinates
def get_pixels (x, y, figx, figy, dpi):
	factor = 100 # resolution ?
	x0, y0 = ax.transData.inverted().transform((0, 0))
	xm, ym = ax.transData.inverted().transform((figx * factor, figy * factor))
	xp = int(figx * factor * factor * (x-x0) / (xm-x0) / dpi)
	yp = int(figy * factor * factor * (y-y0) / (ym-y0) / dpi)
	return (xp, yp)

# function that return datapoints coordinates from the pixels
def get_datapoints (xp, yp, figx, figy, dpi):
	factor = 100 # resolution ?
	x0, y0 = ax.transData.inverted().transform((0, 0))
	xm, ym = ax.transData.inverted().transform((figx * factor, figy * factor))
	x = int(xp * (xm-x0) * dpi / figx / factor / factor + x0)
	y = int(yp * (ym-y0) * dpi / figy / factor / factor + y0)
	return (x, y)

# reading "stations.txt" file 
try:
	stations, colors = np.loadtxt (stations_file, comments='#', usecols=(0,1), dtype=str, unpack=True)
	lats, lons = np.loadtxt (stations_file, comments='#', usecols=(2,3), dtype=float, unpack=True)
except:
	print "Error reading %s file!" % (stations_file)
	sys.exit(1)


# Here we have to get the info about stations involved in the observation
# and then use criterion to compress all lists above to have info only
# for stations used in the observation
try:
	f = open(chisq_file, 'r')
	# ignoring also comments
	infofile = [ff for ff in f.read().splitlines() if re.search("^#", ff) is None]
	f.close()

except:
	print "Error reading %s file!" % (chisq_file)
	sys.exit(1)

pngs=[p.split("file=")[1].split()[0] for p in infofile]
used_ears=[u.split("obs=")[1].split()[0].split("_")[-2] for u in infofile]
chisq=[c.split("S/N=")[1].split()[0] for c in infofile]
station_info={} # dictionary with station name as a key, and value is the tuple with (filename, chi-squared value)
for i in np.arange(len(used_ears)):
	station_info[used_ears[i]] = ('%s' % (pngs[i]), '%.3f' % (chisq[i] != "" and float(chisq[i]) or 0))

antenna=used_ears[0][5:8]                                                         # either LBA or HBA
pulsar=pngs[0].split("_")[-1].split(".pfd")[0]                                    # pulsar name
#obsid="_".join("_".join(pngs[0].split("/")[-1].split("_")[1:]).split("_")[0:-3])  # obsid
obsid="_".join("_".join(pngs[0].split("/")[-1].split("_")[1:]).split("_")[0:-3]).split("_SAP")[0]  # obsid
used_stations=np.unique([u[0:5] for u in used_ears])

#sift_criterium = [s in used_stations for s in stations]
#stations=stations.compress(sift_criterium)
#colors=colors.compress(sift_criterium)
#lats=lats.compress(sift_criterium)
#lons=lons.compress(sift_criterium)

# some modifications
colors=[re.sub('G', 'green', c) for c in colors]
colors=[re.sub('Y', 'yellow', c) for c in colors]
colors=[re.sub('R', 'red', c) for c in colors]
for i in np.arange(len(stations)):
	if lats[i] > 100000:
		lat = 52.91663+(lons[i]-548752.0)*metersToLat
		lons[i] = 6.87018+(lats[i]-254745.0)*metersToLong
		lats[i] = lat

# also printing info about S/N to stdout
print "#"
print "# Station\tS/N"
print "#========================"

# loop over different zooms to make plots for
# international, remote, etc. stations
for z in zooms.keys():

	if z == 'intl':
		# criterium for Intl stations
		criterion = [(s[0] != 'C') & (s[0] != 'R') for s in stations]
	elif z == 'remote':
		# criterium for RS stations
		criterion = [(s[0] == 'R') for s in stations]
	elif z == 'core':
		# criterium for RS stations
		criterion = [(s[0] == 'C') for s in stations]
	else:
		print "Error: Undefined zoom!"
		sys.exit(1)	

	fig = plt.figure(figsize=(figx, figy), dpi=dpi)
	ax = fig.add_axes([0, 0, 1, 1])

	# make specific map
	map=Basemap(llcrnrlon=zooms[z][0],llcrnrlat=zooms[z][1],urcrnrlon=zooms[z][2],urcrnrlat=zooms[z][3], \
        	    projection='%s' % (zooms[z][4]), lat_0=zooms[z][5], lon_0=zooms[z][6], \
		    resolution ='%s' % (zooms[z][7]),area_thresh=zooms[z][8], ax=ax)

	# draw coastlines, country boundaries, fill continents.
	#map.bluemarble()
	map.fillcontinents(color='%s' % (continents_color), lake_color='lightblue')
	if z != 'core':
		map.drawcoastlines()
		map.drawcountries(color='brown')
		map.drawmapboundary(fill_color='lightblue')
	if z == 'remote':
		map.drawrivers(color='blue')

	# draw lat/lon grid lines every 30 degrees.
	if z == 'intl':
		map.drawmeridians(np.arange(0, 360, 30))
		map.drawparallels(np.arange(-90, 90, 30))

	x, y = map (lons.compress(criterion), lats.compress(criterion))
	cols = np.array(colors).compress(criterion)
	ss=stations.compress(criterion)

	# plotting markers of all stations for this zoom
	for i in np.arange(len(x)):
		ax.plot(x[i], y[i], color='%s' % cols[i], marker='o', markersize=10)
	map.plot(x, y, '.', color='black')

	# making necessary profile plots and overlapping them
	is_at_least_one_image = False # will be True if we have at least one plot for the current zoom
	for i in np.arange(len(x)):
		try:
			if cols[i] != 'red':
				label=ss[i]
				if antenna == 'HBA' and z == 'core':
					# checking if both HBA ears were used
					if '%s%s0' % (ss[i], antenna) in station_info and '%s%s1' % (ss[i], antenna) in station_info:
						for ear in ['0', '1']:
							os.system("convert %s -flatten -background white -crop %dx%d+25+16 -background %s -pointsize 14 label:'%s' +swap -gravity Center -append  %s.%c.map.png 2>/dev/null" % \
								(station_info['%s%s%c' % (ss[i], antenna, ear)][0], imw_orig/2, imh_orig-16, (ear == '0' and 'Khaki' or 'Plum'), station_info['%s%s%c' % (ss[i], antenna, ear)][1], ss[i], ear))
							# also printing S/N to stdout
							print "%s%s%c\t%s" % (ss[i], antenna, ear, station_info['%s%s%c' % (ss[i], antenna, ear)][1])
						os.system("convert %s.0.map.png %s.1.map.png +append -scale %dx%d-0 temp%s.map.png 2>/dev/null" % \
							(ss[i], ss[i], imw, imh, ss[i]))
						os.system("rm -f %s.[01].map.png" % (ss[i]))
					else: # only one HBA ear was used
						if '%s%s0' % (ss[i], antenna) in station_info: statkey = '%s%s0' % (ss[i], antenna)
						else: statkey = '%s%s1' % (ss[i], antenna)
						label='%s  %s%c' % (ss[i], antenna, statkey[-1])
						os.system("convert %s -flatten -background white -crop %dx%d+0+16 -background Khaki -pointsize 14 label:'%s' +swap -gravity Center -append -scale %dx%d-0 temp%s.map.png 2>/dev/null" % \
							(station_info[statkey][0], imw_orig, imh_orig-16, station_info[statkey][1], imw, imh, ss[i]))
						# also printing S/N to stdout
						print "%s\t%s" % (statkey, station_info[statkey][1])
				else: # for LBA or remore/international stations
					os.system("convert %s -flatten -background white -crop %dx%d+0+16 -background Khaki -pointsize 14 label:'%s' +swap -gravity Center -append -scale %dx%d-0 temp%s.map.png 2>/dev/null" % \
						(station_info['%s%s' % (ss[i], antenna)][0], imw_orig, imh_orig-16, station_info['%s%s' % (ss[i], antenna)][1], imw, imh, ss[i]))
					# also printing S/N to stdout
					print "%s%s\t%s" % (ss[i], antenna, station_info['%s%s' % (ss[i], antenna)][1])
				# adding the station name label
				os.system("convert temp%s.map.png -background LightBlue -pointsize 10 label:'%s' +swap -gravity Center -append  %s.map.png>/dev/null ; rm -f temp%s.map.png" % \
					(ss[i], label, ss[i], ss[i]))

				im = plt.imread('%s.map.png' % (ss[i]))	
				if offsets[ss[i]][2] == 0:
					xp, yp = get_pixels (x[i], y[i], figx, figy, dpi)
					fig.figimage(im, xp+imw*offsets[ss[i]][0], yp+imh*offsets[ss[i]][1], origin="upper", zorder=10)
				else:
					fig.figimage(im, offsets[ss[i]][0], offsets[ss[i]][1], origin="upper", zorder=10)
					# getting coords for the text label in data coords (rather than pixels)
					#xpt, ypt = get_datapoints (offsets[ss[i]][0], offsets[ss[i]][1]+imh+10, figx, figy, dpi)
					plt.text(xpt, ypt, ss[i], ha='left', va='center', color='black', fontsize=8)
				os.system("rm -f %s.map.png" % (ss[i]))
				is_at_least_one_image = True
		except: pass

	for name, xpt, ypt in zip(ss, x, y):
		plt.text(xpt+text_offsets[name][0], ypt+text_offsets[name][1], name, \
        		 ha='%s' % (text_offsets[name][2]), va='center', color='black', fontsize=8)

	# make plot
	if is_at_least_one_image:
		plt.text(0.45, 0.97, "%s   %s   %s" % (obsid, antenna, pulsar), ha='center', va='center', transform=ax.transAxes, bbox=dict(facecolor="%s" % (continents_color,), alpha=0.5), fontsize=20)
		plt.savefig("%s_%s_%s_%s_status.png" % (obsid, antenna, pulsar, z))
