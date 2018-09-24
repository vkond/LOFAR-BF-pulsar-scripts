#!/usr/bin/env python

import os, sys
import time
from datetime import datetime
import commands
from numpy import arange as arange
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.dates import YearLocator, MonthLocator, DateFormatter
from optparse import OptionParser

######################
#                     #
#      Utilities      #
#                     #
#######################
def version_string():
    r'''
    Returns version number as a string.

    **Example**

    >>> version_string()
    '0.90'
    '''
    return '0.90'

def get_stations():
    full_stations_list = []
    st_list = []
    cs_list = []
    rs_list = []
    int_list = []
    f = file('stations.txt', "r")
    for i,line in enumerate(f):
        full_stations_list.append(line.strip()[:-1])
        if line.strip()[:-1] in ['cs002','cs003','cs004','cs005','cs006','cs007']:
            st_list.append(line.strip()[:-1])
        elif line.strip()[:-1] not in ['cs002','cs003','cs004','cs005','cs006','cs007'] and line.strip()[0] == 'c':
            cs_list.append(line.strip()[:-1])
        elif line.strip()[0] == 'r':
            rs_list.append(line.strip()[:-1])
        else:
            int_list.append(line.strip()[:-1])
    f.close
    return (full_stations_list,st_list,cs_list,rs_list,int_list)

def create_history(databasefile):
    # Get a station list
    ant_dbase = {}
    rcu_dbase = {}
    full_stations_list = []
    f = file('stations.txt', "r")
    for i,line in enumerate(f):
        full_stations_list.append(line.strip()[:-1])
    f.close
    for station in full_stations_list:
        station = station.upper()
        ant_dbase[station[0:5]] = {'LBA': {},
                                    'HBA': {}
                                   }
        rcu_dbase[station[0:5]] = {}
        if 'RS' not in station and 'CS' not in station:
            for i in range(193):
                rcu_dbase[station[0:5]]['{0}'.format(i)] = {}
            for i in range(97):
                ant_dbase[station[0:5]]['LBA']['{0:0>3}'.format(i)] = {}
                ant_dbase[station[0:5]]['HBA']['{0:0>2}'.format(i)] = {}
        else:
            for i in range(97):
                ant_dbase[station[0:5]]['LBA']['{0:0>3}'.format(i)] = {}
                rcu_dbase[station[0:5]]['{0}'.format(i)] = {}
            for i in range(48):
                ant_dbase[station[0:5]]['HBA']['{0:0>2}'.format(i)] = {}

    #f = file('DATA/hardware_states_latest.txt', "r")
    f = file(databasefile, "r")

    for line in f:
        if 'LOFAR.PIC' in line:
            if 'RCU' in line:
                cols = line.split('|')
                subcols_0 = cols[0].split('.')
                rcu_dbase[subcols_0[3]][subcols_0[7][3:]][cols[2].strip()] = cols[1].strip()
            else:
                cols = line.split('|')
                subcols_0 = cols[0].split('.')
                ant_dbase[subcols_0[3]][subcols_0[4][0:3]][subcols_0[4][3:]][cols[2].strip()] = cols[1].strip()
    f.close
    return ant_dbase, rcu_dbase

def create_mapping():
    modes = ['LBA-INNER','LBA-OUTER','HBA-DUAL','HBA-DUAL-INNER']
#    modes = ['LBA-INNER','LBA-OUTER','HBA-DUAL','HBA-DUAL-INNER', 'HBA-ZERO', 'HBA-ONE', 'HBA-JOINED']
    station_sets = ['CS','RS','INT']
    mapping_ant = {}
    mapping_rcu = {}
    for st_set in station_sets:
        mapping_ant[st_set] = {}
        mapping_rcu[st_set] = {}
        if st_set in ['CS']:
            mapping_ant[st_set][modes[0]] = [i for i in range(48)]
            mapping_rcu[st_set][modes[0]] = [[2*i,2*i+1] for i in mapping_ant[st_set][modes[0]]]
            mapping_ant[st_set][modes[1]] = [i+48 for i in range(48)]
            mapping_rcu[st_set][modes[1]] = [[2*(i-48)+1,2*(i-48)] for i in mapping_ant[st_set][modes[1]]]
            mapping_ant[st_set][modes[2]] = [i for i in range(48)]
            mapping_rcu[st_set][modes[2]] = [[2*i,2*i+1] for i in mapping_ant[st_set][modes[2]]]
            mapping_ant[st_set][modes[3]] = [i for i in range(48)]
            mapping_rcu[st_set][modes[3]] = [[2*i,2*i+1] for i in mapping_ant[st_set][modes[3]]]
        if st_set in ['RS']:
            mapping_ant[st_set][modes[0]] = [i for i in range(48)]
            mapping_rcu[st_set][modes[0]] = [[2*i,2*i+1] for i in mapping_ant[st_set][modes[0]]]
            mapping_ant[st_set][modes[1]] = [i+48 for i in range(48)]
            mapping_rcu[st_set][modes[1]] = [[2*(i-48)+1,2*(i-48)] for i in mapping_ant[st_set][modes[1]]]
            mapping_ant[st_set][modes[2]] = [i for i in range(48)]
            mapping_rcu[st_set][modes[2]] = [[2*i,2*i+1] for i in mapping_ant[st_set][modes[2]]]
            mapping_ant[st_set][modes[3]] = [5,6,10,11,12,13,17,18,19,20,21,22,25,26,27,28,29,30,34,35,36,37,41,42]
            mapping_rcu[st_set][modes[3]] = [[2*i,2*i+1] for i in mapping_ant[st_set][modes[3]]]
        if st_set in ['INT']:
            mapping_ant[st_set][modes[0]] = [i for i in range(96)]
            mapping_rcu[st_set][modes[0]] = [[2*i,2*i+1] for i in mapping_ant[st_set][modes[0]]]
            mapping_ant[st_set][modes[1]] = [i for i in range(96)]
            mapping_rcu[st_set][modes[1]] = [[2*i,2*i+1] for i in mapping_ant[st_set][modes[1]]]
            mapping_ant[st_set][modes[2]] = [i for i in range(96)]
            mapping_rcu[st_set][modes[2]] = [[2*i,2*i+1] for i in mapping_ant[st_set][modes[2]]]
            mapping_ant[st_set][modes[3]] = [i for i in range(96)]
            mapping_rcu[st_set][modes[3]] = [[2*i,2*i+1] for i in mapping_ant[st_set][modes[3]]]
    return mapping_ant, mapping_rcu, modes, station_sets

def sorted_lists(x,y):
    x, y = zip(*sorted(zip(x, y)))
    x, y = (list(t) for t in zip(*sorted(zip(x, y))))
    return x,y

def parse_command_line(argv):
    mode_to_check = 'HBA-DUAL'
    station_list = 'st'
    now = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    parser  = OptionParser(usage   = 'python %prog [options]',
                           version = '%prog '+version_string())

    parser.add_option('--test', dest    = 'test',
                help    = 'Use in test mode',
                default = False, action = 'store_true')

    parser.add_option('-s', dest    = 'station_list',
                help    = 'Use stationlist instead of all stations',
                metavar = 'STATION_LIST', default = station_list)

    parser.add_option('-m', dest    = 'mode_to_check',
                help    = 'Check for mode <HBA> or <LBA>',
                metavar = 'MODE_TO_CHECK', default = mode_to_check)

    parser.add_option('-d', dest    = 'set_date',
                help    = 'Use specified date instead of {0}'.format(now),
                      metavar = 'SET_DATE', default = now)

    parser.add_option('-i', dest    = 'input_file',
                help    = 'Use specified date instead of date set or {0}'.format(now),
                      metavar = 'INPUT_FILE', default = 'date')

    parser.add_option('--db', dest    = 'databasefile',
                help    = 'Use specified databasefile instead of DATA/hardware_states_latest.txt',
                      metavar = 'DATABASEFILE_FILE', default = 'DATA/hardware_states_latest.txt')

    parser.add_option('-f', dest    = 'file_name',
                help    = 'Append results to filename',
                      metavar = 'FILE_NAME')

    (options, args) = parser.parse_args(argv[1:])

    return options


def my_main(argv):
#if __name__ == '__main__':
    r'''
    Main routine for my_program.

    **Parameters**

    argv : list of strings
        The command line arguments passed to the script. Basically the
        contents of ``sys.argv``.

    **Returns**
    An integer signifying success if 0, and an error otherwise.

    '''

    argv=sys.argv
    options = parse_command_line(argv)
    if options.file_name is not None:
        if os.path.isfile(options.file_name) is True:
            file_name = file(options.file_name,'a')
        else:
            file_name = file(options.file_name,'w')
    mode_to_check = options.mode_to_check
    full_stations_list,st_list,cs_list,rs_list,int_list = get_stations()
    if options.station_list == 'st':
        stations = st_list
    elif options.station_list == 'core':
        stations = st_list+cs_list
    elif options.station_list == 'remote':
        stations = rs_list
    elif options.station_list == 'int':
        stations = int_list
    elif options.station_list == 'all':
        stations = full_stations_list
    elif ',' in options.station_list:
        stations = [i.upper() for i in options.station_list]
    else:
        stations = [options.station_list.upper(),]
    if options.input_file != 'date':
        import h5py
        assert 'h5' in options.input_file
        inputfile=h5py.File(options.input_file,'r')
        proc_date = datetime.strptime(inputfile.attrs['OBSERVATION_START_UTC'],'%Y-%m-%dT%H:%M:%S.000000000Z')
        inputfile.close()
    else:
        proc_date = datetime.strptime(options.set_date,'%Y-%m-%d %H:%M:%S')

    ant_dbase,rcu_dbase = create_history(options.databasefile)

#    for key in sorted(ant_dbase.keys()):
#        if '001' in key:
#            for s1key in sorted(ant_dbase[key].keys()):
#                for s2key in sorted(ant_dbase[key][s1key].keys()):
#                    for s3key in sorted(ant_dbase[key][s1key][s2key].keys()):
#                        try:
#                            date_tag = datetime.strptime(s3key,'%Y-%m-%d %H:%M:%S.%f')
#                        except(ValueError):
#                            date_tag = datetime.strptime(s3key,'%Y-%m-%d %H:%M:%S')
#                        datetime.date(int(date_tag[0:4]),int(date_tag[4:6]),int(date_tag[6:8]))
#                        ant_dbase[key][s1key][s2key][s3key]
#                        print '{0} {1} {2} {3} {4} {5}'.format(key, s1key, s2key, s3key,ant_dbase[key][s1key][s2key][s3key], date_tag)

#    for key in sorted(rcu_dbase.keys()):
#        if '001' in key:
#            for s1key in sorted(rcu_dbase[key].keys()):
#                for s2key in sorted(rcu_dbase[key][s1key].keys()):
#                    try:
#                        date_tag = datetime.strptime(s2key,'%Y-%m-%d %H:%M:%S.%f')
#                    except(ValueError):
#                        date_tag = datetime.strptime(s2key,'%Y-%m-%d %H:%M:%S')
                        #datetime.date(int(datetag[0:4]),int(datetag[4:6]),int(datetag[6:8]))
                        #rcu_dbase[key][s1key][s2key][s3key]
                    #print '{0} {1} {2} {3} {4}'.format(key, s1key, s2key,rcu_dbase[key][s1key][s2key], date_tag)

    modes = ['HBA','LBA-INNER','LBA-OUTER']
    mode = mode_to_check
    print 75*'#'
    print 'Checking for mode: {0} on date: {1} for stations:'.format(mode_to_check, proc_date)
    print stations
    print 75*'#'

    mapping_ant,mapping_rcu,modes,station_sets = create_mapping()

    for station in stations:
        station = station.upper()
        station_type = station[0:2]
        if station_type not in ['CS','RS']: station_type = 'INT'
#        f_pvss = file('DATA/latest_PVSS/{0}_latest_PVSS.dict'.format(station), "w")
        antennas_off = []
        rcus_off = []
        for element in mapping_ant[station_type][mode]:
#            if element != 8: continue
            state = 'on'
            if mode[0:3]=='LBA':
                elementname='{0:0>3}'.format(element)
            else:
                elementname='{0:0>2}'.format(element)

            for i in sorted(ant_dbase[station][mode[0:3]][elementname]):
                try:
                    use_date = datetime.strptime('{0}'.format(i),'%Y-%m-%d %H:%M:%S.%f')
                except(ValueError):
                    use_date = datetime.strptime('{0}'.format(i),'%Y-%m-%d %H:%M:%S')
                use_value = ant_dbase[station][mode[0:3]][elementname][i]
#                print station,element,i,ant_dbase[station][mode[0:3]][elementname][i]
                if use_date < proc_date:
                    if int(use_value) > 10:
#                        print 'date {0} is valid and value {1} suggests OFF'.format(i,use_value)
                        state = 'off'
                    elif int(use_value) <= 10:
#                        print 'date {0} is valid and value {1} suggests ON'.format(i,use_value)
                        state = 'on'
                elif use_date > proc_date:
                    if state == 'off':
#                        print 'OFF state detected last, adding to ANTENNAS OFF'
			if '{0:0>2}'.format(element) not in antennas_off:
                            antennas_off.append('{0:0>2}'.format(element))
                    continue
                if i == sorted(ant_dbase[station][mode[0:3]][elementname])[-1]:
                    if state == 'off':
#                        print 'OFF state detected last, adding to ANTENNAS OFF'
			if '{0:0>2}'.format(element) not in antennas_off:
                            antennas_off.append('{0:0>2}'.format(element))
#                print station,element,i,ant_dbase[station][mode[0:3]]['{0:0>2}'.format(element)][i],state
        if station_type == 'CS' and (mode_to_check == 'HBA-DUAL' or mode_to_check == 'HBA_DUAL_INNER'):# or mode_to_check == 'LBA-INNER' or mode_to_check == 'LBA-OUTER'):
            prstr1 = ['{0} {1}'.format(i,mapping_rcu[station_type][mode_to_check][int(i)]) for i in [k for k in antennas_off if int(k) <= 23]]
            prstr2 = ['{0} {1}'.format(i,mapping_rcu[station_type][mode_to_check][int(i)]) for i in [k for k in antennas_off if int(k) > 23]]
            string2print = '{0}{1}0 -- Antennas OFF: {2} namely: {3}\n'.format(station,mode[0:3],len([k for k in antennas_off if int(k) <= 23]),prstr1)
            string2printf = '{0} {1}{2}0 {3}\n'.format(proc_date,station,mode[0:3],len([k for k in antennas_off if int(k) <= 23]))
            string2print = string2print + '{0}{1}1 -- Antennas OFF: {2} namely: {3}'.format(station,mode[0:3],len([k for k in antennas_off if int(k) > 23]),prstr2)
            string2printf = string2printf + '{0} {1}{2}1 {3}'.format(proc_date,station,mode[0:3],len([k for k in antennas_off if int(k) > 23]))
        elif station_type == 'RS' or mode_to_check == 'LBA-INNER':# or mode_to_check == 'LBA-OUTER':
            prstr1 = ['{0} {1}'.format(i,mapping_rcu[station_type][mode_to_check][int(i)]) for i in [k for k in antennas_off]]
            string2print = '{0}{1}0 -- Antennas OFF: {2}'.format(station,mode[0:3],prstr1)
            string2printf = '{0} {1}{2}1 {3}'.format(proc_date,station,mode[0:3],len([k for k in antennas_off]))
        elif mode_to_check == 'LBA-OUTER':
            prstr1 = ['{0} {1}'.format(i,mapping_rcu[station_type][mode_to_check][int(i)-48]) for i in [k for k in antennas_off]]
            string2print = '{0}{1}0 -- Antennas OFF: {2}'.format(station,mode[0:3],prstr1)
            string2printf = '{0} {1}{2}1 {3}'.format(proc_date,station,mode[0:3],len([k for k in antennas_off]))
        print string2print
        #if station_type == 'CS': # and (mode_to_check == 'HBA-DUAL' or mode_to_check == 'HBA_DUAL_INNER'):
        #    prstr1 = ['{0} {1}'.format(i,mapping_rcu[station_type][mode_to_check][int(i)]) for i in [k for k in antennas_off if int(k) <= 23]]
        #    prstr2 = ['{0} {1}'.format(i,mapping_rcu[station_type][mode_to_check][int(i)]) for i in [k for k in antennas_off if int(k) > 23]]
        #    string2print = '{0}{1}0 -- Antennas OFF: {2} namely: {3}\n'.format(station,mode[0:3],len([k for k in antennas_off if int(k) <= 23]),prstr1)
        #    string2printf = '{0} {1}{2}0 {3}\n'.format(proc_date,station,mode[0:3],len([k for k in antennas_off if int(k) <= 23]))
        #    string2print = string2print + '{0}{1}1 -- Antennas OFF: {2} namely: {3}'.format(station,mode[0:3],len([k for k in antennas_off if int(k) > 23]),prstr2)
        #    string2printf = string2printf + '{0} {1}{2}1 {3}'.format(proc_date,station,mode[0:3],len([k for k in antennas_off if int(k) > 23]))
        #elif station_type == 'RS':
        #    prstr1 = ['{0} {1}'.format(i,mapping_rcu[station_type][mode_to_check][int(i)]) for i in [k for k in antennas_off]]
        #    string2print = '{0}{1}0 -- Antennas OFF: {2}'.format(station,mode[0:3],prstr1)
        #print string2print
        if options.file_name is None:
            continue
        else:
            print >>file_name,string2printf
    if options.file_name is None:
        print 'Nothing stored to file'
    exit(0)



    total_modes = ['HBA','LBA','LBA-INNER','LBA-OUTER']
    total_hardware = [2592,4320,1824,1824]
    x_tot = {total_modes[0]: [],
            total_modes[1]: [],
            total_modes[2]: [],
            total_modes[3]: []}
    y_tot = {total_modes[0]: [],
            total_modes[1]: [],
            total_modes[2]: [],
            total_modes[3]: []}
    y_diff_tot = {total_modes[0]: [],
            total_modes[1]: [],
            total_modes[2]: [],
            total_modes[3]: []}
    for i,total_mode in enumerate(total_modes):
        y_diff_tot[total_mode].append(0)
        x_tot[total_mode].append(datetime.strptime('2010-10-01 00:00:00','%Y-%m-%d %H:%M:%S'))
        y_tot[total_mode].append(total_hardware[i])


    mapping_ant,mapping_rcu,modes,station_sets = create_mapping()
    for station in st_list+cs_list+rs_list+int_list:
#    for station in ['cs032','cs021','rs508']:
        # Open a file to save the latest PVSS status in (only broken element are written
        # per mode, all stations in a separate file)
        f_pvss = file('DATA/latest_PVSS/{0}_latest_PVSS.dict'.format(station), "w")
#    for station in ['CS001']: #['CS001','CS006','RS210']:#stations:
        station = station.upper()
        # open figure at station level
        fig = plt.figure(figsize=(12,6))
#        modes = ['LBA-INNER','LBA-OUTER','HBA-DUAL','HBA-DUAL-INNER']
        for mode in modes:
#            print 'ant mapping (len={0} x1)'.format(len(mapping_ant[station[0:2]][mode])), mapping_ant[station[0:2].upper()][mode]
#            print 'rcu mapping (len={0} x2)'.format(len(mapping_rcu[station[0:2].upper()][mode])), mapping_rcu[station[0:2].upper()][mode]
            x = []
            y = []
            y_diff = []
#            print station, mode, x, y,

######### START MODE ########## INT INNER/OUTER the same
            if mode == 'LBA-OUTER' and station[2:3] == '6':
                print station,mode
                station_latest_broken = {station: {mode:{}}}
                x.append(datetime.strptime('2010-10-01 00:00:00','%Y-%m-%d %H:%M:%S'))
                y.append(len(mapping_ant['INT'][mode]))
                y_diff.append(0)
                rect = 0.06,0.29,0.9,0.19
                ax = ax = fig.add_axes(rect)
#                print 'x0:',x,'y0:',y,'y_diff0:',y_diff
                for element in mapping_ant['INT'][mode]:
#                    print mapping_ant[station[0:2]][mode]
#                    print '#: ', element, ' -- total possible: ', y[0]
                    tmp_broken = False
                    date_entries = {}
#                    print ant_dbase[station]['LBA']['{0:0>3}'.format(element)]
                    for i in ant_dbase[station]['LBA']['{0:0>3}'.format(element)]:
                        date_entries[i] = 'ANT'
                        # Make sure the mapping from RCU to element is done correctly
                    for i in rcu_dbase[station]['{0}'.format(2*(element)+1)]:
                        date_entries[i] = 'RCU'
                    for i in rcu_dbase[station]['{0}'.format(2*(element))]:
                        date_entries[i] = 'RCU'
                    for date_entry in sorted(date_entries):
                        if date_entries[date_entry] == 'ANT':
                            if int(ant_dbase[station]['LBA']['{0:0>3}'.format(element)][date_entry]) > 10:
                                broken = True
                            else: broken = False
                        if date_entries[date_entry] == 'RCU':
                            try:
                                if int(rcu_dbase[station]['{0}'.format(2*(element)+1)][date_entry]) > 10:
                                    broken = True
                            except(KeyError):
                                if int(rcu_dbase[station]['{0}'.format(2*(element))][date_entry]) > 10:
                                    broken = True
                            else: broken = False
                        try:
                            x.append(datetime.strptime(date_entry,'%Y-%m-%d %H:%M:%S.%f'))
                            x_tot[total_modes[1]].append(datetime.strptime(date_entry,'%Y-%m-%d %H:%M:%S.%f'))
                        except(ValueError):
                            x.append(datetime.strptime(date_entry,'%Y-%m-%d %H:%M:%S'))
                            x_tot[total_modes[1]].append(datetime.strptime(date_entry,'%Y-%m-%d %H:%M:%S'))
                        if broken == tmp_broken: tmp = 0
                        if broken == True & broken != tmp_broken: tmp = -1
                        if broken == False & broken != tmp_broken: tmp = 1
                        y_diff.append(tmp)
                        y_diff_tot[total_modes[1]].append(tmp)
                        tmp_broken = broken
                    if broken == True and date_entry == sorted(date_entries.keys())[-1]:
                        station_latest_broken[station][mode][element] = True
                x, y_diff = sorted_lists(x,y_diff)
                for i,iy_diff in enumerate(y_diff[1:]):
                    y.append(y[i]+iy_diff)
                print len(x),len(y_diff),len(y)
                x.append(datetime.strptime('2015-02-10 00:00:00','%Y-%m-%d %H:%M:%S'))
                y.append(y[-1])
                xplot, yplot = create_step_vars(x,y)
                title = ''#{0} :-: {1}'.format(station,'LBA OUTER')
                create_plot(ax,xplot,yplot,title,'orange','','# antennas',82)
                create_legend(ax,'LBA{0}'.format(''),'orange',0)
                x = [x[0]]
                y = [y[0]]
                y_diff = [0]
                del xplot,yplot
                print >>f_pvss, station_latest_broken
                print '========================================='
                        ######### END MODE ##########

######### START MODE ########## INT HBA modes the same
            if mode == 'HBA-DUAL-INNER' and station[2:3] == '6':
                station_latest_broken = {station: {mode:{}}}
                x.append(datetime.strptime('2010-10-01 00:00:00','%Y-%m-%d %H:%M:%S'))
                y.append(len(mapping_ant['INT'][mode]))
                y_diff.append(0)
                rect = 0.06,0.77,0.9,0.19
                ax = ax = fig.add_axes(rect)
                print x,y,y_diff
                for element in mapping_ant['INT'][mode]:
                    tmp_broken = False
                    date_entries = {}
                    for i in ant_dbase[station]['HBA']['{0:0>2}'.format(element)]:
                        date_entries[i] = 'ANT'
                    # Make sure the mapping from RCU to element is done correctly
                    for i in rcu_dbase[station]['{0}'.format(element*2)]:
                        date_entries[i] = 'RCU'
                    for i in rcu_dbase[station]['{0}'.format(1+element*2)]:
                        date_entries[i] = 'RCU'
                    for date_entry in sorted(date_entries):
                        if date_entries[date_entry] == 'ANT':
                            if int(ant_dbase[station]['HBA']['{0:0>2}'.format(element)][date_entry]) > 10:
                                broken = True
                            else: broken = False
                        if date_entries[date_entry] == 'RCU':
                            try:
                                if int(rcu_dbase[station]['{0}'.format(element*2)][date_entry]) > 10:
                                    broken = True
                            except(KeyError):
                                if int(rcu_dbase[station]['{0}'.format(1+element*2)][date_entry]) > 10:
                                    broken = True
                            else: broken = False
                        try:
                            x.append(datetime.strptime(date_entry,'%Y-%m-%d %H:%M:%S.%f'))
                            x_tot[total_modes[0]].append(datetime.strptime(date_entry,'%Y-%m-%d %H:%M:%S.%f'))
                        except(ValueError):
                            x.append(datetime.strptime(date_entry,'%Y-%m-%d %H:%M:%S'))
                            x_tot[total_modes[0]].append(datetime.strptime(date_entry,'%Y-%m-%d %H:%M:%S'))
                        if broken == tmp_broken: tmp = 0
                        if broken == True & broken != tmp_broken: tmp = -1
                        if broken == False & broken != tmp_broken: tmp = 1
                        y_diff.append(tmp)
                        y_diff_tot[total_modes[0]].append(tmp)
#                        print date_entry, ant_dbase[station]['HBA']['{0:0>2}'.format(element)][date_entry], 'prev:',tmp_broken,'now:',broken, 'diff:',tmp
                        tmp_broken = broken
                        if broken == True and date_entry == sorted(date_entries.keys())[-1]:
                            station_latest_broken[station][mode][element] = True
                x, y_diff = sorted_lists(x,y_diff)
                for i,iy_diff in enumerate(y_diff[1:]):
                    y.append(y[i]+iy_diff)
                x.append(datetime.strptime('2015-02-10 00:00:00','%Y-%m-%d %H:%M:%S'))
                y.append(y[-1])
                xplot, yplot = create_step_vars(x,y)
                title = '{0}'.format(station) #{0} :-: {1}'.format(station,'HBA-INNER')
                create_plot(ax,xplot,yplot,title,'red','','# tiles',82)
                create_legend(ax,'HBA'.format(''),'red',0)
                x = [x[0]]
                y = [y[0]]
                y_diff = [0]
                del xplot,yplot
                print >>f_pvss, station_latest_broken
                print '========================================='
######### END MODE ##########



######### START MODE ##########
            if mode == 'LBA-OUTER' and station[2:3] != '6':
                print station,mode,mapping_ant[station[0:2]][mode]
                station_latest_broken = {station: {mode:{}}}
                x.append(datetime.strptime('2010-10-01 00:00:00','%Y-%m-%d %H:%M:%S'))
                y.append(len(mapping_ant[station[0:2]][mode]))
                y_diff.append(0)
                rect = 0.06,0.05,0.9,0.19
                ax = ax = fig.add_axes(rect)
#                print 'x0:',x,'y0:',y,'y_diff0:',y_diff
                for element in mapping_ant[station[0:2]][mode]:
#                    print mapping_ant[station[0:2]][mode]
#                    print '#: ', element, ' -- total possible: ', y[0]
                    tmp_broken = False
                    date_entries = {}
                    for i in ant_dbase[station]['LBA']['{0:0>3}'.format(element)]:
                        date_entries[i] = 'ANT'
                        # Make sure the mapping from RCU to element is done correctly
                    for i in rcu_dbase[station]['{0}'.format(2*(element-48)+1)]:
                        date_entries[i] = 'RCU'
                    for i in rcu_dbase[station]['{0}'.format(2*(element-48))]:
                        date_entries[i] = 'RCU'
                    for date_entry in sorted(date_entries):
                        if date_entries[date_entry] == 'ANT':
                            if int(ant_dbase[station]['LBA']['{0:0>3}'.format(element)][date_entry]) > 10:
                                broken = True
                            else: broken = False
                        if date_entries[date_entry] == 'RCU':
                            try:
                                if int(rcu_dbase[station]['{0}'.format(2*(element-48)+1)][date_entry]) > 10:
                                    broken = True
                            except(KeyError):
                                if int(rcu_dbase[station]['{0}'.format(2*(element-48))][date_entry]) > 10:
                                    broken = True
                            else: broken = False
                        try:
                            x.append(datetime.strptime(date_entry,'%Y-%m-%d %H:%M:%S.%f'))
                            x_tot[total_modes[3]].append(datetime.strptime(date_entry,'%Y-%m-%d %H:%M:%S.%f'))
                            x_tot[total_modes[1]].append(datetime.strptime(date_entry,'%Y-%m-%d %H:%M:%S.%f'))
                        except(ValueError):
                            x.append(datetime.strptime(date_entry,'%Y-%m-%d %H:%M:%S'))
                            x_tot[total_modes[3]].append(datetime.strptime(date_entry,'%Y-%m-%d %H:%M:%S'))
                            x_tot[total_modes[1]].append(datetime.strptime(date_entry,'%Y-%m-%d %H:%M:%S'))
                        if broken == tmp_broken: tmp = 0
                        if broken == True & broken != tmp_broken: tmp = -1
                        if broken == False & broken != tmp_broken: tmp = 1
                        y_diff.append(tmp)
                        y_diff_tot[total_modes[3]].append(tmp)
                        y_diff_tot[total_modes[1]].append(tmp)
                        tmp_broken = broken
                        if broken == True and date_entry == sorted(date_entries.keys())[-1]:
                            station_latest_broken[station][mode][element] = True
                            print 'LATEST:{3} {0}, {1}, {2}'.format(station,mode,element,date_entry)
                x, y_diff = sorted_lists(x,y_diff)
                for i,iy_diff in enumerate(y_diff[1:]):
                    y.append(y[i]+iy_diff)
                print len(x),len(y_diff),len(y)
                x.append(datetime.strptime('2015-02-10 00:00:00','%Y-%m-%d %H:%M:%S'))
                y.append(y[-1])
                xplot, yplot = create_step_vars(x,y)
                title = ''#{0} :-: {1}'.format(station,'LBA OUTER')
                create_plot(ax,xplot,yplot,title,'orange','','# antennas',35)
                create_legend(ax,'LBA{0}'.format('OUTER'),'orange',0)
                x = [x[0]]
                y = [y[0]]
                y_diff = [0]
                del xplot,yplot
                print >>f_pvss, station_latest_broken
                print '========================================='
                        ######### END MODE ##########



######### START MODE ##########
            if mode == 'LBA-INNER' and station[2:3] != '6':
                print station
                station_latest_broken = {station: {mode:{}}}
                x.append(datetime.strptime('2010-10-01 00:00:00','%Y-%m-%d %H:%M:%S'))
                y.append(len(mapping_ant[station[0:2]][mode]))
                y_diff.append(0)
                rect = 0.06,0.29,0.9,0.19
                ax = ax = fig.add_axes(rect)
#                print 'x0:',x,'y0:',y,'y_diff0:',y_diff
                for element in mapping_ant[station[0:2]][mode]:
#                    print mapping_ant[station[0:2]][mode]
#                    print '#: ', element, ' -- total possible: ', y[0]
                    tmp_broken = False
                    date_entries = {}
                    for i in ant_dbase[station]['LBA']['{0:0>3}'.format(element)]:
                        date_entries[i] = 'ANT'
                        # Make sure the mapping from RCU to element is done correctly
                    for i in rcu_dbase[station]['{0}'.format(element*2)]:
                        date_entries[i] = 'RCU'
                    for i in rcu_dbase[station]['{0}'.format(1+element*2)]:
                        date_entries[i] = 'RCU'
                    for date_entry in sorted(date_entries):
                        if date_entries[date_entry] == 'ANT':
                            if int(ant_dbase[station]['LBA']['{0:0>3}'.format(element)][date_entry]) > 10:
                                broken = True
                            else: broken = False
                        if date_entries[date_entry] == 'RCU':
                            try:
                                if int(rcu_dbase[station]['{0}'.format(element*2)][date_entry]) > 10:
                                    broken = True
                            except(KeyError):
                                if int(rcu_dbase[station]['{0}'.format(1+element*2)][date_entry]) > 10:
                                    broken = True
                            else: broken = False
                        try:
                            x.append(datetime.strptime(date_entry,'%Y-%m-%d %H:%M:%S.%f'))
                            x_tot[total_modes[2]].append(datetime.strptime(date_entry,'%Y-%m-%d %H:%M:%S.%f'))
                            x_tot[total_modes[1]].append(datetime.strptime(date_entry,'%Y-%m-%d %H:%M:%S.%f'))
                        except(ValueError):
                            x.append(datetime.strptime(date_entry,'%Y-%m-%d %H:%M:%S'))
                            x_tot[total_modes[2]].append(datetime.strptime(date_entry,'%Y-%m-%d %H:%M:%S'))
                            x_tot[total_modes[1]].append(datetime.strptime(date_entry,'%Y-%m-%d %H:%M:%S'))
                        if broken == tmp_broken: tmp = 0
                        if broken == True & broken != tmp_broken: tmp = -1
                        if broken == False & broken != tmp_broken: tmp = 1
                        y_diff.append(tmp)
                        y_diff_tot[total_modes[2]].append(tmp)
                        y_diff_tot[total_modes[1]].append(tmp)
                        tmp_broken = broken
                        if broken == True and date_entry == sorted(date_entries.keys())[-1]:
                            station_latest_broken[station][mode][element] = True
                            print 'LATEST:{3} {0}, {1}, {2}'.format(station,mode,element,date_entry)
                x, y_diff = sorted_lists(x,y_diff)
                for i,iy_diff in enumerate(y_diff[1:]):
                    y.append(y[i]+iy_diff)
#                print len(x),len(y_diff),len(y)
                x.append(datetime.strptime('2015-02-10 00:00:00','%Y-%m-%d %H:%M:%S'))
                y.append(y[-1])
                xplot, yplot = create_step_vars(x,y)
                title = ''#{0} :-: {1}'.format(station,'LBA INNER')
                create_plot(ax,xplot,yplot,title,'green','','# antennas',35)
                create_legend(ax,'LBA{0}'.format('INNER'),'green',0)
                x = [x[0]]
                y = [y[0]]
                y_diff = [0]
                del xplot,yplot
                print >>f_pvss, station_latest_broken
                print '========================================='
                        ######### END MODE ##########


######### START MODE ##########
            if mode == 'HBA-DUAL-INNER' and station[0:2] == 'CS':
                station_latest_broken = {station: {mode:{}}}
                x.append(datetime.strptime('2010-10-01 00:00:00','%Y-%m-%d %H:%M:%S'))
                y.append(len(mapping_ant[station[0:2]][mode])/2)
                y_diff.append(0)
                for hbasubmode in ['0','1']:
                    rect = 0.06,0.77-int(hbasubmode)*0.24,0.9,0.19
                    ax = ax = fig.add_axes(rect)
                    print x,y,y_diff
                    for element in mapping_ant[station[0:2]][mode][0+int(hbasubmode)*24:24+int(hbasubmode)*24]:
                        tmp_broken = False
                        date_entries = {}
                        for i in ant_dbase[station]['HBA']['{0:0>2}'.format(element)]:
                            date_entries[i] = 'ANT'
                        # Make sure the mapping from RCU to element is done correctly
                        for i in rcu_dbase[station]['{0}'.format(element*2)]:
                            date_entries[i] = 'RCU'
                        for i in rcu_dbase[station]['{0}'.format(1+element*2)]:
                            date_entries[i] = 'RCU'
                        for date_entry in sorted(date_entries):
                            if date_entries[date_entry] == 'ANT':
                                if int(ant_dbase[station]['HBA']['{0:0>2}'.format(element)][date_entry]) > 10:
                                    broken = True
                                else: broken = False
                            if date_entries[date_entry] == 'RCU':
                                try:
                                    if int(rcu_dbase[station]['{0}'.format(element*2)][date_entry]) > 10:
                                        broken = True
                                except(KeyError):
                                    if int(rcu_dbase[station]['{0}'.format(1+element*2)][date_entry]) > 10:
                                        broken = True
                                else: broken = False
                            try:
                                x.append(datetime.strptime(date_entry,'%Y-%m-%d %H:%M:%S.%f'))
                                x_tot[total_modes[0]].append(datetime.strptime(date_entry,'%Y-%m-%d %H:%M:%S.%f'))
                            except(ValueError):
                                x.append(datetime.strptime(date_entry,'%Y-%m-%d %H:%M:%S'))
                                x_tot[total_modes[0]].append(datetime.strptime(date_entry,'%Y-%m-%d %H:%M:%S'))
                            if broken == tmp_broken: tmp = 0
                            if broken == True & broken != tmp_broken: tmp = -1
                            if broken == False & broken != tmp_broken: tmp = 1
                            y_diff.append(tmp)
                            y_diff_tot[total_modes[0]].append(tmp)
                            tmp_broken = broken
                            if broken == True and date_entry == sorted(date_entries.keys())[-1]:
                                station_latest_broken[station][mode][element] = True
                    x, y_diff = sorted_lists(x,y_diff)
                    for i,iy_diff in enumerate(y_diff[1:]):
                        y.append(y[i]+iy_diff)
                    x.append(datetime.strptime('2015-02-10 00:00:00','%Y-%m-%d %H:%M:%S'))
                    y.append(y[-1])
                    xplot, yplot = create_step_vars(x,y)
                    if hbasubmode == '0': title = '{0}'.format(station) #{0} :-: {1}'.format(station,'HBA [0/1/DUAL(-INNER)]')
                    else: title = ''
                    colors = ['blue','red']
                    create_plot(ax,xplot,yplot,title,colors[int(hbasubmode)],'','# tiles',12)
                    create_legend(ax,'HBA{0}'.format(hbasubmode),colors[int(hbasubmode)],int(hbasubmode))
                    x = [x[0]]
                    y = [y[0]]
                    y_diff = [0]
                    del xplot,yplot
                    print >>f_pvss, station_latest_broken
                    print '========================================='
######### END MODE ##########

######### START MODE ##########
            if mode == 'HBA-DUAL-INNER' and station[0:2] == 'RS':
                station_latest_broken = {station: {mode:{}}}
                x.append(datetime.strptime('2010-10-01 00:00:00','%Y-%m-%d %H:%M:%S'))
                y.append(len(mapping_ant[station[0:2]][mode]))
                y_diff.append(0)
                rect = 0.06,0.53,0.9,0.19
                ax = ax = fig.add_axes(rect)
                print x,y,y_diff
                for element in mapping_ant[station[0:2]][mode]:
                    tmp_broken = False
                    date_entries = {}
                    for i in ant_dbase[station]['HBA']['{0:0>2}'.format(element)]:
                        date_entries[i] = 'ANT'
                    # Make sure the mapping from RCU to element is done correctly
                    for i in rcu_dbase[station]['{0}'.format(element*2)]:
                        date_entries[i] = 'RCU'
                    for i in rcu_dbase[station]['{0}'.format(1+element*2)]:
                        date_entries[i] = 'RCU'
                    for date_entry in sorted(date_entries):
                        if date_entries[date_entry] == 'ANT':
                            if int(ant_dbase[station]['HBA']['{0:0>2}'.format(element)][date_entry]) > 10:
                                broken = True
                            else: broken = False
                        if date_entries[date_entry] == 'RCU':
                            try:
                                if int(rcu_dbase[station]['{0}'.format(element*2)][date_entry]) > 10:
                                    broken = True
                            except(KeyError):
                                if int(rcu_dbase[station]['{0}'.format(1+element*2)][date_entry]) > 10:
                                    broken = True
                            else: broken = False
                        try:
                            x.append(datetime.strptime(date_entry,'%Y-%m-%d %H:%M:%S.%f'))
#                            x_tot[mode].append(datetime.strptime(date_entry,'%Y-%m-%d %H:%M:%S.%f'))
                        except(ValueError):
                            x.append(datetime.strptime(date_entry,'%Y-%m-%d %H:%M:%S'))
#                            x_tot[mode].append(datetime.strptime(date_entry,'%Y-%m-%d %H:%M:%S'))
                        if broken == tmp_broken: tmp = 0
                        if broken == True & broken != tmp_broken: tmp = -1
                        if broken == False & broken != tmp_broken: tmp = 1
                        y_diff.append(tmp)
#                        y_diff_tot[mode].append(tmp)
#                           print date_entry, ant_dbase[station]['HBA']['{0:0>2}'.format(element)][date_entry], 'prev:',tmp_broken,'now:',broken, 'diff:',tmp
                        tmp_broken = broken
                        if broken == True and date_entry == sorted(date_entries.keys())[-1]:
                            station_latest_broken[station][mode][element] = True
                x, y_diff = sorted_lists(x,y_diff)
                for i,iy_diff in enumerate(y_diff[1:]):
                    y.append(y[i]+iy_diff)
                x.append(datetime.strptime('2015-02-10 00:00:00','%Y-%m-%d %H:%M:%S'))
                y.append(y[-1])
                xplot, yplot = create_step_vars(x,y)
                title = '' #{0} :-: {1}'.format(station,'HBA-INNER')
                create_plot(ax,xplot,yplot,title,'blue','','# tiles',12)
                create_legend(ax,'HBA{0}'.format('INNER'),'blue',0)
                x = [x[0]]
                y = [y[0]]
                y_diff = [0]
                del xplot,yplot
                print >>f_pvss, station_latest_broken
                print '========================================='
######### END MODE ##########

######### START MODE ##########
            if mode == 'HBA-DUAL' and station[0:2] == 'RS':
                station_latest_broken = {station: {mode:{}}}
                x.append(datetime.strptime('2010-10-01 00:00:00','%Y-%m-%d %H:%M:%S'))
                y.append(len(mapping_ant[station[0:2]][mode]))
                y_diff.append(0)
                rect = 0.06,0.77,0.9,0.19
                ax = ax = fig.add_axes(rect)
                print x,y,y_diff
                for element in mapping_ant[station[0:2]][mode]:
                    tmp_broken = False
                    date_entries = {}
                    for i in ant_dbase[station]['HBA']['{0:0>2}'.format(element)]:
                        date_entries[i] = 'ANT'
                    # Make sure the mapping from RCU to element is done correctly
                    for i in rcu_dbase[station]['{0}'.format(element*2)]:
                        date_entries[i] = 'RCU'
                    for i in rcu_dbase[station]['{0}'.format(1+element*2)]:
                        date_entries[i] = 'RCU'
                    for date_entry in sorted(date_entries):
                        if date_entries[date_entry] == 'ANT':
                            if int(ant_dbase[station]['HBA']['{0:0>2}'.format(element)][date_entry]) > 10:
                                broken = True
                            else: broken = False
                        if date_entries[date_entry] == 'RCU':
                            try:
                                if int(rcu_dbase[station]['{0}'.format(element*2)][date_entry]) > 10:
                                    broken = True
                            except(KeyError):
                                if int(rcu_dbase[station]['{0}'.format(1+element*2)][date_entry]) > 10:
                                    broken = True
                            else: broken = False
                        try:
                            x.append(datetime.strptime(date_entry,'%Y-%m-%d %H:%M:%S.%f'))
                            x_tot[total_modes[0]].append(datetime.strptime(date_entry,'%Y-%m-%d %H:%M:%S.%f'))
                        except(ValueError):
                            x.append(datetime.strptime(date_entry,'%Y-%m-%d %H:%M:%S'))
                            x_tot[total_modes[0]].append(datetime.strptime(date_entry,'%Y-%m-%d %H:%M:%S'))
                        if broken == tmp_broken: tmp = 0
                        if broken == True & broken != tmp_broken: tmp = -1
                        if broken == False & broken != tmp_broken: tmp = 1
                        y_diff.append(tmp)
                        y_diff_tot[total_modes[0]].append(tmp)
#                           print date_entry, ant_dbase[station]['HBA']['{0:0>2}'.format(element)][date_entry], 'prev:',tmp_broken,'now:',broken, 'diff:',tmp
                        tmp_broken = broken
                        if broken == True and date_entry == sorted(date_entries.keys())[-1]:
                            station_latest_broken[station][mode][element] = True
                x, y_diff = sorted_lists(x,y_diff)
                for i,iy_diff in enumerate(y_diff[1:]):
                    y.append(y[i]+iy_diff)
                x.append(datetime.strptime('2015-02-10 00:00:00','%Y-%m-%d %H:%M:%S'))
                y.append(y[-1])
                xplot, yplot = create_step_vars(x,y)
                title = '{0}'.format(station) #{0} :-: {1}'.format(station,'HBA-INNER')
                create_plot(ax,xplot,yplot,title,'red','','# tiles',35)
                create_legend(ax,'HBA{0}'.format(''),'red',0)
                x = [x[0]]
                y = [y[0]]
                y_diff = [0]
                del xplot,yplot
                print >>f_pvss, station_latest_broken
                print '========================================='
######### END MODE ##########


        fig.savefig('PNG/{0}.png'.format(station))
        print station, mode, 'saved figure'
#        fig.close()


    fig_tot = plt.figure(figsize=(12,6))
    rects = [[0.06,0.77,0.9,0.19],[0.06,0.53,0.9,0.19],[0.06,0.29,0.9,0.19],[0.06,0.05,0.9,0.19]]
    total_hardware = [2592.,4320.,1824.,1824.]
    minscales = [2292./total_hardware[0],4120./total_hardware[1],1724./total_hardware[2],1724./total_hardware[3]]
    for j,mode in enumerate(total_modes):
        rect = rects[j]
        print rect,'-1:',len([a for a in y_diff_tot[mode] if a == -1]), \
                   '+1:',len([a for a in y_diff_tot[mode] if a == 1])
        ax = fig_tot.add_axes(rect)
        x_tot[mode], y_diff_tot[mode] = sorted_lists(x_tot[mode],y_diff_tot[mode])
        for i,iy_diff_tot in enumerate(y_diff_tot[mode][1:]):
#            print x_tot[mode][i],y_diff_tot[mode][i],y_tot[mode][i]+iy_diff
            y_tot[mode].append(y_tot[mode][i]+iy_diff_tot)
        x_tot[mode].append(datetime.strptime('2015-02-10 00:00:00','%Y-%m-%d %H:%M:%S'))
        y_tot[mode].append(y_tot[mode][-1])
        xplot, yplot_tmp = create_step_vars(x_tot[mode],y_tot[mode])
        yplot = [yy/float(total_hardware[j]) for yy in yplot_tmp]
        print mode,len(xplot),len(yplot),min(yplot),max(yplot)
        if j == 0:
            title = '{0}'.format('TOTALS') #{0} :-: {1}'.format(station,'HBA-INNER')
        else:
            title = ''
        create_totals_plot(ax,xplot,yplot,title,'b','','# tiles',0.88)#minscales[j])
        create_legend(ax,'{0}'.format(mode),'b',0)
#        ax.xaxis.set_major_locator(MaxNLocator(4))
        ax.yaxis.set_major_locator(MaxNLocator(4))
        del xplot,yplot
    fig_tot.savefig('PNG/{0}.png'.format('000'))

    return exit_status

#### END MAIN BODY ####




########################
#       M A I N        #
########################

if __name__ == '__main__':
    sys.exit(my_main(sys.argv))
