#!/usr/bin/perl -w
#
# Calculate the local sidereal time

# Vlad Kondratiev (c)
#
use Getopt::Long;

# Longitude
#$lambda = -79.8398333333; # GBT  79 deg 50' 23.40" W
#$lambda = 6.601882; # WSRT
$lambda = 6.869883; # LOFAR

$UTC = "";
# if 1 then given input time is LST rather than UTC
$is_lst = 0;

# if there are no parameters in command line
# help is called and program normally exit
#&help($0) && exit if $#ARGV < 0;

# Parse command line
GetOptions ( "t=s" => \$UTC,   # - time
             "lon=f" => \$lambda, # - longitude (degrees)
             "lst" => \$is_lst,   # flag to switch input time to be LST (works only if -t option is used)
             "h" => \$help);   # - to print the help

if ($help) { &help($0); exit 0; }

if ($lambda eq "") {
 &error ("Longitude is not defined!");
} 

if ($UTC eq "") {
 $UTC = &get_curtime();
 $is_lst = 0;
}

# Given time is UTC and should be converted to LST
if ($is_lst == 0) {

# also checking if UTC is given in another format, namely YYYY-MM-DDTHH:MM:SS[.SSS]
# if it's in different format, then redo it to usual format
if ($UTC =~ /^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}\.?\d*$/) {
 $UTC = &time_reformat($UTC);
}

# running "jd.pl" first to get JD for a given UTC time
$jd = `jd.pl \"$UTC\" | grep JD | grep -v MJD | awk '{print \$2}' -`;
chomp $jd;

# get the sidereal time
$sid = &sidereal($jd, $lambda);
$sid += 0.01/3600.;  # adding fake 1/100 sec to avoid, e.g. 23:59:59.99
if ($sid >= 24.) { $sid -= 24; }

# converting LST to readable format
$LST = &hours2str($sid);
} # $is_lst == 0



# Given time is LST and should be converted to UTC
if ($is_lst != 0) {

# from the command line option the input time was assigned to $UTC, re-organzing it here
# leave only date for UTC and assign time to LST depending on given format
$utc = $UTC;
if ($UTC =~ /^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}\.?\d*$/) {
 $utc =~ s/^(\d{4}-\d{2}-\d{2}T)\d{2}:\d{2}:\d{2}\.?\d*$/$1/;   # format to use in the output
 $UTC = &time_reformat($UTC);
} else { 
 $utc =~ s/^(\d{1,2}\.\d{1,2}\.[-+]?\d+)\s+\d{2}:\d{2}:\d{2}\.?\d*$/$1/;
}
 $LST = $UTC; $LST =~ s/^\d{1,2}\.\d{1,2}\.[-+]?\d+\s+(\d{2}:\d{2}:\d{2}\.?\d*)$/$1/;
 $UTC =~ s/^(\d{1,2}\.\d{1,2}\.[-+]?\d+)\s+\d{2}:\d{2}:\d{2}\.?\d*$/$1/;

# coefficient to convert interval of sidereal time to interval of solar time
$sol_in_sid = 0.99726956632908432554;

# converting LST to hours
$lst = &time2hours($LST);

# running "jd.pl" to get reference JD for a given date 
$ref_jd = `jd.pl \"$UTC 12:00:00\" | grep JD | grep -v MJD | awk '{print \$2}' -`;
chomp $ref_jd;
# get the sidereal time for reference JD
$ref_gmst = &sidereal($ref_jd, $lambda);

$time = ($lst - $ref_gmst);
if ($time >= 24.) { $time -= 24; }
if ($time < 0.) { $time += 24; }
$time *= $sol_in_sid;
$time += 12.;
# checking if $time (UTC) is more than duration of sidereal day. In this case
# we subtract it to give the smaller UTC time for the same LST
if ($time >= 24.*$sol_in_sid) { $time -= 24.*$sol_in_sid; }

# converting UTC to readable format
$time = &hours2str($time);
$UTC = "$utc$time";
} # 


# output
printf ("UTC: %s\n", &time2str($UTC));
printf ("LST: %s\n", &time2str($LST));





# help
sub help {
 local ($prg) = @_;
 $prg = `basename $prg`;
 chomp $prg;
 print "$prg: calculates LST for a given UTC time and given longitude\n";
 print "Usage: $prg [options]\n";
 print "        -t   TIME   - UTC time in format \"DD.MM.YYYY hh:mm:ss.sss\" or \"YYYY-MM-DDThh:mm:ss.sss\"\n";
 print "                      if no time is pointed than current UTC time will be used\n";
 print "        -lon LAMBDA - longitude in degrees (default - LOFAR)\n";
 print "                      Western longitudes are negative!\n";
 print "        -lst        - input time is LST and it is to be converted to UTC; only works if -t option is used\n";
 print "                      Note: if there are two UTC epochs within the same date with the same LST,\n";
 print "                      the earliest one will be given\n";
 print "        -h          - print this help\n";
}


# calculates the sidereal time given the JD and longitude
# from  http://www.usno.navy.mil/USNO/astronomical-applications/astronomical-information-center/approx-sider-time
sub sidereal {
 my ($jd, $lambda) = @_;
 $part = $jd - int ($jd);
 if ($part >= 0.5) { $jd0 = int($jd) + 0.5; $H = ($part - 0.5) * 24.; }
  else { $jd0 = int($jd) - 0.5; $H = ($part + 0.5) * 24.; }

 $D  = $jd  - 2451545.0;
 $D0 = $jd0 - 2451545.0; 
 $T = $D / 36525; 

 $gmst = 6.697374558 + 0.06570982441908 * $D0 + 1.00273790935 * $H + 0.000026 * $T * $T + $lambda/15.;
 $gmst = $gmst - (int($gmst/24.)) * 24.;

 return $gmst;
}

# get current UTC time
sub get_curtime {
 $sec = `date +%S`; chomp $sec;
 $min = `date +%M`; chomp $min;
 $hour = `date +%H`; chomp $hour;
 $day = `date +%d`; chomp $day;
 $month = `date +%m`; chomp $month;
 $year = `date +%Y`; chomp $year;
 $timezone = `date +%z`; chomp $timezone;
 $minadd = $timezone; $minadd =~ s/^[+-]?\d\d(\d\d)$/$1/;
 $houradd = $timezone; $houradd =~ s/^[+-]?(\d\d)\d\d$/$1/;
 if ($timezone =~ /^-/) {
  $min += $minadd;
  if ($min >= 60) { $min -= 60; $hour += 1; }
  $hour += $houradd;
  if ($hour >= 24) { $hour -= 24; $day += 1; }
 } else {
  $min -= $minadd;
  if ($min < 0) { $min += 60; $hour -= 1; }
  $hour -= $houradd;
  if ($hour < 0) { $hour += 24; $day -= 1; }
 }
 $min = sprintf ("%02d", $min);
 $hour = sprintf ("%02d", $hour);
 $day = sprintf ("%02d", $day);
 return "$day.$month.$year $hour:$min:$sec";
}

# converts the date/time from format YYYY-MM-DDTHH:MM:SS.SSS to "DD.MM.YYYY HH:MM:SS.SSS"
sub time_reformat {
 my ($tgiven) = @_;
 $year = $tgiven; $year =~ s/^(\d{4})-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}\.?\d*$/$1/;
 $month = $tgiven; $month =~ s/^\d{4}-(\d{2})-\d{2}T\d{2}:\d{2}:\d{2}\.?\d*$/$1/;
 $day = $tgiven; $day =~ s/^\d{4}-\d{2}-(\d{2})T\d{2}:\d{2}:\d{2}\.?\d*$/$1/;
 $time = $tgiven; $time =~ s/^\d{4}-\d{2}-\d{2}T(\d{2}:\d{2}:\d{2}\.?\d*)$/$1/;
 return "$day.$month.$year $time";
}

# gets time format string and leaves only 2 decimal digits of seconds
sub time2str {
 my ($time) = @_;
 $tmp = $time;
 $tmp =~ s/^(.*\d\d\:\d\d\:\d\d\.?\d?\d?).*$/$1/;
 return $tmp;
}

# to convert time string to hours
sub time2hours {
 my ($timestr) = @_;
 $h = $timestr; $h =~ s/^(\d{2}):\d{2}:\d{2}\.?\d*$/$1/;
 $m = $timestr; $m =~ s/^\d{2}:(\d{2}):\d{2}\.?\d*$/$1/;
 $s = $timestr; $s =~ s/^\d{2}:\d{2}:(\d{2}\.?\d*)$/$1/;
 return $h + ($m + $s/60.)/60.;
}

# converting time in hours to readable format
sub hours2str {
 my ($time) = @_;
 $h = sprintf ("%02d", int ($time));
 $m = sprintf ("%02d", int (($time - $h) * 60.));
 $s = ($time - $h - $m/60.) * 3600.;
 if ($s < 10.) { $s = "0" . $s; }
 return "$h:$m:$s";
}
