#!/usr/bin/perl 
#
# script calculates the azimuth and zenith angle (elevation)
# of the source with pointed RA and DEC for pointed time
# the script also has functions to calculate siderial time
# (based on my sid.pl) for specific time or at "now" time
#
# Vlad Kondratiev (c)
#
use Math::Trig;
use Getopt::Long;

$pi = 3.14159265359;
$rad = $pi/180.;

# Latitude
#$phi = 38.4331194444; # GBT  38 deg 25' 59.23" N
#$phi = 52.915604; # WSRT
$phi = 52.915111; # LOFAR
# Longitude
#$lambda = -79.8398333333; # GBT  79 deg 50' 23.40" W
#$lambda = 6.601882; # WSRT
$lambda = 6.869883; # LOFAR

$ra = "";
$dec = "";
$UTC = "";

# if there are no parameters in command line
# help is called and program normally exit
&help($0) && exit if $#ARGV < 0;

# Parse command line
GetOptions ( "t=s" => \$UTC,   # - time
             "ra=s" => \$ra,   # - RA of the source "hh:mm:ss.sss"
             "dec=s" => \$dec, # - DEC of the source "[+|-]dd:mm:ss.ssss"
             "lat=f" => \$phi, # - latitude (degrees)
             "lon=f" => \$lambda, # - longitude (degrees)
             "h" => \$help);   # - to print the help

if ($help) { &help($0); exit 0; }

if ($UTC eq "") {
 $UTC = &get_curtime();
}

# also checking if UTC is given in another format, namely YYYY-MM-DDTHH:MM:SS[.SSS]
# if it's in different format, then redo it to usual format
if ($UTC =~ /^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}\.?\d*$/) {
 $UTC = &time_reformat($UTC);
}

if ($ra eq "" || $dec eq "") {
 &error ("No RA or DEC commanded!");
} 

if ($phi eq "") {
 &error ("Latitude is not defined!");
} 

if ($lambda eq "") {
 &error ("Longitude is not defined!");
} 

$LST = `sid.pl -t \"$UTC\" -lon $lambda | grep LST | awk '{print \$2}' -`;
chomp $LST;
$HA = &get_HA ($LST, $ra);
$dec1 = &dec2rad ($dec);

# in degrees
$ZA = &get_ZA ($dec1, $phi, $HA);
$ZA = sprintf ("%.1f", $ZA);
# in degrees
$EL = 90. - $ZA;
$EL = sprintf ("%.1f", $EL);
# in degrees (from South clockwise)
$AZ = &get_AZ ($dec1, $phi, $HA);
$AZ = sprintf ("%.1f", $AZ);

print "\n";
print "Source: RA = $ra  DEC = $dec\n";
print "Site: LAT = $phi deg  LON = $lambda deg\n";
print "\n";
printf ("UTC: %s\n", &time2str($UTC));
printf ("LST: %s\n", &time2str($LST));
$hapres = $HA / $rad / 15.; 
$hapres_label = "h";
if (abs($hapres) < 1) { 
 $hapres *= 60.; 
 $hapres_label = "min"; 
 if (abs($hapres) < 1) {
  $hapres *= 60.; 
  $hapres_label = "sec"; 
 }
}
printf ("HA = %.2f %s\n", $hapres, $hapres_label);
print "EL = $EL deg (ZA = $ZA deg)\n";
print "AZ = $AZ deg\n";
print "\n";



# help
sub help {
 local ($prg) = @_;
 $prg = `basename $prg`;
 chomp $prg;
 print "$prg: calculates AZ, EL(ZA), HA  for a given source and time at given observatory\n";
 print "Usage: $prg [options]\n";
 print "        -t   TIME   - UTC time in format \"DD.MM.YYYY hh:mm:ss.sss\" or \"YYYY-MM-DDThh:mm:ss.sss\"\n";
 print "                      if no time is pointed than current UTC time will be used\n";
 print "        -ra  RA     - right assention of the source where RA is in \"hh:mm:ss.sss\"\n";
 print "        -dec DEC    - declination of the source where DEC is in \"[+\-]dd:mm:ss.sss\"\n";
 print "        -lat PHI    - latitude in degrees (default - LOFAR)\n";
 print "        -lon LAMBDA - longitude in degrees (default - LOFAR)\n";
 print "                      Western longitudes are negative!\n";
 print "        -h          - print this help\n";
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

# unique subroutine to output errors
# one input parameter (string): the error message
sub error {
 my ($error) = @_;
 print "Error: $error\n";
 exit 1;
}

# transform time format into radians
sub time2rad {
 my ($time) = @_;
 $hh = $time; $hh =~ s/^(\d\d)\:.*$/$1/;
 $mm = $time; $mm =~ s/^\d\d\:(\d\d)\:.*$/$1/;
 $ss = $time; $ss =~ s/^\d\d\:\d\d\:(\d\d.*)$/$1/;
 $t = `echo \"scale=20\n($hh + $mm / 60. + $ss / 3600.) * 15. * $rad\" | bc -l`;
 chomp $t;
 return $t;
}

# transform dec format into radians
sub dec2rad {
 my ($d) = @_;
 $deg = $d; $deg =~ s/^([+-]?\d\d)\:.*$/$1/;
 $deg *= 1;
 $mm = $d; $mm =~ s/^[+-]?\d\d\:(\d\d)\:.*$/$1/;
 $ss = $d; $ss =~ s/^[+-]?\d\d\:\d\d\:(\d\d.*)$/$1/;
 if ($deg < 0) {
  $newd = `echo \"scale=20\n($deg - $mm / 60. - $ss / 3600.) * $rad\" | bc -l`;
 } else {
  $newd = `echo \"scale=20\n($deg + $mm / 60. + $ss / 3600.) * $rad\" | bc -l`;
 }
 chomp $newd;
 return $newd;
}

# get hour angle (in rad)
sub get_HA {
 my ($s, $alpha) = @_;
 $srad = &time2rad ($s);
 $alpharad = &time2rad ($alpha);
 $HA = `echo \"scale=20\n$srad - $alpharad\" | bc -l`;
 chomp $HA;
 if ($HA / $rad / 15. < -12.) { $HA += 2.*$pi; }
 if ($HA / $rad / 15. >= 12.) { $HA -= 2.*$pi; }
 return $HA;
}

# get zenith angle (in degrees)
sub get_ZA { 
 my ($delta, $shirota, $hour_angle) = @_;
 
 $zen = acos ( sin ($delta) * sin ($shirota * $rad) + cos ($delta) * cos ($shirota * $rad) * cos ($hour_angle) ); 
 $zen /= $rad; 
 return $zen;
}

# get azimuth (in degrees, from North clockwise)
sub get_AZ {
 my ($delta, $shirota, $hour_angle) = @_;

 $azim = atan2 ( cos ($delta) * sin ($hour_angle), -sin ($delta) * cos ($shirota * $rad) + cos ($delta) * sin ($shirota * $rad) * cos ($hour_angle) );
 $azim /= $rad;
 $azim += 180.;
 while ($azim >= 360.) { $azim -= 360.; }
 return $azim; 
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
 $tmp =~ s/^(\d\d\:\d\d\:\d\d\.?\d?\d?).*$/$1/;
 return $tmp;
}
