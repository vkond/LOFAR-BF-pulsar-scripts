#!/usr/bin/perl 
#
# script calculates the azimuth and hour angle (and LST)
# (actually the min and max values of HA and LST)
# of the source with pointed RA and DEC for pointed elevation (zenith angle)
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


$is_never_set = 0;    # flag, if 1 then source never sets
$is_never_rise = 0;   # flag, if 1 then source never rises
$is_never_above = 0;  # flag, if 1 then source is never above this EL (no matter circumpolar or not)
$ra = "";
$dec = "";
$EL = 30.;

# if there are no parameters in command line
# help is called and program normally exit
&help($0) && exit if $#ARGV < 0;

# Parse command line
GetOptions ( "el=f" => \$EL,   # - elevation (in deg)
             "ra=s" => \$ra,   # - RA of the source "hh:mm:ss.sss"
             "dec=s" => \$dec, # - DEC of the source "[+|-]dd:mm:ss.ssss"
             "lat=f" => \$phi, # - latitude (degrees)
             "lon=f" => \$lambda, # - longitude (degrees)
             "h" => \$help);   # - to print the help

if ($help) { &help($0); exit 0; }

if ($ra eq "" || $dec eq "") {
 &error ("No RA or DEC commanded!");
} 

if ($phi eq "") {
 &error ("Latitude is not defined!");
} 

if ($lambda eq "") {
 &error ("Longitude is not defined!");
} 

if ($EL eq "") {
 &error ("Elevation is not defined!");
}

$dec1 = &dec2rad ($dec);


# zenith angle
$ZA = 90. - $EL;
$ZA = sprintf ("%.1f", $ZA);
$EL = sprintf ("%.1f", $EL);

print "\n";
print "Source: RA = $ra  DEC = $dec\n";
print "Site: LAT = $phi deg  LON = $lambda deg\n";
print "\n";
print "EL = $EL deg (ZA = $ZA deg)\n";

# checking if the source is circumpolar or not (never set or rise) at this latitude and elevation
# I'm using 0.001 as a tolerance to compare correctly two float values
if (($phi >= 0. && $dec1/$rad >= 90. - $phi + $EL - 0.001) || ($phi < 0. && $dec1/$rad <= -(90. + $phi + $EL - 0.001))) {  # never sets
  $is_never_set = 1;
  # calculating the minimum EL of the source
  if ($phi >= 0.) { $mEL = sprintf ("%.1f", $dec1/$rad + $phi - 90.); }
   else { $mEL = sprintf ("%.1f", -1. * ($dec1/$rad + $phi) - 90.); }
}
if (($phi >= 0. && $dec1/$rad < -(90. - $phi) + $EL + 0.001) || ($phi < 0. && $dec1/$rad > 90. + $phi - $EL - 0.001)) { # never rise
  $is_never_rise = 1;
  # calculating the maximum EL of the source
  if ($phi >= 0.) { $mEL = sprintf ("%.1f", $dec1/$rad - $phi + 90.); }
   else { $mEL = sprintf ("%.1f", $phi - $dec1/$rad + 90.); }
}

# getting the EL at transit to compare with given EL 
$EL_transit = sprintf ("%.1f", 90. - &get_ZA ($dec1, $phi, 0.0));
if ($EL >= $EL_transit - 0.001) {
 $is_never_above = 1;
}

# Calculating the HA, AZ range and LST set and rise for sources
# that are not circumpolar
if ($is_never_set == 0 && $is_never_rise == 0 && $is_never_above == 0) {
 # hour angle (absolute value) in rad
 # actually there are 2 values +- $HA
 $HA = &get_HA ($ZA, $dec1, $phi);
 $HAmin = -$HA;
 $HAmax = $HA;
 # in degrees (from South clockwise)
 $AZmin = &get_AZ ($dec1, $phi, $HAmin);
 $AZmax = &get_AZ ($dec1, $phi, $HAmax);


 $alphah = &time2hour($ra);
 $sidmin = `echo \"scale=20\n$alphah - ($HA / $rad / 15.)\" | bc -l`;
 $sidmax = `echo \"scale=20\n$alphah + ($HA / $rad / 15.)\" | bc -l`;
 if ($sidmin < 0.) { $sidmin += 24.; }
 if ($sidmax >= 24.) { $sidmax -= 24.; }
 $LSTmin = &hour2time ($sidmin);
 $LSTmax = &hour2time ($sidmax);

 $is_always = ($phi >= 0. ? "   [always North]" : "   [always South]");
 printf ("AZ = [ %.1f ; %.1f ] deg%s\n", $AZmin, $AZmax, abs($dec1/$rad)>=abs($phi) ? $is_always : "");
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
 printf ("EL at transit = %.1f deg\n", $EL_transit);
 printf ("HA =  +/- %.2f %s  [Duration = %.2f %s]\n", $hapres, $hapres_label, 2 * $hapres, $hapres_label);
 printf ("LST = [rise: %s] [set: %s]\n", $LSTmin, $LSTmax);
 print "\n";
}

# if source never sets
if ($is_never_set == 1) {
 if (abs($dec1/$rad)>=abs($phi)) {
  printf ("AZ - always %s\n", $phi >= 0. ? "North" : "South");
 }
 print "Circumpolar source (never sets) at this LAT and EL\n";
 print "The minimum EL = $mEL deg\n";
}
# if source never rises
if ($is_never_rise == 1) {
 print "Circumpolar source (never rises) at this LAT and EL\n";
 print "The maximum EL = $mEL deg\n";
}
# if source is never above given EL
if ($is_never_above == 1 && $is_never_rise == 0) {
 print "The source is never above the given EL\n";
 print "Tha maximum EL = $EL_transit deg\n";
}



# help
sub help {
 local ($prg) = @_;
 $prg = `basename $prg`;
 chomp $prg;
 print "$prg - calculates AZ and LST of rise and set for a given source at given observatory\n";
 print "Usage: $prg [options]\n";
 print "        -el  EL     - elevation in degrees, default = 30 deg\n";
 print "        -ra  RA     - right assention of the source where RA is in \"hh:mm:ss.sss\"\n";
 print "        -dec DEC    - declination of the source where DEC is in \"[+\-]dd:mm:ss.sss\"\n";
 print "        -lat PHI    - latitude in degrees (default - LOFAR)\n";
 print "        -lon LAMBDA - longitude in degrees (default - LOFAR)\n";
 print "                      Western longitudes are negative!\n";
 print "        -h          - print this help\n";
}

# unique subroutine to output errors
# one input parameter (string): the error message
sub error {
 my ($error) = @_;
 print "Error: $error\n";
 exit 1;
}

# transform dec format into radians
sub dec2rad {
 my ($d) = @_;
 $deg = $d; $deg =~ s/^([+-]?\d\d)\:.*$/$1/;
 $deg *= 1;
 $mm = $d; $mm =~ s/^[+-]?\d\d\:(\d\d)\:.*$/$1/;
 $ss = $d; $ss =~ s/^[+-]?\d\d\:\d\d\:(\d\d.*)$/$1/;
 if ($deg < 0) {
  $newd = `echo \"($deg - $mm / 60. - $ss / 3600.) * $rad\" | bc -l`;
 } else {
  $newd = `echo \"($deg + $mm / 60. + $ss / 3600.) * $rad\" | bc -l`;
 }
 chomp $newd;
 return $newd;
}

# get hour angle (absolute value) in radians
# actually there are 2 values +- 
sub get_HA {
 my ($zen, $delta, $shirota) = @_;

 $ha = acos ( (cos ($zen * $rad) - sin ($delta) * sin ($shirota * $rad)) / (cos ($delta) * cos ($shirota * $rad)) );
 return $ha;
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

# transform time format into hours
sub time2hour {
 my ($time) = @_;
 $hh = $time; $hh =~ s/^(\d\d)\:.*$/$1/;
 $mm = $time; $mm =~ s/^\d\d\:(\d\d)\:.*$/$1/;
 $ss = $time; $ss =~ s/^\d\d\:\d\d\:(\d\d.*)$/$1/;
 $t = `echo \"scale=20\n($hh + $mm / 60. + $ss / 3600.)\" | bc -l`;
 chomp $t;
 return $t;
}

# transform hours into time format
sub hour2time {
 my ($hours) = @_;
 $hh = int ($hours);
 $mm = int ( 60. * ($hours - $hh) );
 $ss = (($hours - $hh) * 60. - $mm) * 60.;
 $time = sprintf ("%02d:%02d:%02d", $hh, $mm, $ss);
 return $time;
}
