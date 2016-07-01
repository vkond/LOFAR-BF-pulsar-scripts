#!/usr/bin/perl
#
# Calculate JD, MJD given usual date, and vice versa
#
# Vlad Kondratiev (c)
#
use Math::BigFloat; 

@mon = (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31);

$ref_jd = -0.5;         # 01.01.-4712 00:00:00 (i.e. 4713 BC)
$ref_year = -4712; 
$N = 0;
$dg = 0;
@grig = (1582, 10, 4, 15);

if ($#ARGV == -1) { &help($0); exit; }

#
# If date is set
#
if ($ARGV[0] =~ /^\d{1,2}\.\d{1,2}\.[-+]?\d+\s+\d{1,2}:\d{1,2}:\d+\.?\d*$/) {

$year = $month = $day = $hour = $min = $sec = $ARGV[0];

$day   =~ s/^(\d{1,2})\..*/$1/;
$month =~ s/^\d{1,2}\.(\d{1,2})\..*/$1/;
$year  =~ s/^\d{1,2}\.\d{1,2}\.([-+]?\d+)\s+.*/$1/;
$hour  =~ s/.*\s+(\d{1,2}):.*/$1/;
$min   =~ s/.*\s+\d{1,2}:(\d{1,2}):.*/$1/;
$sec   =~ s/.*\s+\d{1,2}:\d{1,2}:(\d+\.?\d*)$/$1/;

if ($year == 0) { &error("year", $year, ""); }
if ($year < 0) { $year += 1; }
 $stol = int ($year/100);
 $god  = $year % 100;
  if ($god % 4 == 0) { 
    if ($year >= $grig[0] && $stol % 4 != 0 && $god == 0) 
     { $mon[1] = 28; } else { $mon[1] = 29; } 
   } else { $mon[1] = 28; } 
if ($month < 1 || $month > 12) { &error("month", $month, ""); }
if ($day < 1 || $day > $mon[$month-1]) { &error("day", $day, ""); }
if ($hour > 23) { &error("hour", $hour, ""); }
if ($min > 59) { &error("minute", $min, ""); }
if ($sec >= 60.) { &error("second", $sec, ""); }
if ($year == $grig[0] && $month == $grig[1] && ($day > $grig[2] && $day < $grig[3])) {
  &error("day", $day, "\nBecause of Grigorian reform you write bad date!"); }


if ($year >= $ref_year) {

 for ($i=$ref_year;$i<$year;$i++) {
   $stol = int ($i/100);
   $god  = $i % 100;
    if ($i == $grig[0]) { $N += 355; next; }
    if ($god % 4 == 0) { if ($i >= $grig[0] && $stol % 4 != 0 && $god == 0) 
                          { $N += 365; } else { $N += 366; } 
                        } else { $N += 365; } 
 }
} else { for ($i=$year;$i<$ref_year;$i++) {
           $stol = int ($i/100);
           $god  = $i % 100;
           if ($god % 4 == 0) {  $N += 366; }  else { $N += 365; } 
          }
         $N *= -1;
        }

 for ($i=1;$i<$month;$i++) { $dg += $mon[$i-1]; } 
 $dg += $day - 1;
 if ($year == $grig[0]) {
   if ($month > $grig[1] || ($month == $grig[1] && $day > $grig[2])) { $dg -= 10; } }

 $remainder = Math::BigFloat->new((($sec/60 + $min)/60 + $hour)/24); 

 $dg  += $remainder;
 $jd  = $ref_jd + $dg + $N;
 $mjd = $jd - 2400000.5;

 $dg += 1;

} 
#
# If JD or MJD is set
#
elsif ($ARGV[0] =~ /^[JjMm][-+]?\d+\.?\d*$/) {

  $jd = $ARGV[0];
  $jd =~ s/^[JjMm]([-+]?\d+\.?\d*)$/$1/;

  $jd  = Math::BigFloat->new($jd);
  $mjd = Math::BigFloat->new($jd);
  if ($ARGV[0] =~ /^[Mm].*/) { $jd += 2400000.5; }
   else { $mjd -= 2400000.5; }
  $diff = int ($jd - $ref_jd);
  $ost = abs ($jd - $ref_jd - $diff);

  if ($diff !~ /^-.*/) { $year = $ref_year;

   for (;$diff>=0;$diff-=$N,$year++) {
    $stol = int ($year/100);
    $god  = $year % 100;
     if ($year == $grig[0]) { $N = 355; next; }
     if ($god % 4 == 0) { if ($year >= $grig[0] && $stol % 4 != 0 && $god == 0) 
                           { $N = 365; } else { $N = 366; } 
                         } else { $N = 365; } 
   }
    $diff += $N; $year--; $month = 1;
      if ($god % 4 == 0) { if ($year >= $grig[0] && $stol % 4 != 0 && $god == 0) 
                             { $mon[1] = 28; } else { $mon[1] = 29; } 
                            } else { $mon[1] = 28; } 
   for (;$diff>=0;$diff-=$N,$month++) { $N = $mon[$month-1]; }
    $diff += $N; $month--;
    $day = $diff + 1;
     if ($year == $grig[0] && $month >= $grig[1]) { 
       if ($month > $grig[1] || ($month == $grig[1] && $day > $grig[2])) {
          $day += 10;
          if ($day > $mon[$month-1]) { $day -= $mon[$month-1]; $month++; }
        } 
      }
} else { $year = $ref_year-1;

    for (;$diff<0;$diff+=$N,$year--) {
    $stol = int ($year/100);
    $god  = $year % 100;
     if ($god % 4 == 0) { $N = 366; } else { $N = 365; }
   }
    $diff -= $N; $year++; $month = 12;
     if ($god % 4 == 0) {  $mon[1] = 29; } else { $mon[1] = 28; }
   for (;$diff<0;$diff+=$N,$month--) { $N = $mon[$month-1]; $dg -= $N; }
    $diff -= $N; $month++;
    $day = int ($mon[$month-1] + $diff - $ost + 1);
    if ($day < 1) { $month--; if ($month < 1) { $year--; $month = 12; } 
                    $day += $mon[$month-1]; 
                   }
    $ost = 1 - $ost;
}

    $dg = 0;
    $stol = int ($year/100);
    $god  = $year % 100;
    if ($god % 4 == 0) { $mon[1] = 29; } else { $mon[1] = 28; } 
    for($i=1;$i<$month;$i++) { $dg += $mon[$i-1]; } 
    $dg += $day + $ost;

    $hour = int($ost*24);
    $min  = int (($ost*24-$hour)*60);
    $sec  = (($ost*24-$hour)*60-$min)*60;
}
#
# If day in the year is set
#
 elsif ($ARGV[0] =~ /^[Dd]\d{1,3}\.?\d*\s+[-+]?\d+$/) {

    $dg = $ost = $year = $ARGV[0];
    $dg   =~ s/^[Dd](\d{1,3}).*/$1/;
    $ost  =~ s/^[Dd]\d{1,3}(\.?\d*).*/$1/;
    $year =~ s/^[Dd]\d{1,3}\.?\d*\s+([-+]?\d+)$/$1/;

    if ($ost =~ /^\.?\d+/) { $ost = Math::BigFloat->new($ost); }
     else { $ost = 0; }

    if ($dg <= 0) { &error("day of year", $dg, ""); }

    if ($year == 0) { &error("year", $year, ""); }
    if ($year < 0) { $year += 1; }
     $stol = int ($year/100);
     $god  = $year % 100;
      if ($god % 4 == 0) {
        if ($year >= $grig[0] && $stol % 4 != 0 && $god == 0)
         { $mon[1] = 28; } else { $mon[1] = 29; }
       } else { $mon[1] = 28; }
    if ($year == $grig[0] && $dg > 355) {
     &error("day of year", $dg, "\nBecause of Grigorian reform you write bad day of year!"); }
    if ($dg > 366) { &error("day of year", $dg, ""); }
    if ($mon[1] == 28 && $dg > 365) { &error("day of year", $dg, "\nYear $year is not leap year!"); } 
 
    $hour = int ($ost*24);
    $min  = int (($ost*24-$hour)*60);
    $sec  = (($ost*24-$hour)*60-$min)*60;


    $diff = $dg;
    for ($i=0;$i<12;$i++) {
     if ($diff > $mon[$i]) { $diff -= $mon[$i]; } 
      else { $month = $i+1; $day = $diff; last; }
    }

    if ($year >= $ref_year) {

     for ($i=$ref_year;$i<$year;$i++) {
       $stol = int ($i/100);
       $god  = $i % 100;
        if ($i == $grig[0]) { $N += 355; next; }
        if ($god % 4 == 0) { if ($i >= $grig[0] && $stol % 4 != 0 && $god == 0)
                              { $N += 365; } else { $N += 366; }
                            } else { $N += 365; }
     }
    } else { for ($i=$year;$i<$ref_year;$i++) {
               $stol = int ($i/100);
               $god  = $i % 100;
               if ($god % 4 == 0) {  $N += 366; }  else { $N += 365; }
              }
             $N *= -1;
            }

     $N -= 1;
     $jd = $ref_jd + $dg + $N + $ost;
     $mjd = $jd - 2400000.5;
     $dg += $ost;

} 
#
# argument is wrong
#
else { if ($ARGV[0] =~ /^[Dd]\d{1,3}\.?\d*\s*([-+]?\s*)$/) {
        &error("year", $1, "");
       } 
       if ($ARGV[0] =~ /^[Dd](.*)\s+[-+]?\d+$/) {
        &error("day of year", $1, "");
       }
       if ($ARGV[0] =~ /^[Dd](.*)\s*\S*$/) {
        &error("day of year", $1, "");
       }
       if ($ARGV[0] =~ /^[Dd](.*)$/) { 
        &error("day of year", $1, "");
       }
       if ($ARGV[0] =~ /^[Jj](.*)$/) {
        &error("julian day", $1, ""); 
       }
       if ($ARGV[0] =~ /^[Mm](.*)$/) {
        &error("modified julian day", $1, ""); 
       }
        &error("date", $ARGV[0], ""); 
      }



# preparing for the output

if ($year <= 0) { $year--; }
if ($day !~ /^\d{2}$/) { $day = "0$day"; }
if ($month !~ /^\d{2}$/) { $month = "0$month"; }
if ($hour !~ /^\d{2}$/) { $hour = "0$hour"; }
if ($min !~ /^\d{2}$/) { $min  = "0$min"; }
if ($sec =~ /^\d{0,2}\.\d*/) { 
  if ($sec =~ /^\d{0,2}\.\d+/) { $s = $c = $sec;
                                 $s =~ s/^(\d{0,2})\.\d+/$1/; 
                                 $c =~ s/^\d{0,2}\.(\d+)/$1/;
                                 if ($s !~ /^\d{2}$/) { $s  = "0$s"; }
                                 if ($s !~ /^\d{2}$/) { $s  = "0$s"; }
                                 $sec = $s . "." . $c;
                                }
  if ($sec =~ /^\d{1,2}\.$/) { $sec =~ s/^(\d{1,2})\.$/$1/;
                               if ($sec !~ /^\d{2}$/) { $sec  = "0$s"; }
                              }
 }
if ($sec =~ /^\d{1}$/) { $sec  = "0$sec"; }



# print output

print "[DATE]: $day.$month.$year $hour:$min:$sec\n";
print "[JD]:   $jd\n";
print "[MJD]:  $mjd\n";
print "Day:    $dg\n";









# error message
sub error {
 my ($str, $par, $msg) = @_;
 print ("Error: Bad value of $str: $par !$msg\n");
 exit 1;
}


# help
sub help {
 my ($prg) = @_;
 $prg = `basename $prg`;
 chomp $prg;
 print "\nThis program calculates Julian (modified) dates\n";
 print "given usual date, or vice versa\n\n";
 print "Simple usage: $prg \"dd.mm.yyyy hh:mm:ss[.ddd]\"\n";
 print "              $prg Jddddddd[.ddd]\n";
 print "              $prg Mddddd[.ddd]\n";
 print "              $prg Dddd[.ddd] yyyy\n\n";
 print "Expert: $prg \"dd.mm.[-|+]year hh:mm:seconds[.[part_of_seconds]]\"\n";
 print "        $prg j|J[-|+]julian_day[.[part_of_julian_day]]\n";
 print "        $prg m|M[-|+]modified_julian_day[.[part_of_modified_julian_day]]\n";
 print "        $prg \"d|Dday_from_year_begin[.[part_of_day]] [-|+]year\"\n";
}
