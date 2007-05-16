#!/usr/bin/perl

# $Id$

# The suffixes of files to be measured.
@code_sufx = (
              "hpp",
              "h",
              "cpp",
              "c"
              );

# The directories, relative to the root directory, of the files to be measured.
@code_dirs = (
              "steps",
              "steps/csteps", 
              "steps/events", 
              "steps/math", 
              "steps/model", 
              "steps/parser",
              "steps/rng",
              "steps/sim", 
              "steps/stepsml",
              "steps/system", 
              "stepsim"
              );

# First, we parse the optional root path option.
die "$0 can have only 1 argument: the root path of the STEPS source dode\n"
if $#ARGV > 1;
$rootpath = "./";
if ($#ARGV == 1)
{
	$rootpath = $ARGV[0];
}

# Next, create a list of all files to measure.
$line = "cncc";
foreach (@code_dirs)
{
	$linebuf = $rootpath . $_ . "/*.";
     foreach (@code_sufx)
     {
         $line .= " " . $linebuf . $_;
     }
}

# Measure, output and exit.
exec ($line) or die "Error executing cncc command\n";
print "You might want to delete any remaining ci2/ci3/hi2/hi3 files...\n";

# END
