#!/usr/bin/perl
## setup_e - set up GCM modelE run
## This script is not supposed to be executed by users but rather
## is a part of Makefile functionality. It is run as a part of
##           gmake setup RUN=RunID
## Inside the Makefile it is called (from .../decks): setup_e RunID
## It looks for $HOME/.modelErc file and extracts all necessary options from
## there. If such file is not present the default options which are
## specified below will be used. Those default options are adjusted for
## working environment of ra.giss.nasa.gov.
##
## The options setup_e is looking for are:
## CMRUNDIR - directory to which all run directories will be linked.
## EXECDIR - path to directory with modelE scripts and with some executables.
## SAVEDISK - a directory (big) where all run directories will be created.
## GCMSEARCHPATH - directory to search for gcm input files.
## MAILTO - email address of the user. (default `whoami`)
## UMASK - the value of 'umask' to be used for model runs.

## default settings
$CMRUNDIR="/u/cmrun";
$EXECDIR="/u/exec";
$SAVEDISK="/raid1";
$GCMSEARCHPATH="/u/cmrun";
$MAILTO="";
$UMASK=002;

## if $HOME/.modelErc is present get settings from there

$modelerc = (getpwuid($>))[7]."/.modelErc";

if ( -f $modelerc ) {
    print "Using settings from ~/.modelErc\n";
    open MODELERC, $modelerc or die "can't open $modelerc";
    while(<MODELERC>) {
	$CMRUNDIR = $1 if /^ *CMRUNDIR *= *(\S+)/;
	$EXECDIR = $1 if /^ *EXECDIR *= *(\S+)/;
	$SAVEDISK = $1 if /^ *SAVEDISK *= *(\S+)/;
	$GCMSEARCHPATH = $1 if /^ *GCMSEARCHPATH *= *(\S+)/;
	$MAILTO = $1 if /^ *MAILTO *= *(\S+)/;
	$UMASK = oct "$1" if /^ *UMASK *= *(\S+)/;
    }
    close MODELERC;
} else {
    print "$HOME/.modelErc is not present. Using default settings.\n";
}

if ( $MAILTO =~ /^\s*$/ ) { $MAILTO = `whoami`; chop $MAILTO; }
print "CMRUNDIR = $CMRUNDIR\n";
print "EXECDIR = $EXECDIR\n";
print "GCMSEARCHPATH = $GCMSEARCHPATH\n";
print "SAVEDISK = $SAVEDISK\n";
print "MAILTO = $MAILTO\n";
printf "UMASK = %03lo\n", $UMASK;

if ( $#ARGV != 0 ) { 
    print "This script is not supposed to be run outside of Makefile\n";
    print "Inside the Makefile it is called: setup_e.pl runID\n";
    print "$#ARGV\n";
    exit 1; 
    }

$runID = shift;
$rfile = "$runID.R";
$RunDir = "$SAVEDISK/$runID";
$lockfile = "$RunDir/lock";
$CMRUN = $CMRUNDIR;
$DeckDir = `pwd`; chop $DeckDir;
$umask_inv = $UMASK ^ 0777;
$umask_str = sprintf "%lo", $UMASK;

## check if this run is already running
if ( -f $lockfile ) {
    print "            **********************                \n";
    print "$1 seems to be already running in $SAVEDISK/$runID\n";
    print "If you think it is an error, then most probably this\n";
    print "task was interrupted in an unusual way. Please check.\n";
    print "Then remove the lock file:\n";
    print "$lockfile\n";
    print "and restart the setup.\n";
    print "            **********************                \n";
    exit 1;
}

## check if the rundeck is present
if ( ! -s $rfile ) {
   print "File does not exist: $rfile\n";
   exit 1;
}

print "setting up run $runID\n";
print "output files will be saved in $SAVEDISK/$runID\n";

## Create the run directory if necessary
if ( ! -d $RunDir ) {
    mkdir $RunDir, 0777 & $umask_inv 
	or die "Can't create $RunDir. Aborting setup.\n";
}

## Check that link is not already correct (added by gavin)
if ( -e $runID ) {
    if ( `ls -l $runID` !~ /-> *$RunDir$/ ) {
	print "./$runID exists and is pointing to something else:\n";
	print "Please check. Aborting setup.\n";
	exit 1;
    }
} else {
    symlink $RunDir, $runID or die "Can't create link $runID in local dir";
}

## Also link to $CMRUN if different from $SAVEDISK (added by gavin)
if ( $CMRUN ne $SAVEDISK ) {
    if ( -e "$CMRUN/$runID" ) {
	die "Can't create a link $CMRUN/$runID" if ! -l "$CMRUN/$runID";
	if (`ls -l $CMRUN/$runID` !~ /-> *$RunDir$/ ) {
	    unlink "$CMRUN/$runID" or die "Can't rm old link $CMRUN/$runID";
	    symlink "$RunDir", "$CMRUN/$runID" or die "Can't create link";
	}
    } else {
        symlink "$RunDir", "$CMRUN/$runID" or die "Can't create link";
    }
}

## Make sure that we really can write to dir/link $runID
if ( ! ( -d $runID && -r $runID && -w $runID && -x $runID ) ) {
    print "Couldn't create the link $runID\n";
    print "or run directory has wrong permissions. Aborting setup.\n";
    exit 1;
}

## Switching to Run directory
chdir "$runID" or die "Can't chdir to $runID : $!\n";

open PRT, ">$runID".".PRT" or die "can't open ${runID}.PRT for writing\n";
print PRT "0Run $runID\n";

## Move executable from "$runID"_bin directory to $runID
if ( ! -f "$DeckDir/$runID"."_bin/$runID.exe" ) {
    print "$DeckDir/$runID"."_bin/$runID.exe","not found\n";
    exit 1;
}
`mv -f "$DeckDir/$runID"_bin/$runID.exe $runID.exe`;
chmod 0777 & $umask_inv, "$runID.exe";

open RFILE, "$DeckDir/$rfile" or die "can't open $DeckDir/$rfile";

## Check signature
$_ = <RFILE>;
if ( ! /^\s*$rfile/ ) {
    print "inconsistent Naming: $rfile is not $_\n" ;
    exit 1;
}
print PRT " $_";

## Read rfile until "Data input files:" is encountered
while (<RFILE>) {
    print PRT " $_";
    if ( /Data *input *files/ ) { last; }
}

## Use the subsequent information to create the link and unlink files
while (<RFILE>) {
    if ( /Label *and *Namelist/ ) { last; }
    print PRT " $_";
    s/!.*//;
    push @data_files, /([\w.+_-]+\s*=\s*[\w.\/+_-]+)/g;
}

open RUNIDLN, ">$runID"."ln" or die "can't open ${runID}ln for writing\n";
open RUNIDULN, ">$runID"."uln" or die "can't open ${runID}ln for writing\n";

foreach $_ ( @data_files ) {
    ($name, $dest) = split /\s*=\s*/;
    if ( ! -e "$GCMSEARCHPATH/$dest" ) {
	print "$dest not found in $GCMSEARCHPATH\n";
	exit 1;
    }
    if ( $name !~ /^(AIC|OIC|GIC)$/ ) {
	print RUNIDLN "ln -s $GCMSEARCHPATH/$dest $name\n";
	print RUNIDULN "rm $name\n";
    } else {
	`ln -s $GCMSEARCHPATH/$dest $name` ; 
	print "using $GCMSEARCHPATH/$dest for IC only\n";
    }
}

## Architecture-dependent settings
$uname = `uname`;
if ( $uname =~ /IRIX64/ ) {
    print RUNIDLN "export PAGESIZE_DATA=64 PAGESIZE_STACK=64\n";
} elsif ( $uname =~ /AIX/ ) {
    print RUNIDLN "export NAMELIST=OLD\n";
}

close RUNIDLN;
close RUNIDULN;
chmod 0777 & $umask_inv, "${runID}ln", "${runID}uln";

open I, ">I" or die "can't open 'I' for writing\n";
## Check signature
$_ = <RFILE>;
if ( ! /^\s*$runID/ ) {
    print "inconsistent Naming: $rfile is not start of $_\n" ;
    exit 1;
}
print I;
$_ = <RFILE>;
print I;
while (<RFILE>) {
    s/!.*//g;
    s/^\s+//;
    print I " $_";
}
close I;
close PRT;
chmod 0666 & $umask_inv, "I", "${runID}.PRT";

## Run the 1st hour then clean up
print "starting the execution \n";
print "current dir is ", `pwd`;

## Switching to background
if ($pid = fork) {
    print "Starting 1st hour in the background.\n\n";
    exit 0;
}

if ( ! defined $pid ) {
    print "Fork failed. Continuing in foreground.\n";
}

## Running the model
print <<`EOC`;
    umask $umask_str
    touch lock
    rm -f error_message 2> /dev/null
    ./"$runID"ln
    ./"$runID".exe < I >> ${runID}.PRT
    rc=\$?
    rm -f AIC GIC OIC
    ./"$runID"uln
    rm -f lock
    exit \$rc
EOC

$rcode = $? >> 8;
if ( $rcode != 13 && $rcode != 12 ) {
    print " Problem encountered while running hour 1 :\n"; 
    $error_message = `cat error_message`; chop $error_message;
    print " >>> $error_message <<<\n";
    exit 4 ;
} else {
    print "1st hour completed sucessfully\n";
}

## Create executable script file RUNID
open RUNID, ">$runID" or die "can't open $runID for writing\n";
print RUNID <<EOF;
\#!/bin/sh
    umask $umask_str
    if [ -f lock ] ; then
      echo 'lock file present - aborting' ; exit 1 ; fi
    touch lock
    rm -f error_message 2> /dev/null
    PRTFILE=${runID}.PRT
    if [ \$\# -ge 1 ] && [ $1='-q' ]; then
      PRTFILE='/dev/null'; fi
    ./${runID}ln
    ./${runID}.exe < I > \$PRTFILE
    rc=\$?
    ./${runID}uln
    rm -f lock
    exec $EXECDIR/runpmE $runID \$rc $MAILTO
EOF
close RUNID;
chmod 0777 & $umask_inv, $runID;
## end of RUNID script
`ln -sf $runID E`;

## overwrite RUNID.I omitting ISTART=.. line.
`umask $umask_str; grep -v "ISTART=" I > I.tmp; mv -f I.tmp I`;

