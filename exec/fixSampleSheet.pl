#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;


##############################
# By Matt Cannon
# Date: 
# Last modified: 
# Title: .pl
# Purpose: 
##############################

##############################
# Options
##############################


my $verbose;
my $help;
my $sampleInfo;
my $sampleSheet;

# i = integer, s = string
GetOptions ("verbose"           => \$verbose,
            "help"              => \$help,
            "sampleInfo=s"      => \$sampleInfo,
            "sampleSheet=s"     => \$sampleSheet
      ) or pod2usage(0) && exit;

pod2usage(1) && exit if ($help);


##############################
# Global variables
##############################
my %sampleInfoHash;
my $printing = "TRUE";
my $sampleNum = 0;

##############################
# Code
##############################

##############################
### Build hash of arrays of sample info using column header as key

open my $siFH, "$sampleInfo" or die "Could not open sample info input\n$!";

my $siHeader = <$siFH>;
chomp $siHeader;
my @siHeaderArray = split "\t", $siHeader;

while (my $input = <$siFH>){
    chomp $input;
    my @lineArray = split "\t", $input;

    for (my $i = 0; $i < scalar(@siHeaderArray); $i++) {
        push @{ $sampleInfoHash{$siHeaderArray[$i]} }, $lineArray[$i];
    }

    $sampleNum++;
}

##############################
### Stuff
### More stuff

open my $ssFH, "$sampleSheet" or die "Could not open sample sheet input\n$!";
while (my $input = <$ssFH>){
    chomp $input;
    if ($printing eq "TRUE") {
        print $input, "\n";
    }

    if ($input =~ /^Lane/) {
        my @columnNames = split ",", $input;

        for (my $i = 0; $i < $sampleNum; $i++) {
            for (my $j = 0; $j < scalar(@columnNames); $j++) {
                if(exists($sampleInfoHash{$columnNames[$j]})) {
                    print $sampleInfoHash{$columnNames[$j]}[$i];
                }

                # Don't want a comma after the last column
                print "," if ($j < scalar(@columnNames) - 1);
            }
            print "\n";
        }
        $printing = "FALSE";
    }
}


#open R1OUTFILE, ">", $outputNameStub . "_R1.fastq";

##############################
# POD
##############################

#=pod
    
=head SYNOPSIS

Summary:    
    
    xxxxxx.pl - generates a consensus for a specified gene in a specified taxa
    
Usage:

    perl xxxxxx.pl [options] 


=head OPTIONS

Options:

    --verbose
    --help

=cut
