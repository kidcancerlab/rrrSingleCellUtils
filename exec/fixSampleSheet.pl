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
my $debug;
my $sampleInfo;
my $sampleSheet;

# i = integer, s = string
GetOptions ("verbose"           => \$verbose,
            "help"              => \$help,
            "debug"             => \$debug,
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
my $printBuffer;

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
    print STDERR $sampleNum, "\n" if ($debug);
}
close $siFH;

##############################
### Read in the original sample sheet, drop the excess header and replace with our info

open my $ssFH, "$sampleSheet" or die "Could not open sample sheet input\n$!";
while (my $input = <$ssFH>){
    chomp $input;

    print STDERR $input, "\n" if ($debug);

    # Keep only the data header
    # If the header is kept, the parsed reads have varying lengths
    if ($input =~ /^\[Data\]/) {
        $printBuffer .= $input . "\n";
    }

    # Once we get to the "Lane" line of the file, print out headers and all data
    # from the sample info sheet
    if ($input =~ /^Lane/) {
        $printBuffer .= $input . "\n";
        my @columnNames = split ",", $input;

        for (my $i = 0; $i < $sampleNum; $i++) {
            for (my $j = 0; $j < scalar(@columnNames); $j++) {
                if (exists($sampleInfoHash{$columnNames[$j]})) {
                    $printBuffer .= $sampleInfoHash{$columnNames[$j]}[$i];
                }

                # Don't want a comma after the last column
                $printBuffer .= "," if ($j < scalar(@columnNames) - 1);
            }
            $printBuffer .= "\n";
        }
    }
}

close $ssFH;

# Wait to print until the end so we can overwrite the original file
print $printBuffer;

##############################
# POD
##############################

#=pod
    
=head SYNOPSIS

Summary:    
    
    xxxxxx.pl - 
    
Usage:

    perl xxxxxx.pl [options] 


=head OPTIONS

Options:

    --verbose
    --help

=cut
