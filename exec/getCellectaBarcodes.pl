#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;

##############################
# By Matt Cannon
# Date: 4/22/21
# Last modified:
# Title: getCellectaBarcodes.pl
# Purpose: Pull Cellecta barcode information out of bam files
##############################

##############################
# Options
##############################

my $verbose;
my $help;
my $sam;
#my $outputStub = "";
my $printBuffer = 1000;

# i = integer, s = string
GetOptions ("verbose"           => \$verbose,
            "help"              => \$help,
            "sam=s"             => \$sam,
            "printBuffer=i"     => \$printBuffer
      ) or pod2usage(0) && exit;

#                  "outputStub=s"      => \$outputStub, # kicking this out as I'm planning on pulling it straight into R


pod2usage(1) && exit if ($help);

##############################
# Global variables
##############################
my $lineCounter = 0;
my $foundCounter = 0;
my @printArray = ("cid\tlt_14\tlt_30"); # put header in array to begin

##############################
# Code
##############################
## cid.fq is in CB:Z:
## lineage tracing barcode is in the read

##############################
### Stuff
### More stuff

#open my $ltOutFH, '>', $outputStub . "lt.fq";

# Let me pull in bam files
$sam =~ s/(.*\.bam)\s*$/samtools view < $1|/;

##############################
### Go through sam/bam file and pull out lineage tracing info for each line

open my $samFH, "$sam" or die "Could not open sam input\n";
while (my $input = <$samFH>){
    chomp $input;
    # Skip headers
    if($input !~ /^@/) {
      if ($input =~ /CGCA([ACGT]{14})TGGT([ACGT]{30})TGGT.+CB:Z:([ATGC]{16})/) {
          my $lt_14 = $1;
          my $lt_30 = $2;
          my $cid = $3;

          printOut($cid, $lt_14, $lt_30);

          $foundCounter++;
      }
      if($verbose) {
          $lineCounter++;
          if($lineCounter % 1000000 == 0) {
              print STDERR commify($lineCounter), " sam entries processed. ", commify($foundCounter), " barcodes found                                  \r"
          }
      }
    }
}

##############################
### Print out last few lines
if(scalar(@printArray) > 0) {
    print join("\n", @printArray), "\n";
}

##############################
### Subfunctions

sub printOut {
    my $cid = $_[0];
    my $lt_14 = $_[1];
    my $lt_30 = $_[2];

    push @printArray, $cid . "\t". $lt_14 . "\t" . $lt_30;

    if(scalar(@printArray) >= $printBuffer) {
        print join("\n", @printArray), "\n";
        @printArray = ();
    }
}

sub commify { # function stolen from web
    my $text = reverse $_[0];
    $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
    return scalar reverse $text;
}

##############################
# POD
##############################

#=pod

=head SYNOPSIS

Summary:

    getCellectaBarcodes.pl - pull out and print

Usage:

    perl xxxxxx.pl [options]


=head OPTIONS

Options:

    --verbose
    --help

=cut
