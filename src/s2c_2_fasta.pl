#!/usr/bin/perl -w
################################################################
#                                                              #
#  Program:                                                    #
#     s2c_2_fasta.pl (version 0.1A)
#                                                              #
#  Copyright (C) 2010                                          #
#  Jianhong Ou @umassmed                                       #
#  All Rights Reserved.                                        #
#                                                              #
#  Author:                                                     #
#     Jianhong Ou
#                                                              #
#  Email:                                                      #
#     jianhong.ou@umassmed.edu                                 #
#                                                              #
#  Date Created:                                               #
#     Feb.25, 2013
#                                                              #
################################################################
use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;
use Carp;
use File::Spec::Functions;
#use Data::Dumper::Simple;
#use IO::All;

my $inFile;
my $outFile;
my $prefix  = "SEQ_";
my $nofilterN = 0;
my $nofilterL =0;
my $verbose = 0;
my $man     = 0;
my $help    = 0;

=head1 NAME

s2c_2_fasta.pl

=head1 SYNOPSIS

    perl s2c_2_fasta.pl --in <file> [options]

    Required:
    --in          <file>   s2c input file
    --out         <file>   fasta output file
    
    Optional:
    --prefix      <string> prefix of sequence name, default "SEQ_"
    --notfilterN           filter "N" or not, default filter.
    --notfilterLow         filter low complexity reads or not, default filter.
    --help                 Print help message
    --man                  Print manual page  

=head1 DESCRIPTION

    Convert a sequence-to-count (s2c) formated file to 
    restricted FASTA format for use with programs such as
    miRdeep, mirTools

=cut

my $argv = GetOptions    # get command line options
  (
	"in=s"  => \$inFile,
	"out=s" => \$outFile,
	"prefix=s" => \$prefix,
	notfilterN => \$nofilterN,
    notfilterLow => \$nofilterL,
	verbose => \$verbose,
	help    => sub { pod2usage(1) },
	man     => sub { pod2usage(-verbose => 2) },
  )
  or pod2usage(2);
pod2usage("--in is a required parameter\n")  unless (defined $inFile);
pod2usage("--out is a required parameter\n") unless (defined $outFile);

# convert the s2c data to a fasta file for BioProspector

open my $S2C, "<", $inFile or die "Couldn't open \"$inFile\": $!";
open my $FASTA, ">", $outFile or die "Couldn't open \"$outFile\": $!";
my $Nscnt = 0;
my $Ncnt =0;
my $Lscnt = 0;
my $Lcnt =0;
while (<$S2C>)
{
	if (/^([ACGTN]+)\t(\d+(?:\.\d+)?)/)
	{
		my $seq    = $1;
		my $counts = $2;
        if(!$nofilterL){
            if (!/(?=.*A)(?=.*C)(?=.*G)(?=.*T)/) {
                chomp;
                print STDERR "$_ low complexity \n" ;
                $Lscnt +=1;
                $Lcnt +=$counts;
                next;
            }
        }
		if($nofilterN){
			print $FASTA ">$prefix$._x$counts\n";
			print $FASTA "$seq\n";
		}else{
			if(/N/){
				chomp;
				print STDERR "$_ contain N \n" ;
				$Nscnt +=1;
				$Ncnt +=$counts;
			}else{
				print $FASTA ">$prefix$._x$counts\n";
				print $FASTA "$seq\n";
			}
		}

	}
	else
	{
		print STDERR "line $. of file ($inFile) is not in valid s2c format\n";
	}
}
close $S2C;
close $FASTA;
print "total No. of N containing sequence: $Ncnt (unique count: $Nscnt)\n" if(!$nofilterN);
print "total No. of low complexity sequence: $Lcnt (unique count: $Lscnt)\n" if(!$nofilterL);
