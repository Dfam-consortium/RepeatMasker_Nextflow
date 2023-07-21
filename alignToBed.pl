#!/usr/bin/env perl
##---------------------------------------------------------------------------##
##  File:
##      @(#) alignToBed.pl
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##      A simple script to convert RepeatMasker *.out files to 
##      BED format.  
##
#******************************************************************************
#*  This software is provided ``AS IS'' and any express or implied            *
#*  warranties, including, but not limited to, the implied warranties of      *
#*  merchantability and fitness for a particular purpose, are disclaimed.     *
#*  In no event shall the authors or the Institute for Systems Biology        *
#*  liable for any direct, indirect, incidental, special, exemplary, or       *
#*  consequential damages (including, but not limited to, procurement of      *
#*  substitute goods or services; loss of use, data, or profits; or           *
#*  business interruption) however caused and on any theory of liability,     *
#*  whether in contract, strict liability, or tort (including negligence      *
#*  or otherwise) arising in any way out of the use of this software, even    *
#*  if advised of the possibility of such damage.                             *
#*                                                                            *
#******************************************************************************
#
# ChangeLog
#
#     $Log$ 
#
###############################################################################
#
# To Do:
#
=head1 NAME

alignToBed.pl - Convert RepeatMasker *.align to BED format

=head1 SYNOPSIS

  alignToBed.pl [-version] [-noSimple] [-fullAlign] 
                *.align

=head1 DESCRIPTION

A simple script to convert RepeatMasker's default alignment
format ( *.align ) to BED format.  For more details on BED
see: http://genome.ucsc.edu/FAQ/FAQformat.html#format1

The *.align fields are converted as:

=begin text

    BED Field               RepeatMasker Field
    =========               ==================
    chrom                   Query Sequence Identifier
    chromStart              Query Start - 1
    chromEnd                Query End
    name                    The entire *.out line 
 

=end text

The options are:

=over 4

=item -version

Displays the version of the program

=item -noSimple

Filter out repeats with "Simple_repeat" or "Low_complexity" 
repeat classes.

=item -fullAlign

Encode the full alignment data in the BED 'name' field.  Line
terminations are encoded using the '$' character.

=back

=head1 SEE ALSO

=head1 COPYRIGHT

Copyright 2014 Robert Hubley, Institute for Systems Biology

=head1 LICENSE

This code may be used in accordance with the Open Source License v. 3.0
http://opensource.org/licenses/OSL-3.0

=head1 AUTHOR

Robert Hubley <rhubley@systemsbiology.org>

=cut

#
# Module Dependence
#
use strict;
use Getopt::Long;
use Data::Dumper;

my $Version = "1.0";

#
# Magic numbers/constants here
#  ie. my $PI = 3.14159;
#   
my $DEBUG = 0;

#
# Option processing
#  e.g.
#   -t: Single letter binary option
#   -t=s: String parameters
#   -t=i: Number paramters
#
my @getopt_args = (
    '-version', # print out the version and exit
    '-noSimple',
    '-fullAlign'
);

my %options = ();
Getopt::Long::config("noignorecase", "bundling_override");
unless (GetOptions(\%options, @getopt_args)) {
    usage();
}

sub usage {
  print "$0 - $Version\n\n";
  exec "pod2text $0";
  exit;
}

if ($options{'version'}) {
  print "$Version\n";
  exit;
}


#
# ARGV Processing
#
#if ( !@ARGV  ) {
#  usage();
#}

my $align_str = "";
my $summ_str = "";
my $bed_id = "";
my $bed_start = 0;
my $bed_end = 0;
while (<>)
{
  if ( /^\s*(\d+\s+\d+\.\d+\s+\d+\.\d+\s+\d+\.\d+.*)/ )
  {
    $align_str = "";
    next if ( $options{'noSimple'} &&
              /Simple_repeat|Low_complexity/ );
    my @fields = split(/\s+/, $1);
    if ( $fields[8] ne "C" && $fields[8] ne "+" )
    {
      splice( @fields, 8, 0, "+");
    }
    $summ_str = join(" ",@fields);
    # Zero Based Half Open Coords
    $bed_id = $fields[4];
    $bed_start = ($fields[5] - 1 );
    $bed_end = $fields[6];
  }
  $align_str .= $_;
  # Gap_init rate = 0.01 (3 / 223), avg. gap size = 1.67 (5 / 3)
  if ( /^Gap_init/ ) {
    if ( $options{'fullAlign'} ) {
      $align_str =~ s/\n/\$/g;
      print "$bed_id\t$bed_start\t$bed_end\t$align_str\n";
    }else {
      print "$bed_id\t$bed_start\t$bed_end\t$summ_str\n";
    }
  }
}
