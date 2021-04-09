#!/usr/bin/perl
##---------------------------------------------------------------------------##
##  File:
##      @(#) bedToAlign.pl
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##      A simple script to convert alignment BED format to RepeatMasker
##      *.align files.
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

bedToAlign.pl - Convert BED to *.align

=head1 SYNOPSIS

  bedToAlign.pl [-version]

=head1 DESCRIPTION

A simple script for reversing the ./alignToBed -fullAlign *.align
conversion to a BED file.

The options are:

=over 4

=item -version

Displays the version of the program

=back

=head1 SEE ALSO

=head1 COPYRIGHT

Copyright 2021 Robert Hubley, Institute for Systems Biology

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

while (<>)
{
  if ( /^\S+\t\d+\t\d+\t(\S.*)/ ) 
  {
    my $align_str = $1;
    $align_str =~ s/\$/\n/g;
    print "$align_str\n";
  }
}
