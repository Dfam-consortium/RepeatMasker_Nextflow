#!/usr/bin/env perl
##---------------------------------------------------------------------------##
##  File:
##      @(#) renumberIDs.pl
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##      Convert batch-uniq identifiers back to numeric only 
##      identifiers. 
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

renumberIDs.pl - Convert batch uinique identifiers to numeric only

=head1 SYNOPSIS

  renumberIDs.pl [-translation file] *.out or *.align

=head1 DESCRIPTION

A simple script to convert unique alpha-numeric identifiers
to numeric ascending ones.  If a translation file is not 
provided the ID '1' will be the first used.

The options are:

=over 4

=item -translation file

A simple tab delimited file in the following format:

oldID\tnewID

=back

=head1 SEE ALSO

=head1 COPYRIGHT

Copyright 2023 Robert Hubley, Institute for Systems Biology

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
    '-translation=s'
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

my $useTranslation = 0;
my %translation = ();
if ( $options{'translation'} ) {
  open IN,"<$options{'translation'}" or die;
  while ( <IN> ) {
    my @flds = split(/\t/);
    $translation{$flds[0]} = $flds[1];
  }
  close IN;
  $useTranslation = 1;
}

my %lookup = ();
my $newID;
my $cIdx = 1;
while (<>)
{
  # OUT:
  #   271 20.2  5.1  0.0 CM016569.1            3257     3335  (12585) C rnd-5_family-191 LTR/ERV1            (0)     261     179 b1_12       *
  # ALIGN:
  #  19 21.50 1.11 9.64 CM016569.1 2534 2623 (13297) (AAAT)n#Simple_repeat 1 83 (0) b1_10 m_b1s252i1
  if ( /^\s*(\d+\s+\d+\.\d+\s+\d+\.\d+\s+\d+\.\d+.*)/ )
  {
    my @fields = split(/\s+/, $1);
    my $f_idx = $#fields;
    $f_idx-- if ( $fields[$f_idx] eq "*" || $fields[$f_idx] =~ /\[?[mc]_b\d+s/ );
    if ( $useTranslation ) {
      if ( exists $translation{$fields[$f_idx]} ) {
        $newID = $translation{$fields[$f_idx]};
      }else {
        die "Could not find translation for ID: $fields[$f_idx]\n    Found in line: $_";
      }
    }else {
      if ( exists $lookup{$fields[$f_idx]} ) {
        $newID = $lookup{$fields[$f_idx]};
      }else {
        $newID = $cIdx++;
        $lookup{$fields[$f_idx]} = $newID;
      }
    }
    s/ $fields[$f_idx] / $newID /;
  }
  print;
}

open OUT,">translation-out.tsv" or die;
foreach my $key ( keys(%lookup) ) {
  print OUT "$key\t" . $lookup{$key} . "\n";
}
close OUT;

