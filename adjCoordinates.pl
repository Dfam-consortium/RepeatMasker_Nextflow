#!/usr/bin/env perl
use strict;
use Data::Dumper;
#use lib "/lustre/work/daray/software/RepeatMasker";
use lib "/usr/local/RepeatMasker";
use SearchResult;
use SearchResultCollection;
use CrossmatchSearchEngine;

my $batchBED = $ARGV[0];
my $resultFile = $ARGV[1];

my %seqRanges = ();
open IN,"<$batchBED" or die;
while (<IN>){
  my @flds = split();
  # seqID start end name seqID_len
  # e.g chr1 1 100 seq-1 200
  $seqRanges{$flds[3]} = [$flds[0],$flds[1],$flds[4]];
}
close IN;


my $resultCollection =
    CrossmatchSearchEngine::parseOutput( searchOutput => $resultFile );

open OUT,">$resultFile.adjusted" or die;
for ( my $i = 0 ; $i < $resultCollection->size() ; $i++ ) {
  my $result = $resultCollection->get( $i );
  my $qID    = $result->getQueryName();
  my $qBeg   = $result->getQueryStart();
  my $qEnd   = $result->getQueryEnd();
  my $qRem   = $result->getQueryRemaining();
  my $qSeq   = $result->getQueryString();

  if ( $seqRanges{$qID} ) {
    $result->setQueryName($seqRanges{$qID}->[0]);
    $qBeg += $seqRanges{$qID}->[1];
    $result->setQueryStart($qBeg);
    $qEnd += $seqRanges{$qID}->[1];
    $result->setQueryEnd($qEnd);
    $result->setQueryRemaining($seqRanges{$qID}->[2] - $qEnd);
  }else {
    die "Could not find $qID in BED file!\n";
  }
 
  if ( $qSeq ) { 
    print OUT "" . $result->toStringFormatted( SearchResult::AlignWithQuerySeq );
  }else {
    print OUT "" . $result->toStringFormatted( SearchResult::OutFileFormat );
  }
}
close OUT;

#open IN,"<$resultFile" or die;
#while ( <IN> ) 
#{
#  # RepeatMasker *.out format
#  if ( /^\s*(\d+\s+\d+\.\d+\s+\d+\.\d+.*)/)
#  {
#    my $line = $1;
#    my @flds = split(/\s+/,$line);
##    my $seqID = $flds[4];
#    my $seqStart = $flds[5];
#    my $seqEnd = $flds[6];
#    if ( $seqRanges{$seqID} ) {
#      $flds[4] = $seqRanges{$seqID}->[0];
#      $flds[5] += $seqRanges{$seqID}->[1];
#      $flds[6] += $seqRanges{$seqID}->[1];
#      $flds[7] = "(" . ($seqRanges{$seqID}->[2] - $flds[5]) . ")";
#    }else {
#      die "Could not find $seqID in BED file! Line = $_\n";
#    }
#    print OUT "" . join(" ",@flds) . "\n";
#  }else{
#    print OUT;
#  }
#}
#close IN;
##close OUT;

