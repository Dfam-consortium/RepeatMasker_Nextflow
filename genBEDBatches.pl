#!/usr/bin/env perl
use strict;
use Data::Dumper;

# Where this script can find twoBitInfo 
my $ucscBinDir = "/opt/ucsc_tools";
if ( -d $ENV{'UCSCTOOLSDIR'} && -s $ENV{'UCSCTOOLSDIR'} . "/twoBitInfo" ) {
  $ucscBinDir = $ENV{'UCSCTOOLSDIR'};
  $ucscBinDir .= "/" if ( $ucscBinDir !~ /^.*\/$/ );
}

my $seqFile = $ARGV[0];
my $maxSize = $ARGV[1];


my $batches = getBatches( $seqFile, $maxSize );

my $batchNum = 1;
my $seqIdx = 1;
foreach my $batch ( @{$batches} ) {
    my $fileName = "batch-$batchNum.bed";
    open OUT,">$fileName" or die;
    foreach my $seq ( @{$batch} )
    {
      # seqID start end name seqID_len
      print OUT "".$seq->[0] . "\t" . 
                   $seq->[1] . "\t" . 
                   $seq->[2]  . "\t" .
                   "seq-$seqIdx" . "\t" .
                   $seq->[3] . "\n";
      $seqIdx++;
    }
    close OUT;
    print "batch-$batchNum.bed\n";
    $batchNum++;
}


sub getBatches {
  my $seqFile = shift;
  my $maxSeqSize = shift;
  my $numBatches = shift;
  
  # 
  # Tabulate original sizes
  #
  my $totalSize = 0;
  my %seqSizes = ();
  my @batches = ();
  my $batchSize = 0;
  my @seqList = ();

  if ( $seqFile =~ /.*\.2bit$/ ) {
    open IN, "$ucscBinDir" . "twoBitInfo $seqFile stdout|"
      or die "Could not run $ucscBinDir" . "twoBitInfo on $seqFile!\n";
    while ( <IN> )
    {
      if ( /^(\S+)\s+(\d+)/ )
      {
        $seqSizes{$1} = $2;
      }
    }
    close IN;
  }else {
    if ( $seqFile =~ /.*\.gz$/ )
    {
      open IN,"gunzip -c $seqFile|" or die;
    }else {
      open IN,"<$seqFile" or die;
    }
    my $seq;
    my $id;
    while (<IN>){
      if ( /^>(\S+)/ ){
        my $tmpId = $1;
        if ( $seq ne "" ) 
        {
          $seqSizes{$id} = length($seq);
        }
        $seq = "";
        $id = $tmpId;
        next;
      }
      s/[\n\r\s]+//g;
      $seq .= $_;
    }
    if ( $seq )
    {
      $seqSizes{$id} = length($seq);
    }
    $seq = "";
    $id = "";
    close IN;
  }

  if ( defined $numBatches && $numBatches > 0 )
  {
    $maxSeqSize = int($totalSize / $numBatches) + 1;
  }

  foreach my $seq ( keys(%seqSizes ) )
  {
    my $seqSize = $seqSizes{$seq};
    if ( ($batchSize + $seqSize) <= $maxSeqSize )
    {
      # Add it and move on 
      push @seqList, [ $seq, 0, $seqSize, $seqSizes{$seq} ];
      $batchSize += $seqSize;
      next;
    }

    # First chunk
    my $firstChunkSize = $maxSeqSize - $batchSize;
    push @seqList, [$seq, 0, $firstChunkSize, $seqSizes{$seq}];
    push @batches, [ @seqList ];
    $batchSize = 0;
    @seqList = ();


    # Middle chunks
    my $start = $firstChunkSize;
    for ( my $i = 0; $i < int(($seqSize-$firstChunkSize)/$maxSeqSize); $i++ )
    {
      push @seqList, [$seq, $start, $start + $maxSeqSize, $seqSizes{$seq}];
      push @batches, [ @seqList ];
      $batchSize = 0;
      @seqList = ();
      $start += $maxSeqSize;
    }
        
    # Last chunk
    if ( $start < $seqSize ) 
    {
       push @seqList, [$seq, $start, $seqSize, $seqSizes{$seq}];
       $batchSize += $seqSize - $start;
    }
  }
  if ( @seqList )
  {
    push @batches, [ @seqList ];
  }
  return \@batches;
} # getBatches
