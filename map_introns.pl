#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

# page coords = 612 x 792

my $SCALE  = 1;
my $HEADER = 1;
my $SCALEBAR_WIDTH = 100;

parseArgs();

my $SCALEBAR_LABEL = $SCALEBAR_WIDTH;

my $INFILE = shift @ARGV;

$SCALEBAR_WIDTH *= $SCALE;

my $y = 750;
my $left_margin = 40;
my $space_between_sequences = 20;
my $with_data_color = "black";
my $no_data_color   = "gray";
my @file; #structure --> ( [$id, [@introns], [@sequence]], [$id, [@introns], [@sequence]], [$id, [@introns], [@sequence]] )

my $outfile = $INFILE;
$outfile =~ s/\..*$//;
$outfile .= ".ps";

my $pdf = $outfile;
$pdf =~ s/\..*$/\.pdf/;

my $header = $INFILE;
$header =~ s/\..*$//;

parseFasta( \@file );

# includes left margin and space to accommodate the longest name; all features are drawn relative to this left-most point
my $feature_offset = $left_margin + 160; 

# print postscript header
open( OUT, ">$outfile" ) || die "Can't open $outfile: $!\n";
print OUT "%!PS\n";
close OUT;

# print filename header
if( $HEADER ){
  open( OUT, ">$outfile" ) || die "Can't open $outfile: $!\n";
  print OUT "/Arial findfont\n";
  print OUT "12 scalefont\n";
  print OUT "setfont\n";
  print OUT "0 setgray\n";
  
  print OUT $left_margin, " ", $y, " moveto\n";
  print OUT "($header) show\n";
  $y -= 40;
}

# iterate over all sequences in the file -- the real business starts here
for( my $i = 0; $i < @file; $i++ ){
  # print $file[$i][1][0] and exit;
  my( $counter, $j, $seq, $id, @missing_coords, @data_coords, @ranges, @introns, @edits );


  # NOTE: CHANGE '[n]' to '[n-]' TO COLOR MISSING DATA AND GAPS THE SAME WAY; COULD ADD AN ELSIF
  # TO GIVE GAPS A DIFFERENT COLOR
  for( $j = 0; $j < @{$file[$i][2]}; $j++ ){
    if( $file[$i][2][$j] =~ /[n]/i ){ #find and store coordinates for regions with missing data
      push @missing_coords, $j+1;
    }else{                             #find and store coordinates for regions with data
      push @data_coords, $j+1;
    }
  }
  
  # get ranges of parts of sequence *with* data
  if( @data_coords ){
    @data_coords = sort {$a <=> $b} @data_coords;
    getRanges( \@data_coords, \@ranges, $with_data_color );
  }

  # get ranges for gaps and missing data
  if( @missing_coords ){
    @missing_coords = sort {$a <=> $b} @missing_coords;
    getRanges( \@missing_coords, \@ranges, $no_data_color );
  }
  
  #write taxon name
  writeName( $left_margin, $y, $file[$i][0] );
  
  #draw line segments (gray if missing data, black if there is data)
  @ranges = sort {$a->[0] <=> $b->[0]} @ranges;
  draw_gene( $feature_offset, $y, \@ranges, $outfile );
  # print "$left_margin\n$y\n$ranges[0][0]\n"

  #map introns
  mapIntrons( $feature_offset, $y, $file[$i][1], $file[$i][2], $outfile );
  
  $y -= $space_between_sequences;

}

#draw scale bar
scaleBar( $y, $outfile );

#covert to pdf
system( "ps2pdf $outfile $pdf" );

exit;
######################################################SUBROUTINES######################################################
sub parseFasta{
  my $file = shift;
  my $tmp;

  #read file into string
  open( FASTA, $INFILE ) || die "Couldn't open $INFILE: $!\n";
  while( <FASTA> ){
    $tmp .= $_;
  }
  close FASTA;

  # parse string into array, each element is an array => ($header,[@introns],[@sequence])
  foreach my $sequence ( split />/, $tmp ){
    $sequence =~ /\S+/ or next;
    my( $in_sequence, $seq, @introns, $id );
    foreach( split /\n/, $sequence ){
      if( $in_sequence ){
	$seq .= $_;
      }else{
	@introns = split( /\s+/, $_ );
	$id = shift @introns;
	$in_sequence = 1;
      }
    }
    $seq =~ s/[0-9\s]//g;
    my @chars = split( //, $seq );
    my $R = [$id, [@introns], [@chars]];
    push @$file, $R;
  }
  
  $tmp = undef;
  
}
########################################################################################################################
sub getRanges{
  my( $coords, $coord_ranges, $color ) = @_;
  my $range_cnt = 0;
  my $maxgap    = 1;
  my ( @ranges, $i );

  $ranges[0][0] = $coords->[0];
  for ($i = 1; $i < @$coords; $i++) {
    if ($coords->[$i] - $coords->[$i - 1] > $maxgap) {
      $ranges[$range_cnt][1] = $coords->[$i - 1];
      $range_cnt++;
      $ranges[$range_cnt][0] = $coords->[$i];
    }
  }
  $ranges[$range_cnt][1] = $coords->[$#$coords];
  
  for ($i = 0; $i < @ranges; $i++) {
    my $R = [$ranges[$i][0], $ranges[$i][1], $color];   
    push @$coord_ranges, $R;
    # print $color . "\n";
    # printf "%d\t%d\t%d\t%d\n",  $i + 1, $ranges[$i][1] - $ranges[$i][0] + 1, $ranges[$i][0], $ranges[$i][1];
  }

}
########################################################################################################################
sub draw_gene {
  my( $left_margin, $y, $coord_ranges, $outfile ) = @_;
  my( $k, $line_color );
  
  open( OUT, ">>$outfile" ) || die "Can't open $outfile: $!\n";

  print OUT "1 setlinewidth\n";

  foreach( $k = 0; $k < @$coord_ranges; $k++ ){
    $$coord_ranges[$k][2] eq "black" and $line_color = "0 setgray\n";
#   $$coord_ranges[$k][2] eq "gray"  and $line_color = "0.6 setgray\n";
    $$coord_ranges[$k][2] eq "gray"  and $line_color = "1 0 0 setrgbcolor\n";
    if( $k == 0 ){
      print OUT "newpath\n";
      print OUT "$left_margin $y moveto\n", $left_margin + (${$coord_ranges}[$k][1] * $SCALE), " $y lineto\n";
      print OUT $line_color;
      print OUT "stroke\n";
    }else{
      print OUT "newpath\n";
      print OUT $left_margin + (${$coord_ranges}[$k-1][1] * $SCALE), " $y moveto\n", $left_margin + (${$coord_ranges}[$k][1] * $SCALE), " $y lineto\n";
      print OUT $line_color;
      print OUT "stroke\n";
    }

    #print "${$coord_ranges}[$k][2] ${$coord_ranges}[$k][0] ${$coord_ranges}[$k][1]\n";
  }
  close OUT; 
}
########################################################################################################################
sub mapIntrons{
  my( $left_margin, $y, $introns, $sequence, $outfile ) = @_;
  my $intron_offset = $SCALE * 0.5; #put intron marker midway between exon coordinates
  my $missing;

  open( OUT, ">>$outfile" ) || die "Can't open $outfile: $!\n";

  print OUT "0 setgray\n";
  print OUT "0.5 setlinewidth\n";

  for( my $i = 0; $i < @$introns; $i++ ){
    if( $introns->[$i] =~ /m/ ){
      $missing = 1;
      $introns->[$i] =~ s/m//;
    }
    
    # print $introns->[$i] . " ";
    #adjust intron positions for alignment gaps
    for( my $j = 0; $j < $introns->[$i]; $j++ ){
      $sequence->[$j] eq "-" and $introns->[$i]+=1;
    }
    # print $introns->[$i] . "\n";
    
    print OUT "newpath\n";
    print OUT $left_margin + ($SCALE * ($introns->[$i] + 1) + $intron_offset), " ", $y+1, " moveto\n";
    print OUT $left_margin + ($SCALE * ($introns->[$i] + 1) + $intron_offset) - 7, " ", $y+14, " lineto\n";
    print OUT $left_margin + ($SCALE * ($introns->[$i] + 1) + $intron_offset) + 7, " ", $y+14, " lineto\n";
    print OUT $left_margin + ($SCALE * ($introns->[$i] + 1) + $intron_offset), " ", $y+1, " lineto\n";
    print OUT "closepath\n";
    if( $missing ){
      print OUT "stroke\n";
    }else{
      print OUT "gsave\n";
      print OUT "fill\n";
      print OUT "grestore\n";
      print OUT "stroke\n";
    }
    $missing = undef;
  }
  close OUT;
}
########################################################################################################################
sub scaleBar{
  my( $y, $outfile ) = @_;
  
  $y -= 10;

  open( OUT, ">>$outfile" ) || die "Can't open $outfile: $!\n";

  print OUT "0 setgray\n";
  print OUT "1 setlinewidth\n";
  print OUT "newpath\n";
  print OUT "$left_margin $y moveto\n", $left_margin + $SCALEBAR_WIDTH, " $y lineto\n";
  print OUT "stroke\n";

  print OUT "/Arial findfont\n";
  print OUT "10 scalefont\n";
  print OUT "setfont\n";

  print OUT "newpath\n";
  print OUT $left_margin + $SCALEBAR_WIDTH + 6, " ", $y-3, " moveto\n";
  print OUT "($SCALEBAR_LABEL nt) show\n";
  
  close OUT;
  
}
########################################################################################################################
sub writeName{
  my( $left_margin, $y, $name ) = @_;

  $name =~ s/_/ /g;

  open( OUT, ">>$outfile" ) || die "Can't open $outfile: $!\n";

  print OUT "0 setgray\n";

  print OUT "/Arial findfont\n";
  print OUT "10 scalefont\n";
  print OUT "setfont\n";
  print OUT $left_margin, " ", $y-3, " moveto\n";
  print OUT "($name) show\n";
  
  close OUT;
  
}
########################################################################################################################
sub parseArgs{

  my $usage = "\nUsage: $0 [options] cDNA.fasta
   
   options
          --scale - proportion to scale width of gene line (default: none)
          --bar_width - width of scale bar (default: 100 nt)
          --header    - print header at top of page (filename minus file extension) (default: yes)

   NOTE: FASTA header can have list of intron coordinates (e.g., \">Citrullus 500 600 700m\"); introns are mapped as filled triangles and missing introns (e.g. 700m) are mapped as open triangles; it is assumed that intron locations do not consider any gap characters that might exist in the alignment; that is, intron locations are adjusted based on the number of gap characters (\"-\") that precede them


\n\n";
	
  my $result = GetOptions
	(
	 'file=s'       => \$INFILE,
	 'scale=s'      => \$SCALE,
	 'bar_width=s'  => \$SCALEBAR_WIDTH,
         'header!'      => \$HEADER,
	);
  
  $ARGV[0] or die $usage;

}
#######################################################################################################################
