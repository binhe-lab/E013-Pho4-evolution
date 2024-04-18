#!/usr/bin/perl

use warnings;
use strict;

package analyze_agilent;

use lib './';
use Statistics::Regression;

require Exporter;

our @ISA = qw(Exporter);
our @EXPORT = qw(
              median_neighborhood
		      cy3_regression
		      alexa488_normalize
		      alexa635_normalize
		      median
		      convert_to_letters
		     );

#####################################################################
#### For each spot on a 4x44K or 8x15K or 8x60K array, takes the neighborhood
####   of spots as the block of radius R (up, down, left,
####   right by R rows and columns), and computes the median
####   intensity for the block.
#### 4x44k array has 45220 spots: 170 rows x 266 columns.
#### 8x15K array has 16032 spots: 167 rows x 96 columns.
#### 8x60K array has 62976 spots: 328 rows x 192  columns.
#### 4x180K array has 180880 spots: 340 rows x 532  columns.
#### 2x400K array has~420288 spots: 796 rows x 528  columns.
#### Adjusts for rounded corners and margins (11/20/06).  If
####   the spot is too close to the corner, go in further
####   (diagonally) so the block size is constant.  If the
####   spot is too close to the margin, go in further
####   (laterally) so the block size is constant.
####
#### M. Berger
######################################################################


###$final_index = median_neighborhood(\@data_matrix,$radius,"Alexa488",$array_type,$Alexa488_median,$final_index);

sub median_neighborhood {
	my $data_matrix_Aref = shift;	
	my $radius = shift;
	my $fluor = shift;
	my $array_type = shift;
	my $chip_median = shift;
	my $lastcell = shift;
	my $adjbsi;
	my $flags;

	my $nextcell = $lastcell+1;
	
	if ($fluor eq "Alexa532") {$adjbsi=3; $flags=4;}
	elsif ($fluor eq "Alexa488") {$adjbsi=$lastcell; $flags=6;}
	elsif ($fluor eq "Alexa635") {$adjbsi=$lastcell; $flags=6;}
	else {die "Incorrect usage of median_neighborhood_4x44k";}

    my $toprow; my $topcol;
	my $upperleft; my $upperright; my $lowerleft; my $lowerright;
   	if ($array_type eq "4x44k") {
		$toprow = 170;
		$topcol = 266;
		$upperleft = 20;
		#$upperright = 245;
		$upperright = 246;
		$lowerleft = 150;
		$lowerright = 416;
   	}
    elsif ($array_type eq "8x15k") {
		$toprow = 164;
		$topcol = 96;
		$upperleft = 10;
		$upperright = 86;
		$lowerleft = 154;
		$lowerright = 250;
    }
    elsif ($array_type eq "8x60k") {
		$toprow = 328;
		$topcol = 192;
		$upperleft = 20;
		$upperright = 172;
		$lowerleft = 308;
		$lowerright = 500;
    }
    elsif ($array_type eq "4x180k") {
		$toprow = 340;
		$topcol = 532;
		$upperleft = 20;
		$upperright = 512;
		$lowerleft = 320;
		$lowerright = 852;
   }
    elsif ($array_type eq "2x400k") {
		$toprow = 796;
		$topcol = 528;
		$upperleft = 10;
		$upperright = 518;
		$lowerleft = 786;
		$lowerright = 1314;
   }
   else {die "Array type is not 4x44k or 8x15k or 8x60k or 4x180k or 2x400k for median_neighborhood.pm"}

	for (my $row=1; $row<=$toprow; $row++) {
		for (my $col=1; $col<=$topcol; $col++) {
			my $left_col;
			my $right_col;
			my $top_row;
			my $bottom_row;
			my $adj_col;
			my $adj_row;
			my $counter;

			##### Upper Left Corner
			##### Empty spots defined by the line "Row + Column < $upperleft".
			if ( ($col + $row) < ($upperleft + 2*$radius) ) {
				$counter = 0;
				while (($row + $counter + $col + $counter) < ($upperleft + 2*$radius)) {
					$counter++;
				}
				$adj_col = $col + $counter;
				$adj_row = $row + $counter;
			}
	
			##### Upper Right Corner
			##### Empty spots defined by the line "Column - Row > $upperright".
			elsif ( ($col - $row) > ($upperright - 2*$radius) ) {
				$counter = 0;
				while ((($col - $counter) - ($row + $counter)) > ($upperright - 2*$radius)) {
					$counter++;
				}
				$adj_col = $col - $counter;
				$adj_row = $row + $counter;
			}
	
			##### Lower Left Corner
			##### Empty spots defined by the line "Row - Column > $lowerleft".
			elsif ( ($row - $col) > ($lowerleft - 2*$radius) ) {
				$counter = 0;
				while ((($row - $counter) - ($col + $counter)) > ($lowerleft - 2*$radius)) {
					$counter++;
				}
				$adj_row = $row - $counter;
				$adj_col = $col + $counter;
			}
				
			##### Lower Right Corner
			##### Empty spots defined by the line "$Row + Column > $lowerright".
			elsif ( ($row + $col) > ($lowerright - 2*$radius) ) {
				$counter = 0;
				while (($row - $counter + $col - $counter) > ($lowerright - 2*$radius)) {
					$counter++;
				}
				$adj_row = $row - $counter;
				$adj_col = $col - $counter;
			}
	
			##### All Other Spots
			else {
				$adj_col = $col;
				$adj_row = $row;
			}
	
			##### Column Boundaries
			if (($adj_col > $radius) && ($adj_col <= ($topcol - $radius))) {
				$left_col = ($adj_col - $radius);
				$right_col = ($adj_col + $radius);
			}
			elsif ($adj_col <= $radius) {
				$left_col = 1;
				$right_col = 1 + (2 * $radius);
			}
			elsif ($adj_col > ($topcol - $radius)) {
				$left_col = $topcol - (2 * $radius);
				$right_col = $topcol;
			}
	
			##### Row Boundaries
			if (($adj_row > $radius) && ($adj_row <= ($toprow - $radius))) {
				$top_row = ($adj_row - $radius);
				$bottom_row = ($adj_row + $radius);
			}
			elsif ($adj_row <= $radius) {
				$top_row = 1;
				$bottom_row = 1 + (2 * $radius);
			}
			elsif ($adj_row > ($toprow - $radius)) {
				$top_row = $toprow - (2 * $radius);
				$bottom_row = $toprow;
			}
	
			##### Make array of elements in block
	
			my @block_spots = ();
	
			for (my $block_col = $left_col; $block_col <= $right_col; $block_col++) {
				for (my $block_row = $top_row; $block_row <= $bottom_row; $block_row++) {
					if($$data_matrix_Aref[$block_col][$block_row][$adjbsi]ne"NA" && $$data_matrix_Aref[$block_col][$block_row][$flags]>-100 && ($$data_matrix_Aref[$block_col][$block_row][1]=~"dBr" || $$data_matrix_Aref[$block_col][$block_row][1]=~"Ctrl" || $$data_matrix_Aref[$block_col][$block_row][1]=~"Pho4or" || $$data_matrix_Aref[$block_col][$block_row][0]=~"Cbf1Pho4Tye7")) {
						push (@block_spots, $$data_matrix_Aref[$block_col][$block_row][$adjbsi]);
					}
				}
			}
			my $median;

			##### MFB added 2/10/07; if more than half of spots missing from neighborhood, 

			my $blocksize = (2 * $radius + 1) * (2 * $radius + 1);
			
			if ($#block_spots <= $blocksize/2) {
			    $$data_matrix_Aref[$col][$row][$nextcell] = $chip_median;
			}

			else {
			    $median = median (\@block_spots);
			    $$data_matrix_Aref[$col][$row][$nextcell] = $median;
			    #print "blocksize $blocksize; median $median\n";
			}

			if ($$data_matrix_Aref[$col][$row][$flags]>-100 && $fluor eq "Alexa532" && ($$data_matrix_Aref[$col][$row][1]=~"dBr" || $$data_matrix_Aref[$col][$row][1]=~"Ctrl" || $$data_matrix_Aref[$col][$row][1]=~"Pho4or" || $$data_matrix_Aref[$col][$row][0]=~"Cbf1Pho4Tye7")) {
				$$data_matrix_Aref[$col][$row][$nextcell+1] = $$data_matrix_Aref[$col][$row][$adjbsi]/$$data_matrix_Aref[$col][$row][$nextcell]*$chip_median;
			}
			elsif ($$data_matrix_Aref[$col][$row][$flags]>-100 && (($fluor eq "Alexa488")||($fluor eq "Alexa635")) && $$data_matrix_Aref[$col][$row][$lastcell] ne "NA") {
				$$data_matrix_Aref[$col][$row][$nextcell+1] = $$data_matrix_Aref[$col][$row][$adjbsi]/$$data_matrix_Aref[$col][$row][$nextcell]*$chip_median;
			}
			else {$$data_matrix_Aref[$col][$row][$nextcell+1] = "NA";}
		}
	}
	
	$nextcell + 1;
}

#########################################################################
# Computes the regression coefficients for Cy3 using the (unflagged)
#  combinatorial "dBr" spots.  Appends expected Cy3 to data matrix.
#  Calculates observed/expected Cy3 for all custom-designed deBruijn
#  and control spots.
#    -- adapted from cpan.org, A. Philippakis, M. Berger
#########################################################################

sub cy3_regression{
    my $data_matrix_Aref = shift;
    my $regressionorder = shift;
    my $combinatoriallength=$regressionorder+36;  #MFB changed 2/10/07
    my $array_type = shift;  ###MFB changed 7/5/07 (can be 4x44k or 8x15k)
    my $output = shift;
    my $lastcell = shift;

    my $nextcell = $lastcell+1;

    my $toprow; my $topcol;
    if ($array_type eq "4x44k") {
		$toprow = 170;
		$topcol = 266;
    }
    elsif ($array_type eq "8x15k") {
		$toprow = 164;
		$topcol = 96;
    }
    elsif ($array_type eq "8x60k") {
		$toprow = 328;
		$topcol = 192;
    }
    elsif ($array_type eq "4x180k") {
		$toprow = 340;
		$topcol = 532;
    }
    elsif ($array_type eq "2x400k") {
		$toprow = 796;
		$topcol = 528;
    }
    else {die "Array type is not 4x44k or 8x15k or 8x60k or 4x180k or 2x400k for cy3_regression.pm"}

    my $i; my $j; my $k; my $key;
    my @Xcomponentnames=qw(intercept);
    my %Xcomponentvalues;
    my @theta;
    my %regressors;
    my $predictedvalue;
    my $numbercomponents=(4**$regressionorder)+1;
    my $currentletter;
    my $sequence; my $signal;
    my @tmp;

    for($i=0; $i<4**$regressionorder; $i++){
	$currentletter=convert_to_letters($i,$regressionorder);
	$currentletter="A".$currentletter;
	$Xcomponentvalues{$currentletter}=0;
	push @Xcomponentnames,$currentletter;
    }

    my $reg=Statistics::Regression->new($numbercomponents,"components",[@Xcomponentnames]);

    for (my $row=1; $row<=$toprow; $row++) {
	for (my $col=1; $col<=$topcol; $col++) {
	    if ($$data_matrix_Aref[$col][$row][1]=~"dBr" && $$data_matrix_Aref[$col][$row][$lastcell] ne "NA"  && $$data_matrix_Aref[$col][$row][4] > -100) 
	    {		
			foreach $key (sort keys %Xcomponentvalues){
			    $Xcomponentvalues{$key}=0;
			}
			
			$sequence=substr($$data_matrix_Aref[$col][$row][2],0,$combinatoriallength);
			if ($lastcell > 6) {
			    $signal=$$data_matrix_Aref[$col][$row][$lastcell];
			}
			else {
			    $signal=$$data_matrix_Aref[$col][$row][3];
			}

			@tmp=();
			push @tmp,1;
	
			foreach $key (sort keys %Xcomponentvalues){
			    while($sequence=~m/$key/gi){
				$Xcomponentvalues{$key}++;
				pos($sequence)-=$regressionorder;
			    }
			    push @tmp, $Xcomponentvalues{$key};
			}
			$reg->include($signal , [@tmp] );
			for($i=0; $i<=$#tmp; $i++){
			    push @{$regressors{$sequence}},$tmp[$i];
			}
			push @{$regressors{$sequence}},$signal;
	    }
	}
    }

    @theta=$reg->theta;

    open(OUTPUT,">$output") || die "Couldn't open 'regression' output file\n";
    print OUTPUT "*****************************************************\n";
    print OUTPUT "Regression 'components'\n";
    print OUTPUT "*****************************************************\n";
    for (my $x=0; $x<=$#theta; $x++) {
		print OUTPUT "Theta[$x='$Xcomponentnames[$x]']=\t$theta[$x]\n";
    }
    my $rsquared = sprintf("%.3f", $reg->rsq());
    print OUTPUT "R^2= $rsquared\n";
    print OUTPUT "*****************************************************\n";
    close (OUTPUT);

    for (my $row=1; $row<=$toprow; $row++) {
	for (my $col=1; $col<=$topcol; $col++) {
	    if(($$data_matrix_Aref[$col][$row][1]=~"dBr"||$$data_matrix_Aref[$col][$row][1]=~"Ctrl") && $$data_matrix_Aref[$col][$row][$lastcell]ne"NA" && $$data_matrix_Aref[$col][$row][4]>-100) 
	    {
			foreach $key (sort keys %Xcomponentvalues) {
			    $Xcomponentvalues{$key}=0;
			}
			my $array_sequence = substr($$data_matrix_Aref[$col][$row][2], 0, $combinatoriallength);
			@tmp=();
			push @tmp,1;
	
			foreach $key (sort keys %Xcomponentvalues) {
			    while($array_sequence=~m/$key/gi){
					$Xcomponentvalues{$key}++;
					pos($array_sequence)-=$regressionorder;
			    }
			    push @tmp, $Xcomponentvalues{$key};
			}
			$predictedvalue = 0;
			for($i=0; $i<=$#tmp; $i++) {
			    $predictedvalue += $tmp[$i]*$theta[$i];
			}
			$$data_matrix_Aref[$col][$row][$nextcell] = $predictedvalue;
			if ($lastcell > 6) {
			    $$data_matrix_Aref[$col][$row][$nextcell+1] = $$data_matrix_Aref[$col][$row][$lastcell]/$predictedvalue;
			}
			else {$$data_matrix_Aref[$col][$row][$nextcell+1] = $$data_matrix_Aref[$col][$row][3]/$predictedvalue;}
		}
	    else {
			$$data_matrix_Aref[$col][$row][$nextcell] = "NA";
			$$data_matrix_Aref[$col][$row][$nextcell+1] = "NA";
	    }
	}
    }
    $nextcell + 1;
}


#########################################################################
# Normalizes Alexa488 signal by (Observed/Expected) Cy3 and stores
#   array of all Alexa488 intensities for custom-designed 
# deBruijn (combinatorial) and control probes.
#########################################################################

sub alexa488_normalize{
    my $data_matrix_Aref = shift;
    my $alexa488_list = shift;
    my $array_type = shift;  ###MFB changed 7/5/07 (can be 4x44k or 8x15k)
    my $lastcell = shift;
    my $nextcell = $lastcell + 1;

    my $toprow; my $topcol;
    if ($array_type eq "4x44k") {
		$toprow = 170;
		$topcol = 266;
    }
    elsif ($array_type eq "8x15k") {
		$toprow = 164;
		$topcol = 96;
    }
    elsif ($array_type eq "8x60k") {
		$toprow = 328;
		$topcol = 192;
    }
    elsif ($array_type eq "4x180k") {
		$toprow = 340;
		$topcol = 532;
   }
    elsif ($array_type eq "2x400k") {
		$toprow = 796;
		$topcol = 528;
   }
   else {die "Array type is not 4x44k or 8x15k or 8x60k or 4x180k or 2x400k for alexa488_normalize.pm"}

    for (my $row=1; $row<=$toprow; $row++) {
	for (my $col=1; $col<=$topcol; $col++) {
	    if ($$data_matrix_Aref[$col][$row][$lastcell] ne "NA") {
		if ($$data_matrix_Aref[$col][$row][$lastcell] > 0.5 && $$data_matrix_Aref[$col][$row][$lastcell] < 2 && $$data_matrix_Aref[$col][$row][6] > -100) {
		    $$data_matrix_Aref[$col][$row][$nextcell] = $$data_matrix_Aref[$col][$row][5] / $$data_matrix_Aref[$col][$row][$lastcell];
		    push @{$alexa488_list}, $$data_matrix_Aref[$col][$row][$nextcell];
		}
		else {$$data_matrix_Aref[$col][$row][$nextcell]="NA";}
	    }
	    else {$$data_matrix_Aref[$col][$row][$nextcell]="NA";}
	}
    }
    $nextcell;
}

#########################################################################
# Returns the median of an array of numbers -- written by A. Philippakis
#########################################################################

sub median{
	my $Aref = shift; #ref to input array;
	my $median;
	my $i; my $j; my $center;
	my @tmpA; #will store the sorted $Aref;

	my $oddeven = (($#{$Aref}+1)%2);
	@tmpA = sort {$a <=> $b} @{$Aref};

	if($oddeven){
		$center = $#tmpA/2;
		$median = $tmpA[$center];
	}
	else{
		$center = int($#tmpA/2);
		$median = ($tmpA[$center]+$tmpA[$center+1])/2;
	}
	return($median);
}


#######################################################################################
# Takes a number and a "k" and returns the ACGT-equivalent of that number for length k
#    -- written by A. Philippakis
#######################################################################################

sub convert_to_letters {
    my $number = shift;
    my $k = shift;
    my @letters = ("A", "C", "G", "T");
    my $string = "";
    my $i;
    for ($i = ($k-1)*2; $i >= 0; $i-=2) {
	$string .= $letters[($number >> $i) & 0x3];
    }
    return $string;
}
