#!/usr/bin/perl
use strict;
use warnings;

package prob;
use MIME::Base64;
use Data::Dumper;
use Switch;
use Math::CDF qw(:all);

# define constants that will be required 
# to calulate probabilities and what not
my @peng_bits = (26, 45, 64, 48, 8, 21, 46, 30, 36, 10, 51, 14, 38, 60, 52, 53, 42, 16, 33, 61, 40, 56, 18);
my $prob_success  = 0.6;
my @bin_ids;

=comment
%alosome_freq = ('bit' => [anc freq, atl freq])
=cut

sub autosome_prob{
	@bin_ids = (unpack('B*', decode_base64($_[0])) , unpack('B*', decode_base64($_[1])) );
	check_version();
	my $prob = 0; my $goodData = 1;

	# determine how many to skip over
	my @fBits = (bin2dec(substr($bin_ids[0],0,6)) ,
			bin2dec(substr($bin_ids[1],0,6)) );

	# variables to pass to calculate matches
	my $more_missing = ($fBits[1] > $fBits[0])? 1 : 0; 
	my @bits;   # array that stores the bits that should be checked
	my $skip=0; # the number of bits to skip

	# determine what case we have here
	switch ($fBits[$more_missing]){
		case 0 {
			@bits = 7..66;
		}
		case 59 {
			@bits = @peng_bits;
			push(@bits,65); 
			push(@bits,66);
		}
		case 63 {
			@bits = 7..66;
			$goodData = 0;
		}
		case [60..62]{
			@bits = 7..66;
			$skip = $fBits[$more_missing] - 58;
		}
		else {
			@bits = 7..66;
			@bits = grep {$_ != $fBits[$more_missing]+6 } @bits;
		}
	}

	my $k = cAutosome_match(\@bits,$skip,$more_missing);
	my $n = scalar @bits - $skip;

# FIX THIS
	# modify if odd or even missing
	if($skip % 2){

	}
	# what if only one bit was missed?
	elsif (scalar @bits == 58){
		$n--;
	}
	# even number of missing bits
	else{
		$k -= $skip/2;
		$prob = 1 - pbinom($k,$n,$prob_success);
	}

	return ($prob,$goodData);
}

# calculate probability for sex marker bits
sub allosome_prob{
	if(length $_[0] != 20 || length $_[1] != 20){
		die "Id for allosome probability must be 20 characters\n";
	}

	@bin_ids = (unpack('B*', decode_base64($_[0])) , unpack('B*', decode_base64($_[1])) );
	check_version();

	my $prob = 1;
	
	# determine if same gender
	if(substr($bin_ids[0],64,1) eq substr($bin_ids[1],64,1) ){
		# determine if Y chromosome present
		if(substr($bin_ids[0],64,1)){

		}

		# no Y chromosome present
		else{
			# loop over XMB
			for(my $i=72;$i<119;$i +=2){
				my $marker1 = substr($bin_ids[0],$i,2);
				my $marker2 = substr($bin_ids[1],$i,2);
				
				# ensure markers are same and defined
				next unless $marker1 eq $marker2;
				next unless $marker1 ne "10";
				
				# get frequencies
				my @freq = (0.4,0.6); #@{allosome_freq{$i+1} };

				# multiple probability over
				if($marker1 eq "01"){
					$prob *= 2 * ( $freq[0] * $freq[1] ) **2;
				}
				else{
					$prob *= $freq[substr($marker1,0,1)] ** 2;
				}

			}
		}
	}

	# not same gender, more difficult
	else{

	}

	return $prob;
}

sub mayRelated{
	if(length $_[0] != 20 || length $_[1] != 20){
		die "Id must be 20 characters\n";
	}
	@bin_ids = (unpack('B*', decode_base64($_[0])) , unpack('B*', decode_base64($_[1])) );
	check_version();
	my $bitcount = 1;

	# determine if not same gender
	if( substr($bin_ids[0],64,1) ne substr($bin_ids[1],64,1) ){
		# loop over right handed bits, ensure same marker
		for(my $i=73; $i<120;$i+=2){
			last unless $bitcount != 9;

			$bitcount++;
			next if(substr($bin_ids[0],$i,1) eq substr($bin_ids[1],$i,1) );

			return 0;
		}
	}

	return 1;
}


# subroutine to loop through bits and compare matches
sub cAutosome_match{
	# store input variables 
	my @bits = @{$_[0]};
	my $skip = $_[1];
	my $moreMiss = $_[2];

	my $matches = 0;

	foreach my $bit (@bits){
		$bit = $bit-1;

		if($skip != 0 && !substr($bin_ids[$moreMiss],$bit,1)){
			$skip--;
			next;
		}

		$matches += (substr($bin_ids[0],$bit,1) eq substr($bin_ids[1],$bit,1));
	}

	return $matches;
}


# check version block to ensure they are compatible
sub check_version{
	my $version_id1 = substr($bin_ids[0],66,6);
	my $version_id2 = substr($bin_ids[0],66,6);

	if($version_id1 ne $version_id2){
		die "Genome IDs are incompatible versions\n";
	}
}

#convert binary to decimal
sub bin2dec {
	unpack("N", pack("B32", substr("0" x 32 . shift, -32)));
}


package main;

my $id1 = "A2U2Kg1WQkJBqqqqqqqq";
my $id2 = "/AAAAAAAAABBvqurhquu";

my $h = prob::allosome_prob($id2,$id2);
my ($p,$j) = prob::autosome_prob("A2U2Kg1WQkJB","A3U2KK3WR/qB");

my $v = prob::mayRelated($id1,$id2);
print "$v\n";
