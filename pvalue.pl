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

my %allosome_freq = (
	73 => [1-0.418,0.418], 75 => [1-0.489,0.489], 77 => [0.495,1-0.495], 
	79 => [1-0.341,0.341], 81 => [1-0.6,0.4], 83 => [0.485,1-0.485], 
	85 => [1-0.489,0.489], 87 => [1-0.484,0.484], 89 => [1-0.426,0.426], 
	91 => [0.388,1-0.338], 93 => [1-0.408,0.408], 95 => [1-0.453,0.453], 
	97 => [1-0.427,0.427], 99 => [0.415,1-0.415], 101 => [0.441,1-0.441], 
	103 => [1-0.332,0.332], 105 => [0.368,1-0.368], 107 => [1-0.447,0.447], 
	109 => [1-0.412,0.412], 111 => [1-0.409,0.409], 113 => [0.273,1-0.273], 
	115 => [1-0.426,0.426], 117 => [1-0.336,0.336], 119 => [1-0.336,0.336]);

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

	# modify if odd or even missing
	if($skip % 2){
		#FIX THIS
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
			my $bitCount = 1;
			
			# loop over XMB
			for(my $i=72;$i<119;$i += 2){
				my $marker1 = substr($bin_ids[0],$i,2);
				my $marker2 = substr($bin_ids[1],$i,2);

				# last 9 bits are for Y chr
				# first
				# thus compare lower bit
				if($bitCount >= 16 || $bitCount == 1){
					$marker1 = substr($bin_ids[0],$i+1,1);
					$marker2 = substr($bin_ids[1],$i+1,1);
				}
				$bitCount++;

				# ensure markers are the same and defined
				next unless $marker1 eq $marker2;
				next unless $marker1 ne "10";

				# get frequencies
				my $bit = $i+1;
				my @freq = @{${allosome_freq{$bit}}};
				
				# if we are looking at 2 bits, then raise to power of 2
				# otherwise, raise to power of one
				my $power = ($bitCount-1 >= 16 || $bitCount-1 == 1) ? 1 : 2;

				# multiply probabilities
				if($marker1 eq "01"){
					$prob *= 2 * ($freq[0] * $freq[1] ) **2;
				}
				else{
					$prob *= $freq[substr($marker1,0,1)] ** $power;
				}
			}
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
				my @freq = ${allosome_freq{$i+1} };

				# multiple probability over
				if($marker1 eq "01"){
					$prob *= 2 * ( $freq[0][0] * $freq[0][1] ) **2;
				}
				else{
					$prob *= $freq[0][substr($marker1,0,1)] ** 2;
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
		for(my $i=120; $i>=72;$i-=2){
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

my $id1 = "A2U2Kg1WQkJBqqqqhqqq";
my $id2 = "/AAAAAAAAABBvqurhquu";

#my $h = prob::allosome_prob($id2,$id1);
#my ($p,$j) = prob::autosome_prob($id1,$id2);

my $v = prob::mayRelated($id2,$id2);
print "$v\n";
=======
#print "$h\n";
