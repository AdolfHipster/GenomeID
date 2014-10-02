#!/usr/bin/perl
use strict;
use warnings;
use MIME::Base64;
use Data::Dumper;
use Switch;
use Math::CDF qw(:all);

# define constants that will be required 
# to calulate probabilities and what not
my @peng_bits = (26, 45, 64, 48, 8, 21, 46, 30, 36, 10, 51, 14, 38, 60, 52, 53, 42, 16, 33, 61, 40, 56, 18);
my $prob_success  = 0.6;

my $id1 = "A3U2KK3WR/qB";
my $id2 = "B2U2Kg1WQkJB";

my @bin_ids = (unpack('B*', decode_base64($id1)) , unpack('B*', decode_base64($id2)) );
check_version();

# determine if 12 or 20 character id given
if(length $id1 == 12 || $id2 == 12){
   my $prob = 0; my $goodData = 1;

   # force to first 72 bits 
   $bin_ids[0] = substr($bin_ids[0],0,72);
   $bin_ids[1] = substr($bin_ids[1],0,72);

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
   
   print "$k \t $n\n";
   print "$prob \t $goodData\n";
}

# 20 character id given
else{

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
