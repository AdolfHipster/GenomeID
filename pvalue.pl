#!/usr/bin/perl
use strict;
use warnings;
use MIME::Base64;
use Data::Dumper;
use Switch;

my $id1 = "A3U2KK3WR/qB";
my $id2 = "A2U2Kg1WQkJB";

my @bin_ids = (unpack('B*', decode_base64($id1)) , unpack('B*', decode_base64($id2)) );

# determine if 12 or 20 character id given
if(length $id1 == 12 || $id2 == 12){
   # force to first 72 bits 
   $bin_ids[0] = substr($bin_ids[0],0,72);
   $bin_ids[1] = substr($bin_ids[1],0,72);

   # determine how many to skip over
   my @fBits = (bin2dec(substr($bin_ids[0],0,6)) ,
	 bin2dec(substr($bin_ids[1],0,6)) );

   my $more_missing = ($fBits[1] > $fBits[0])? 1 : 0; 
   my $k = 0;

   # determine what case we have here
   switch ($fBits[$more_missing]){
      case 0 {
	 print "all here\n";
	 # determine number of shared elements
	 for(my $i =5; $i <72;$i++){
	    $k += (substr($bin_ids[0],$i,1) eq substr($bin_ids[1],$i,1) );
	 }
      }
      case 59 {
	 print "only peng\n";
      }
      case [1..58]{
	 print "one marker missing\n";
      }
      case [60..62]{
	 print "less than 4 markers missing\n";
      }
      else {
	 print "more than 4 missing\n";
      }
   }

   print "$k\n";
   print "$bin_ids[0] \n$bin_ids[1]\n";
}

# 20 character id given
else{

}

#convert binary to decimal
sub bin2dec {
   unpack("N", pack("B32", substr("0" x 32 . shift, -32)));
}
