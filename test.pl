 #!/usr/bin/perl
use strict;
use warnings;
use Bio::DB::Sam;
use MIME::Base64;

package genomeID;
use MIME::Base64; 
use Data::Dumper;
use Tie::IxHash;


# define input variables
my $VERSION = 1;
my $file; my $type; my $noise; my $sex;
my $prefix; my $ucn; my $ref; my $baq;
my $guess_hg;

# define static constants
my %bit_loc_static; my %pengelly_static; my %autoSomes_static;
my %allosomes_static; my %ref_allel_static;

my %bit_loc; my %pengelly; my %autoSomes; my %allosomes;
my %ref_allel;

# define bit holders
my $MADIB; my @AMB; my @SMB; 
my $GVB; my @XMB;

# some precompiled regexes
my $homo1 = qr/^(1[\/|]1)[^\/]?/;
my $homo2 = qr/^0[\/|]0[^\/]?/;
my $misR = qr/^(\.[\/|]\.)/;

sub regen{
# define marker relative bit locations
%bit_loc_static = ( 'hg19' => ['1:45973928'=> 1, '1:67861520'=> 2, '1:158582646'=> 3, '1:179520506'=> 4,
					'1:209968684'=> 5, '1:228431095'=> 6, '2:75115108'=> 7, '2:169789016'=> 8,
					'2:215820013'=> 9, '2:227896976'=> 10, '2:49381585'=> 11, '3:4403767'=> 12,
					'3:45989044'=> 13, '3:148727133'=> 14, '4:5749904'=> 15, '4:86915848'=> 16,
					'5:13719022'=> 17, '5:40981689'=> 18, '5:55155402'=> 19, '5:82834630'=> 20,
					'5:138456815'=> 21, '5:171849471'=> 22, '7:34009946'=> 23, '7:48450157'=> 24,
					'7:100804140'=> 25, '7:151254175'=> 26, '8:94935937'=> 27, '9:27202870'=> 28,
					'9:77415284'=> 29, '9:100190780'=> 30, '10:69926097'=> 31, '10:100219314'=> 32,
					'11:6629665'=> 33, '11:16133413'=> 34, '11:30255185'=> 35, '12:993930'=> 36,
					'12:52200742'=> 37, '13:25466955'=> 38, '13:39433606'=> 39, '14:50769717'=> 40,
					'14:76045858'=> 41, '16:70303580'=> 42, '17:7192091'=> 43, '17:42449789'=> 44,
					'17:71197748'=> 45, '18:21413869'=> 46, '19:10267077'=> 47, '19:33353464'=> 48,
					'19:55441902'=> 49, '20:6100088'=> 50, '20:19970705'=> 51, '20:35865054'=> 52,
					'20:52786219'=> 53, '21:44323590'=> 54, '22:21141300'=> 55, '22:37469591'=> 56,
					'6:56471402'=> 57, '6:146755140'=> 58,
					'X:4066743'=>1,	 'X:5616964'=>13,  'X:11882557'=>25, 'X:18151389'=>37,
			   		'X:25779585'=>3, 'X:34091679'=>15, 'X:35632777'=>27, 'X:40871567'=>39,
			   		'X:42848162'=>5, 'X:94465580'=>17, 'X:116521416'=>29,'X:139920048'=>41,
			   		'X:47686005'=>7, 'X:95114611'=>19, 'X:119826608'=>31,'X:144091597'=>43,
			   		'X:67200648'=>9, 'X:104268001'=>21,'X:120864696'=>33,'X:145069663'=>45,
			   		'X:78397143'=>11,'X:106598707'=>23,'X:126325138'=>35,'X:151150133'=>47]);
 
# declar pengelly marker locations
%pengelly_static = ( 'hg19' => ['1:45973928'=> 0, '1:67861520'=> 1, '1:158582646'=> 0, '1:179520506'=> 1,
					'1:209968684'=> 0, '1:228431095'=> 0, '2:75115108'=> 0, '2:169789016'=> 1,
					'2:215820013'=> 0, '2:227896976'=> 1, '2:49381585'=> 0, '3:4403767'=> 1,
					'3:45989044'=> 0, '3:148727133'=> 0, '4:5749904'=> 1, '4:86915848'=> 0,
					'5:13719022'=> 0, '5:40981689'=> 0, '5:55155402'=> 0, '5:82834630'=> 1,
					'5:138456815'=> 0, '5:171849471'=> 0, '7:34009946'=> 0, '7:48450157'=> 1,
					'7:100804140'=> 0, '7:151254175'=> 0, '8:94935937'=> 1, '9:27202870'=> 0,
					'9:77415284'=> 0, '9:100190780'=> 1, '10:69926097'=> 0, '10:100219314'=> 1,
					'11:6629665'=> 0, '11:16133413'=> 1, '11:30255185'=> 0, '12:993930'=> 1,
					'12:52200742'=> 0, '13:25466955'=> 0, '13:39433606'=> 1, '14:50769717'=> 1,
					'14:76045858'=> 0, '16:70303580'=> 1, '17:7192091'=> 0, '17:42449789'=> 0,
					'17:71197748'=> 1, '18:21413869'=> 1, '19:10267077'=> 1, '19:33353464'=> 0,
					'19:55441902'=> 0, '20:6100088'=> 1, '20:19970705'=> 0, '20:35865054'=> 0,
					'20:52786219'=> 0, '21:44323590'=> 1, '22:21141300'=> 1, '22:37469591'=> 0,
					'6:56471402'=> 0, '6:146755140'=> 1]);

# declare autosomal markers
%autoSomes_static = ('hg19' => ['1:45973928'=> 'A:G', '1:67861520'=> 'A:C', '1:158582646'=> 'C:T', '1:179520506'=> 'A:G',
					'1:209968684'=> 'A:C', '1:228431095'=> 'C:T', '2:75115108'=> 'A:G', '2:169789016'=> 'C:T',
					'2:215820013'=> 'A:G', '2:227896976'=> 'C:T', '2:49381585'=> 'A:G', '3:4403767'=> 'C:T',
					'3:45989044'=> 'G:T', '3:148727133'=> 'A:G', '4:5749904'=> 'T:C', '4:86915848'=> 'C:T',
					'5:13719022'=> 'A:C', '5:40981689'=> 'A:C', '5:55155402'=> 'C:T', '5:82834630'=> 'C:T',
					'5:138456815'=> 'C:T', '5:171849471'=> 'A:G', '7:34009946'=> 'C:T', '7:48450157'=> 'C:T',
					'7:100804140'=> 'C:T', '7:151254175'=> 'C:T', '8:94935937'=> 'C:T', '9:27202870'=> 'C:T',
					'9:77415284'=> 'A:C', '9:100190780'=> 'C:T', '10:69926097'=> 'A:G', '10:100219314'=> 'A:G',
					'11:6629665'=> 'C:T', '11:16133413'=> 'A:G', '11:30255185'=> 'C:T', '12:993930'=> 'C:T',
					'12:52200742'=> 'A:C', '13:25466955'=> 'C:T', '13:39433606'=> 'A:G', '14:50769717'=> 'A:G',
					'14:76045858'=> 'A:G', '16:70303580'=> 'C:T', '17:7192091'=> 'C:T', '17:42449789'=> 'C:T',
					'17:71197748'=> 'A:G', '18:21413869'=> 'C:T', '19:10267077'=> 'A:G', '19:33353464'=> 'A:G',
					'19:55441902'=> 'C:T', '20:6100088'=> 'C:T', '20:19970705'=> 'A:G', '20:35865054'=> 'C:T',
					'20:52786219'=> 'A:G', '21:44323590'=> 'G:T', '22:21141300'=> 'C:T', '22:37469591'=> 'A:G',
					'6:56471402'=> 'A:G', '6:146755140'=> 'A:G']);

# declare sex chromosome markers
%allosomes_static = ( 'hg19' => ['X:4066743'=>'G:A','X:5616964'=>'G:A','X:11882557'=>'T:C','X:18151389'=>'G:A',
			   'X:25779585'=>'T:C','X:34091679'=>'T:C','X:35632777'=>'C:T','X:40871567'=>'G:A',
			   'X:42848162'=>'C:T','X:94465580'=>'T:G','X:116521416'=>'G:T','X:139920048'=>'T:C',
			   'X:47686005'=>'C:T','X:95114611'=>'T:C','X:119826608'=>'T:G','X:144091597'=>'G:C',
			   'X:67200648'=>'T:C','X:104268001'=>'T:C','X:120864696'=>'C:A','X:145069663'=>'G:A',
			   'X:78397143'=>'T:C','X:106598707'=>'A:G','X:126325138'=>'G:T','X:151150133'=>'G:T'],

            'hg38' => ['X:4148702'=>'G:A','X:5698923'=>'G:A','X:11864438'=>'T:C','X:18133269'=>'G:A',
            'X:25761468'=>'T:C','X:34073562'=>'T:C','X:35614660'=>'C:T','X:41012314'=>'G:A',
            'X:42988913'=>'C:T','X:95210581'=>'T:G','X:117387453'=>'G:T','X:140837883'=>'T:C',
            'X:47826606'=>'C:T','X:95859612'=>'T:C','X:120692753'=>'T:G','X:145010077'=>'G:C',
            'X:67980806'=>'T:C','X:105023319'=>'T:C','X:121730843'=>'C:A','X:145988145'=>'G:A',
            'X:79141646'=>'T:C','X:107355477'=>'A:G','X:127191155'=>'G:T','X:151981661'=>'G:T']);


# declare reference sex chromosome markers
%ref_allel_static = ( 
   'hg19' => ['X:4066743'=>'G','X:5616964'=>'A','X:11882557'=>'C','X:18151389'=>'A',
			   'X:25779585'=>'T','X:34091679'=>'T','X:35632777'=>'C','X:40871567'=>'G',
			   'X:42848162'=>'C','X:94465580'=>'T','X:116521416'=>'G','X:139920048'=>'C',
			   'X:47686005'=>'T','X:95114611'=>'C','X:119826608'=>'G','X:144091597'=>'G',
			   'X:67200648'=>'T','X:104268001'=>'C','X:120864696'=>'A','X:145069663'=>'G',
			   'X:78397143'=>'C','X:106598707'=>'A','X:126325138'=>'G','X:151150133'=>'G',

            # hg19 ref allele at the hg38 positions
            'X:4148702'=>'G','X:5698923'=>'A','X:11864438'=>'A','X:18133269'=>'C',
            'X:25761468'=>'A','X:34073562'=>'A','X:35614660'=>'G','X:41012314'=>'G',
            'X:42988913'=>'A','X:95210581'=>'C','X:117387453'=>'C','X:140837883'=>'G',
            'X:47826606'=>'','X:95859612'=>'C','X:120692753'=>'G','X:145010077'=>'T',
            'X:67980806'=>'T','X:105023319'=>'C','X:121730843'=>'T','X:145988145'=>'G',
            'X:79141646'=>'C','X:107355477'=>'T','X:127191155'=>'T','X:151981661'=>'G'],

   'hg38' => ['X:4148702'=>'G','X:5698923'=>'A','X:11864438'=>'C','X:18133269'=>'A',
            'X:25761468'=>'T','X:34073562'=>'T','X:35614660'=>'C','X:41012314'=>'G',
            'X:42988913'=>'C','X:95210581'=>'T','X:117387453'=>'G','X:140837883'=>'C',
            'X:47826606'=>'T','X:95859612'=>'C','X:120692753'=>'G','X:145010077'=>'G',
            'X:67980806'=>'T','X:105023319'=>'C','X:121730843'=>'A','X:145988145'=>'G',
            'X:79141646'=>'C','X:107355477'=>'A','X:127191155'=>'G','X:151981661'=>'G',

            # hg38 ref allele at the hg19 positions
            'X:4066743'=>'C','X:5616964'=>'A','X:11882557'=>'A','X:18151389'=>'G',
            'X:25779585'=>'T','X:34091679'=>'A','X:35632777'=>'T','X:40871567'=>'T',
            'X:42848162'=>'G','X:94465580'=>'A','X:116521416'=>'T','X:139920048'=>'C',
            'X:47686005'=>'C','X:95114611'=>'T','X:119826608'=>'G','X:144091597'=>'A',
            'X:67200648'=>'C','X:104268001'=>'A','X:120864696'=>'G','X:145069663'=>'A',
            'X:78397143'=>'A','X:106598707'=>'T','X:126325138'=>'A','X:151150133'=>'T']);

$prefix = "";
$guess_hg  =0;

# define bit holders
$MADIB = ""; @AMB = ("0") x 58; @SMB = ("0") x 2;
$GVB = "";  @XMB = ("0") x 48;
}


sub generate_id{
   # requires input as hash
   my (%input) = @_;

   regen();

   # determine hg version
   if(exists $input{'hg'}){
      
      # determine if hg version is supported
      if(! exists $bit_loc_static{$input{'hg'}}){
	 die "Version: $input{'hg'} is not supported\n";
      }

      (%pengelly) = @{$pengelly_static{$input{'hg'}}}; 
      (%bit_loc) = @{$bit_loc_static{$input{'hg'}}};
      (%autoSomes) = @{$autoSomes_static{$input{'hg'}}};
      (%allosomes) = @{$allosomes_static{$input{'hg'}}};
      (%ref_allel) = @{$ref_allel_static{$input{'hg'}}};
   }

   # otherwise guess it by looking at marker alignment
   else{
      $guess_hg = 1;

      # order the hash locations
      tie my %autoSomes, 'Tie::IxHash';
      tie my %allosomes, 'Tie::IxHash';

      # start with a guess of hg19
      my $hg = 'hg19';
      (%pengelly) = @{$pengelly_static{ $hg }};
      (%bit_loc) = @{$bit_loc_static{ $hg }};
      (%autoSomes) = @{$autoSomes_static{ $hg }};
      (%allosomes) = @{$allosomes_static{ $hg }};
      (%ref_allel) = @{$ref_allel_static{ $hg }};
   }

   # store input flags
   $baq = (exists $input{'baq'})? $input{'baq'} : 30; 
   $noise = (exists $input{'noise'})? $input{'noise'} : 0.05;
   $sex = $input{'sex'}; $file = $input{'file'}; 
   $ucn = (exists $input{'ucn'})? $input{'ucn'} : 0;
   $ref = (exists $input{'ref'})? $input{'ref'} : 0;

   # for the XMB, we assume everything is
   # undefined, later we will fix this 
   # if necessary
   for(my $i=0;$i< scalar @XMB; $i +=2){
      $XMB[$i] = 1
   }
	 	
   # determine if BAM or VCF
   if ($input{'type'} eq 'bam'){
      # open file using sam tools
      $file = Bio::DB::Sam->new(-bam =>$file);
	       
      # determine if prefix is required
      my @chrs = $file->seq_ids;
      if($chrs[0] =~ /^([a-z]+)\d+/){
         $prefix = $1;
      }
      bam();
   }
   elsif($input{'type'} eq 'vcf'){
      vcf();
   }
   elsif($input{'type'} eq 'tbi'){
      tbi();
   }
   else{
      die "No supported file type was specified\n";
   }

   # generate version block
   $GVB = 0 x (6-length(sprintf("%b",$VERSION))) . sprintf("%b",$VERSION);
      
   # generate SIB
   # if ucn or non informative
   # sex marker bit, then this is 1
   $SMB[1] = ($ucn || !$SMB[0] )? 1: "0";
      	
   my $XMB = ($sex)? join("",@XMB) : "";
   my $AMB = join("",@AMB);
   my $SMB = join("",@SMB);

   my $genString = "$MADIB$AMB$SMB$GVB$XMB";
   print "$genString\n";
    
   # print base 64 code
   #return encode_base64 pack 'B*', $genString;
}

=comment

Tabix Parser

=cut
sub tbi{
   my $missing_markers =0; my $only_peng = 1; my $bit_mis = 0;
   
   $prefix = "chr";

   # loop through autoSomal markers
   foreach my $key(keys %autoSomes){
      # grab marker location to query vcf file 
      my ($chr, $pos) = split(":",$key);
      my $isPengelly = $pengelly{$key};
      my $bit  = $bit_loc{$key} -1;
      my $mis = 0; my $zyg = 0;

      # query tabix file
      # need to figure out optimal way to determine
      # if prefix is used
      my $data = `tabix $file -b 2 -e 2 $chr:$pos-$pos`;
      if($data eq ""){
         $data = `tabix $file -b 2 -e 2 $prefix$chr:$pos-$pos`;
      }

      my @rows = split(/\n/,$data);
      my @col = split(/\t/,$data);
     
      # sometimes query may give more 
      # than one row (i dunno why)
      if(scalar @rows > 1){
         foreach my $row (@rows){
            @col = split(/\t/,$row);
            if ($col[1] eq $pos){
               $data = $row;
               last;
            }
         }
      }

      # determine zygosity
      if(! $data eq ""){
         $zyg = ($col[9] =~ $homo1 or $col[9] =~ $homo2)? 0:1;
         $mis = ($col[9] =~ $misR)? 1:0;
      }

      # no read was found
      if ($data eq "" || $mis){
         $AMB[$bit] = 0;
         $missing_markers++;
         $bit_mis = $bit;

         if($isPengelly){$only_peng =0;}
      }
      
      # read was found
      else{
         # append amb
         $AMB[$bit] = $zyg;
	 
         # if this is a non pengelly
         # then not only pengelly's are found
         if(!$isPengelly){$only_peng =0;}
      }
   }

   # generate MADIB
   genMADIB('miss_count'=>$missing_markers,'oPeng'=>$only_peng, 'bMis'=>$bit_mis);
   
   if(!$sex){
      return;
   }

   # loop through allosomal markers
   foreach my $key (keys %allosomes){
      # grab marker location to query vcf file
      my ($chr,$pos) = split(':', $key);
      my ($anc,$alt) = split(":", $allosomes{$key});
      my $bit = $bit_loc{$key} - 1;

      # query tabix file
      my $data = `tabix $file -b 2 -e 2 $chr:$pos-$pos`;
      if($data eq ""){
         $data = `tabix $file -b 2 -e 2 $prefix$chr:$pos-$pos`;
      }
      my @rows = split(/\n/,$data);

      # determine if there was results,
      # if not, assign to reference if flag
      # was specified
      if($data eq ""){
         if($ref){
            my $ref_al = $ref_allel{$key};

            # homozygous ancestral
            if($ref_al eq $anc && !$SMB[0]){
               $XMB[$bit] = 0;
            }

            # homozygous alternate
            else{
               $XMB[$bit+1] = 1;
            }
         }
         next;
      }

      # ensure that only one row was returned in from
      # the tabix query
      if(scalar @rows > 1){
         foreach my $row (@rows){
            my @col = split(/\t/,$row);
            
            if($col[1] eq $pos){
               $data = $row;
               last;
            }
         }
      }

      # build the xmb bits
      xmb_builder( 'key' => $key, 'data' => $data );
   }
}

=comment

BAM Parser

=cut
sub bam{
   my $missing_markers = 0; my $only_peng = 1; my $bit_mis =0;

   # loop through autoSomal markers
   foreach my $key (keys %autoSomes){  
      # grab marker location to query bam file
      my ($chr,$pos) = split(":",$key);
      my ($fir,$las) = split(":",$autoSomes{$key});
      my $isPengelly = $pengelly{$key};	
      my $bit = $bit_loc{$key} -1;
	        
      # get zygosity 
      my $zyg = bam_zygosity('chr'=>"$prefix$chr",'loc'=>$pos,'al1'=>$fir,'al2'=>$las);
		
      # if no reads
      if($zyg == -1){
         # append to AMB
         $AMB[$bit] = "0"; 
         $missing_markers++;
         $bit_mis = $bit;

         # if this is a missing pengelly, 
         # then not only pengelly's are read
         if($isPengelly){ $only_peng =0;}
      }
      
      # otherwise a read was found
      else{
         $AMB[$bit] = $zyg;

         # if this is a read non pengelly,
         # then not only pengelly's are read
         if(!$isPengelly){$only_peng = 0;}
      } 
   }
   genMADIB('miss_count'=>$missing_markers,'oPeng'=>$only_peng,'bMis'=>$bit_mis);

   $SMB[0] = (grep $_ eq  "$prefix"."Y", $file->seq_ids);
   
   # generate additional sex markers
   if($sex){
         # loop through allosomes
         foreach my $key (keys %allosomes){
         # grab marker information used to query bam file
         my ($chr,$pos) = split(":",$key);
         my ($anc,$alt) = split(":",$allosomes{$key});
         my $bit = $bit_loc{$key} -1;

         # determine zygosity
         my $zyg = bam_zygosity('chr'=>"$prefix$chr",'loc'=>$pos,'anc'=>$anc,'alt'=>$alt,'smb'=>$SMB[0]);
      
         # modify both bits
         if($chr eq 'X'){
            my ($bit1,$bit2) = split(":",$zyg);
            $XMB[$bit] = $bit1;
            $XMB[$bit+1] = $bit2;
         }
         
         # Y chromosome modify left bit
         else{

         }

      }
   }
}


=comment
   
VCF Parser

=cut
sub vcf{
   # define useful constants
   my $missing_markers = 0;
   my $only_peng  = 1; my $bit_mis = 0;
   my $found = 0; my $found_peng = 0; 
	
   open my $vcf, '<', $file or die "$!\n";
   while(<$vcf>){
      # make sure its a chromosome
      next if /^#/;
		
      # get required information
      my @fields = split(/\t/,$_);
      my $chr = $fields[0]; $chr =~ s/^chr//;
      my $loc = $fields[1];
      my $key = "$chr:$loc";

      # determine if this is marker
      next unless exists $bit_loc{$key};

      my $bit = $bit_loc{$key}-1;
      my $mis = ($fields[9] =~ $misR)? 1:0;
	        
      # determine if sex chromosome marker
      if(exists $allosomes{$key} && $sex){
         delete $ref_allel{$key};

         xmb_builder( 'key' => $key,'data' => $_ );

         next;
      }

      # if its not a sex chromosome marker, then
      # it must be an autosome
      next unless exists $autoSomes{$key};

      my $zyg = ($fields[9] =~ $homo1 or $fields[9] =~ $homo2)? 0:1;

      #fill in the correct amb bit
      $AMB[$bit] = $zyg * (!$mis);
      my $isPengelly = $pengelly{$key};
      $found_peng += $isPengelly;	

      # determine if missed read
      if($mis){
         $missing_markers++;
         $bit_mis = $bit;

         if($isPengelly){$only_peng = 0;}

      }elsif(!$isPengelly){$only_peng = 0;}
		
      $found += (!$mis);
   }

   # if there was a flag to match
   # sex information to the reference genome, then
   # loop through and do so
   if($ref){
      foreach my $key (keys %ref_allel){
         my $bit = $bit_loc{$key} - 1;
         my ($anc,$alt) = split(":",$allosomes{$key});
   	    
         # determine if ref is ancestral or alternate
         # and update the bit accordingly
   	    
         # homozygous ancestral
         # check if y is present, if so
         # no modification to the array (already set to 0)
         if($ref_allel{$key} eq $anc && !$SMB[0]){
            $XMB[$bit] = 0;
         }

         # homozygous alternate
         else{
            $XMB[$bit+1] = 1;
         }
      }
   }

   # deal with missing markers
   $missing_markers += 58 - $found;
   $only_peng = ($only_peng && $found_peng == 24);
	
   genMADIB('miss_count'=>$missing_markers,'oPeng'=>$only_peng,'bMis'=>$bit_mis);	
}


=comment

   UTILITY/HELPER SUBROUTINES 
      
      - not directly involved in the process
      - but do help do stuff

   genMADIB
      generates the MADIB for all forms of input files

      @param miss_count    number of missing bits
      @param oPeng         1|0 only the Pengelly markers read
      @param bMis          the last location of a missed read
      @return              void
   
   xmb_builder
      For vcf and tabix files, this will build the 
      xmb marker bits

      @param key        a string of the form <chr>:<position>
      @param data       VCF row for a specific location
      @return           void

   bam_zygosity
      For bam files, this will automate the quering of 
      the bam file and determine if the specific diploid
      is homozygous or not

      @param chr        the chromosome location
      @param loc        the location on the chromosome
      @param smb        1|0, determine sex or not
      @param anc        ancestral allele
      @param alt        alternate allele
      @return           returns the zygosity of the snp
=cut

sub genMADIB{
   my (%input) = @_;
   
   # if no missing markers, MADIB = 0
   if($input{'miss_count'} == 0){
      $MADIB = "000000";
   }
   # if greater than 4, set to 63
   elsif($input{'miss_count'} > 4){
      $MADIB = sprintf("%b",63);
   }
   # only Pengelly, set to 59
   elsif($input{'oPeng'}){
      $MADIB = sprintf("%b",59);
   }
   # only one missing marker
   elsif($input{'miss_count'} == 1){
      $MADIB = "0" x (6-length(sprintf("%b",$input{'bMis'}))) . sprintf("%b",$input{'bMis'});
   }
   # else just show the number of missing markers
   else{
      $MADIB = "0" x (6-length(sprintf("%b",$input{'miss_count'}))) . sprintf("%b",$input{'miss_count'});
   }
}

sub xmb_builder{
   my (%input) = @_;

   # start to define useful constants
   # inorder to build the XMB given
   # a vcf line and allosome key
   my $key = $input{'key'};
   my ($chr,$pos) = split(':',$input{'key'});
   my $row = $input{'data'};
   my @col = split(/\t/,$row);
   my $mis = 0;
   my $isY = 0;

   # get information that directly is involved
   # with the value/location of the bit
   my ($anc,$alt) = split(":", $allosomes{$key});
   my $bit = $bit_loc{$key} - 1;

   # determine if it is a Y chromosome
   if($chr eq 'Y'){
      $isY = 1;
      $SMB[0] = 1;
   }

   # determine if multiple alternates were listed
   if(! $col[4] =~ /^\s*\w{1}\s*/){
      $ucn =1;

      # determine if alternate is listed
      if( index($col[4],$alt) != -1){
         $col[4] = $alt;
      }
      else{ $col[4] = "P"; }
   }

   $mis = ($col[9] =~ $misR)? 1:0;

   # homozygous to ancestral allele
   if($col[9] =~ $homo2 && $col[3] eq $anc){
      if($isY){

      }
      else{ $XMB[$bit] = 0;}
   }

   # homozygous to alternate allele
   elsif($col[9] =~ $homo1 && $col[4] eq $alt){
      if($isY){

      }
      else{ $XMB[$bit+1] = 1;}
   }

   # heterozygous
   elsif(!$mis && $col[3] eq $anc && $col[4] eq $alt){
      if($isY){

      }
      elsif($SMB[0]){
         $XMB[$bit+1] = 1;
      }
      else{
         $XMB[$bit +1] = 1; $XMB[$bit] = 0;
      }
   }
}

sub bam_zygosity{
   my (%input) = @_;
   
   # get aligned reads at given position
   my @pairs = $file->get_features_by_location(-type=>'read_pair',-seq_id=>$input{'chr'},
                  -start=>$input{'loc'},-end=>$input{'loc'});
   # counters
   my $reads = 0;
   my %alleles = ('G' => 0, 'C'=> 0, 'T'=> 0, 'A' => 0);
   
   # loop through pairs
   for my $pair(@pairs){
      my ($f,$s) = $pair->get_SeqFeatures;
      my $start = $input{'loc'} - $f->start;
      
      # filter if score is not good enough
      my $score = ($f->qscore)[$start];
      next unless $score >= $baq;
      
      # increment the allele counters
      my $bp = substr($f->query->dna,$start,1);
      $alleles{$bp}++;
      
      $reads++; #increment reads, accounts for noise
   }

   # determine the 2 most populated alleles
   # and compare to see if homo or hetero
   my $al_1 = (sort {$alleles{$a} <=> $alleles{$b}} keys %alleles)[(keys %alleles)-1];
   my $al_2 = (sort {$alleles{$a} <=> $alleles{$b}} keys %alleles)[(keys %alleles)-2];
   my $al1 = $alleles{$al_1};
   my $al2 = $alleles{$al_2};

   # ensure there were reads
   if ( $al1 ==0 && $al2 == 0  ){
      return (exists $input{'smb'})? "1:0" :-1;
   }

   # ensure that it was greater than noise
   my $zyg_prob = $al1/($al1 + $al2) * 100;
   if($zyg_prob**$reads < $noise){
      return (exists $input{'smb'})? "1:0": 0;
   }
      
   # handle if this is sex chromosome bits
   if(exists $input{'smb'} ){
      # if chr is Y chromosome
      if($input{'chr'} eq 'Y'){

      }
    
      # determine if homozygous
      if($al1 >= 0.90 * ($al1 + $al2)){
         # if homo, determine for which allele
         if($al_1 eq $input{'anc'}){
            return "0:0";
         }
         elsif($al_1 eq $input{'alt'}){
            return "1:1";
         }
         return "1:0";
      }
      
      # then must be some form of 
      # hetero zygous
      else{
         # ensure its hetero between the 
         # ancestral and alternate alleles
         if( $alleles{$input{'anc'}} == 0  or
            $alleles{$input{'alt'}} == 0 ){
               return "1:0";
         }
         
         # heterozygous
         return "0:1";
      }
       
   }

   # return zygosity
   # ensuring that $al1 is over 90% of reads
   return ( $al1  >= 0.90 * ($al1+$al2) )? 0:1; 
}


package main;

my $vcfFile = "/export/home/yusuf/geneomeID/sample2.bam";
   #$vcfFile = "/export/home/yusuf/geneomeID/FC08-NGS051.mapreads.diBayes.chrX.vcf.gz";

my $genID = genomeID::generate_id('type'=>'bam','file'=>$vcfFile,'sex'=>0,'hg'=>'hg19','ref'=>1);










