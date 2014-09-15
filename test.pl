#!/usr/bin/perl
use strict;
use warnings;
use Bio::DB::Sam;
use MIME::Base64;

package genomeID;
use MIME::Base64; 
use Data::Dumper;

my $VERSION = 1;
# declar pengelly marker locations
my %bit_loc = ( 'hg19' => ['1:45973928'=> 1, '1:67861520'=> 2, '1:158582646'=> 3, '1:179520506'=> 4,
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
					'6:56471402'=> 57, '6:146755140'=> 58]);

my %pengelly = ( 'hg19' => ['1:45973928'=> 0, '1:67861520'=> 1, '1:158582646'=> 0, '1:179520506'=> 1,
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


# declare marker locations 
my %autoSomes = ('1' => [[45973928,'A','G'], [67861520,'A','C'], [158582646,'C','T'], 
			 [179520506,'A','G'], [209968684,'A','C'], [228431095,'C','T']],

		 '2' => [[49381585,'A','G'], [75115108,'A','G'], [169789016,'A','G'],
			 [215820013,'A','G'], [227896976,'C','T']],

		 '3' => [[4403767,'C','T'], [45989044,'G','T'], [148727133,'A','G']],

		 '4' => [[5749904,'A','G'], [86915848,'C','T']],

		 '5' => [[13719022,'A','C'], [40981689,'A','C'],[55155402,'C','T'],
			 [82834630,'A','G'],[138456815,'A','G'],[171849471,'A','G']],

		 '6' => [[56471402,'A','G'],[146755140,'A','G']],

		 '7' => [[34009946,'C','T'],[48450157,'C','T'],[100804140,'C','T'],[151254175,'C','T']],

		 '8' => [[94935937,'C','T']],

		 '9' => [[27202870,'C','T'],[77415284,'A','C'],[100190780,'C','T']],

		 '10'=> [[69926097,'A','G'],[100219314,'A','G']],

		 '11'=> [[6629665,'C','T'],[16133413,'A','G'],[30255185,'C','T']],

		 '12'=> [[993930,'C','T'],[52200742,'A','C']],

		 '13'=> [[25466955,'C','T'],[39433606,'A','G']],

		 '14'=> [[50769717,'A','G'],[76045858,'A','G']],

		 '16'=> [[70303580,'C','T']],

		 '17'=> [[7192091,'C','T'],[42449789,'C','T'],[71197748,'A','G']],

		 '18'=> [[21413869,'C','T']],

		 '19'=> [[10267077,'A','G'],[33353464,'A','G'],[55441902,'C','T']],

		 '20'=> [[6100088,'C','T'],[19970705,'A','G'],[35865054,'C','T'],[52786219,'A','G']],

		 '21'=> [[44323590,'G','T']],

		 '22'=> [[21141300,'C','T'],[37469591,'A','G']]
		);

# declare strings making the ID
my $MADIB = ""; my $AMB = ""; my $SMB = "";

# declare some useful global vars
my $file; my $type; my $baq; my $noise; my $sex;

sub generate_id{
	# requires input as hash
	my (%input) = @_;
	
	# determine hg version	
	(%pengelly) = @{$pengelly{$input{'hg'}}};
	(%bit_loc) = @{$bit_loc{$input{'hg'}}};
	
	# store input flags
	$baq = $input{'baq'}; $noise = $input{'noise'};
	$sex = $input{'sex'}; $file = $input{'file'}; 
	
	# determine if BAM or VCF
	if ($input{'type'} eq 'bam'){
		# open file using sam tools
		$file = Bio::DB::Sam->new(-bam =>$file);		
	
		# generate MADIB and AMB
		bam_AMB();
		
		# generate sex and version block  markerBits
		genSMB((exists $input{'ucn'})? 1:0);
	}
	elsif($input{'type'} eq 'vcf'){
		# generate MADIB and AMB
		vcf_AMB();
	}
	else{
		die "No supported file type was specified\n";
	}
	
	# print base 64 code
	return encode_base64 pack 'B*', "$MADIB$AMB$SMB";
}

sub genSMB{
	# if Y chromosome is present
	$SMB = (grep $_ eq 'Y', $file->seq_ids)	. 
		(($sex )? $_[0] : (!(grep $_ eq 'Y', $file->seq_ids))? "1":"0");
	
	# append version block
	$SMB = "$SMB" . 0 x (6-length(sprintf("%b",$VERSION))) . sprintf("%b",$VERSION);  
}

sub bam_AMB{
	my (%input) = @_;
	my $missing_markers = 0; my $only_peng = 1;my $bit=0; my $bit_mis =0;
	
	# loop through autoSomal markers
	foreach my $chr (keys %autoSomes){
		
		# loop through each marker per chromosome
		foreach (@{$autoSomes{$chr}}){
			# grabbing data from the ugly autosomal structure
		 	my $pos = @{$_}[0]; my $f = ${$_}[1]; 
			my $s = ${$_}[2]; $bit++; my $zyg = 0;
			my $isPengelly = $pengelly{'$chr:$pos'};
				
			# get zygosity 
			$zyg = bam_zygosity('chr'=>$chr,'loc'=>$pos,'al1'=>$f,'al2'=>$s);
			
			# if no reads
			if($zyg == -1){
				# append to and AMB
				$AMB = "$AMB" . "0";
				$missing_markers++;
				$bit_mis = $bit;				

				# determine if its pengelly unkown
				if($isPengelly){
					$only_peng =0;
				}
			}else{
				# append as usual
				$AMB = "$AMB$zyg";
				
				# determine if only pengelly
				if(!$isPengelly){
					$only_peng = 0;
				}
			}
		}
	}
	genMADIB('miss_count'=>$missing_markers,'oPeng'=>$only_peng,'bMis'=>$bit_mis);
}

sub vcf_AMB{
	my $chr = qr/^\d{0,}[XY]?\t\d+/s;
	my $homo = qr/^(1[\/|]1)[^\/]/;
	my $misR = qr/^(\.[\/|]\.)/;

	# define useful constants
	my $missing_markers = 0;
	my $only_peng  = 1; my $bit_mis = 0;
	my @amb = ("0") x 58;
	
	open my $vcf, '<', $file or die "$!\n";
	while(<$vcf>){
		# make sure its a chromosome
		next unless $_ =~ $chr;
		
		# get required information
		my @fields = split(/\t/,$_);
		my $chr = $fields[0];
		my $loc = $fields[1];
		my $key = "$chr:$loc";

		# determine if this is marker
		next unless exists $pengelly{$key};
		
		my $zyg = ($fields[9] =~ $homo)? "0" : "1";
		my $mis = ($fields[9] =~ $misR)? "1":"0";
		my $bit = $bit_loc{$key};
		
		#fill in the correct amb bit
		$amb[$bit-1] = $zyg * (!$mis);
		my $isPengelly = $pengelly{$key};		
		
		# determine if missed read
		if($mis){
			$missing_markers++;
			$bit_mis = $bit;

			if($isPengelly){
				$only_peng = 0;
			}
		}elsif(!$isPengelly){
			$only_peng = 0;
		}
	}
	$AMB = join("",@amb);
	genMADIB('miss_count'=>$missing_markers,'oPeng'=>$only_peng,'bMis'=>$bit_mis);	
}

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

sub bam_zygosity{
	my (%input) = @_;
	
	# get aligned reads at given position
	my @pairs = $file->get_features_by_location(-type=>'read_pair',-seq_id=>$input{'chr'},
						-start=>$input{'loc'},-end=>$input{'loc'});
	# counters
	my $al1=0; my $al2 = 0; my $reads = 0;
	# loop through pairs
	for my $pair(@pairs){
		my ($f,$s) = $pair->get_SeqFeatures;
		my $start = $input{'loc'} - $f->start;
		
		# filter if score is not good enough
		my $score = ($f->qscore)[$start];
		next unless $score >= $baq;
		
		# increment the allele counters
		my $bp = substr($f->query->dna,$start,1);
		if($bp eq $input{'al1'}){
			$al1++;
		}
		elsif($bp eq $input{'al2'}){
			$al2++;
		}
			
		$reads++; #increment reads, accounts for noise
	}
	
	# ensure there were reads
	if ( $al1 ==0 && $al2 == 0  ){
		return -1;
	}
	
	# ensure that it was greater than noise
	my $zyg_prob =($al1,$al2)[$al2>$al1]/($al1 +$al2)*100;
	if($zyg_prob**$reads < $noise){
		return 0;
	}		
	
	# return zygosity, must be within 5% to be homo	
	return (abs ($al1-$al2)/($al1 + $al2) > 0.10)? 1:0;
}












#	Ev+vv9f/d92B
#	000100101111111110101111101111111101011111111111011101111101110110000001



package main;

my $vcfFile = "/export/home/yusuf/geneomeID/HG00157.1000g.vcf";
   #$vcfFile = "/export/home/yusuf/geneomeID/HG00157.mapped.ILLUMINA.bwa.GBR.exome.20120522.bam";

my $genID = genomeID::generate_id('type'=>'vcf','file'=>$vcfFile,'sex'=>0,'baq'=>30,'noise'=>0.05,'hg'=>'hg19');
print "$genID\n";











