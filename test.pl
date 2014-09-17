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
my %autoSomes = ('hg19' => ['1:45973928'=> 'A:G', '1:67861520'=> 'A:C', '1:158582646'=> 'C:T', '1:179520506'=> 'A:G',
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

# declare some useful global vars
my $file; my $type; my $baq; my $noise; my $sex;
my $MADIB = ""; my $AMB = ""; my $SMB = "";

sub generate_id{
	# requires input as hash
	my (%input) = @_;
	
	# determine hg version
	if(exists $input{'hg'}){
		(%pengelly) = @{$pengelly{$input{'hg'}}};
		(%bit_loc) = @{$bit_loc{$input{'hg'}}};
		(%autoSomes) = @{$autoSomes{$input{'hg'}}};
	}
	else{
		#somehow guess it
	}
	
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
	#return encode_base64 pack 'B*', "$MADIB$AMB$SMB";
	return "$MADIB$AMB$SMB\n";
}

sub bam_AMB{
	my (%input) = @_;
	my $missing_markers = 0; my $only_peng = 1; my $bit_mis =0;
	my @amb = ("0") x 58;	

	# loop through autoSomal markers
	foreach my $key (keys %autoSomes){  
		# grab marker location to query bam file
		my ($chr,$pos) = split(":",$key);
		my ($fir,$las) = split(":",$autoSomes{$key});
		my $isPengelly = $pengelly{$key};		
		
		# get zygosity
		my $zyg = bam_zygosity('chr'=>$chr,'loc'=>$pos,'al1'=>$fir,'al2'=>$las);
		
		# if no reads
		if($zyg == -1){
			# append to AMB
			$amb[$bit_loc{$key}-1] = "0"; 
			$missing_markers++;
			$bit_mis = $bit_loc{$key}-1;

			# if this is a missing pengelly, 
			# then not only pengelly's are read
			if($isPengelly){ $only_peng =0;}
		}
		else{
			# otherwise, append as usual
			$amb[$bit_loc{$key}-1] = $zyg;

			# if this is a read non pengelly,
			# then not only pengelly's are read
			if(!$isPengelly){$only_peng = 0;}
		} 
	}
	$AMB = join("",@amb); 
	genMADIB('miss_count'=>$missing_markers,'oPeng'=>$only_peng,'bMis'=>$bit_mis);
}

sub vcf_AMB{
	my $homo1 = qr/^(1[\/|]1)[^\/]/;
	my $homo2 = qr/^(0[\/|]0)[^\/]/;
	my $misR = qr/^(\.[\/|]\.)/;

	# define useful constants
	my $missing_markers = 0;
	my $only_peng  = 1; my $bit_mis = 0;
	my @amb = ("0") x 58; 
	my $found = 0; my $found_peng = 0;
	
	open my $vcf, '<', $file or die "$!\n";
	while(<$vcf>){
		# make sure its a chromosome
		next if /^#/;
		
		# get required information
		my @fields = split(/\t/,$_);
		my $chr = $fields[0];
		my $loc = $fields[1];
		my $key = "$chr:$loc";

		# determine if this is marker
		next unless exists $pengelly{$key};
		
		my $zyg = ($fields[9] =~ $homo1 or $fields[9] =~ $homo2)? 0:1;
		my $mis = ($fields[9] =~ $misR)? 1:0;
		my $bit = $bit_loc{$key};
		
		#fill in the correct amb bit
		$amb[$bit-1] = $zyg * (!$mis);
		my $isPengelly = $pengelly{$key};
		$found_peng += $isPengelly;	
		
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
		$found += (!$mis);
	}
	$AMB = join("",@amb);
	
	# deal with missing markers
	$missing_markers += 58 - $found;
	$only_peng = ($only_peng && $found_peng == 24);
	
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
	# determine the majority alleles
	# and compare to see if homo or hetero
	my $al_1 = (sort {$alleles{$a} <=> $alleles{$b}} keys %alleles)[(keys %alleles)-1];
	my $al_2 = (sort {$alleles{$a} <=> $alleles{$b}} keys %alleles)[(keys %alleles)-2];
	my $al1 = $alleles{$al_1};
	my $al2 = $alleles{$al_2};

	# ensure there were reads
	if ( $al1 ==0 && $al2 == 0  ){
		return -1;
	}
	
	# ensure that it was greater than noise
	my $zyg_prob = $al1/($al1 +$al2)*100;
	if($zyg_prob**$reads < $noise){
		return 0;
	}		

	return ( $al1  >= 0.90 * $reads )? 0:1;	
}

sub genSMB{
	# if Y chromosome is present
	$SMB = (grep $_ eq 'Y', $file->seq_ids)	. 
		(($sex )? $_[0] : (!(grep $_ eq 'Y', $file->seq_ids))? "1":"0");
	
	# append version block
	$SMB = "$SMB" . 0 x (6-length(sprintf("%b",$VERSION))) . sprintf("%b",$VERSION);  
}









package main;

my $vcfFile = "/export/home/yusuf/geneomeID/HG00157.1000g.vcf";
   #$vcfFile = "/export/home/yusuf/geneomeID/HG00157.mapped.ILLUMINA.bwa.GBR.exome.20120522.bam";

my $genID = genomeID::generate_id('type'=>'vcf','file'=>$vcfFile,'sex'=>0,'baq'=>30,'noise'=>0.05,'hg'=>'hg19');
print "$genID\n";











