#!/usr/bin/perl
use strict;
use warnings;
use Bio::DB::Sam;
use MIME::Base64;

package genomeID;
use MIME::Base64;

my $VERSION = 1;

# declare marker locations
my %autoSomes = ('1' => [[45973928,'A','G'],
			 [67861520,'A','C'],
			 [158582646,'C','T'],
			 [179520506,'A','G'],
			 [209968684,'A','C'],
			 [228431095,'C','T']],
		 '2' => [[49381585,'A','G'],
			 [75115108,'A','G'],
			 [169789016,'A','G'],
			 [215820013,'A','G'],
			 [227896976,'C','T']],
		 '3' => [[4403767,'C','T'],
			 [45989044,'G','T'],
			 [148727133,'A','G']
			],
		 '4' => [[5749904,'A','G'],
			 [86915848,'C','T']]
		 
		);

# declare strings making the ID
my $MADIB = "";
my $AMB = "";
my $SMB = "";

# declare some useful global vars
my $file; my $bam; my $baq; my $noise; my $sex;

sub generate_id{
	# requires input as hash
	my (%input) = @_;
	
	# store input flags
	$baq = $input{'baq'}; $noise = $input{'noise'};
	$sex = $input{'sex'}; 
	
	# determine if BAM or VCF
	if (exists $input{'bam'}){
		$bam=1; $file = $input{'bam'};
	
		# generate MADIB and AMB
		genMADIB_AMB();
		
		# generate sex and version block  markerBits
		genSMB((exists $input{'ucn'})? 1:0);			
	}
	elsif(exists $input{'vcf'}){
		$bam=0; $file = $input{'vcf'};
	}
	else{
		die "No BAM or VCF file was specified\n";
	}
	#my $j = encode_base64('10');
	#print "$j";
	print "$MADIB $AMB $SMB\n";
}

sub genSMB{
	# if Y chromosome is present
	$SMB = (grep $_ eq 'Y', $file->seq_ids)	. 
		(($sex )? $_[0] : (!(grep $_ eq 'Y', $file->seq_ids))? "1":"0");
	
	# append version block
	$SMB = "$SMB" . 0 x (6-length(sprintf("%b",$VERSION))) . sprintf("%b",$VERSION);  
}

sub genMADIB_AMB{
	# open file using sam tools
	$file = Bio::DB::Sam->new(-bam =>$file);
	
	# loop through autoSomal markers
	foreach my $chr (keys %autoSomes){
		
		foreach (@{$autoSomes{$chr}}){
			# grabbing data from this ugly data structure
		 	my $pos = @{$_}[0]; my $f = ${$_}[1];
			my $s = ${$_}[2];
			
			# get zygosity
			my $zyg = zygosity('chr'=>$chr,'loc'=>$pos,'al1'=>$f,'al2'=>$s);
			
			#append to MADIB and AMB	
			$MADIB = "$MADIB" . ( ($zyg==-1)? "0" : "1" );				
			$AMB = "$AMB" . (($zyg == -1)? "0" : $zyg);
		}
	}
}

sub zygosity{
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
	return (abs ($al1-$al2)/($al1 + $al2) > 0.05)? 1:0;
}
















package main;

my $bamFile = "/export/home/yusuf/geneomeID/HG00157.mapped.ILLUMINA.bwa.GBR.exome.20120522.bam";
my %test_inputs = ('bam'=>$bamFile,'sex'=>0,'baq'=>30,'noise'=>0.05);

genomeID::generate_id(%test_inputs);

















