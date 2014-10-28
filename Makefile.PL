use 5.006000;
use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    NAME              => 'Achri::GenomeID',
    VERSION_FROM      => 'lib/Achri/GenomeID.pm', # finds $VERSION
    PREREQ_PM         => {
	 	'Bio::DB::Sam' => 0.01,
		'MIME::Base64' => 0.01,
		'List::MoreUtils' => 0.01
	 }, 
    ($] >= 5.005 ?     ## Add these new keywords supported since 5.005
      (ABSTRACT_FROM  => 'lib/Achri/GenomeID.pm', # retrieve abstract from module
       AUTHOR         => 'A. U. Thor <yusuf@localdomain>') : ()),
);
