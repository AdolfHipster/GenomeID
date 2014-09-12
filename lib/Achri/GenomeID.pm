package Achri::GenomeID;

use 5.006000;
use strict;
use warnings;
require Exporter;
require Bio::DB::Sam;

use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

$VERSION = 1.00;
@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw(generate_id);
%EXPORT_TAGS ( DEFAULT => [qw(&generate_id)] );

sub generate_id {
	# hash 

}



1;
__END__
# Below is stub documentation for your module. You'd better edit it!

=head1 NAME

Achri::GenomeID - Perl extension for blah blah blah

=head1 SYNOPSIS

  use Achri::GenomeID;
  blah blah blah

=head1 DESCRIPTION

Stub documentation for Achri::GenomeID, created by h2xs. It looks like the
author of the extension was negligent enough to leave the stub
unedited.

Blah blah blah.

=head2 EXPORT

None by default.



=head1 SEE ALSO

Mention other useful documentation such as the documentation of
related modules or operating system documentation (such as man pages
in UNIX), or any relevant external documentation such as RFCs or
standards.

If you have a mailing list set up for your module, mention it here.

If you have a web site set up for your module, mention it here.

=head1 AUTHOR

A. U. Thor, E<lt>yusuf@localdomainE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2014 by A. U. Thor

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.10.1 or,
at your option, any later version of Perl 5 you may have available.


=cut
