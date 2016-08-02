#!/usr/bin/perl
use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::Tools::GFF;
use Bio::SeqIO;
use Data::Dumper;
use Switch;

# get the command line arguments
my $man         = 0;
my $help        = 0;

pod2usage("$0: No arguments specified.") if ( @ARGV == 0 );
my %command_args;
%{ $command_args{arguments} } = (@ARGV);
GetOptions(
	'help|?'      	=> \$help,
	man       		=> \$man
) or pod2usage(2);
pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

pod2usage(1) unless $#ARGV+1 == 1;

## Input & output
my $infile = $ARGV[0];

# specify input via -fh or -file
my $gffio = Bio::Tools::GFF->new(-file => $infile, -gff_version => 2);

# loop over the input stream
print STDERR "reading file\n";
my %entry;
while(my $feature = $gffio->next_feature()) {
	my ($transcriptID) = $feature->get_tag_values('transcript_id');

	switch ($feature->primary_tag) {
		case 'transcript' {
			$entry{$transcriptID}{'start'} 				= int($feature->location->start);
			$entry{$transcriptID}{'end'}				= int($feature->location->end);
			$entry{$transcriptID}{'chromosome'}			= $feature->seq_id;
			$entry{$transcriptID}{'strand'}				= ($feature->location->strand > 0) ? '+' : '-';
		}
		case 'exon' {
			unless ( $entry{$transcriptID} ) {
				$entry{$transcriptID}{'chromosome'}	= $feature->seq_id;
				$entry{$transcriptID}{'strand'}		= ($feature->location->strand > 0) ? '+' : '-';
				$entry{$transcriptID}{'start'} 		= int($feature->location->start);
				$entry{$transcriptID}{'end'}		= int($feature->location->end);
			} else {
				$entry{$transcriptID}{'start'} 		= int($feature->location->start) if int($feature->location->start) < $entry{$transcriptID}{'start'} ;
				$entry{$transcriptID}{'end'}		= int($feature->location->end)   if int($feature->location->end) > $entry{$transcriptID}{'end'};
			}
			$entry{$transcriptID}{'tmp_exons'} = [] unless( $entry{$transcriptID}{'tmp_exons'} );
			push(@{$entry{$transcriptID}{'tmp_exons'}}, {
				start	=> int($feature->location->start), 
				end 	=> int($feature->location->end)});
		}
	}

}
$gffio->close();

print STDERR "sorting exons and generating BED file\n";

while(my ($transcriptID, $t) = each(%entry)) {
	$entry{$transcriptID}{'exons'} = [];

	# sort exons
	foreach my $exon ( sort { $a->{start} <=> $b->{start} } @{$entry{$transcriptID}{'tmp_exons'}} ) {
		push(@{$entry{$transcriptID}{'exons'}}, $exon);
	}

	delete($entry{$transcriptID}{'tmp_exons'});
	
	my $blockCount	= 0;
	my $blockSizes	= '';
	my $blockStarts	= '';
	
	foreach my $e ( @{$entry{$transcriptID}{'exons'}} ){
		$blockCount++;
		$blockSizes 	.= ( $e->{end} - $e->{start} + 1 ).",";
		$blockStarts	.= ( $e->{start} - $t->{start} ).",";
	}

	# transcripts with no exons corrupt the file
	next unless ($blockCount > 0);

	# chromosome must have a trailing "chr"
	my $chr = $t->{chromosome};
	$chr = "chr".$chr unless ($chr =~ /^chr/);

	#chr1	529838	532878	TCONS_00000125	0.0	+	529838	532878	0,0,0	2	298,200,	0,2840,
	print STDOUT 
		$chr			."\t".
		($t->{start}-1)	."\t".
		$t->{end} 		."\t".
		$transcriptID 	."\t".
		"0"			."\t".
		$t->{strand} 	."\t".
		($t->{start}-1) ."\t".
		$t->{end} 		."\t".
		"0,0,0"			."\t".
		$blockCount		."\t".
		$blockSizes		."\t".
		$blockStarts	."\n";

}


__END__

=head1 NAME

gtf_to_bed.pl - Convert GTF to BED files

=head1 SYNOPSIS

gtf_to_bed.pl INPUT > OUTPUT

INPUT is a gtf file
OUTPUT is a bed file

Use -? to see options

=head1 OPTIONS

=over 8

=item B<-? --help>	Print a brief help message and exits.

=item B<-m --man>		Displays the manual page.

=item B<-i --in>		Name of the input file.

=back

=head1 DESCRIPTION

Save GTF information to MongoDB

=cut