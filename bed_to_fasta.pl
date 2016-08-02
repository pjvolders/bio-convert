#!/usr/bin/perl 
use strict;
use IPC::Open2;
use Getopt::Long;
use Pod::Usage;
use Bio::SeqIO;
use Bio::Seq;

my $man = 0;
my $help = 0;

## Parse options and print usage if there is a syntax error,
## or if usage was explicitly requested.
GetOptions(
	'help|?' => \$help, 
	man => \$man
) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;

pod2usage(1) unless $#ARGV+1 == 2;

## Fastahack code
my $cmd = 'fastahack -c /opt/genomes/9606/GRCh37/fastahack/GRCh37.fa';
my $pid = open2(\*DATAOUT, \*DATAIN, $cmd);

## Input & output
my $infile = $ARGV[0];
my $outfile = $ARGV[1];

# open output file
open(OUTPUT, '>', $outfile) or die;
my $fasta_file = Bio::SeqIO->new(-fh => \*OUTPUT, -format => 'fasta' );


open INFILE, $infile or die $!;
while (<INFILE>) {
	chomp $_;
	my @gen = split(/\t/, $_);
	
	my $chr = $gen[0];
	my $start = int($gen[1])+1; #WHYYYYYYY?????!!!
	my $end = int($gen[2]);
	my $name = $gen[3];
	my $strand = $gen[5];

	my $seq;

	# BED6 vs BED12
	if ( scalar(@gen) == 12 && $gen[9] > 1 ) {
		# entry had exon information => determine sequence based on exon positions

		my @blockSizes = split(/,/, $gen[10]);
		my @blockStarts = split(/,/, $gen[11]);
		#print($gen[10]."\t".$gen[11]."\n");
		
		my @exons;
		my $exonRank = 1;
		
		for( my $i = 0; $i < @blockStarts; $i++ ) {
			my $exonStart = int($start+$blockStarts[$i]);
			my $exonEnd = int($start+$blockStarts[$i]+$blockSizes[$i]-1); #WHY 2
			
			my $exon_seq = getSequence($chr, $exonStart, $exonEnd);
			#getSequence('chr12',54356092,54357908);
			
			if ( $strand eq '-' ) {
				$exon_seq = RevComp($exon_seq);
				$seq = $exon_seq.$seq;
			} else {
				$seq = $seq.$exon_seq;
			}
		}
	} else {
		$seq = getSequence($chr, $start, $end);
		if ( $strand eq '-' ) {
			$seq = RevComp($seq);
		}
	}
	
	my $seq_obj = Bio::Seq->new(-seq => uc($seq),
								-display_id => $name,
								-alphabet => "dna" );
	$fasta_file->write_seq($seq_obj);
	
}

close INFILE;



sub getSequence($$$) {
	my ($chr, $start, $end) = @_;
	$chr =~ s/chr//g;
	print DATAIN "$chr:$start..$end"."\n";
	my $out = <DATAOUT>;
	chomp $out;
	#print "SEQ: $out\n";
	return $out;
}


# maak het reverse compliment met een subroutine
# gebruik: $b = RevComp($a);
sub RevComp {
   my ($seq) = @_;
   $seq = reverse($seq);
   $seq =~ tr/ATGC/TACG/;
   return $seq;
}
__END__

=head1 NAME

sample - Using GetOpt::Long and Pod::Usage

=head1 SYNOPSIS

bed_to_fasta.pl INPUTFILE OUTPUT

 Options:
   -help            brief help message
   -man             full documentation

=head1 OPTIONS

=over 8

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> will read the given input file(s) and do something
useful with the contents thereof.

=cut