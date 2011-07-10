use strict;
use G;
use Cache::KyotoTycoon;
use Storable;

my %files;
$files{'cds_fasta'} = "$ENV{'HOME'}/.glang/restauro_g/tmp/cds_fasta.fasta";
$files{'sprot_fasta'} = "$ENV{'HOME'}/.glang/restauro_g/uniprot/sprot_fasta.fasta";
$files{'trembl_fasta'} = "$ENV{'HOME'}/.glang/restauro_g/uniprot/trembl_fasta.fasta";
$files{'blast_output'} = "$ENV{'HOME'}/.glang/restauro_g/tmp/blast_output.fasta";


open my $CDS_FASTA, '>', $files{'cds_fasta'} or die;
for my $cds ( $gb->feature('CDS') ) {
	my $translate = $gb->{$cds}->{'translate'} or translate $gb->get_geneseq($cds);
	print $CDS_FASTA '>'.$cds."\n";
}

close $CDS_FASTA;

my $last_query;
open my $BLAST_OUTPUT, $files{'blast_output'} or die;
while(<$BLAST_OUTPUT>){
	chomp;
	my @lines = split /\t/;
	if ($line[0] ne $last_query or eof()){
		my @args = sort {$values->{$a}->{'e_value'} <=> $values->{$b}->{'e_value'}
			or $values->{$a}->{'identity'} <=> $values->{$b}->{'identity'}}
	}
	$last_query = $line[0];
}

sub store_uniprot_to_DMB {

}

sub store_uniprot_to_kt {
	my $kt = Cache::KyotoTycoon->new(host => '127.0.0.1', port => 1978);
	$kt->set('foo' => 'bar');
	$kt->get('foo'); # => 'bar'
}
