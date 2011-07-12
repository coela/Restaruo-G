use strict;
use G;
use SWISS::Entry;
use SWISS::IDs;
use Storable qw(nstore store_fd nstore_fd freeze thaw dclone);
use TokyoCabinet;
use Data::Dumper;

my @databases = ('sprot');

my %files;
$files{'cds_fasta'} = "$ENV{'HOME'}/.restauro_g/tmp/cds_fasta.fasta";
$files{'cds_fasta'} = "./hoge.fasta";
$files{'trembl_fasta'} = "$ENV{'HOME'}/.glang/restauro_g/uniprot/trembl_fasta.fasta";
$files{'blast_output'} = "hoge.out";
$files{'dbm_uniprot'} = "./foo.db";
$files{'sprot'} = "./uniprot_sprot_bacteria.dat";
$files{'sprot_fasta'} = "./uniprot_sprot_bacteria.fasta";

my $tc = TokyoCabinet::HDB->new();
if(!$tc->open("uniprot.tch", $tc->OWRITER | $tc->OCREAT)){
	my $ecode = $tc->ecode();
	printf STDERR ("open error: %s\n", $tc->errmsg($ecode));
}

my $tc_length = TokyoCabinet::HDB->new();
if(!$tc_length->open("uniprot_length.tch", $tc_length->OWRITER | $tc_length->OCREAT)){
	my $ecode = $tc_length->ecode();
	printf STDERR ("open error: %s\n", $tc_length->errmsg($ecode));
}

###
my $gb = load ($ARGV[0]);
$gb->disable_pseudogenes();

foreach my $database (@databases){
	create_fasta($gb);
	perform_blast_search($gb, $database);
}
###

sub perform_blast_search {
	my $gb = shift;
	my $database = shift;
#system("blastp -query $files{'cds_fasta'} -db $files{"${database}_fasta"} -outfmt 6 -out $files{'blast_output'} -max_target_seqs 5 -num_threads 4");
	my $last_query;
	my $blast;
	open my $BLAST_OUTPUT, $files{'blast_output'} or die;
	while(<$BLAST_OUTPUT>){
		chomp;
		my @line = split /\t/;
		my $subject_seq_length = $tc_length->get($line[1]);
		my $identity = $line[3]/$subject_seq_length;
		my $evalue = $line[10];

		if ($evalue <= 1e-70 and $identity >= 0.98){
			$blast->{$line[1]}->{'level'} = 1;
		}
		elsif($evalue <=1e-50 and $identity >= 0.95){
			$blast->{$line[1]}->{'level'} = 2;
		}
		elsif($evalue <=1e-30 and $identity >= 0.90){
			$blast->{$line[1]}->{'level'} = 3;
		}
		elsif($evalue <=1e-10 and $identity >= 0.80){
			$blast->{$line[1]}->{'level'} = 4;
		}
		else{
			$blast->{$line[1]}->{'level'} = 5;
		}

		$blast->{$line[1]}->{'e_value'} = $evalue;
		$blast->{$line[1]}->{'identity'} = $identity;

		if ($line[0] ne $last_query or eof){
			my @args = sort { $blast->{$a}->{'level'} <=> $blast->{$b}->{'level'}	
				or $blast->{$a}->{'e_value'} <=> $blast->{$b}->{'e_value'}
				or $blast->{$b}->{'identity'} <=> $blast->{$a}->{'identity'}
			} keys %$blast;
			my $tophit = shift @args;
			my $uniprot = thaw $tc->get($tophit);

			$gb->{$last_query}->{on} = 0 if $blast->{$tophit}->{"level"} <= 2;
			print $last_query, " ",$uniprot->ID," ",$blast->{$tophit}->{e_value},"\n";
			$blast = {};
		}
		$last_query = $line[0];
	}
}

sub create_fasta {
	my $gb = shift;
	open my $CDS_FASTA, '>', $files{'cds_fasta'} or die;
	for my $cds ( $gb->feature('CDS') ) {
		my $translate = $gb->{$cds}->{'translation'} or translate $gb->get_geneseq($cds);
		print $CDS_FASTA '>'.$cds."\n";
		print $CDS_FASTA $translate."\n";
	}
	close $CDS_FASTA;
}

sub store_uniprot_to_kt {
	open my $UNIPROT, $files{'sprot'} or die;	
	$/ = "\/\/\n";
	while (<$UNIPROT>){
		my $entry = SWISS::Entry->fromText($_);
		print $entry->ID, "\n";
		$tc->put($entry->ID, freeze $entry);
		my $aa_length = length $entry->SQ;
		$tc_length->put($entry->ID, $aa_length);
	}	
	$/ = "/n";
}

sub create_fasta_from_uniprot_flatfile {
	open my $UNIPROT, $files{'sprot'} or die;
	open my $UNIPROT_FASTA, '>', $files{'sprot_fasta'} or die;
	$/ = "\/\/\n";
	while (<$UNIPROT>){
		my $entry = SWISS::Entry->fromText($_);
		print $UNIPROT_FASTA '>'.$entry->ID."\n";
		print $UNIPROT_FASTA $entry->SQ."\n";
	}  
	$/ = "/n";
}

1;
