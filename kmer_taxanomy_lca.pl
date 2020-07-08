#!/usr/bin/perl -w
#use strict;


#taxa.list must include k,p,c,o,f,g,s
#the seq id of seq.fna must be consistent with IDs in taxa.list
@ARGV >= 2 || die "Usage: perl $0 taxa.list seq.fna 31\n";
my ($taxalist, $seqfile, $ksize) = @ARGV;

$ksize ||= 31;

my %ID_taxa = hash_build($taxalist);

my %kmer_taxa;
open SEQ, $seqfile || die $!;
$/ = '>';
<SEQ>;
while(<SEQ>){
	chomp;
	my @temp = split /\n/;
	my $tempid = shift @temp;
	my $id = (split /\s+/, $tempid)[0];
	my $seq = join("", @temp);
	my @kmers = &kmer($seq, $ksize);
#	print $id, "\t", $ID_taxa{$id}, "\t", join(";", @kmers), "\n";
	for my $k (@kmers){
		push @{$kmer_taxa{$k}}, $ID_taxa{$id};
	}
}
$/ = "\n";
close SEQ;
my $w = '#' x $ksize ;

open OUT1, '>all_kmer_taxa.list' || die $!;
open OUT2, '>lca_kmer_taxa.list' || die $!;
for my $k (sort {$a cmp $b } keys %kmer_taxa){
	print OUT1 $k, "\t", join("\n$w\t", @{$kmer_taxa{$k}}), "\n\n";
	print OUT2 $k, "\t", lca_taxanomy(\@{$kmer_taxa{$k}}), "\n";
}
close OUT1;
close OUT2;

sub kmer{
	my ($seq, $ksize) = @_;
	$seq = uc($seq);
	my $len = length($seq);
	my @kmer_pool;
	for my $i (0..$len-$ksize){
		my $kmer = substr($seq, $i, $ksize);
		push @kmer_pool, $kmer;
	}
	return @kmer_pool;
}


sub hash_build{
	my ($taxa_file) = @_;
	open TAXA, $taxa_file;
	my %idinfo;
	while(<TAXA>){
		chomp;
		my @temp = split /\t+/;
		$idinfo{$temp[0]} = $temp[1];
	}
	close TAXA;
	return %idinfo;
}

sub lca_taxanomy{
	#taxon must contains k,p,c,o,f,g,s
	my ($taxa) = @_;
	my @Mtaxa = @{$taxa};
	my @units = qw(k p c o f g s);
	my %hash;
	for my $taxon (@Mtaxa){
		my @temp = split /;/, $taxon;
		for my $i (0..$#units){
			push @{$hash{$units[$i]}}, $temp[$i];
		}
	}
	my @taxa_signs;
	for my $u (@units){
		my @uniq = &uniq( @{$hash{$u}} );
		$hash{$u} = \@uniq;
		push @taxa_signs, scalar(@uniq);
	}
	my $lca_taxa;
	if(scalar(@{$hash{'k'}}) != 1){
		return 'NO';
	}else{
		my $end_pos = search_top_same_chars(\@taxa_signs);
		my @lca_tmp = map { ${$hash{$_}}[0] } @units[0..$end_pos];
		
		return join(";", @lca_tmp);
	}

}

sub uniq{
	my %count;
	return grep { ++$count{$_} < 2 } @_;
}


sub search_top_same_chars{
	my ($chars, $i) = @_;
	$i ||= 0;
	my @CHARS = @{$chars};
	if(scalar(@CHARS) == 1){
		 return $i;
	}elsif($CHARS[0] ne $CHARS[1]){
		 return $i; 
	}else{
		$i++;
		shift @CHARS;
		return search_top_same_chars(\@CHARS, $i)
	}
}
