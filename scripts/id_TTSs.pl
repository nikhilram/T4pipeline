#!/usr/bin/perl -w 
use strict;
use warnings;
use List::Util qw( min max );
 
#Modules

use genes_scripts;
use common_scripts;
use end3_scripts;
use sequencing_data_class;
use Getopt::Long qw(GetOptions);


my $basic_usage = qq^

############################################################################################################
#               
#                                               Id TTSs
#
############################################################################################################
#
# Required:
#
# -files <string>       : 3' end coverage files (.3p)
#
##
############################################################################################################
#
# Optional :
#
# -fasta <string>       : Genome fasta file (optional after setting default)
# -annot <string>       : Genome annotation file (optional after setting default)
# -max_3utr_len <int>   : Maximum length downstream of translational stop (default: 150)
# -min_cov <int>        : Minimum first base coverage of 3'-end reads (default: 10)
# -help                 : Print basic usage
# 
# Typical usage:
#
#        perl id_TTS.pl -files <.s2> -fasta <.fasta> -annot <.genes>
############################################################################################################
 
^;

if(!@ARGV){print "$basic_usage\n";
           exit;}

GetOptions(
        'files=s{1,}' => \ my @file,
        'fasta:s' => \ my $contig_in,
        'annot:s' => \ my $genesfile_in,
        'max_3utr_len:i' => \ my $max_3utr_len_in,
        'min_cov:i' => \ my $min_cov_in,
        'help|h' => \ my $help) or die("$basic_usage\n");


if (defined $help) {
print("$basic_usage\n");
exit;
}

my $min_cov = $min_cov_in || 10;
my $max_3UTR = $max_3utr_len_in || 150;
my ($min_avg_cov,$min_reps,$max_5utr_len) = (10,1,500);
my $contig = $contig_in || "/cluster/home/nmr07003/Termseq/v3/rrna/NC_003028.v3.17.fasta";
my $genome_str = common_scripts->get_genome_sequence($contig);
my $genome_length = length($genome_str);
my $genesfile = $genesfile_in || "/cluster/home/nmr07003/Termseq/v3/NC_003028.v3.17.ncrna.genes";
my @genes_array = @{genes_scripts->get_genes($genesfile,$genome_length)};
my %pillar_hash;
my %secondary_pillar_hash;
foreach my $n (@file) {
                print "Reading file $n\n";
                my ($strain,$condition,$exp_type) = $n =~ /^.+(N.+)\.(.+)\.(.+)\.3p/;
                open (P3, $n) or die "can't open 3P file in 3'end $n\n";
                while (<P3>) {
                                my @line = split(/\t/); my ($pos,$strand,$coverage) = @line[1,2,3];
				if(!$pillar_hash{$pos}){
                                $pillar_hash{$pos} = end3_scripts->new_3end($pos,$strand,$coverage,$condition,$strain);
				} elsif($pillar_hash{$pos} and $pillar_hash{$pos}{'_st'} ne $strand) {
                                $secondary_pillar_hash{$pos} = end3_scripts->new_3end($pos,$strand,$coverage,$condition,$strain);
				} else {next;}
                }
                close P3;
}
print "Extracted strand specific 5' end coverage from the above files.\n";

print "Mapping coverages to 3'-UTRs of annotated features.....\n ";
my $out = 'TTS.' . $min_cov . '.' .  $min_reps . '.' . $max_3UTR;
my $out_bedgraph = $out . '.bedGraph';
open(TTS,'>',$out);
open(BED,'>',$out_bedgraph);
print TTS "Locus\tFrom\tTo\tStrand\tTTS\tCoverage\t3'UTR_length\t40nt_upstream\t20nt_downstream\tFold\tMFE\tNo._of_Us\tstem_loop\n";
print BED "track type=bedGraph color=17,48,66 visibility=full yLineOnOff=on autoScale=on yLineMark='0.0' alwaysZero=on graphType=bar maxHeightPixels=128:75:11 windowingFunction=maximum smoothingWindow=off transformFunch LOG\n";
my @pillar_positions =  sort {$a<=>$b} keys %pillar_hash;
my @pillar_positions_secondary =  sort {$a<=>$b} keys %secondary_pillar_hash;

foreach my $gene_obj (@genes_array) {
        my ($loc,$fr,$to,$strand) = $gene_obj->get_basic_gene_data;
        my ($fr_3UTR,$to_3UTR,$span_3UTR) = (0,0,0);
        if ($strand eq '+') {
                ($fr_3UTR,$to_3UTR) = ($to,$to + $max_3UTR);
        }
        else {
                ($fr_3UTR,$to_3UTR) = ($fr - $max_3UTR,$fr);
        }
        my %hits;
        my @positions;
        my $max_pos;
        my $max_cov;
        my $len;
        my @pillars_in_bounds = grep {$_ >= $fr_3UTR and $_ <= $to_3UTR and $pillar_hash{$_}{'_st'} ne $strand} @pillar_positions;
        my @pillars_in_bounds_sec = grep {$_ >= $fr_3UTR and $_ <= $to_3UTR and $secondary_pillar_hash{$_}{'_st'} ne $strand} @pillar_positions_secondary;
	if(@pillars_in_bounds){
	foreach my $pill_pos (@pillars_in_bounds) {
		my $pos_cov = $pillar_hash{$pill_pos}{_cov};
		my $position = $pillar_hash{$pill_pos}{_pos};
		if ($pos_cov >= $min_cov){
		$hits{$position} = $pos_cov;
		} else {
		next;
		}
	  }
	}
	elsif(@pillars_in_bounds_sec){
        foreach my $pill_pos (@pillars_in_bounds_sec) {
                my $pos_cov = $secondary_pillar_hash{$pill_pos}{_cov};
                my $position = $secondary_pillar_hash{$pill_pos}{_pos};
                if ($pos_cov >= $min_cov){
                $hits{$position} = $pos_cov;
                } else {
                next;
                }
          }
        }
	else {next;}
	if (%hits){
	@positions = sort {$hits{$a} <=> $hits{$b}} keys %hits;
	$max_pos = $positions[-1];
	$max_cov = $hits{$max_pos};
		if ($strand eq '+'){
		$len = $max_pos - $to;
		} else {
		$len = $fr - $max_pos;
		}
	my $term_sequence = &get_seq($genome_str,$max_pos,$strand,40);
        my $term_seq_dnstream = &get_seq($genome_str,$max_pos,$strand,-20);
        chomp (my @rna_fold_output = `echo $term_sequence | RNAfold `);### fold sequence using RNAfold
        $rna_fold_output[1] =~ /^(.+)\s+\((.+)\)/;
        my ($fold,$energy) = ($1,$2);
	$max_pos_1 = $max_pos + 1;
	my $num_of_Us = &count_U_stretch($term_sequence);
        my $stem_loop = &is_stem_loop($fold);
       	print TTS "$loc\t$fr\t$to\t$strand\t$max_pos\t$max_cov\t$len\t$term_sequence\t$term_seq_dnstream\t$fold\t$energy\t$num_of_Us\t$stem_loop\n";
	print BED "NC_003028\t$max_pos\t$max_pos_1\t$max_cov\n";
	undef %hits;
	undef @positions;
	} else {next;}
}


system("sort -k2n $out_bedgraph -o $out_bedgraph");

print "Output files written.\n";

sub get_seq {
        my ($genome_str,$pos,$strand,$span) = @_;
        my $seq;
        if ($strand eq '+')     {
                if ($span > 0) {
                        $seq = substr($genome_str,$pos-1-$span,$span);
                }
                else {
                        $span *= -1;
                        $seq = substr($genome_str,$pos-1,$span);
                }
        }
        else {
                if ($span > 0) {
                        $seq = substr($genome_str,$pos-1,$span);
                }
                else {
                        $span *= -1;
                        $seq = substr($genome_str,$pos-1-$span,$span);
                }
                $seq = reverse $seq; $seq =~ tr/ATGC/TACG/;
        };
        $seq =~ tr/T/U/;
        return $seq;
}

sub is_stem_loop {
        my ($fold) = @_;
        $fold =~ m/(\(+\.{0,2}\(+\.{0,2}\(+)(\.+)(\)+\.{0,2}\)+\.{0,2}\)+).{0,10}$/;
        if ($1 and $2 and $3 ) {
                my $stem_diff = abs(length($1) - length($3));
                my $loop_stem_diff_1 = length($2) / length($1);
                my $loop_stem_diff_2 = length($2) / length($3);
                if ($stem_diff <= 4 and $loop_stem_diff_1 < 2 and $loop_stem_diff_2 < 2) {
                        my $stem_loop = $1 . $2 . $3;
                        return $stem_loop;
                }
        }
        return 'no_stemloop';
}

sub count_U_stretch {
        my ($seq) = @_;
        my $tail = substr($seq,length($seq)-8,8);
        my @nucs = split(//,$tail);
        my $number_of_U_residues = 0;
        foreach my $nuc (@nucs) {
                $number_of_U_residues++ if ($nuc eq 'U');
        }
        return $number_of_U_residues;
}


