#!/usr/bin/perl -w
#use strict;
use warnings;

use common_scripts;
use genes_scripts;
use TSS_scripts;
use end3_scripts;
use sequencing_data_class;
use Getopt::Long qw(GetOptions);


my $basic_usage = qq^

############################################################################################################
#               
#       Id possible secondary TSSs [Second highest ratio atleast 10 nts away from primary TSS]
#
############################################################################################################
#
# Required:
#
# -files <string>       : 5' end coverage files (.s2)
#
##
############################################################################################################
#
# Optional :
#
# -fasta <string>       : Genome fasta file (optional after setting default)
# -annot <string>       : Genome annotation file (optional after setting default)
# -max_5utr_len <int>   : Maximum length upstream of translational start (default: 500)
# -min_cov <int>        : Minimum TAP/Processed coverage (default: 2)
# -min_ratio <int>      : Minimum TAP/noTAP|Processsed/unprocessed ratio (default: 1)
# -help                 : Print basic usage
# 
# Typical usage:
#
#        perl id_possible_secondary_TSS.pl -files <.s2> -fasta <.fasta> -annot <.genes>
############################################################################################################
 
^;

if(!@ARGV){print "$basic_usage\n";
           exit;}

GetOptions(
        'files=s{1,}' => \@file,
        'fasta:s' => \$contig_in,
        'annot:s' => \$genesfile_in,
        'max_5utr_len:i' => \$max_5utr_len_in,
        'min_cov:i' => \$min_cov,
        'min_ratio:i' => \$min_rat,
        'help|h' => \$help) or die("$basic_usage\n");


if (defined $help) {
print("$basic_usage\n");
exit;
}

######## Read genome fasta and annotations #########
my $contig = $contig_in || "/cluster/home/nmr07003/Termseq/v3/rrna/NC_003028.v3.17.fasta"; ##### Edit to default location of fasta file
my $genome_str = common_scripts->get_genome_sequence($contig);
my $genome_length = length($genome_str);
my $genesfile = $genesfile_in || "/cluster/home/nmr07003/Termseq/v3/NC_003028.v3.17.ncrna.genes"; ##### Edit to default location of annotation file
my $genes_array_ref = genes_scripts->get_genes($genesfile,$genome_length);
foreach my $gene_obj (@{$genes_array_ref}) { 
        $gene_obj->save_sequencing_data_obj( sequencing_data_class->new_seq_data_obj() );
}

my $max_5utr_len = $max_5utr_len_in || 500;
my $is_use_intergenic =0;
my $TSS_min_cov = $min_cov || 2;
my $min_tapnotap = $min_rat || 1;
my %end5_position_hash;
my %options =  ( TSS_min_cov => $TSS_min_cov, '5UTR_max_len' => $max_5utr_len, min_TAPnoTAP_ratio => $min_tapnotap,use_intergenic => 1, min_r5utr_len => 70
                );

######## Read 5'-end single nucleotide coverage for processed and unprocessed files and calculate ratios #########
foreach my $n (@file) {
#		print "Reading 5'end single nucleotide coverage file $n\n";
                my ($strain,$condition,$exp_type) = $n =~ /(^.+)\.(.+)\.(.+)\.s2/;
                open (S2, $n) or die "can't open S2 file in 5'end $n\n";
                while (<S2>) {
                        my @line = split(/\t/); my ($pos,$coverage,$strand) = @line[1,2,3]; #[2,3,4];                   
                        if (!$end5_position_hash{$pos}) {
                                $end5_position_hash{$pos} = TSS_scripts->new_5end($pos,$strand,$exp_type,$coverage);
                                $end5_position_hash{$pos}->set_exp_coverage($exp_type,$coverage);
                                }
                        else {
                                $end5_position_hash{$pos}->set_exp_coverage($exp_type,$coverage);
                        }
                }
                close S2;
        }
        while ( my ($pos,$pill_obj) = each %end5_position_hash ) {
               my ($strand,$tap,$notap,$ratio) = ($pill_obj->get_st,$pill_obj->get_TAP,$pill_obj->get_noTAP,$pill_obj->get_ratio);
                if ($notap == 0) {
                        $ratio = $tap;
                }
                else {
                        $ratio = $tap / $notap;
                }
                $pill_obj->set_ratio($ratio);
                my $a = $pill_obj->get_ratio;
                if ($a == 0){next;}
                else {
        }

}

my @genes_array = @{$genes_array_ref};
my @pillar_positions = sort {$a <=> $b} keys %end5_position_hash;
my $gene_idx = 0;
        PILLARS: for (my $pillar_idx=0;$pillar_idx<= $#pillar_positions;$pillar_idx++) { 
                my $pill_obj = $end5_position_hash{$pillar_positions[$pillar_idx]};
                my ($pill_pos,$pill_st,$TAP,$noTAP,$ratio) = $pill_obj->get_basic_pill_data;
                my $gene_obj = $genes_array[$gene_idx];
                my $next_gene_obj;
                $next_gene_obj = $genes_array[$gene_idx+1];
                my $return_message = TSS_scripts->map_pill_to_gene('5end',$pill_obj,$gene_obj,$max_5utr_len,$next_gene_obj) if ($gene_obj);
		if ($return_message){
                if ($return_message eq 'next_gene') {
                        $gene_idx++;
                        redo PILLARS;
                }
                elsif ($return_message eq 'next_pill') {
                        next PILLARS;
                }
                elsif ($return_message eq 'compare2nextGene' and $next_gene_obj) {
                        TSS_scripts->map_pill_to_gene('5end',$pill_obj,$next_gene_obj,$max_5utr_len);
                        next PILLARS;
                }
		} else {next;}
        }

print "Locus\tStart\tEnd\tStrand\tTSS\tTAP_cov\tRatio\tDistance_from_Start\tSec_TSS\tSec_TAP_cov\tSec_Ratio\tSec_distance_from_start\n";
foreach my $gene_obj (@{$genes_array_ref}) {
        ($gene_obj->get_sequencing_data_obj)->determine_TSS_and_5utr($gene_obj,\%options);
	my $l =$gene_obj->{_locus};
	my $fro = $gene_obj->{_fr};
	my $t = $gene_obj->{_to};
	my $stra = $gene_obj->{_st};
	my $pos = $gene_obj->{_seq_data}->{_5end}->{_sense}->{_5utr}=>;
	foreach my $n ($pos){
		foreach my $key (sort keys $n){
			my $posi = $n->{$key}->{_pos};
			my $b = $n->{$key}->{_ratio};
			my $tap_cov = $n->{$key}->{_TAP};
			my $notap_cov = $n->{$key}->{_noTAP};
			if ($tap_cov >= $options{'TSS_min_cov'} and $b >= $options{'min_TAPnoTAP_ratio'}){
				$hits{$posi}=$tap_cov;
				$ratios{$tap_cov}=$b;
			}

}
}

	if(%hits){
        @positions = sort {$hits{$a} <=> $hits{$b}} keys %hits;
        $max_pos = $positions[-1];
	$max_cov = $hits{$max_pos};
	$max_ratio = $ratios{$max_cov};
	foreach $p (@positions){
	 if($p >= $max_pos + 10 or $p <= $max_pos - 10){
		$pos_cov = $hits{$p};
		$pos_ratio = $ratios{$pos_cov};
		$sec_tss{$p}=$pos_cov;
		$sec_Ratio{$pos_cov}=$pos_ratio;
	  } else {
		next;
		}	
	}
	if(%sec_tss){
	@sec_positions = sort {$sec_tss{$a} <=> $sec_tss{$b}} keys %sec_tss;
	$secondary_tss = $sec_positions[-1];
	$secondary_cov = $sec_tss{$secondary_tss};
	$secondary_ratio = $sec_Ratio{$secondary_cov};
	if ($stra eq '+'){
                        $primary_dist = abs($fro - $max_pos);
			$secondary_dist = abs($fro - $secondary_tss);
                } else {
                        $primary_dist = abs($t - $max_pos);
			$secondary_dist = abs($t - $secondary_tss);
        }

	print "$l\t$fro\t$t\t$stra\t$max_pos\t$max_cov\t$max_ratio\t$primary_dist\t$secondary_tss\t$secondary_cov\t$secondary_ratio\t$secondary_dist\n" if ($primary_dist <= $max_5utr_len && $secondary_dist <= $max_5utr_len);
	} else {
	next;
	}
	undef @positions;
        undef %hits;
        undef %ratios;
	undef @sec_positions;
	undef %sec_tss;
}
}
