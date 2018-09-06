#!/usr/bin/perl -w 
use strict;

###############################################   calculate_riboactivity_from_RNAseq.pl    ###########################################################
######	Use RNA-seq coverage to calcaulte the ratio of terminated vs. readthrough transcripts across all regulators
######	Takes the putative regulator table generated by the map_putative_regulatory_5UTRs.pl script
######	experiment_file specifies where to find RNA-seq coverage maps for calculation
######	run example: perl calculate_riboactivity_from_RNAseq.pl regulator_table experiment_file
######################################################################################################################################################
use lib '/cluster/home/nmr07003/Termseq/aad9822';
#modules
use common_scripts;

#command line inputs
my ($candidate_table,$experiment_file,$shift_correction) = @ARGV;
$shift_correction = 0 if (!$shift_correction);
#collect genes from candidate table 
open(CAND,$candidate_table) or die "can't open candidate file $candidate_table\n";
my @gene_candidates; 
chomp (my $line_header = <CAND>);
while (<CAND>) { 
	chomp; 
	my $line = $_;
	my $gene;
	my @line = split(/\t/);
	$gene->{_line} = $line;
	$gene->{_fr} = $line[1];
	$gene->{_to} = $line[2];
	$gene->{_st} = $line[3];
	$gene->{_tss_pos} = $line[6];
	$gene->{_term_pos} = $line[11];
	push @gene_candidates,$gene;
	}

### Collect data
open(EXPS,$experiment_file) or die "can't open experiment file $experiment_file\n";
my %exp_hash; 
my @experiments;
my @samples;
while (<EXPS>) { 
	chomp;
	my ($type,$project_dir,$run_dir,$sample_name,$exp_name) = split(/\t/);
        next unless ($type eq 'WT');
#	if ($_ =~ /\//) {$path = $_;}
#	else {$path = "/chome/nmr07003/Termseq/output/$run_dir.coverage";}
#	my ($type,$project_dir,$run_dir,$sample_name,$exp_name) = split(/\t/); 
#	next unless ($type eq 'WT');
	my $path = "/chome/nmr07003/Termseq/Karens/T4antibiotics/coverage/$run_dir.coverage"; 
	open (SAMPLE,$path) or die "can't open Coverage file for sample $run_dir at $path\n";
	push @experiments,$exp_name;
	push @samples,$sample_name;
	$exp_hash{$exp_name} = 1;
	my @genome_coverage;
	my $pos = 0; 
	while (<SAMPLE>) { 
		my ($contig,$merged,$sense,$antisense) = split(/\t/);
		$genome_coverage[$pos] = $merged;
		$pos++;
	}
	foreach my $gene (@gene_candidates) {
		my (@reg_array,@gene_array);
		my $readthrough = 'NA';
		if ($gene->{_st} eq '+') { 
			@reg_array = @genome_coverage[ ($gene->{_tss_pos}-1 + $shift_correction)..($gene->{_term_pos}-1) ];
			@gene_array = @genome_coverage[ ($gene->{_fr}-1)..($gene->{_to}-1) ];
		}
		else { 
			@reg_array = @genome_coverage[ ($gene->{_term_pos}-1)..($gene->{_tss_pos}-1 - $shift_correction) ];
			@gene_array = @genome_coverage[ ($gene->{_fr}-1)..($gene->{_to}-1) ];
		}
		my $reg_med_cov = common_scripts->average(\@reg_array);
		my $gene_med_cov = common_scripts->average(\@gene_array);		
		$readthrough = ($gene_med_cov+0.01) / ($reg_med_cov+0.01);
		$readthrough = 1 if ($readthrough > 1); #limit to 100%
		$gene->{_readthrough}->{$exp_name}->{$sample_name} = $readthrough;
		$gene->{_gene_cov}->{$exp_name}->{$sample_name} = $gene_med_cov;
		$gene->{_reg_cov}->{$exp_name}->{$sample_name} = $reg_med_cov;
	}
}

# calculate mean readthrough and stdv
foreach my $gene (@gene_candidates) {  
	my @exps =  keys %{$gene ->  { _readthrough }}; 
	foreach my $exp_name (@exps) { 
		my @reathrough_per_sample_in_exp; 
		my @samples = keys %{$gene->{_readthrough}->{$exp_name}}; #($gene->{_readthrough}->{$exp});
		foreach my $sample_name (@samples) { 
			my $sample_readthrough = $gene->{_readthrough}->{$exp_name}->{$sample_name};
			push @reathrough_per_sample_in_exp,$sample_readthrough;
		}
		my $exp_avg = common_scripts->average(\@reathrough_per_sample_in_exp);
		my $exp_stdv = common_scripts->std_dev($exp_avg,\@reathrough_per_sample_in_exp);
		$gene->{_averages}->{$exp_name} = $exp_avg;
		$gene->{_stdvs}->{$exp_name} = $exp_stdv;
	}
}

#report readthrough calculations
my $gene_field_headers =$line_header; 
my $mean_per_exp_header;
my $stdv_per_exp_header;
foreach my $exp_name (keys %exp_hash) { 
	$mean_per_exp_header  .= "\t${exp_name}_mean";
	$stdv_per_exp_header  .= "\t${exp_name}_stdv";
}
my $readthrough_per_sample_header;
for (my $i = 0; $i<=$#experiments; $i++) { 
	my $exp_name = $experiments[$i];
	my $sample_name = $samples[$i];
	my $output_name = $exp_name . "-" . $sample_name;
	$readthrough_per_sample_header .= "\t$output_name";
}
my $raw_cov_gene_and_reg_header;
for (my $i = 0; $i<=$#experiments; $i++) { 
	my $exp_name = $experiments[$i];
	my $sample_name = $samples[$i];
	my $output_name_gene = $exp_name . "-" . $sample_name . "_gene";
	my $output_name_reg = $exp_name . "-" . $sample_name . "_r5utr";
	my $output_readthrough = $exp_name . "-" . $sample_name . "_readthrough";
	$raw_cov_gene_and_reg_header .= "\t$output_name_reg\t$output_name_gene\t$output_readthrough";
}
print $gene_field_headers;
print $mean_per_exp_header;
print $stdv_per_exp_header;
print $readthrough_per_sample_header;
print $raw_cov_gene_and_reg_header;
print "\n";

foreach my $gene (@gene_candidates) { 
	my $line = $gene->{_line};
	print "$line";	
	foreach my $exp_name (keys %exp_hash) {
		my $print_val = $gene->{_averages}->{$exp_name};
		print "\t$print_val";
	}
	foreach my $exp_name (keys %exp_hash) {
		my $print_val = $gene->{_stdvs}->{$exp_name};
		print "\t$print_val";
	}
	for (my $i = 0; $i<=$#experiments; $i++) { 
		my $exp_name = $experiments[$i];
		my $sample_name = $samples[$i];
		my $sample_readthrough = $gene->{_readthrough}->{$exp_name}->{$sample_name}; 
		print "\t$sample_readthrough";
	}
	for (my $i = 0; $i<=$#experiments; $i++) { 
		my $exp_name = $experiments[$i];
		my $sample_name = $samples[$i];
		my $gene_cov = $gene->{_gene_cov}->{$exp_name}->{$sample_name};
		my $reg_cov = $gene->{_reg_cov}->{$exp_name}->{$sample_name};
		my $readthrough = int 100* $gene_cov / ($reg_cov+1);
		$readthrough = 100 if ($readthrough > 100);
		print "\t$reg_cov\t$gene_cov\t$readthrough";
	}	
	print "\n";
}
