#!/bin/env perl -w 
use strict;

#######################################   calculate_stemloop_polyU_from_gene_table.pl  ##########################################################
######	Uses the mapped primary 3' ends table and measures the polyU length and whether there is a stem-loop fold
######	run example: perl calculate_stemloop_polyU_from_gene_table.pl termseq_genes_table
#################################################################################################################################################

#command line inputs
my ($term_table) = @ARGV;

open(TERM,$term_table) or die "missing term table $term_table\n";
chomp (my $header = <TERM>);
#$header .= "\tFold\t" . "\tMFE\t" . "\t#ofUs\t" . "stem_loop\n";
$header .= "\tFold\t" . "\tMFE\t" . "\t#ofUs\n";
print $header;
while (<TERM>) { 
	chomp; 
	my @line = split(/\t/); 
	my $table_line = $_;
#	my ($up_seq,$fold) = @line[7,9]; #for all terms  @line[6,8];
	my $term_sequence = $line[15];
	my $strand = $line[3];
	$term_sequence =~ tr/T/U/;
	chomp (my @rna_fold_output = `echo $term_sequence | RNAfold `);### fold sequence using RNAfold
        $rna_fold_output[1] =~ /^(.+)\s+\((.+)\)/;
        my ($fold,$energy) = ($1,$2);
	if ($strand eq "+"){
        my $num_of_Us = &count_U_stretch($term_sequence);
	my $stem_loop = &is_stem_loop($fold);
#       $table_line .= "\t$num_of_Us\t" . "$stem_loop\n";
        $table_line .= "\t$fold\t$energy\t$num_of_Us\n";
        print $table_line;
	} else {
	my $reverse = reverse $term_sequence;
	my $reverse_fold = reverse $fold;
	my $num_of_Us = &count_U_stretch($reverse);
        my $stem_loop = &is_stem_loop($reverse_fold);
#	$table_line .= "\t$num_of_Us\t" . "$stem_loop\n";
	$table_line .= "\t$fold\t$energy\t$num_of_Us\n";
	print $table_line;
}
}




###########################################################################################
############################# Additional Functions ########################################
###########################################################################################

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
	my ($seq,$strand) = @_;
	my $tail = substr($seq,length($seq)-8,8); 
	my @nucs = split(//,$tail); 
	my $number_of_U_residues = 0;
	foreach my $nuc (@nucs) { 
		$number_of_U_residues++ if ($nuc eq 'U');
	}
	return $number_of_U_residues;
}



