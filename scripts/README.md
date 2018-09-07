Perl scripts and modules used to identify the transcription start and termination sites.


[id_TSSs.pl](https://github.com/nikhilram/T4pipeline/blob/master/scripts/id_TSSs.pl) was used to identify the transcription start sites.    

```
############################################################################################################
#		
#						Id TSSs
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
#        perl id_TSS.pl -files <.s2> -fasta <.fasta> -annot <.genes>
############################################################################################################
 ```
