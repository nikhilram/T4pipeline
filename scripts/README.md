Perl scripts and modules used to identify the transcription start and termination sites.


## TSS prediction -
[id_TSSs.pl](https://github.com/nikhilram/T4pipeline/blob/master/scripts/id_TSSs.pl) was used to identify the transcription start sites.    

Example -
```
id_TSSs.pl -files *TAP.s2 -fasta NC_003028.v3.17.fasta -annot NC_003028.v3.17.ncrna.genes
```
```id_TSSs.pl -help``` for basic usage. 

This will generate a table with the TSS predictions for each feature in the annotation. If there is a predicted TSS upstream of a feature, the the row includes the TSS position, whether the TSS is intergenic or within the coding region, the processe/unprocessed ratio, processed coverage, and the unprocessed coverage. If there is not predicted TSS with the defined cutoff, the row reads Unknown_TSS.



## TTS prediction -
[id_TTSs.pl](https://github.com/nikhilram/T4pipeline/blob/master/scripts/id_TTSs.pl) was used to identify the transcription termination sites downstream of the annotated genes

Example -
```
id_TTSs.pl -files *.3p -fasta NC_003028.v3.17.fasta -annot NC_003028.v3.17.ncrna.genes
```
```id_TTSs.pl -help``` for basic usage. 

This generates two output files. 1. A sorted bedGraph file of the predicted TTSs and their coverages for easy visualization. 2. A table with all the predicted TTSs, their coverages, 3'UTR-lengths, 40 nucleotides upstream and 20 nucleotides downstream of the TTS, dot-bracket notation of the secondary struture identified in the 40 nucleotides upstream, number of uridines immediately upstream, and the presence of a classical stem-loop structure.
