Perl scripts and modules used to identify the transcription start and termination sites.


## [id_TSSs.pl](https://github.com/nikhilram/T4pipeline/blob/master/scripts/id_TSSs.pl) ## was used to identify the transcription start sites.    

Example -
```
id_TSSs.pl -files *TAP.s2 -fasta NC_003028.v3.17.fasta -annot NC_003028.v3.17.ncrna.genes
```
```id_TSSs.pl -help``` for basic usage. 



## [id_TTSs.pl](https://github.com/nikhilram/T4pipeline/blob/master/scripts/id_TTSs.pl) was used to identify the transcription termination sites downstream of the annotated genes

Example -
```
id_TTSs.pl -files *.3p -fasta NC_003028.v3.17.fasta -annot NC_003028.v3.17.ncrna.genes
```
```id_TTSs.pl -help``` for basic usage. 
