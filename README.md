Process to concatenate transcriptomes from Trinity from several species, clean them to keep one transcript per locus (longest protein) using evigene, transdecoder and DIAMOND. Clean "cDNA" from genoimc anotations (NCBI and Mikado de novo annotations).
Run Orthofinder to obtain Orthogroups.
Obtain and align Single-copy orthogroups with all speices, using MAFFT.
Trimm alignments with trimAL.
Run BPP from trimmed, aligned and concatenated Single-copy orthogroups with all speices in Phyllip format.
