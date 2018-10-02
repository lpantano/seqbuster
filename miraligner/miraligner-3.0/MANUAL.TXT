========================
WHAT DO YOU NEED
========================

This tool comes with the version 16.0 of miRBase repository located in DB folder.
If you want to upgrade the files, you can replace them, by files from miRBase repository: 
hairpin.fa and miRNA.str they can be found at: ftp://mirbase.org/pub/mirbase/CURRENT/

========================
QUICK START
========================

java -jar miraligner/miraligner.jar 1 3 3 hsa miraligner-1.0/test/test.fa miraligner-1.0/DB 

Parameter 1: mismatches allowed (0/1)
Parameter 2: nucleotides allowed for trimming (max 3)
Parameter 3: nucleotides allowed for addition (max 3)
Parameter 4: species (3 letter, human=>hsa)
Parameter 5:fasta file
Parameter 6:folder where miRBase files are(one copy at miraligner-1.0/DB folder)
========================
INPUT FORMAT
========================

A fasta file reads:
>seq
CACCGCTGTCGGGGAACCGCGCCAATTT

or tabular file with count information:
CACCGCTGTCGGGGAACCGCGCCAATTT 45

========================
OUTPUT FORMAT
========================

Track file *.mirna.opt: information about the process


Header of the *.mirna.ann file:
*tr5: nucleotides at 5' end different from the annonated sequence in miRBase:

	precursor 	  => cctgtggttagctggttgcatatcc
	annotated miRNA   => ||TGTGGTTAGCTGGTTGCATAT
	sequence tr5: qCC => CCTGTGGTTAGCTGGTTGCATAT
	sequence tr5: tTG =>     TGGTTAGCTGGTTGCATAT

*tr3: nucleotides at 3' end different from the annotated sequence in miRBase:

	precursor         => cctgtggttagctggttgcatatcc
	annotated miRNA   =>   TGTGGTTAGCTGGTTGCATAT||
        sequence tr3: qCC =>   TGTGGTTAGCTGGTTGCATATCC
        sequence tr3: tAT =>   TGTGGTTAGCTGGTTGCAT

*addition: nucleotides at 3' end added
	
	precursor         => cctgtggttagctggttgcatatcc
        annotated miRNA   =>   TGTGGTTAGCTGGTTGCATAT..
        sequence add: qTT =>   TGTGGTTAGCTGGTTGCATATTT

*s5: offset nucleotides at the begining of the annotated miRNAs
	
	precursor         => agcctgtggttagctggttgcatatcc
	annotated miRNA   => ||||TGTGGTTAGCTGGTTGCATAT
	s5                => AGCCTGTG

*s3:offset nucleotides at the ending of the annotated miRNAs

	precursor         =>  cctgtggttagctggttgcatatccgc
	annotated miRNA   =>    TGTGGTTAGCTGGTTGCATAT||||
        s3                =>			 ATATCCGC

*type: mapped on precursor or miRNA sequences

*Example:

seq			miRNA		start	end	mism	tr5	tr3	add	s5	s3	type
TGGCTCAGTTCAGCAGGACC    hsa-mir-24-2    50      67      0       qCC     0       0       0       0       precursor
ACTGCCCTAAGTGCTCCTTCTG  hsa-miR-18a*    47      68      0       0       0       tG      ATCTACTG        CTGGCA  miRNA
 
