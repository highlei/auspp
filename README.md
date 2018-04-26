LICENSE
    auspp

    Copyright (c) 2018- Lei Gao

    This program is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    This program is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
    Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program. If not, see <http://www.gnu.org/licenses/>.

SYNOPSIS
    auapp : a universal short-read pre-processing package/program

AUTHOR
    Lei Gao, Shenzhen University, leigao@szu.edu.cn or highlei@hotmail.com

VERSION
    1.0 : April 22, 2018

Easy INSTALL for Linux/Unix OS platform and Perl in /usr/bin/perl
     For administrator/root
	perl MAKEFILE.pl -i /usr/local/bin/
     For other user, put all scripts to your PATH (check your PATH by "echo $PATH"), e.g.
	perl MAKEFILE.pl -i ~/bin/

     Then you can call auspp directly, e.g. test the example:
	auspp -M sRNAexample -e example/

INSTALL -- if "easy INSTALL" doesn't work, please check the following:
  Dependencies - Linux/Unix OS platform
     auspp was developed on linux (Ubuntu 14.04.5 LTS), and hasn't been
     tested on other OS platform.

  Dependencies - Required Perl in /usr/bin/perl
    auspp is a perl program, so it needs perl installed on your system. It
    was developed on perl version 5.14, and hasn't been tested on other
    versions (but there is no reason to suspect problems with other perl
    5.x versions). auspp will not compile. Getopt::Std, FileHandle, strict
    and Cwd 'abs_path' are pre-loaded into most (all?) Perl distros. auspp
    expects to find perl in /usr/bin/perl .. if not, edit line 1 
    accordingly for all the perl scripts in bin/

  Dependencies - PATH executables
	one of the following short read aligners:
            bowtie  and bowtie-build
	or  bowtie2 and bowtie2-build
        or  soap    and 2bwt-builder
	or  hisat2  and hisat2-build
	or  bwa

	samtools (Version 1.0 +, http://www.htslib.org/)

	blast+, if you'd like to use the step 4 (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/)

	bam2wig, if you'd like to use the step 7 (https://github.com/MikeAxtell/bam2wig )

    All of the above must be executable from your PATH. Depending on the
    mode of the auspp run (see below), only a subset of these programs
    may be required for a given run.

  Installation
    Except for the above dependencies, there is no "real" installation. Put
    all the scripts in bin/ into your working directory,  then you can call 
    auspp with

            ./auspp

    For convenience, you can add it to your PATH. e.g.

            sudo mv bin/* /usr/local/bin/


############### run ##################
auspp -i Col.fastq -x Col -D referenceIndex -M mRNA
auspp -i Col.fastq -x Col -D referenceIndex -M smallRNA
auspp -i Col.fastq -x Col -G referenceGENOME -M chip

############### USAGE ##################
auspp includes auspp and 7 modules (other 7 perl scripts searchseq, searchLineACList, trim_adaptor, soap2sam_gl, 
collapseFasta, fastq2fasta, blast_m8). For convenience, you can check their USAGE by "-h", e.g.
            trim_adaptor -h
The 7 modules can be used individually for other purpose.

==========================| auspp  start |=================================
Now = 2018-04-22 08:53:27

Version :   1.0
Author  :   Lei Gao   <highlei@hotmail.com> or <leigao@szu.edu.cn>

Usage:   auspp -i fastq_file -x sampleID -M Modes {-D index | -G genome} [options]
Usage:   auspp -i fastq_file -x sampleID -M degradome [options]

   -i <str>   input the fastq file (Could be gzip'ed (extension: .gz)). eg: Col.fastq or Col_r1.fastq
   -I <str>   input the other mate if paired-end sequencing. eg: Col_r2.fastq
   -x <str>   input the sampleID for -i library. eg: Col
   -D <str>   reference sequence index: soap index or bwa or bowtie(2) or hisat2 index.
   -G <str>   reference sequence/genome in fasta format. (Required when -P soap and step 7.)
   -M <str>   Modes: presets for supported SEQ:
            smallRNA   same as   -P soap -s 124567 -L "20-25;21;22;23;24;All";
            mRNA       same as   -P hisat2 -s 1367 -L All;
            ribo       same as   -P soap -s 12467 -L All;
            chip       same as   -P bowtie -s 167 -L All;
            snp        same as   -P baw mem -s 167 -L All;
            pseudo     same as   -P soap -s 1267 -L All;
            nucleosome same as   -P bowtie2 -s 167 -L All;
            degradome  same as   -s 12
            sRNAexample will run example

   Customized settings by user:
   -s <str>   running step (eg: 1-7 or 1367):
            1 quality control,2 trim,3 collapse,4 filter,5 length,6 mapping,7 GenomeBrowser
   -P <str>   align program: soap or "bwa aln" or "bwa mem" or bowtie(2) or hasat2 or tophat2. [soap]
   -L <str>   the read lenth range. eg: "20-25;21;22;23;24;All" [All]

   Required by special step:
   -R <str>   r/t/sn/snoRNA or repeats or other database in fasta format for filter;
            must be makeblastdb by blast+. Required when step 4 activated
   -a <str>   adaptor sequence. Required when step 2 activated eg: TGGAATTCTCGGG or AAAAAAAAAAA
   -A <str>   adaptor sequence for mate if paired-end. eg: TGGAATTCTCGGG or AAAAAAAAAAA
   -f <str>   gtf file for hisat2 or gff File for tophat2 if have. eg: TAIR10_GFF3_genes.gff.gtf or TAIR10_GFF3_genes_transposons.gff

   Defult settings are recommended:
   -d <str>   new fold for store files. [fasta,trim,filter,map2gnm]
   -S <str>   the parameter set for soap or bwa samse or bowtie(2) or hasat2 or tophat2.
   -T <str>   parameter settings for trim_adaptor ["-l 9 -m 18"]
   -p <str>   the path for all perl scripts. [./]
   -Q <str>   quality control. ["-q 20 -c 5"]
   -C <str>   copy number filter. e.g. "-c 5,10" to discard reads with copy>10 or copy<5. ["-c 1,"]

   -h   display this help

Example:
auspp -i Col.fastq -x Col -M smallRNA -D tair10.Chr.fa.index
auspp -i Col.fastq -x Col -M RNA -G tair10.Chr.fa
auspp -i Col.fastq -x Col -M degradome
==========================| auspp  end   |=================================


############### Modules' descriptions ##################
1. searchseq
   extract seqences from a fasta file according to the sequence name.
2. searchLineACList
   extract line from a txt file
3. trim_adaptor
   trim the adaptor sequence
4. soap2sam_gl
   transfer soap result to sam format
5. collapseFasta
   cluster the same sequences in fasta file into one
6. fastq2fasta
   transfer fastq file to fasta format
7. blast_m8
   filter the blast results (-outfmt 6)

