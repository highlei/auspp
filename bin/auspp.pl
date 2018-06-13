#!/usr/bin/perl
# Copyright (c)  2014-
# Program:			auspp # a universal short-read pre-processing package/program
# Author:			Gaolei <highlei@hotmail.com or leigao@szu.edu.cn>
# Program Date:		2014.05.01
# Modifier:			Gaolei <highlei@hotmail.com or leigao@szu.edu.cn>
# Last Modified:	2018.03.23
# Description:	pre-processing raw data from NGS and aligning to the reference genome
#**************************
# Version: 1.0
#**************************
# e-mail:highlei@hotmail.com

my $version="1.0";
print STDERR ("\n==========================| $0  start |=================================\n");

my $start = time();
my $Time_Start = sub_format_datetime(localtime(time())); 
print STDERR "Now = $Time_Start\n\n";
my $path=$0;


use Getopt::Std;
getopts("hi:I:o:l:x:P:D:S:L:R:p:s:d:E:a:A:r:Q:C:M:T:F:G:f:u:e:");
my $infile	= $opt_i; # for one library, you can just give the fastq file
my $jnfile	= $opt_I; # for paird-end library, the other mate
#my $inlist	= $opt_I; # for two or more libraries, you should give txt file which include the "fastq file\tprefix" with one row one library.
my $lisfile	= $opt_l;
my $prefix	= (defined $opt_x) ? $opt_x : "";#auspp_$start e.g. Col_r1";
my $outfile	= (defined $opt_o) ? $opt_o : "";#$prefix.".result";
my $soap	= (defined $opt_P) ? $opt_P : "soap"; # or bowtie or bowtie2 or tophat
my $gffFile = (defined $opt_f) ? $opt_f : "";#"~/Ath/data/TAIR10_GFF3_genes_transposons.gff"
my $soapdb	= (defined $opt_D) ? $opt_D : "";#($soap =~/soap/ ? "~/Ath/data/tair10.Chr.fa.index" : ($soap =~/bwa/ ? "~/Ath/data/tair10.Chr.fa" : ($soap =~/hisat2/ ? "~/Ath/data/tair10.Chr_tranTE" : "~/Ath/data/tair10.Chr")));
my $soapset	= (defined $opt_S) ? $opt_S : ($soap =~/soap/ ? "-M 0 -r 2 -v 0" : ($soap =~/bowtie2/ ? "-t -p 3" : 
											($soap=~/bowtie/ ? "-t -k 1 -v 2 --best" : ($soap=~/bwa/ ? "-t 3 -n 100000" : ($soap=~/hisat2/ ? "-p 4" : ($gffFile ne "" ? "-p 3 -G $gffFile" : "-p 3") )) ) )); # --strata
my $soaplen	= (defined $opt_L) ? $opt_L : "";#20-25;21;22;23;24;All";
my $repdb	= (defined $opt_R) ? $opt_R : "";#"~/Ath/data/ath.repeats.fa";#
my $prodir	= (defined $opt_p) ? $opt_p : ($path=~/^(.+\/)[^\/]+pl/? $1 : "./");#"~/pro"; # where all the scripts
my $step	= (defined $opt_s) ? $opt_s : 0; # to do from step x to step y.
my $newdir	= (defined $opt_d) ? $opt_d : "fasta,trim,filter,map2gnm"; #
my $logfile	= (defined $opt_E) ? $opt_E : ""; # output the log file of perl scripts
my $adaptor	= (defined $opt_a) ? $opt_a : "AGATCGGAAGA";#TGGAATTCTCGGG";
my $adaptor2= (defined $opt_A) ? $opt_A : "AGATCGGAAGA";
my $repRM	= (defined $opt_r) ? $opt_r : 0; # default: remove repeats and r/t/sn/snoRNA; 1: donot remove
my $qc		= (defined $opt_Q) ? $opt_Q : "-q 20 -c 5"; # default setting for quality control
my $copyFilter= (defined $opt_C) ? $opt_C : "-c 1,"; # default setting for copy number filter
my $mode	= (defined $opt_M) ? $opt_M : ""; # smallRNA;mRNA;ribo;igv
my $trimSet	= (defined $opt_T) ? $opt_T : "-l 9 -m 18"; # parameter setting for trim_adaptor
my $middleF	= (defined $opt_F) ? $opt_F : 0; # 0: ; 1: skip the step before the 1st step.
my $genome	= (defined $opt_G) ? $opt_G : "";#"~/Ath/data/tair10.Chr.fa"; # reference sequence/genome
my $cleanup	= (defined $opt_u) ? $opt_u : 1;#0: Clean up all intermediate files, 1: retain 
my $example	= (defined $opt_e) ? $opt_e : "./";#where is the example/?(It locates in the directory of source code, e.g. yoursoft/auspp/example)

# fastq2fasta, trim_adaptor, collapseFasta, blast_m8, searchseq, searchLineACList, soap2sam_gl #do_ssearch
# soap or bwa or bowtie or bowtie2 or tophat2, blast+(blastn), bam2wig #ssearch36 (fasta36)
#smallRNA (= -s 2-6 -L 20-25;21;22;24;all);mRNA (= -s 36);ribo (= -s 236);chip (= -s 6); nucleosome (= -s 6)";
my $op1st="";
if ($soap eq "bwa") {
	print STDERR "\nPlease make sure \"bwa aln\" or \"bwa mem\"?\n";
	sub_end_program();
}
if ($mode=~/smallRNA/i) {
	$step   =(defined $opt_s) ? $opt_s : "124567";
	$soaplen=(defined $opt_L) ? $opt_L : "20-25;21;22;23;24;All";
	$soap   =(defined $opt_P) ? $opt_P : "soap";
	$soapset	= (defined $opt_S) ? $opt_S : "-M 0 -r 2 -v 0"; 
} elsif ($mode=~/sRNAexample/i) {
	$step   =(defined $opt_s) ? $opt_s : "124567";
	$soaplen=(defined $opt_L) ? $opt_L : "20-25;21;22;23;24;All";
	$soap   =(defined $opt_P) ? $opt_P : "soap";
	$soapset	= (defined $opt_S) ? $opt_S : "-M 0 -r 2 -v 0"; 
	if (-e "$example/reads/sRNA.fq.gz" && -e "$example/reference/tair10.2M.fa" && -e "$example/lib/ath.RNA.fa") {
		$infile	= "$example/reads/sRNA.fq.gz";
		$genome = "$example/reference/tair10.2M.fa";
		$repdb	= "$example/lib/ath.RNA.fa";
	} else {
		print STDERR "\nPlease tell me where is the example/?(It locates in the directory of source code, e.g. yoursoft/auspp/example/)";
		$example = <STDIN>;
		if (-e "$example/reads/sRNA.fq.gz" && -e "$example/reference/tair10.2M.fa" && -e "$example/lib/ath.RNA.fa") {
			$infile	= "$example/reads/sRNA.fq.gz";
			$genome = "$example/reference/tair10.2M.fa";
			$repdb	= "$example/lib/ath.RNA.fa";
		} else {
			print STDERR "I cann't find the example files in $example\n";
			usage();
		}
	}
	$prefix	= "test.sRNA";
} elsif ($mode=~/RNA/i) {
	$step   =(defined $opt_s) ? $opt_s : "1367";
	$soaplen=(defined $opt_L) ? $opt_L : "All";
	$soap   =(defined $opt_P) ? $opt_P : "hisat2";
	$soapset	= (defined $opt_S) ? $opt_S : "-p 4"; # --strata
} elsif ($mode=~/pseudo/i) {
	$step   =(defined $opt_s) ? $opt_s : "1267";
	$soaplen=(defined $opt_L) ? $opt_L : "All";
	$soap   =(defined $opt_P) ? $opt_P : "soap";
	$soapset	= (defined $opt_S) ? $opt_S : "-M 0 -r 2 -v 0"; 
} elsif ($mode=~/ribo/i) {
	$step   =(defined $opt_s) ? $opt_s : "12467";
	$soaplen=(defined $opt_L) ? $opt_L : "All";
	$soap   =(defined $opt_P) ? $opt_P : "hisat2";
	$soapset	= (defined $opt_S) ? $opt_S : "-p 4"; # --strata
} elsif ($mode=~/chip/i) {
	$step   =(defined $opt_s) ? $opt_s : "167";
	$soaplen="All";
	$soap   =(defined $opt_P) ? $opt_P : "bowtie";
	$soapset	= (defined $opt_S) ? $opt_S : "-t -k 1 -v 2 --best"; 
} elsif ($mode=~/nucleosome|MNase/i) {
	$step   =(defined $opt_s) ? $opt_s : "167";
	$soaplen=(defined $opt_L) ? $opt_L : "All";
	$soap   =(defined $opt_P) ? $opt_P : "bowtie2";
	$soapset	= (defined $opt_S) ? $opt_S : "-t -p 3"; # --strata
} elsif ($mode=~/snp/i) {
	$step   =(defined $opt_s) ? $opt_s : "167";
	$soaplen=(defined $opt_L) ? $opt_L : "All";
	$soap   =(defined $opt_P) ? $opt_P : "bwa mem";
	$soapset	= (defined $opt_S) ? $opt_S : "-t 3 -n 100000"; # --strata
} elsif ($mode=~/degr/i) {
	$step   =(defined $opt_s) ? $opt_s : "12";
	$soaplen=(defined $opt_L) ? $opt_L : "";
	$soap   =(defined $opt_P) ? $opt_P : "";
	$soapset	= (defined $opt_S) ? $opt_S : "";
	$newdir		= (defined $opt_d) ? $opt_d : "fasta,trim";
} elsif ($mode=~/igv/i) {
	$step = ($step == 0) ? "267" : $step;
	$trimSet="-m 12";
#	$soaplen=(defined $opt_L) ? $opt_L : "All";
	$soap   =(defined $opt_P) ? $opt_P : "bwa aln";
	$soapset = (defined $opt_S) ? $opt_S : "-n 100000";
#	$soapdb	= "/mnt/disk2/Ath/data/tair10.Chr.fa";
	$soapdb	= (defined $opt_D) ? $opt_D : "~/Ath/data/tair10.Chr.fa";

} elsif ($mode=~/igvL/i) {
	$step = ($step == 0) ? "367" : $step;
#	$trimSet="-m 12";
#	$soaplen=(defined $opt_L) ? $opt_L : "All";
	$soapset = (defined $opt_S) ? $opt_S : "-n 100000";
	$soap    = (defined $opt_P) ? $opt_P : "bwa mem";
#	$soapdb	= "/mnt/disk2/Ath/data/tair10.Chr.fa";
	$soapdb	= (defined $opt_D) ? $opt_D : "~/Ath/data/tair10.Chr.fa";
} elsif ($mode ne "") {
	print STDERR "Please ennter a correct mode:\t$mode\n";
	print "\t-M <str>	Current recommended preset for some SEQ: 
	\t\t\tsmallRNA\tsame as\t-P soap -s 124567 -L \"20-25;21;22;23;24;All\";
	\t\t\tmRNA\t\tsame as\t-P hisat2 -s 1367 -L All;
	\t\t\tribo\t\tsame as\t-P soap -s 12467 -L All;
	\t\t\tchip\t\tsame as\t-P bowtie -s 167 -L All; 
	\t\t\tsnp\t\tsame as\t-P baw mem -s 167 -L All;
	\t\t\tpseudo\t\tsame as\t-P soap -s 1267 -L All;
	\t\t\tnucleosome\tsame as\t-P bowtie2 -s 167 -L All;
	\t\t\tdegradome\tsame as\t-s 12;
	\t\t\tsRNAexample\twill run example";
	print " \n";
	sub_end_program();
} 
$op1st .= "Current Settings:\n\tmode=$mode\t$soap\t$soapset\t$soapdb\n";


if ($opt_h || ($infile eq "" && $prefix eq "" && $mode !~/^sRNAexample/i)) {# && $inlist eq ""
	usage();
}
use FileHandle;
use strict;
use Cwd 'abs_path';

my $head = 100;
my $tail = 100;

my ($i,$j,$k,$m,$n,$k1,$k2,$k3,$file,$line,$count,$flag,$block,$a,$b,$end);
my (@buf,@tmp,@array,@ndir,@nlen,@rnum,@rstep,@rdir);
my $ii=0;
my $jj=0;
my $kk=0;
my %snp=();
my %gnm=();
my %cro=();
my $key="";
my $k4 ="";



#===========================================================================================================
#====================                  main
#===========================================================================================================
$rstep[0]="";
for ($k1 = 1; $k1 <= 7 ;$k1++) {
	$rstep[$k1]=0;
}
if ($step == 0) {
	$step = "1-7";
}
if ($step =~/(\d+)-(\d+)/) {
#	@buf = split(/-/,$step);
	for ($k1 = $1; $k1 <= $2 ;$k1++) {
		$rstep[$k1]=$k1;
	}
} elsif ($step =~/(\d+)\D+(\d+)/) {
	@buf = split(/\D+/,$step);
	for ($k1 = 0; $k1 < @buf ;$k1++) {
		$rstep[$buf[$k1]]=$buf[$k1];
	}
} elsif ($step=~/(\d+)/) {
	@buf = split(//,$step);
	for ($k1 = 0; $k1 < @buf ;$k1++) {
		$rstep[$buf[$k1]]=$buf[$k1];
	}
#	$rstep[$1] = $1;
}
if ($rstep[7]!=0 && $soap eq "soap" && $genome eq "") {
	print STDERR "\nPlease provide reference sequence by setting the parameter -G!\n";
	sub_end_program();
}
if ($rstep[4]!=0 && $repdb eq "") {
	print STDERR "\nPlease provide r/t/sn/snoRNA or repeats or other database in fasta format for filter by setting the parameter -R! ; must be formatdb by blast+!\n";
	sub_end_program();
} elsif ($rstep[4]!=0) {
	if ((-e "$repdb.nhr") && (-e "$repdb.nin") && (-e "$repdb.nsq")) {
	} else {
		$i=readpipe("makeblastdb -in $repdb -dbtype nucl");
	}
	if ((-e "$repdb.nhr") && (-e "$repdb.nin") && (-e "$repdb.nsq")) {
	} else {
		print STDERR "\nPlease provide r/t/sn/snoRNA or repeats or other database in fasta format for filter by setting the parameter -R! ; must be formated by makeblastdb of blast+!\n";
		print STDERR "e.g.: makeblastdb -in $repdb -dbtype nucl\n";
		sub_end_program();
	}
}

if ($newdir ne "") {
#	$newdir .= ",$outfile";
	$i = "";
	@ndir = split(",",$newdir); #"fasta,trim,filter,soap";
	@rdir = split(",",$newdir); #"fasta,trim,filter,soap";
	for ($k1 = 0; $k1 < @ndir ;$k1++) {
		if (-e $ndir[$k1]) {
			$rdir[$k1] = 0;
		} else {
			mkdir("$ndir[$k1]"); $rdir[$k1] = 1;
			$i.=$ndir[$k1].", ";
		}
	}
#	if ($i ne "") {
		$op1st .= "\tMake new fold: $i\n";
#	}
} else {
#	$newdir .= "$outfile";
	$i = "";
	$ndir[0] = $newdir; #"fasta,trim,filter,soap";
	$rdir[0] = $newdir; #"fasta,trim,filter,soap";
	for ($k1 = 0; $k1 < @ndir ;$k1++) {
		if (-e $ndir[$k1]) {
			$rdir[$k1] = 0;
		} else {
			mkdir("$ndir[$k1]"); $rdir[$k1] = 1;
			$i.=$ndir[$k1].", ";
		}
	}
#	if ($i ne "") {
		$op1st .= "\tMake new fold: $i\n";
#	}
}

if ($soaplen =~/[\;|\,]/ ) {
	@nlen = split(/\;|\,/,$soaplen);
	$op1st .= "\tLength range: @nlen\n";
} elsif ($soaplen ne "") {
	$nlen[0] = $soaplen;
	$op1st .= "\tLength range: @nlen\n";
} else {
	$nlen[0] = "All";
	$op1st .= "\tLength range: @nlen\n";
}

print STDERR $op1st;
#$i = substr(join(",",@rstep),1);
print STDERR "\tStep range [123456]: @rstep\n";
#sub_end_program();

print STDERR "\nRunning:";
if ($infile =~/\.gz$/) {
	$i = $infile; $i =~s/[\\|\/|\s]/_/g; $i=~s/\.gz$//;
	$i = ".auspptmp.".$i; print STDERR "\ngzip -d output = $i\n";
	system("gzip -d -c $infile > $i");
	$infile = $i; $ii = 1;
}
if ($jnfile =~/\.gz$/) {
	$i = $jnfile; $i =~s/[\\|\/|\s]/_/g; $i=~s/\.gz$//;
	$i = ".auspptmp.".$i; print STDERR "\ngzip -d output = $i\n";
	system("gzip -d -c $jnfile > $i");
	$jnfile = $i; $jj = 1;
}

if ($rstep[6]!=0 && $soapdb eq "" && $genome ne "") {
	$j = abs_path("$genome");
	$j =~/(.+\/)([^\/]+)$/;
#	if (-W "$1") {
#	} else {
#		print STDERR "\nSince the directory of $genome isn't writable, the $genome file will be copy to the current directory!\n";
		system("cp $genome .$2.auspptmp.fa");
		$k4 = $genome;
		$genome = ".".$2.".auspptmp.fa";
#	}
	$soapdb = $genome;
	if ($soap =~/soap/) {
		$i=readpipe("2bwt-builder $genome"); $soapdb = $genome.".index";
	} elsif ($soap =~/hisat2/) {
		if ($gffFile ne "") {
			$i=readpipe("extract_splice_sites.py $gffFile > $gffFile.auspp.ss");
			$i=readpipe("extract_exons.py $gffFile > $gffFile.auspp.exon");
			$i=readpipe("hisat2-build -p 4 $genome --ss $gffFile.auspp.ss --exon $gffFile.auspp.exon $genome");
		} else {
			$i=readpipe("hisat2-build $genome $genome");
		}
	} elsif ($soap =~/bowtie2/) {
		$i=readpipe("bowtie2-build $genome $genome");
	} elsif ($soap =~/bowtie/) {
		$i=readpipe("bowtie-build $genome $genome");
	} elsif ($soap =~/bwa/) {
		$i=readpipe("baw index $genome");
	} elsif ($soap =~/tophat2/) {
		$i=readpipe("bowtie2-build $genome $genome");
	} 
}

($k,$j) = sub_raw2align($infile,$jnfile,$prefix,0);

if ($k4 ne "") {
	$genome.="*";
	system("rm $genome");
}
if ($ii == 1) {
	system("rm $infile");
}
if ($jj == 1) {
	system("rm $jnfile");
}
print STDERR "Info:\t$k\nSamples:\t$prefix\n";

if ($rstep[1]!=0) {
	print STDERR "Total\t$rnum[0][1][0]\nQualityControl\t$rnum[0][1][1]\n";
} else {
	print STDERR "Total\t$rnum[0][1][0]\n";
}
if ($rstep[2]!=0) {
	print STDERR "beforeTrim\t$rnum[0][2][0]\nafterTrim\t$rnum[0][2][1]\n";
}
if ($rstep[3]!=0) {
	print STDERR "Collapse\t$rnum[0][3][1]\nunique\t$rnum[0][3][0]\n";
}
if ($rstep[4]!=0) {
	print STDERR "Filter\t$rnum[0][4][0]\nunique\t$rnum[0][4][1]\n";
}
if ($rstep[5]!=0) {
	for ($k1=0; $k1 < @nlen ;$k1++) {
		if ($nlen[$k1]=~/(\d+)/) {
			print STDERR "Length_$nlen[$k1]\t$rnum[0][5][$k1][0]\nunique\t$rnum[0][5][$k1][1]\n";
		}
	}
}
if ($rstep[6]!=0) {
	for ($k1=0; $k1 < @nlen ;$k1++) {
		if ($nlen[$k1]=~/(\w+)/) {
			print STDERR "Map_$nlen[$k1]\t$rnum[0][6][$k1][0]\nunique\t$rnum[0][6][$k1][1]\n";
		}
	}
}
#for ($k1 = 0; $k1 < $k ;$k1++) {
#	print STDERR $rnum[0][$k1],"\n";
#}
#if ($outfile ne "") {
#	if ($rdir[-1] == 1) {
#		system("mv $ndir[-1] $outfile")
#	}
#}
if (@ndir > 1 && $cleanup == 0) {
	for ($k1 = 0; $k1 < @ndir-1 ;$k1++) {
		if ($rdir[$k1] == 1) {
			system("rm -fr $ndir[$k1]");
		}
	}
}


if ($logfile ne "") {
	open(LOGF, ">$logfile") || die("Can not open logfile: $logfile\n");
	print LOGF "\nLog: \n$j\n";
	close(LOGF);
}

sub_end_program();


#############################################################################################################
####################################                                         ################################
####################################              "main end"                 ################################
####################################                                         ################################
#############################################################################################################
sub usage
{
#	print "Program :\t$0\n";
   print "Version :   $version\n";
   print "Author  :   Lei Gao   <highlei\@hotmail.com> or <leigao\@szu.edu.cn>\n";
   print "\nUsage:   $0 -i fastq_file -x sampleID -M Modes {-D index | -G genome} [options]";
   print "\nUsage:   $0 -i fastq_file -x sampleID -M degradome [options]\n\n";
   print "   -i <str>   input the fastq file (Could be gzip'ed (extension: .gz)).";
   print " eg: Col.fastq or Col_r1.fastq\n";
   print "   -I <str>   input the other mate if paired-end sequencing.";
   print " eg: Col_r2.fastq\n";
#   print "   -I <str>   input the file include two or more libraries.";
#   print " eg: library.txt\n";
   print "   -x <str>   input the sampleID for -i library.";
   print " eg: Col\n";
   print "   -D <str>   reference sequence index: soap index or bwa or bowtie(2) or hisat2 index.";
   print " \n";#tair10.Chr.fa.index
   print "   -G <str>   reference sequence/genome in fasta format. (Required when -P soap and step 7.)";
   print " \n";#tair10.Chr.fa
   print "   -M <str>   Modes: presets for supported SEQ: 
            smallRNA   same as   -P soap -s 124567 -L \"20-25;21;22;23;24;All\";
            mRNA       same as   -P hisat2 -s 1367 -L All;
            ribo       same as   -P soap -s 12467 -L All;
            chip       same as   -P bowtie -s 167 -L All; 
            snp        same as   -P baw mem -s 167 -L All;
            pseudo     same as   -P soap -s 1267 -L All;
            nucleosome same as   -P bowtie2 -s 167 -L All;
            degradome  same as   -s 12;
            sRNAexample will run example";
   print " \n";

   print "\n   Customized settings by user:\n";
   print "   -s <str>   running step (eg: 1-7 or 1367):
            1 quality control,2 trim,3 collapse,4 filter,5 length,6 mapping,7 GenomeBrowser\n";
#   print "";
   print "   -P <str>   align program: soap or \"bwa aln\" or \"bwa mem\" or bowtie(2) or hasat2 or tophat2.";
   print " [$soap]\n";
   print "   -L <str>   the read lenth range.";
   print " eg: \"20-25;21;22;23;24;All\" [All]\n";

   print "\n   Required by special step:\n";
   print "   -R <str>   r/t/sn/snoRNA or repeats or other database in fasta format for filter;
            must be makeblastdb by blast+. Required when step 4 activated";
   print " \n";
   print "   -a <str>   adaptor sequence. Required when step 2 activated. e.g.TGGAATTCTCGGG.";
   print " [$adaptor]\n";
   print "   -A <str>   adaptor sequence for mate if paired-end. e.g. AAAAAAAAAAA.";
   print " [$adaptor]\n";
   print "   -f <str>   gtf file for hisat2 or gff File for tophat2 if have.";
   print " eg: TAIR10_GFF3_genes.gff.gtf\n";# or TAIR10_GFF3_genes_transposons.gff\n";
   print "   -e <str>   where is the example/?(It locates in the directory of source code). Required when mode=sRNAexample.";
   print " eg: yoursoft/auspp/example/\n";# or TAIR10_GFF3_genes_transposons.gff\n";
#   print "   -o <str>   name of directory for output (will be created if doesn't exist).";
#   print " [$outfile]\n";
#   print "   -u <int>   0: clean up all intermediate files. 1: retain";
#   print " [$cleanup]\n";
#   print "   -r   <int>   0: remove repeats and r/t/sn/snoRNA; 1: donot remove.";
#   print " [$repRM]\n";


   print "\n   Defult settings are recommended:\n";
   print "   -d <str>   name of directory for store intermediate/output files.";
   print " [$newdir]\n";
   print "              Normally, fasta/: output of 1st step; trim/: output of 2nd and 3rd step; filter/: output of 4th and 5th; map2gnm/: output of 6th step; map2gnm/bam2wigM/: output of 7th step\n";
   print "   -S <str>   the parameter set for soap or bwa samse or bowtie(2) or hasat2 or tophat2.";
   print " \n";
   print "   -T <str>   parameter settings for trim_adaptor";
   print " [\"$trimSet\"]\n";
#   print "   -p <str>   the path for all perl scripts.";
 #  print " [$prodir]\n";
   print "   -Q <str>   quality control.";
   print " [\"$qc\"]\n";
   print "   -C <str>   copy number filter. e.g. \"-c 5,10\" to discard reads with copy>10 or copy<5.";
   print " [\"$copyFilter\"]\n";
   print "\n   -h   display this help\n";
   print "\nExample:\n";
   print "$0 -i Col.fastq -x Col -M smallRNA -D tair10.Chr.fa.index\n";
   print "$0 -i Col.fastq -x Col -M RNA -G tair10.Chr.fa\n";
   print "$0 -i Col.fastq -x Col -M degradome\n";
#   print "$0 -I library.txt -P $soap -D $soapdb -S \"$soapset\" -L $soaplen -R $repdb -p $prodir\n";
   print "\nDocumentation:\n\tperldoc auspp\n\n";
   print ("==========================| $0  end   |=================================\n\n");

    exit(0);
}
############################################################################################################
######################                  sub_format_datetime
############################################################################################################

sub sub_format_datetime 
{
    my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

############################################################################################################
######################                  sub_end_program
############################################################################################################
sub sub_end_program
{
	print STDERR ("\n............................................................\n");
	my $Time_End = sub_format_datetime(localtime(time()));
	print STDERR "Running from [$Time_Start] to [$Time_End]\n";
	$end = time();
	printf STDERR ("Total execute time : %.2f s\n",$end-$start);
	print STDERR ("==========================| $0  end   |=================================\n\n");
	exit(0);

}

############################################################################################################
######################                  sub_raw2align
############################################################################################################

sub sub_raw2align {
	my ($rin,$rjn,$rfix,$nlib) =@_;
	my @rerr;
	my @drawlen;
	my ($rr,$r1,$r2,$r3,$rs,$r4,$r5,$rt1,$rt2,$r11,$r21); $rr=""; $r11=0;# $r4=0; # $r4 minimum step
	my @lendist;
	for ($rs =1; $rs <= 6 ;$rs++) {
		if ($rstep[$rs] != 0) {
			$r4=$rs; last;
		}
	}
	for ($rs =1; $rs <= 6 ;$rs++) {
		if ($rstep[$rs] != 0) {
			$r5=$rs;
		}
	}
# step 1: quality control
	$rs = 1; print STDERR "\n"; #middleF=$middleF\n
	if ($rstep[$rs] != 0) {
		print STDERR sub_format_datetime(localtime(time())),"\tStep $rs in Lib $nlib, $rfix: Quality control\n";
		$r1=readpipe("fastq2fasta -f 0 -p \"\" -F 1 $qc -i $rin -o $ndir[0]/$rfix.fasta 2>&1");
		$rr.=$r1; $r1=~/OK\!\s(\d+),(\d+),(\d+)/;	$rnum[$nlib][$rs][0] = $1;$rnum[$nlib][$rs][1] = $2;
		if ($rjn ne "") {
			$r1=readpipe("fastq2fasta -f 0 -p \"\" -F 1 $qc -i $rjn -o $ndir[0]/$rfix.r2.fasta 2>&1");
			$rr.=$r1; $r1=~/OK\!\s(\d+),(\d+),(\d+)/;	$rnum[$nlib][$rs][0] .= ",".$1;$rnum[$nlib][$rs][1] .= ",".$2;
		}
	#} elsif (-e "$ndir[0]/$rfix.fasta") {
	#	print STDERR sub_format_datetime(localtime(time())),"\tStep $rs in Lib $nlib, $rfix: Convert fastq to fasta\n";
	#	$r1=readpipe("fastq2fasta -f 0 -p \"\" -F 0 -i $rin -o $ndir[0]/$rfix.fasta 2>&1");
	#	$rr.=$r1; $r1=~/OK\!\s(\d+),(\d+),(\d+)/;	$rnum[$nlib][$rs][0] = $1;$rnum[$nlib][$rs][1] = $2;
	} elsif (($r4 < $rs || !(-e "$ndir[0]/$rfix.fasta")) && $rs < $r5) {
	#	print STDERR sub_format_datetime(localtime(time())),"\tStep $rs in Lib $nlib, $rfix: Convert fastq to fasta\n";
		$r1=readpipe("fastq2fasta -f 0 -p \"\" -F 0 -i $rin -o $ndir[0]/$rfix.fasta 2>&1");
		$rr.=$r1; $r1=~/OK\!\s(\d+),(\d+),(\d+)/;	$rnum[$nlib][$rs][0] = $1;$rnum[$nlib][$rs][1] = $2;
		if ($rjn ne "") {
			$r1=readpipe("fastq2fasta -f 0 -p \"\" -F 0 -i $rjn -o $ndir[0]/$rfix.r2.fasta 2>&1");
			$rr.=$r1; $r1=~/OK\!\s(\d+),(\d+),(\d+)/;	$rnum[$nlib][$rs][0] .= ",".$1;$rnum[$nlib][$rs][1] .= ",".$2;
		}
	}
# step 2: trim adaptor
	$rs=2; 
	if ($rstep[$rs] != 0) {
		print STDERR sub_format_datetime(localtime(time())),"\tStep $rs in Lib $nlib, $rfix: Trim adaptor\n";
		$r1=readpipe("trim_adaptor $trimSet -a $adaptor -I $ndir[0]/$rfix.fasta -o $ndir[1]/ 2>&1");
		$rr.=$r1; $r1=~/Total:\s(\d+).+Good:\s(\d+)/;	$rnum[$nlib][$rs][0] = $1;$rnum[$nlib][$rs][1] = $2;
		if ($rjn ne "") {
			$r1=readpipe("trim_adaptor $trimSet -a $adaptor2 -I $ndir[0]/$rfix.r2.fasta -o $ndir[1]/ 2>&1");
			$rr.=$r1; $r1=~/Total:\s(\d+).+Good:\s(\d+)/;	$rnum[$nlib][$rs][0] .= ",".$1;$rnum[$nlib][$rs][1] .= ",".$2;
		}
	#} elsif (-e "$ndir[1]/$rfix.trim.fasta") {
	} elsif ($r4 > $rs && $middleF == 1) {
	} elsif (($r4 < $rs || !(-e "$ndir[1]/$rfix.trim.fasta")) && $rs < $r5) {
#		print STDERR sub_format_datetime(localtime(time())),"\tStep $rs in Lib $nlib, $rfix: Trim adaptor\n";
		$r1=readpipe("cp $ndir[0]/$rfix.fasta $ndir[1]/$rfix.trim.fasta");
		$rr.=$r1;# $r1=~/Total:\s(\d+).+Good:\s(\d+)/; $rnum[$nlib][$r4++] = $1;$rnum[$nlib][$r4++] = $2;
		if ($rjn ne "") {
			$r1=readpipe("cp $ndir[0]/$rfix.r2.fasta $ndir[1]/$rfix.r2.trim.fasta");
			$rr.=$r1;# $r1=~/Total:\s(\d+).+Good:\s(\d+)/; $rnum[$nlib][$r4++] = $1;$rnum[$nlib][$r4++] = $2;
		}

	}
	if ($rjn ne "") {
		system("cat $ndir[1]/$rfix.r2.trim.fasta >> $ndir[1]/$rfix.trim.fasta");
	}
# step 3: collapse Fasta
	$rs=3; 
	if ($rstep[$rs] != 0) {
		print STDERR sub_format_datetime(localtime(time())),"\tStep $rs in Lib $nlib, $rfix: Collapse reads\n";
		$r1=readpipe("collapseFasta -i $ndir[1]/$rfix.trim.fasta $copyFilter -0 0 -o $ndir[1]/$rfix.fa 2>&1");
		$rr.=$r1; $r1=~/(\d+),(\d+)\sunique,sequence/;	$rnum[$nlib][$rs][0] = $1;$rnum[$nlib][$rs][1] = $2;
		system("perl -e \'my \%a;my \$i=0;while(<>){if(\$_=~/^>\\S+_(\\d+)x/){\$j=\$1;} else {\$_=~s/[\\r|\\n|\\s]+\$//g;\$k=length(\$_);if(exists(\$a{\$k})){\$a{\$k}+=\$j} else{\$a{\$k}=\$j}} } \$m=0;print \"\\t$rfix\\n\"; foreach \$k (sort keys(\%a)) {\$m+=\$a{\$k};print \"\$k\\t\$a{\$k}\\n\"} print \"Total\\t\$m\\n\"\' $ndir[1]/$rfix.fa > $ndir[1]/$rfix.lenDist");
		$r1=readpipe("perl -e \'my \%a;my \$i=0;while(<>){if(\$_=~/^>\\S+_(\\d+)x/){\$j=\$1;} else {\$_=~s/[\\r|\\n|\\s]+\$//g;\$k=length(\$_);if(exists(\$a{\$k})){\$a{\$k}+=\$j} else{\$a{\$k}=\$j}} } \$m=0; print \"Length\\tAfter_Trim\\n\"; foreach \$k (sort keys(\%a)) {\$m+=\$a{\$k};print \"\$k\\t\$a{\$k}\\n\"} \' $ndir[1]/$rfix.fa");
		$lendist[$r11++]={split(/[\t|\n]/,$r1)};
	#} elsif (-e "$ndir[1]/$rfix.fa") {
	} elsif ($r4 > $rs && $middleF == 1) {
	} elsif (($r4 < $rs || !(-e "$ndir[1]/$rfix.fa")) && $rs < $r5) {
	#	print STDERR sub_format_datetime(localtime(time())),"\tStep $rs in Lib $nlib, $rfix: Cluster reads\n";
		$r1=readpipe("awk \'{if(\$_!~/^>/){print \">\" \$1 \"_1_1x\\n\" \$1}}\' $ndir[1]/$rfix.trim.fasta > $ndir[1]/$rfix.fa");
		$rr.=$r1; #$r1=~/(\d+),(\d+)\sunique,sequence/;$rnum[$nlib][$r4++] = $2;
		system("perl -e \'my \%a;my \$i=0;while(<>){if(\$_=~/^>\\S+/){\$j=1;} else {\$_=~s/[\\r|\\n|\\s]+\$//g;\$k=length(\$_);if(exists(\$a{\$k})){\$a{\$k}+=\$j} else{\$a{\$k}=\$j}} } \$m=0;print \"\\t$rfix\\n\"; foreach \$k (sort keys(\%a)) {\$m+=\$a{\$k};print \"\$k\\t\$a{\$k}\\n\"} print \"Total\\t\$m\\n\"\' $ndir[1]/$rfix.fa > $ndir[1]/$rfix.lenDist");
		$r1=readpipe("perl -e \'my \%a;my \$i=0;while(<>){if(\$_=~/^>\\S+/){\$j=1;} else {\$_=~s/[\\r|\\n|\\s]+\$//g;\$k=length(\$_);if(exists(\$a{\$k})){\$a{\$k}+=\$j} else{\$a{\$k}=\$j}} } \$m=0;print \"Length\\tAfter_Trim\\n\"; foreach \$k (sort keys(\%a)) {\$m+=\$a{\$k};print \"\$k\\t\$a{\$k}\\n\"} \' $ndir[1]/$rfix.fa");
		$lendist[$r11++]={split(/[\t|\n]/,$r1)};print STDERR "\nr11=$r11,r1=$r1 end\n";
	}
# step 4: remove repeats with local database
	$rs=4; 
	if ($rstep[$rs] != 0) {
		print STDERR sub_format_datetime(localtime(time())),"\tStep $rs in Lib $nlib, $rfix: Filter reads similar to repeats or/and r/t/sn/snoRNAA or/and other\n";

		system("blastn -db $repdb -query $ndir[1]/$rfix.fa -task \'blastn\' -num_alignments 1 -outfmt 6 -dust no -num_threads 3 -out $ndir[2]/$rfix.VS.filter.blastn");
#blastn -db /home/lgao/Ath/data/ath.repeats.fa -query trim/col_r1.fa -num_alignments 1 -outfmt 6 -dust no -num_threads 3 -out tmp.t.blastn -task 'blastn'
#		system("blastall -p blastn -d $repdb -i $ndir[1]/$rfix.fa -v 1 -b 1 -m 8 -F F -a 3 -o $ndir[2]/$rfix.VS.filter.blastn");

		@rerr=readpipe("blast_m8 -i $ndir[2]/$rfix.VS.filter.blastn -d $ndir[1]/$rfix.fa -r 0.8 -t 80 -L 10 -m 10 -N 4 -R 2 -o $ndir[2]/$rfix.VS.filter.blastn.lis 2>/dev/null");
		$r1=join("\n",@rerr); $rr.=$r1;
#		@rerr=readpipe("searchseq1.7s -i $ndir[1]/$rfix.fa -l $ndir[2]/$rfix.VS.filter.blastn.lis > $ndir[2]/$rfix.VS.filter.blastn.lis.fasta 2>/dev/null");
	#	$r1=join("\n",@rerr); $rr.=$r1;
#		@rerr=readpipe("do_ssearch -i $ndir[2]/$rfix.VS.filter.blastn.lis.fasta -d /media/disk1/miRBase/R20/plant.mature.r20.fa -r 0 -t 10 -l 1 -L 8 -s 1 -D 2 -o $ndir[2]/$rfix.VS.filter.blastn.lis.VS.plant.mature.r20.match 2>&1");
#		$r1=join("\n",@rerr); $rr.=$r1;
#		@rerr=readpipe("searchLineACList1.2s -i $ndir[2]/$rfix.VS.filter.blastn.lis -l $ndir[2]/$rfix.VS.filter.blastn.lis.VS.plant.mature.r20.match -r n > $ndir[2]/$rfix.VS.filter.blastn.lis.no_miRBase_r20.lis 2>/dev/null");
	#	$r1=join("\n",@rerr); $rr.=$r1;
#		@rerr=readpipe("searchseq1.7s -i $ndir[1]/$rfix.fa -l $ndir[2]/$rfix.VS.filter.blastn.lis.no_miRBase_r20.lis -r n > $ndir[2]/$rfix.filter.fa 2>/dev/null");
	#	$r1=join("\n",@rerr); $rr.=$r1;
		@rerr=readpipe("searchseq -i $ndir[1]/$rfix.fa -l $ndir[2]/$rfix.VS.filter.blastn.lis -r n > $ndir[2]/$rfix.filter.fa 2>/dev/null");
		$r1=join("\n",@rerr); $rr.=$r1;

	#	$r1=readpipe("awk -F \"_\" \'BEGIN{k=0;j=0}{if(\$_~/^>[]/){gsub(\"x\",\"\",\$3);k+=\$3;j++}else{} }END{print k,j}\' $ndir[2]/$rfix.filter.fa");
		$r1=readpipe("perl -e \'\$k=0;\$j=0;while(<>){if(\$_=~/^>\\S+_(\\d+)x/){\$k+=\$1;\$j++} elsif(\$_=~/^>\\S+/){\$k++;\$j++;} } print \"\$k,\$j\"\' $ndir[2]/$rfix.filter.fa");
		$r1=~/(\d+),(\d+)/;	$rnum[$nlib][$rs][0] = $1;$rnum[$nlib][$rs][1] = $2; # remove repeats
		system("perl -e \'my \%a;my \$i=0;while(<>){if(\$_=~/^>\\S+_(\\d+)x/){\$j=\$1;} elsif(\$_=~/^>\\S+/){\$j=1;} else {\$_=~s/[\\r|\\n|\\s]+\$//g;\$k=length(\$_);if(exists(\$a{\$k})){\$a{\$k}+=\$j} else{\$a{\$k}=\$j}} } \$m=0;print \"\\t$rfix\\n\"; foreach \$k (sort keys(\%a)) {\$m+=\$a{\$k};print \"\$k\\t\$a{\$k}\\n\"} print \"Total\\t\$m\\n\"\' $ndir[2]/$rfix.filter.fa > $ndir[2]/$rfix.filter.lenDist");
		$r1=readpipe("perl -e \'my \%a;my \$i=0;while(<>){if(\$_=~/^>\\S+_(\\d+)x/){\$j=\$1;} elsif(\$_=~/^>\\S+/){\$j=1;} else {\$_=~s/[\\r|\\n|\\s]+\$//g;\$k=length(\$_);if(exists(\$a{\$k})){\$a{\$k}+=\$j} else{\$a{\$k}=\$j}} } \$m=0;print \"Length\\tAfter_Filter\\n\"; foreach \$k (sort keys(\%a)) {\$m+=\$a{\$k};print \"\$k\\t\$a{\$k}\\n\"} \' $ndir[2]/$rfix.filter.fa");
		$lendist[$r11++]={split(/[\t|\n]/,$r1)};
	#} elsif (-e "$ndir[2]/$rfix.filter.fa") {
	} elsif ($r4 > $rs && $middleF == 1) {
	} elsif (($r4 < $rs || !(-e "$ndir[2]/$rfix.filter.fa")) && $rs < $r5) {
		@rerr=readpipe("cp $ndir[1]/$rfix.fa $ndir[2]/$rfix.filter.fa");
		system("perl -e \'my \%a;my \$i=0;while(<>){if(\$_=~/^>\\S+_(\\d+)x/){\$j=\$1;} elsif(\$_=~/^>\\S+/){\$j=1;} else {\$_=~s/[\\r|\\n|\\s]+\$//g;\$k=length(\$_);if(exists(\$a{\$k})){\$a{\$k}+=\$j} else{\$a{\$k}=\$j}} } \$m=0;print \"\\t$rfix\\n\"; foreach \$k (sort keys(\%a)) {\$m+=\$a{\$k};print \"\$k\\t\$a{\$k}\\n\"} print \"Total\\t\$m\\n\"\' $ndir[2]/$rfix.filter.fa > $ndir[2]/$rfix.filter.lenDist");
		$r1=readpipe("perl -e \'my \%a;my \$i=0;while(<>){if(\$_=~/^>\\S+_(\\d+)x/){\$j=\$1;} elsif(\$_=~/^>\\S+/){\$j=1;} else {\$_=~s/[\\r|\\n|\\s]+\$//g;\$k=length(\$_);if(exists(\$a{\$k})){\$a{\$k}+=\$j} else{\$a{\$k}=\$j}} } \$m=0;print \"Length\\tAfter_Filter\\n\"; foreach \$k (sort keys(\%a)) {\$m+=\$a{\$k};print \"\$k\\t\$a{\$k}\\n\"} \' $ndir[2]/$rfix.filter.fa");
		$lendist[$r11++]={split(/[\t|\n]/,$r1)};
	}
# step 5: length of 20-25bp
	$rs=5; 
	if ($rstep[$rs] != 0) {
		print STDERR sub_format_datetime(localtime(time())),"\tStep $rs in Lib $nlib, $rfix: Running length select: ";
		for ($r2 = 0; $r2 < @nlen ;$r2++) {
			if ($nlen[$r2]=~/(\d+)\D+(\d+)/) {
				print STDERR "$1-$2;";
				`awk '{if(\$_~/^>/){i=\$1}else{if(length(\$1)>=$1 && length(\$1)<=$2){print i "\\n" \$1} }}' $ndir[2]/$rfix.filter.fa > $ndir[2]/$rfix.filter.$1-$2.fa`;
			#	$r1=readpipe("awk -F \"_\" \'BEGIN{k=0;j=0}{if(\$_~/^>/){gsub(\"x\",\"\",\$3);k+=\$3;j++}}END{print k,j}\' $ndir[2]/$rfix.filter.$1-$2.fa");
				$r1=readpipe("perl -e \'\$k=0;\$j=0;while(<>){if(\$_=~/^>\\S+_(\\d+)x/){\$k+=\$1;\$j++} elsif(\$_=~/^>\\S+/){\$k++;\$j++;} } print \"\$k,\$j\"\' $ndir[2]/$rfix.filter.$1-$2.fa");
				$r1=~/(\d+),(\d+)/;	$rnum[$nlib][$rs][$r2][0] = $1;$rnum[$nlib][$rs][$r2][1] = $2;#$rnum[$nlib][$rs][$r2][3] = $nlen[$r2];

			} elsif ($nlen[$r2]=~/(\d+)/) {
				print STDERR "$1;";
				`awk '{if(\$_~/^>/){i=\$1}else{if(length(\$1)==$1){print i "\\n" \$1} }}' $ndir[2]/$rfix.filter.fa > $ndir[2]/$rfix.filter.$1.fa`;
			#	$r1=readpipe("awk -F \"_\" \'BEGIN{k=0;j=0}{if(\$_~/^>/){gsub(\"x\",\"\",\$3);k+=\$3;j++}}END{print k,j}\' $ndir[2]/$rfix.filter.$1.fa");
				$r1=readpipe("perl -e \'\$k=0;\$j=0;while(<>){if(\$_=~/^>\\S+_(\\d+)x/){\$k+=\$1;\$j++} elsif(\$_=~/^>\\S+/){\$k++;\$j++;} } print \"\$k,\$j\"\' $ndir[2]/$rfix.filter.$1.fa");
				$r1=~/(\d+),(\d+)/;	$rnum[$nlib][$rs][$r2][0] = $1;$rnum[$nlib][$rs][$r2][1] = $2;#$rnum[$nlib][$rs][$r2][3] = $nlen[$r2];
			} elsif ($nlen[$r2] eq "All") {
				print STDERR "$nlen[$r2];";
				`cp $ndir[2]/$rfix.filter.fa $ndir[2]/$rfix.filter.$nlen[$r2].fa`;
			}
			$rr.=$r1;
		}
	#} elsif (-e "$ndir[2]/$rfix.filter.All.fa") {
	} elsif ($r4 > $rs && $middleF == 1) {
	} elsif (($r4 < $rs || !(-e "$ndir[2]/$rfix.filter.All.fa")) && $rs < $r5) {
		`cp $ndir[2]/$rfix.filter.fa $ndir[2]/$rfix.filter.All.fa`;
	}
# step 6: map to genome using soap or bowtie
	$rs=6; print STDERR "\n";
	if ($rstep[$rs] != 0) {

	if ($soap =~/soap/) {
		print STDERR sub_format_datetime(localtime(time())),"\tStep $rs in Lib $nlib, $rfix: Running soap: ";
		for ($r2 = 0; $r2 < @nlen ;$r2++) {
			if ($nlen[$r2]=~/(\d+)\D+(\d+)/) {
				print STDERR "$1-$2;";
				system("$soap $soapset -D $soapdb -a $ndir[2]/$rfix.filter.$1-$2.fa -o $ndir[3]/$rfix.$1-$2.soap 2> $ndir[3]/$rfix.$1-$2.soap.err");
		#		system("perl -lane \'\$F[0]=~/_\\d+_(\\d+)x?/;\$k1=\$1;\@buf=\@F; \$i=shift(\@buf); \$j=join(\"\\t\",\@buf); for(\$k=1;\$k<=\$k1;\$k++) {print \"\$i\\_\$k\\t\$j\"}\' $ndir[3]/$rfix.$1-$2.soap > $ndir[3]/$rfix.$1-$2.1by1.soap");
		#		system("soap2sam_gl -G ~/Ath/data/tair10.Chr.fa $ndir[3]/$rfix.$1-$2.1by1.soap > $ndir[3]/$rfix.$1-$2.sam");
		#		system("samtools view -S -b $ndir[3]/$rfix.$1-$2.sam | samtools sort -T $ndir[3]/tmpp_$rfix.$1-$2 -o $ndir[3]/$rfix.$1-$2.bam -");
		#		system("samtools index $ndir[3]/$rfix.$1-$2.bam");
		#		$r1=readpipe("bam2wig -m -D $ndir[3]/bam2wigM $ndir[3]/$rfix.$1-$2.bam 2>&1");
		#		system("rm $ndir[3]/$rfix.$1-$2.1by1.soap");
		#		system("rm $ndir[3]/$rfix.$1-$2.sam");
			} elsif ($nlen[$r2]=~/(\d+)/) {
				print STDERR "$1;";
				system("$soap $soapset -D $soapdb -a $ndir[2]/$rfix.filter.$1.fa -o $ndir[3]/$rfix.$1.soap 2> $ndir[3]/$rfix.$1.soap.err");
		#		system("perl -lane \'\$F[0]=~/_\\d+_(\\d+)x?/;\$k1=\$1;\@buf=\@F; \$i=shift(\@buf); \$j=join(\"\\t\",\@buf); for(\$k=1;\$k<=\$k1;\$k++) {print \"\$i\\_\$k\\t\$j\"}\' $ndir[3]/$rfix.$1.soap > $ndir[3]/$rfix.$1.1by1.soap");
		#		system("soap2sam_gl -G ~/Ath/data/tair10.Chr.fa $ndir[3]/$rfix.$1.1by1.soap > $ndir[3]/$rfix.$1.sam");
		#		system("samtools view -S -b $ndir[3]/$rfix.$1.sam | samtools sort -T $ndir[3]/tmpp_$rfix.$1 -o $ndir[3]/$rfix.$1.bam -");
		#		system("samtools index $ndir[3]/$rfix.$1.bam");
		#		$r1=readpipe("bam2wig -m -D $ndir[3]/bam2wigM $ndir[3]/$rfix.$1.bam 2>&1");
		#		system("rm $ndir[3]/$rfix.$1.1by1.soap");
		#		system("rm $ndir[3]/$rfix.$1.sam");
			} elsif ($nlen[$r2] =~/(\w+)/) {
				print STDERR "$1;";
				system("$soap $soapset -D $soapdb -a $ndir[2]/$rfix.filter.$1.fa -o $ndir[3]/$rfix.$1.soap 2> $ndir[3]/$rfix.$1.soap.err");
				system("cut -f 1 $ndir[3]/$rfix.$1.soap | sort -u | perl -e \'my \%a;my \$i=0;while(<>){if(\$_=~/^(\\S+)_(\\S+)_(\\d+)x/){\$j=\$3;\$k=length(\$1);if(exists(\$a{\$k})){\$a{\$k}+=\$j} else{\$a{\$k}=\$j}} } \$m=0;print \"\\t$rfix\\n\"; foreach \$k (sort keys(\%a)) {\$m+=\$a{\$k};print \"\$k\\t\$a{\$k}\\n\"} print \"Total\\t\$m\\n\"\' > $ndir[3]/$rfix.$1.soap.lenDist");
				$r1=readpipe("cut -f 1 $ndir[3]/$rfix.$1.soap | sort -u | perl -e \'my \%a;my \$i=0;while(<>){if(\$_=~/^(\\S+)_(\\S+)_(\\d+)x/){\$j=\$3;\$k=length(\$1);if(exists(\$a{\$k})){\$a{\$k}+=\$j} else{\$a{\$k}=\$j}} } \$m=0;print \"Length\\tMap2genome\\n\"; foreach \$k (sort keys(\%a)) {\$m+=\$a{\$k};print \"\$k\\t\$a{\$k}\\n\"} \'");
				$lendist[$r11++]={split(/[\t|\n]/,$r1)};
		#		system("perl -lane \'\$F[0]=~/_\\d+_(\\d+)x?/;\$k1=\$1;\@buf=\@F; \$i=shift(\@buf); \$j=join(\"\\t\",\@buf); for(\$k=1;\$k<=\$k1;\$k++) {print \"\$i\\_\$k\\t\$j\"}\' $ndir[3]/$rfix.$1.soap > $ndir[3]/$rfix.$1.1by1.soap");
		#		system("soap2sam_gl -G ~/Ath/data/tair10.Chr.fa $ndir[3]/$rfix.$1.1by1.soap > $ndir[3]/$rfix.$1.sam");
		#		system("samtools view -S -b $ndir[3]/$rfix.$1.sam | samtools sort -T $ndir[3]/tmpp_$rfix.$1 -o $ndir[3]/$rfix.$1.bam -");
		#		system("samtools index $ndir[3]/$rfix.$1.bam");
		#		$r1=readpipe("bam2wig -m -D $ndir[3]/bam2wigM $ndir[3]/$rfix.$1.bam 2>&1");
		#		system("rm $ndir[3]/$rfix.$1.1by1.soap");
		#		system("rm $ndir[3]/$rfix.$1.sam");
			}
			$rr.=$r1;
		}
	} elsif ($soap=~/hisat2/) {
		print STDERR sub_format_datetime(localtime(time())),"\tStep $rs in Lib $nlib, $rfix: Running hisat2: ";
		for ($r2 = 0; $r2 < @nlen ;$r2++) {
			if ($nlen[$r2]=~/(\d+)\D+(\d+)/) {
				print STDERR "$1-$2;";
				system("$soap $soapset -f -x $soapdb -U $ndir[2]/$rfix.filter.$1-$2.fa -S $ndir[3]/$rfix.$1-$2.sam 2> $ndir[3]/$rfix.$1-$2.sam.err");
			} elsif ($nlen[$r2]=~/(\d+)/) {
				print STDERR "$1;";
				system("$soap $soapset -f -x $soapdb -U $ndir[2]/$rfix.filter.$1.fa -S $ndir[3]/$rfix.$1.sam 2> $ndir[3]/$rfix.$1.sam.err");
			} elsif ($nlen[$r2]=~/(\w+)/) {
				print STDERR "$1;";
				system("$soap $soapset -f -x $soapdb -U $ndir[2]/$rfix.filter.$1.fa -S $ndir[3]/$rfix.$1.sam 2> $ndir[3]/$rfix.$1.sam.err");
				system("awk \'{if(\$1!~/^@/ && \$3 != \"*\"){print \$1}}\' $ndir[3]/$rfix.$1.sam | sort -u | perl -e \'my \%a;my \$i=0;while(<>){if(\$_=~/^(\\S+)_(\\S+)_(\\d+)x/){\$j=\$3;\$k=length(\$1);if(exists(\$a{\$k})){\$a{\$k}+=\$j} else{\$a{\$k}=\$j}} } \$m=0;print \"\\t$rfix\\n\"; foreach \$k (sort keys(\%a)) {\$m+=\$a{\$k};print \"\$k\\t\$a{\$k}\\n\"} print \"Total\\t\$m\\n\"\' > $ndir[3]/$rfix.$1.sam.lenDist");
				$r1=readpipe("awk \'{if(\$1!~/^@/ && \$3 != \"*\"){print \$1}}\' $ndir[3]/$rfix.$1.sam | sort -u | perl -e \'my \%a;my \$i=0;while(<>){if(\$_=~/^(\\S+)_(\\S+)_(\\d+)x/){\$j=\$3;\$k=length(\$1);if(exists(\$a{\$k})){\$a{\$k}+=\$j} else{\$a{\$k}=\$j}} } \$m=0;print \"Length\\tMap2genome\\n\"; foreach \$k (sort keys(\%a)) {\$m+=\$a{\$k};print \"\$k\\t\$a{\$k}\\n\"} \' ");
				$lendist[$r11++]={split(/[\t|\n]/,$r1)};
			}
		}
	} elsif ($soap=~/bowtie2/) {
		print STDERR sub_format_datetime(localtime(time())),"\tStep $rs in Lib $nlib, $rfix: Running bowtie2: ";
		for ($r2 = 0; $r2 < @nlen ;$r2++) {
			if ($nlen[$r2]=~/(\d+)\D+(\d+)/) {
				print STDERR "$1-$2;";
				system("$soap $soapset -f -x $soapdb -U $ndir[2]/$rfix.filter.$1-$2.fa -S $ndir[3]/$rfix.$1-$2.sam 2> $ndir[3]/$rfix.$1-$2.sam.err");
			} elsif ($nlen[$r2]=~/(\d+)/) {
				print STDERR "$1;";
				system("$soap $soapset -f -x $soapdb -U $ndir[2]/$rfix.filter.$1.fa -S $ndir[3]/$rfix.$1.sam 2> $ndir[3]/$rfix.$1.sam.err");
			} elsif ($nlen[$r2]=~/(\w+)/) {
				print STDERR "$1;";
				system("$soap $soapset -f -x $soapdb -U $ndir[2]/$rfix.filter.$1.fa -S $ndir[3]/$rfix.$1.sam 2> $ndir[3]/$rfix.$1.sam.err");
				system("awk \'{if(\$1!~/^@/ && \$3 != \"*\"){print \$1}}\' $ndir[3]/$rfix.$1.sam | sort -u | perl -e \'my \%a;my \$i=0;while(<>){if(\$_=~/^(\\S+)_(\\S+)_(\\d+)x/){\$j=\$3;\$k=length(\$1);if(exists(\$a{\$k})){\$a{\$k}+=\$j} else{\$a{\$k}=\$j}} } \$m=0;print \"\\t$rfix\\n\"; foreach \$k (sort keys(\%a)) {\$m+=\$a{\$k};print \"\$k\\t\$a{\$k}\\n\"} print \"Total\\t\$m\\n\"\' > $ndir[3]/$rfix.$1.sam.lenDist");
				$r1=readpipe("awk \'{if(\$1!~/^@/ && \$3 != \"*\"){print \$1}}\' $ndir[3]/$rfix.$1.sam | sort -u | perl -e \'my \%a;my \$i=0;while(<>){if(\$_=~/^(\\S+)_(\\S+)_(\\d+)x/){\$j=\$3;\$k=length(\$1);if(exists(\$a{\$k})){\$a{\$k}+=\$j} else{\$a{\$k}=\$j}} } \$m=0;print \"Length\\tMap2genome\\n\"; foreach \$k (sort keys(\%a)) {\$m+=\$a{\$k};print \"\$k\\t\$a{\$k}\\n\"} \' ");
				$lendist[$r11++]={split(/[\t|\n]/,$r1)};
			}
		}
	} elsif ($soap=~/bowtie/) {
		print STDERR sub_format_datetime(localtime(time())),"\tStep $rs in Lib $nlib, $rfix: Running bowtie: ";
		for ($r2 = 0; $r2 < @nlen ;$r2++) {
			if ($nlen[$r2]=~/(\d+)\D+(\d+)/) {
				print STDERR "$1-$2;";
				system("$soap $soapset -S -f $soapdb $ndir[2]/$rfix.filter.$1-$2.fa $ndir[3]/$rfix.$1-$2.sam 2> $ndir[3]/$rfix.$1-$2.sam.err");
			} elsif ($nlen[$r2]=~/(\d+)/) {
				print STDERR "$1;";
				system("$soap $soapset -S -f $soapdb $ndir[2]/$rfix.filter.$1.fa $ndir[3]/$rfix.$1.sam 2> $ndir[3]/$rfix.$1.sam.err");
			} elsif ($nlen[$r2]=~/(\w+)/) {
				print STDERR "$1;";
				system("$soap $soapset -S -f $soapdb $ndir[2]/$rfix.filter.$1.fa $ndir[3]/$rfix.$1.sam 2> $ndir[3]/$rfix.$1.sam.err");
				system("awk \'{if(\$1!~/^@/ && \$3 != \"*\"){print \$1}}\' $ndir[3]/$rfix.$1.sam | sort -u | perl -e \'my \%a;my \$i=0;while(<>){if(\$_=~/^(\\S+)_(\\S+)_(\\d+)x/){\$j=\$3;\$k=length(\$1);if(exists(\$a{\$k})){\$a{\$k}+=\$j} else{\$a{\$k}=\$j}} } \$m=0;print \"\\t$rfix\\n\"; foreach \$k (sort keys(\%a)) {\$m+=\$a{\$k};print \"\$k\\t\$a{\$k}\\n\"} print \"Total\\t\$m\\n\"\' > $ndir[3]/$rfix.$1.sam.lenDist");
				$r1=readpipe("awk \'{if(\$1!~/^@/ && \$3 != \"*\"){print \$1}}\' $ndir[3]/$rfix.$1.sam | sort -u | perl -e \'my \%a;my \$i=0;while(<>){if(\$_=~/^(\\S+)_(\\S+)_(\\d+)x/){\$j=\$3;\$k=length(\$1);if(exists(\$a{\$k})){\$a{\$k}+=\$j} else{\$a{\$k}=\$j}} } \$m=0;print \"Length\\tMap2genome\\n\"; foreach \$k (sort keys(\%a)) {\$m+=\$a{\$k};print \"\$k\\t\$a{\$k}\\n\"} \' ");
				$lendist[$r11++]={split(/[\t|\n]/,$r1)};
			}
		}
	} elsif ($soap=~/bwa\s+mem/) {
		print STDERR sub_format_datetime(localtime(time())),"\tStep $rs in Lib $nlib, $rfix: Running bwa mem: ";
		for ($r2 = 0; $r2 < @nlen ;$r2++) {
			if ($nlen[$r2]=~/(\d+)\D+(\d+)/) {
				print STDERR "$1-$2;";
				system("$soap $soapdb $ndir[2]/$rfix.filter.$1-$2.fa > $ndir[3]/$rfix.$1-$2.sam 2> $ndir[3]/$rfix.$1-$2.sam.err");
			} elsif ($nlen[$r2]=~/(\d+)/) {
				print STDERR "$1;";
				system("$soap $ndir[2]/$rfix.filter.$1.fa > $ndir[3]/$rfix.$1.sam 2> $ndir[3]/$rfix.$1.sam.err");
			} elsif ($nlen[$r2]=~/(\w+)/) {
				print STDERR "$1;";
				system("$soap $soapdb $ndir[2]/$rfix.filter.$1.fa > $ndir[3]/$rfix.$1.sam 2> $ndir[3]/$rfix.$1.sam.err");
				system("awk \'{if(\$1!~/^@/ && \$3 != \"*\"){print \$1}}\' $ndir[3]/$rfix.$1.sam | sort -u | perl -e \'my \%a;my \$i=0;while(<>){if(\$_=~/^(\\S+)_(\\S+)_(\\d+)x/){\$j=\$3;\$k=length(\$1);if(exists(\$a{\$k})){\$a{\$k}+=\$j} else{\$a{\$k}=\$j}} } \$m=0;print \"\\t$rfix\\n\"; foreach \$k (sort keys(\%a)) {\$m+=\$a{\$k};print \"\$k\\t\$a{\$k}\\n\"} print \"Total\\t\$m\\n\"\' > $ndir[3]/$rfix.$1.sam.lenDist");
				$r1=readpipe("awk \'{if(\$1!~/^@/ && \$3 != \"*\"){print \$1}}\' $ndir[3]/$rfix.$1.sam | sort -u | perl -e \'my \%a;my \$i=0;while(<>){if(\$_=~/^(\\S+)_(\\S+)_(\\d+)x/){\$j=\$3;\$k=length(\$1);if(exists(\$a{\$k})){\$a{\$k}+=\$j} else{\$a{\$k}=\$j}} } \$m=0;print \"Length\\tMap2genome\\n\"; foreach \$k (sort keys(\%a)) {\$m+=\$a{\$k};print \"\$k\\t\$a{\$k}\\n\"} \' ");
				$lendist[$r11++]={split(/[\t|\n]/,$r1)};
				if ($mode=~/igv/i ) {
					system("samtools view -b -o $ndir[3]/$rfix.$nlen[$r2].nost.bam $ndir[3]/$rfix.$nlen[$r2].sam 2>> $ndir[3]/$rfix.$nlen[$r2].sam.err");
					system("samtools sort -T $ndir[3]/$rfix.$nlen[$r2].nosorted -o $ndir[3]/$rfix.$nlen[$r2].bam $ndir[3]/$rfix.$nlen[$r2].nost.bam 2>> $ndir[3]/$rfix.$nlen[$r2].sam.err");
					system("samtools index $ndir[3]/$rfix.$nlen[$r2].bam 2>> $ndir[3]/$rfix.$nlen[$r2].sam.err");
					system("rm $ndir[3]/$rfix.$nlen[$r2].nost.bam");
				}
			}
		}
	} elsif ($soap=~/bwa\s+aln/) {
		print STDERR sub_format_datetime(localtime(time())),"\tStep $rs in Lib $nlib, $rfix: Running bwa aln: ";
		for ($r2 = 0; $r2 < @nlen ;$r2++) {
			if ($nlen[$r2]=~/(\d+)\D+(\d+)/) {
				print STDERR "$1-$2;";
				system("$soap $soapdb $ndir[2]/$rfix.filter.$1-$2.fa > $ndir[3]/$rfix.$1-$2.sai 2> $ndir[3]/$rfix.$1-$2.sam.err");
				system("bwa samse $soapset $soapdb $ndir[3]/$rfix.$1-$2.sai $ndir[2]/$rfix.filter.$1-$2.fa > $ndir[3]/$rfix.$1-$2.sam 2>> $ndir[3]/$rfix.$1-$2.sam.err");
			} elsif ($nlen[$r2]=~/(\d+)/) {
				print STDERR "$1;";
				system("$soap $soapdb $ndir[2]/$rfix.filter.$1.fa > $ndir[3]/$rfix.$1.sai 2> $ndir[3]/$rfix.$1.sam.err");
				system("bwa samse $soapset $soapdb $ndir[3]/$rfix.$1.sai $ndir[2]/$rfix.filter.$1.fa > $ndir[3]/$rfix.$1.sam 2>> $ndir[3]/$rfix.$1.sam.err");
			} elsif ($nlen[$r2]=~/(\w+)/) {
				print STDERR "$1;";
				system("$soap $soapdb $ndir[2]/$rfix.filter.$1.fa > $ndir[3]/$rfix.$1.sai 2> $ndir[3]/$rfix.$1.sam.err");
				system("bwa samse $soapset $soapdb $ndir[3]/$rfix.$1.sai $ndir[2]/$rfix.filter.$1.fa > $ndir[3]/$rfix.$1.sam 2>> $ndir[3]/$rfix.$1.sam.err");
			#	system("$soap $soapset $soapdb $ndir[2]/$rfix.filter.$1.fa > $ndir[3]/$rfix.$1.sam 2> $ndir[3]/$rfix.$1.sam.err");
				system("awk \'{if(\$1!~/^@/ && \$3 != \"*\"){print \$1}}\' $ndir[3]/$rfix.$1.sam | sort -u | perl -e \'my \%a;my \$i=0;while(<>){if(\$_=~/^(\\S+)_(\\S+)_(\\d+)x/){\$j=\$3;\$k=length(\$1);if(exists(\$a{\$k})){\$a{\$k}+=\$j} else{\$a{\$k}=\$j}} } \$m=0;print \"\\t$rfix\\n\"; foreach \$k (sort keys(\%a)) {\$m+=\$a{\$k};print \"\$k\\t\$a{\$k}\\n\"} print \"Total\\t\$m\\n\"\' > $ndir[3]/$rfix.$1.sam.lenDist");
				$r1=readpipe("awk \'{if(\$1!~/^@/ && \$3 != \"*\"){print \$1}}\' $ndir[3]/$rfix.$1.sam | sort -u | perl -e \'my \%a;my \$i=0;while(<>){if(\$_=~/^(\\S+)_(\\S+)_(\\d+)x/){\$j=\$3;\$k=length(\$1);if(exists(\$a{\$k})){\$a{\$k}+=\$j} else{\$a{\$k}=\$j}} } \$m=0;print \"Length\\tMap2genome\\n\"; foreach \$k (sort keys(\%a)) {\$m+=\$a{\$k};print \"\$k\\t\$a{\$k}\\n\"} \' ");
				$lendist[$r11++]={split(/[\t|\n]/,$r1)};
				if ($mode=~/igv/i ) {
					system("samtools view -b -o $ndir[3]/$rfix.$nlen[$r2].nost.bam $ndir[3]/$rfix.$nlen[$r2].sam 2>> $ndir[3]/$rfix.$nlen[$r2].sam.err");
					system("samtools sort -T $ndir[3]/$rfix.$nlen[$r2].nosorted -o $ndir[3]/$rfix.$nlen[$r2].bam $ndir[3]/$rfix.$nlen[$r2].nost.bam 2>> $ndir[3]/$rfix.$nlen[$r2].sam.err");
					system("samtools index $ndir[3]/$rfix.$nlen[$r2].bam 2>> $ndir[3]/$rfix.$nlen[$r2].sam.err");
					system("rm $ndir[3]/$rfix.$nlen[$r2].nost.bam");
				}
			}
		}
	} else { # tophat
		print STDERR sub_format_datetime(localtime(time())),"\tStep $rs in Lib $nlib, $rfix: Running tophat: ";
		for ($r2 = 0; $r2 < @nlen ;$r2++) {
			if ($nlen[$r2]=~/(\d+)\D+(\d+)/) {
				print STDERR "$1-$2;";
				system("$soap $soapset -o $ndir[3]/$rfix.$1-$2.tophat2 $soapdb $ndir[2]/$rfix.filter.$1-$2.fa  2> $ndir[3]/$rfix.filter.$1-$2.tophat.err");
				system("mv $ndir[3]/$rfix.$1-$2.tophat2/accepted_hits.bam $ndir[3]/$rfix.$1-$2.bam");
				system("samtools index $ndir[3]/$rfix.$1-$2.bam");
				system("samtools view -h -o $ndir[3]/$rfix.$1-$2.sam $ndir[3]/$rfix.$1-$2.bam");
			} elsif ($nlen[$r2]=~/(\d+)/) {
				print STDERR "$1;";
				system("$soap $soapset -o $ndir[3]/$rfix.$1.tophat2 $soapdb $ndir[2]/$rfix.filter.$1.fa 2> $ndir[3]/$rfix.filter.$1.tophat.err");
				system("mv $ndir[3]/$rfix.$1.tophat2/accepted_hits.bam $ndir[3]/$rfix.$1.bam");
				system("samtools index $ndir[3]/$rfix.$1.bam 2>> $ndir[3]/$rfix.filter.$1.tophat.err");
				system("samtools view -h -o $ndir[3]/$rfix.$1.sam $ndir[3]/$rfix.$1.bam 2>> $ndir[3]/$rfix.filter.$1.tophat.err");
			} elsif ($nlen[$r2]=~/(\w+)/) {
				print STDERR "$1;";
				system("$soap $soapset -o $ndir[3]/$rfix.$1.tophat2 $soapdb $ndir[2]/$rfix.filter.$1.fa 2> $ndir[3]/$rfix.filter.$1.tophat.err");
				system("mv $ndir[3]/$rfix.$1.tophat2/accepted_hits.bam $ndir[3]/$rfix.$1.bam");
				system("samtools index $ndir[3]/$rfix.$1.bam 2>> $ndir[3]/$rfix.filter.$1.tophat.err");
				system("samtools view -h -o $ndir[3]/$rfix.$1.sam $ndir[3]/$rfix.$1.bam 2>> $ndir[3]/$rfix.filter.$1.tophat.err");
				system("awk \'{if(\$1!~/^@/ && \$3 != \"*\"){print \$1}}\' $ndir[3]/$rfix.$1.sam | sort -u | perl -e \'my \%a;my \$i=0;while(<>){if(\$_=~/^(\\S+)_(\\S+)_(\\d+)x/){\$j=\$3;\$k=length(\$1);if(exists(\$a{\$k})){\$a{\$k}+=\$j} else{\$a{\$k}=\$j}} } \$m=0;print \"\\t$rfix\\n\"; foreach \$k (sort keys(\%a)) {\$m+=\$a{\$k};print \"\$k\\t\$a{\$k}\\n\"} print \"Total\\t\$m\\n\"\' > $ndir[3]/$rfix.$1.sam.lenDist");
				$r1=readpipe("awk \'{if(\$1!~/^@/ && \$3 != \"*\"){print \$1}}\' $ndir[3]/$rfix.$1.sam | sort -u | perl -e \'my \%a;my \$i=0;while(<>){if(\$_=~/^(\\S+)_(\\S+)_(\\d+)x/){\$j=\$3;\$k=length(\$1);if(exists(\$a{\$k})){\$a{\$k}+=\$j} else{\$a{\$k}=\$j}} } \$m=0;print \"Length\\tMap2genome\\n\"; foreach \$k (sort keys(\%a)) {\$m+=\$a{\$k};print \"\$k\\t\$a{\$k}\\n\"} \' ");
				$lendist[$r11++]={split(/[\t|\n]/,$r1)};
			}
		}
	} 
	
	}
# step 7: yield wig file for genome view
	$rs=7; print STDERR "\n";
	if ($rstep[$rs] != 0) {

	if ($soap =~/soap/) {
		print STDERR sub_format_datetime(localtime(time())),"\tStep $rs in Lib $nlib, $rfix: Converted to wiggle: ";
		for ($r2 = 0; $r2 < @nlen ;$r2++) {
			if ($nlen[$r2]=~/(\d+)\D+(\d+)/) {
				print STDERR "$1-$2;";
		#		system("$soap $soapset -D $soapdb -a $ndir[2]/$rfix.filter.$1-$2.fa -o $ndir[3]/$rfix.$1-$2.soap 2> $ndir[3]/$rfix.$1-$2.soap.err");
				system("perl -lane \'if(\$F[0]=~/_\\d+_(\\d+)x?/){\$k1=\$1;}else{\$k1=1;} \@buf=\@F; \$i=shift(\@buf); \$j=join(\"\\t\",\@buf); for(\$k=1;\$k<=\$k1;\$k++) {print \"\$i\\_\$k\\t\$j\"}\' $ndir[3]/$rfix.$1-$2.soap > $ndir[3]/$rfix.$1-$2.1by1.soap");
				system("soap2sam_gl -G $genome $ndir[3]/$rfix.$1-$2.1by1.soap > $ndir[3]/$rfix.$1-$2.sam");
				system("samtools view -S -b $ndir[3]/$rfix.$1-$2.sam | samtools sort -T $ndir[3]/tmpp_$rfix.$1-$2 -o $ndir[3]/$rfix.$1-$2.bam - >> $ndir[3]/$rfix.$1-$2.soap.err");
				system("samtools index $ndir[3]/$rfix.$1-$2.bam");
				$r1=readpipe("bam2wig -m -D $ndir[3]/bam2wigM $ndir[3]/$rfix.$1-$2.bam 2>&1");
				system("rm $ndir[3]/$rfix.$1-$2.1by1.soap");
		#		system("rm $ndir[3]/$rfix.$1-$2.sam");
			} elsif ($nlen[$r2]=~/(\d+)/) {
				print STDERR "$1;";
		#		system("$soap $soapset -D $soapdb -a $ndir[2]/$rfix.filter.$1.fa -o $ndir[3]/$rfix.$1.soap 2> $ndir[3]/$rfix.$1.soap.err");
				system("perl -lane \'if(\$F[0]=~/_\\d+_(\\d+)x?/){\$k1=\$1;}else{\$k1=1;} \@buf=\@F; \$i=shift(\@buf); \$j=join(\"\\t\",\@buf); for(\$k=1;\$k<=\$k1;\$k++) {print \"\$i\\_\$k\\t\$j\"}\' $ndir[3]/$rfix.$1.soap > $ndir[3]/$rfix.$1.1by1.soap");
				system("soap2sam_gl -G $genome $ndir[3]/$rfix.$1.1by1.soap > $ndir[3]/$rfix.$1.sam");
				system("samtools view -S -b $ndir[3]/$rfix.$1.sam | samtools sort -T $ndir[3]/tmpp_$rfix.$1 -o $ndir[3]/$rfix.$1.bam - >> $ndir[3]/$rfix.$1.soap.err");
				system("samtools index $ndir[3]/$rfix.$1.bam");
				$r1=readpipe("bam2wig -m -D $ndir[3]/bam2wigM $ndir[3]/$rfix.$1.bam 2>&1");
				system("rm $ndir[3]/$rfix.$1.1by1.soap");
		#		system("rm $ndir[3]/$rfix.$1.sam");
			} elsif ($nlen[$r2] =~/(\w+)/) {
				print STDERR "$1;";
		#		system("$soap $soapset -D $soapdb -a $ndir[2]/$rfix.filter.$1.fa -o $ndir[3]/$rfix.$1.soap 2> $ndir[3]/$rfix.$1.soap.err");
		#		system("cut -f 1 $ndir[3]/$rfix.$1.soap | sort -u | perl -e \'my \%a;my \$i=0;while(<>){if(\$_=~/^(\\S+)_(\\S+)_(\\d+)x/){\$j=\$3;\$k=length(\$1);if(exists(\$a{\$k})){\$a{\$k}+=\$j} else{\$a{\$k}=\$j}} } \$m=0;print \"\\t$rfix\\n\"; foreach \$k (sort keys(\%a)) {\$m+=\$a{\$k};print \"\$k\\t\$a{\$k}\\n\"} print \"Total\\t\$m\\n\"\' > $ndir[3]/$rfix.$1.soap.lenDist");
				system("perl -lane \'if(\$F[0]=~/_\\d+_(\\d+)x?/){\$k1=\$1;}else{\$k1=1;} \@buf=\@F; \$i=shift(\@buf); \$j=join(\"\\t\",\@buf); for(\$k=1;\$k<=\$k1;\$k++) {print \"\$i\\_\$k\\t\$j\"}\' $ndir[3]/$rfix.$1.soap > $ndir[3]/$rfix.$1.1by1.soap");
				system("soap2sam_gl -G $genome $ndir[3]/$rfix.$1.1by1.soap > $ndir[3]/$rfix.$1.sam");
				system("samtools view -S -b $ndir[3]/$rfix.$1.sam | samtools sort -T $ndir[3]/tmpp_$rfix.$1 -o $ndir[3]/$rfix.$1.bam - >> $ndir[3]/$rfix.$1.soap.err");
				system("samtools index $ndir[3]/$rfix.$1.bam");
				$r1=readpipe("bam2wig -m -D $ndir[3]/bam2wigM $ndir[3]/$rfix.$1.bam 2>&1");
				system("rm $ndir[3]/$rfix.$1.1by1.soap");
		#		system("rm $ndir[3]/$rfix.$1.sam");
			}
			$rr.=$r1;
		}
	} elsif ($soap=~/tophat/) { # tophat
		print STDERR sub_format_datetime(localtime(time())),"\tStep $rs in Lib $nlib, $rfix: Converted to wiggle: ";
		for ($r2 = 0; $r2 < @nlen ;$r2++) {
			if ($nlen[$r2]=~/(\d+)\D+(\d+)/) {
				print STDERR "$1-$2;";
				$r1=readpipe("bam2wig -m -D $ndir[3]/bam2wigM $ndir[3]/$rfix.$1-$2.bam 2>&1");
			} elsif ($nlen[$r2]=~/(\d+)/) {
				print STDERR "$1;";
				$r1=readpipe("bam2wig -m -D $ndir[3]/bam2wigM $ndir[3]/$rfix.$1.bam 2>&1");
			} elsif ($nlen[$r2]=~/(\w+)/) {
				print STDERR "$1;";
				$r1=readpipe("bam2wig -m -D $ndir[3]/bam2wigM $ndir[3]/$rfix.$1.bam 2>&1");
			}
		}
	} else { #elsif ($soap=~/hisat2/) {
		print STDERR sub_format_datetime(localtime(time())),"\tStep $rs in Lib $nlib, $rfix: Converted to wiggle: ";
		for ($r2 = 0; $r2 < @nlen ;$r2++) {
			if ($nlen[$r2]=~/(\d+)\D+(\d+)/) {
				print STDERR "$1-$2;";
				system("samtools view -S -b $ndir[3]/$rfix.$1-$2.sam | samtools sort -T $ndir[3]/tmpp_$rfix.$1-$2 -o $ndir[3]/$rfix.$1-$2.bam - >> $ndir[3]/$rfix.$1-$2.sam.err");
				system("samtools index $ndir[3]/$rfix.$1-$2.bam");
				$r1=readpipe("bam2wig -m -D $ndir[3]/bam2wigM $ndir[3]/$rfix.$1-$2.bam 2>&1");
			} elsif ($nlen[$r2]=~/(\d+)/) {
				print STDERR "$1;";
				system("samtools view -S -b $ndir[3]/$rfix.$1.sam | samtools sort -T $ndir[3]/tmpp_$rfix.$1 -o $ndir[3]/$rfix.$1.bam - >> $ndir[3]/$rfix.$1.sam.err");
				system("samtools index $ndir[3]/$rfix.$1.bam");
				$r1=readpipe("bam2wig -m -D $ndir[3]/bam2wigM $ndir[3]/$rfix.$1.bam 2>&1");
			} elsif ($nlen[$r2]=~/(\w+)/) {
				print STDERR "$1;";
				system("samtools view -S -b $ndir[3]/$rfix.$1.sam | samtools sort -T $ndir[3]/tmpp_$rfix.$1 -o $ndir[3]/$rfix.$1.bam - >> $ndir[3]/$rfix.$1.sam.err");
				system("samtools index $ndir[3]/$rfix.$1.bam");
				$r1=readpipe("bam2wig -m -D $ndir[3]/bam2wigM $ndir[3]/$rfix.$1.bam 2>&1");
			}
		}
	} 
	
	}
# step 8: info for output
#	$rs=8; 
	print STDERR "\t";
	if ($rstep[$rs] != 0) {

#	print STDERR sub_format_datetime(localtime(time())),"\tStep $rs in Lib $nlib, $rfix: Get the info: ";
#	if ($rstep[1]==1) {
#		$r1=readpipe("awk -F \"_\" \'BEGIN{k=0;j=0}{if(\$_~/^>/){gsub(\"x\",\"\",\$3);k+=\$3;j++}}END{print k,j}\' $ndir[2]/$rfix.filter.fa");
#	#	print STDERR "r1=$r1,"; 
#		$r1=~/(\d+)\s+(\d+)/;$rnum[$nlib][$r4++] = $1;$rnum[$nlib][$r4++] = $2; # remove repeats
#	} else {
#		$rnum[$nlib][$r4++] = 0;$rnum[$nlib][$r4++] = 0; # remove repeats
#	}

#	$r4 = 5;
	if ($soap =~/soap/) {
		for ($r2 = 0; $r2 < @nlen ;$r2++) {
			if ($nlen[$r2]=~/(\d+)\D+(\d+)/) {
				print STDERR "$1-$2;";
				$r3="$1-$2";
			#	$r1=readpipe("awk -F \"_\" \'BEGIN{k=0;j=0}{if(\$_~/^>/){gsub(\"x\",\"\",\$3);k+=\$3;j++}}END{print k,j}\' $ndir[2]/$rfix.filter.$1-$2.fa");
			#	$r1=~/(\d+)\s+(\d+)/;	$rnum[$nlib][$rs][$r2][0] = $1;$rnum[$nlib][$rs][$r2][1] = $2;$rnum[$nlib][$rs][$r2][3] = $nlen[$r2];

			#	$r1=readpipe("awk -F \"_\" \'BEGIN{k=0;j=0}{if(\$_~/^>/){gsub(\"x\",\"\",\$3);k+=\$3;j++}}END{print k,j}\' $ndir[2]/$rfix.filter.$r3.fa");
			#	print STDERR "r1=$r1,"; 
			#	$r1=~/(\d+)\s+(\d+)/;$rnum[$nlib][6][0] = $1;$rnum[$nlib][6][1] = $2; # length
				$r1=readpipe("cut -f 1 $ndir[3]/$rfix.$r3.soap | sort -u | awk -F \"_\" \'BEGIN{k=0;j=0}{gsub(\"x\",\"\",\$3);k+=\$3;j++}END{print k,j}\' ");
			#	print STDERR "r1=$r1,"; 
			#	$r1=~/(\d+)\s+(\d+)/;$rnum[$nlib][6][0] = $1;$rnum[$nlib][6][1] = $2; # length
				$r1=~/(\d+)\s+(\d+)/;$rnum[$nlib][6][$r2][0] = $1;$rnum[$nlib][6][$r2][1] = $2;
			} elsif ($nlen[$r2]=~/(\w+)/) {
				print STDERR "$1;";
				$r3=$1;
			#	$r1=readpipe("awk -F \"_\" \'BEGIN{k=0;j=0}{if(\$_~/^>/){gsub(\"x\",\"\",\$3);k+=\$3;j++}}END{print k,j}\' $ndir[2]/$rfix.filter.$r3.fa");
			#	print STDERR "r1=$r1,"; 
			#	$r1=~/(\d+)\s+(\d+)/;$rnum[$nlib][$r4++] = $1;$rnum[$nlib][$r4++] = $2; # length
				$r1=readpipe("cut -f 1 $ndir[3]/$rfix.$r3.soap | sort -u | awk -F \"_\" \'BEGIN{k=0;j=0}{gsub(\"x\",\"\",\$3);k+=\$3;j++}END{print k,j}\' ");
			#	print STDERR "r1=$r1,"; 
			#	$r1=~/(\d+)\s+(\d+)/;$rnum[$nlib][$r4++] = $1;$rnum[$nlib][$r4++] = $2; # length
				$r1=~/(\d+)\s+(\d+)/;$rnum[$nlib][6][$r2][0] = $1;$rnum[$nlib][6][$r2][1] = $2;
			}
		}
	} elsif ($soap =~/bowtie/) {
		for ($r2 = 0; $r2 < @nlen ;$r2++) {
			if ($nlen[$r2]=~/(\d+)\D+(\d+)/) {
				print STDERR "$1-$2;";
				$r3="$1-$2";
			#	$r1=readpipe("awk -F \"_\" \'BEGIN{k=0;j=0}{if(\$_~/^>/){gsub(\"x\",\"\",\$3);k+=\$3;j++}}END{print k,j}\' $ndir[2]/$rfix.filter.$r3.fa");
			#	print STDERR "r1=$r1,"; 
			#	$r1=~/(\d+)\s+(\d+)/;$rnum[$nlib][$r4++] = $1;$rnum[$nlib][$r4++] = $2; # length
				$r1=readpipe("awk \'{if(\$1!~/^\@/ && \$3~/[a-zA-Z]+/){print \$1}}\' $ndir[3]/$rfix.$r3.sam | sort -u | awk -F \"_\" \'BEGIN{k=0;j=0}{gsub(\"x\",\"\",\$3);k+=\$3;j++}END{print k,j}\' ");
			#	print STDERR "r1=$r1,"; 
			#	$r1=~/(\d+)\s+(\d+)/;$rnum[$nlib][$r4++] = $1;$rnum[$nlib][$r4++] = $2; # length
				$r1=~/(\d+)\s+(\d+)/;$rnum[$nlib][6][$r2][0] = $1;$rnum[$nlib][6][$r2][1] = $2;
			} elsif ($nlen[$r2]=~/(\w+)/) {
				print STDERR "$1;";
				$r3=$1;
			#	$r1=readpipe("awk -F \"_\" \'BEGIN{k=0;j=0}{if(\$_~/^>/){gsub(\"x\",\"\",\$3);k+=\$3;j++}}END{print k,j}\' $ndir[2]/$rfix.filter.$r3.fa");
			#	print STDERR "r1=$r1,"; 
			#	$r1=~/(\d+)\s+(\d+)/;$rnum[$nlib][$r4++] = $1;$rnum[$nlib][$r4++] = $2; # length
				$r1=readpipe("awk \'{if(\$1!~/^\@/ && \$3~/[a-zA-Z]+/){print \$1}}\' $ndir[3]/$rfix.$r3.sam | sort -u | awk -F \"_\" \'BEGIN{k=0;j=0}{gsub(\"x\",\"\",\$3);k+=\$3;j++}END{print k,j}\' ");
			#	print STDERR "r1=$r1,"; 
			#	$r1=~/(\d+)\s+(\d+)/;$rnum[$nlib][$r4++] = $1;$rnum[$nlib][$r4++] = $2; # length
				$r1=~/(\d+)\s+(\d+)/;$rnum[$nlib][6][$r2][0] = $1;$rnum[$nlib][6][$r2][1] = $2;
			}
		}
	} elsif ($soap =~/bwa/) {
		for ($r2 = 0; $r2 < @nlen ;$r2++) {
			if ($nlen[$r2]=~/(\d+)\D+(\d+)/) {
				print STDERR "$1-$2;";
				$r3="$1-$2";
			#	$r1=readpipe("awk -F \"_\" \'BEGIN{k=0;j=0}{if(\$_~/^>/){gsub(\"x\",\"\",\$3);k+=\$3;j++}}END{print k,j}\' $ndir[2]/$rfix.filter.$r3.fa");
			#	print STDERR "r1=$r1,"; 
			#	$r1=~/(\d+)\s+(\d+)/;$rnum[$nlib][$r4++] = $1;$rnum[$nlib][$r4++] = $2; # length
				$r1=readpipe("awk \'{if(\$1!~/^\@/ && \$3~/[a-zA-Z]+/){print \$1}}\' $ndir[3]/$rfix.$r3.sam | sort -u | awk -F \"_\" \'BEGIN{k=0;j=0}{gsub(\"x\",\"\",\$3);k+=\$3;j++}END{print k,j}\' ");
			#	print STDERR "r1=$r1,"; 
			#	$r1=~/(\d+)\s+(\d+)/;$rnum[$nlib][$r4++] = $1;$rnum[$nlib][$r4++] = $2; # length
				$r1=~/(\d+)\s+(\d+)/;$rnum[$nlib][6][$r2][0] = $1;$rnum[$nlib][6][$r2][1] = $2;
			} elsif ($nlen[$r2]=~/(\w+)/) {
				print STDERR "$1;";
				$r3=$1;
			#	$r1=readpipe("awk -F \"_\" \'BEGIN{k=0;j=0}{if(\$_~/^>/){gsub(\"x\",\"\",\$3);k+=\$3;j++}}END{print k,j}\' $ndir[2]/$rfix.filter.$r3.fa");
			#	print STDERR "r1=$r1,"; 
			#	$r1=~/(\d+)\s+(\d+)/;$rnum[$nlib][$r4++] = $1;$rnum[$nlib][$r4++] = $2; # length
				$r1=readpipe("awk \'{if(\$1!~/^\@/ && \$3~/[a-zA-Z]+/){print \$1}}\' $ndir[3]/$rfix.$r3.sam | sort -u | awk -F \"_\" \'BEGIN{k=0;j=0}{gsub(\"x\",\"\",\$3);k+=\$3;j++}END{print k,j}\' ");
			#	print STDERR "r1=$r1,"; 
			#	$r1=~/(\d+)\s+(\d+)/;$rnum[$nlib][$r4++] = $1;$rnum[$nlib][$r4++] = $2; # length
				$r1=~/(\d+)\s+(\d+)/;$rnum[$nlib][6][$r2][0] = $1;$rnum[$nlib][6][$r2][1] = $2;
				if ($mode=~/igv/i ) {
					system("rm $ndir[3]/$rfix.$nlen[$r2].sam");
				}
			}
		}
	} else { # tophat
		for ($r2 = 0; $r2 < @nlen ;$r2++) {
			if ($nlen[$r2]=~/(\d+)\D+(\d+)/) {
				print STDERR "$1-$2;";
				$r3="$1-$2";
			#	$r1=readpipe("awk -F \"_\" \'BEGIN{k=0;j=0}{if(\$_~/^>/){gsub(\"x\",\"\",\$3);k+=\$3;j++}}END{print k,j}\' $ndir[2]/$rfix.filter.$r3.fa");
			#	print STDERR "r1=$r1,"; 
			#	$r1=~/(\d+)\s+(\d+)/;$rnum[$nlib][$r4++] = $1;$rnum[$nlib][$r4++] = $2; # length
				$r1=readpipe("awk \'{if(\$1!~/^\@/ && \$3~/^[a-zA-Z]+/){print \$1}}\' $ndir[3]/$rfix.$r3.sam | sort -u | awk -F \"_\" \'BEGIN{k=0;j=0}{gsub(\"x\",\"\",\$3);k+=\$3;j++}END{print k,j}\' ");
			#	print STDERR "r1=$r1,"; 
			#	$r1=~/(\d+)\s+(\d+)/;$rnum[$nlib][$r4++] = $1;$rnum[$nlib][$r4++] = $2; # length
				$r1=~/(\d+)\s+(\d+)/;$rnum[$nlib][6][$r2][0] = $1;$rnum[$nlib][6][$r2][1] = $2;
			} elsif ($nlen[$r2]=~/(\w+)/) {
				print STDERR "$1;";
				$r3=$1;
			#	$r1=readpipe("awk -F \"_\" \'BEGIN{k=0;j=0}{if(\$_~/^>/){gsub(\"x\",\"\",\$3);k+=\$3;j++}}END{print k,j}\' $ndir[2]/$rfix.filter.$r3.fa");
			#	print STDERR "r1=$r1,"; 
			#	$r1=~/(\d+)\s+(\d+)/;$rnum[$nlib][$r4++] = $1;$rnum[$nlib][$r4++] = $2; # length
				$r1=readpipe("awk \'{if(\$1!~/^\@/ && \$3~/^[a-zA-Z]+/){print \$1}}\' $ndir[3]/$rfix.$r3.sam | sort -u | awk -F \"_\" \'BEGIN{k=0;j=0}{gsub(\"x\",\"\",\$3);k+=\$3;j++}END{print k,j}\' ");
			#	print STDERR "r1=$r1,"; 
			#	$r1=~/(\d+)\s+(\d+)/;$rnum[$nlib][$r4++] = $1;$rnum[$nlib][$r4++] = $2; # length
				$r1=~/(\d+)\s+(\d+)/;$rnum[$nlib][6][$r2][0] = $1;$rnum[$nlib][6][$r2][1] = $2;
			}
		}
	}
	
	}
	################################ draw length distribution ###################################
	@drawlen = sort {$a<=>$b} keys %{$lendist[0]};print "\n====drawlen====\n";print " ",@drawlen; print "\n====drawlen====\n";
	$rs = "";
	for ($r1 = $drawlen[0]; $r1 <= $drawlen[-1] ;$r1++) {
		$rs .= "$r1";
		for ($r21 = 0; $r21 < $r11 ;$r21++) {
			if (exists($lendist[$r21]{$r1})) {
				$rs .= "\t$lendist[$r21]{$r1}";
			} else {
				$rs .= "\t0";
			}
		}
		$rs .= "\n";
	}
	print "\n====rs====\n$rs\n====rs====\n";
	################################ delete the middle file   ###################################
	for ($rs = $r4; $rs <= $r5 ;$rs++) { # delete the middle file from copying directly
		if ($rstep[$rs] == 0 && $rs == 1) {
			system("rm $ndir[0]/$rfix.fasta")
		} elsif ($rstep[$rs] == 0 && $rs == 2) {
			system("rm $ndir[1]/$rfix.trim.fasta")
		} elsif ($rstep[$rs] == 0 && $rs == 3) {
			system("rm $ndir[1]/$rfix.fa")
		} elsif ($rstep[$rs] == 0 && $rs == 4) {
			system("rm $ndir[2]/$rfix.filter.fa")
		} elsif ($rstep[$rs] == 0 && $rs == 5) {
			system("rm $ndir[2]/$rfix.filter.All.fa")
		} elsif ($rstep[$rs] == 0 && $rs == 6) {
	#		system("rm $ndir[0]/$rfix.fasta")
		}
	}
	print STDERR "\n";
	print STDERR sub_format_datetime(localtime(time())),"\tLib $nlib, $rfix finished!\n";
	return $r4,$rr;
}


__END__

=head1 LICENSE

auspp

Copyright (c) 2018- Lei GAO

This program is free software: you can redistribute it and/or modify                                                                                                                
it under the terms of the GNU General Public License as published by                                                                                                                
the Free Software Foundation, either version 3 of the License, or                                                                                                                   
(at your option) any later version.                                                                                                                                                 
                                                                                                                                                                                    
This program is distributed in the hope that it will be useful,                                                                                                                     
but WITHOUT ANY WARRANTY; without even the implied warranty of                                                                                                                      
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                                                                                                       
GNU General Public License for more details.                                                                                                                                        
                                                                                                                                                                                    
You should have received a copy of the GNU General Public License                                                                                                                   
along with this program.  If not, see <http://www.gnu.org/licenses/>. 

=head1 SYNOPSIS

auapp : a universal short-read pre-processing package/program

=head1 AUTHOR

Lei Gao, Shenzhen University, leigao@szu.edu.cn or highlei@hotmail.com

=head1 VERSION

1.0 : April 22, 2018

=head1 Easy INSTALL

     For administrator/root
	  perl MAKEFILE.pl -i /usr/local/bin/
     For other user, put all scripts to your PATH (check your PATH by "echo $PATH"), e.g.
	  perl MAKEFILE.pl -i ~/bin/

     Then you can call auspp directly, e.g. test the example:
	  auspp -M sRNAexample -e example/


=head1 INSTALL -- if "easy INSTALL" doesn't work, please check the following:

=head2 Dependencies - Linux/Unix OS platform
     auspp was developed on linux (Ubuntu 14.04.5 LTS), and hasn't been
     tested on other OS platform.

=head2 Dependencies - Required Perl in /usr/bin/perl
    auspp is a perl program, so it needs perl installed on your system. It
    was developed on perl version 5.14, and hasn't been tested on other
    versions (but there is no reason to suspect problems with other perl
    5.x versions). auspp will not compile. Getopt::Std, FileHandle, strict
    and Cwd 'abs_path' are pre-loaded into most (all?) Perl distros. auspp
    expects to find perl in /usr/bin/perl .. if not, edit line 1 
    accordingly for all the perl scripts in bin/

=head2 Dependencies - PATH executables

	one of the following short read aligners:
            bowtie  and bowtie-build
	or  bowtie2 and bowtie2-build
        or  soap    and 2bwt-builder
	or  hisat2  and hisat2-build
	or  bwa

	samtools (Version 1.0 +, http://www.htslib.org/)

	blast+, if use the step 4 (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/)

	bam2wig, if use the step 7 (https://github.com/MikeAxtell/bam2wig )

    All of the above must be executable from your PATH. Depending on the
    mode of the auspp run (see below), only a subset of these programs
    may be required for a given run.

=head2 Installation

    Except for the above dependencies, there is no "real" installation. Put
    all the scripts in bin/ into your working directory,  then you can call 
    auspp with

            ./auspp

    For convenience, you can add it to your PATH. e.g.

            sudo mv bin/* /usr/local/bin/

=head1 run

 auspp -i Col.fastq -x Col -D referenceIndex -M mRNA
 auspp -i Col.fastq -x Col -D referenceIndex -M smallRNA
 auspp -i Col.fastq -x Col -G referenceGENOME -M chip

=head1 USAGE

auspp includes auspp and 7 modules (other 7 perl scripts searchseq, searchLineACList, trim_adaptor, soap2sam_gl, 
collapseFasta, fastq2fasta, blast_m8). For convenience, you can check their USAGE by "-h", e.g.

			trim_adaptor -h

The 7 modules can be used individually for other purpose.

=head2 USAGE for auspp
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

=head2 Modules' descriptions
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

