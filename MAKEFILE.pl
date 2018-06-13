#!/usr/bin/perl
# Copyright (c)  2018-
# Program:			MAKEFILE for auspp #
# Author:			Gaolei <highlei@hotmail.com or leigao@szu.edu.cn>
# Program Date:		2018.04.19
# Modifier:			Gaolei <highlei@hotmail.com or leigao@szu.edu.cn>
# Last Modified:	2018.04.19
# Description:	makefile for auspp
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
getopts("hi:");
my $infile	= $opt_i; # for one library, you can just give the fastq file
#my $jnfile	= $opt_I; # for paird-end library, the other mate
#my $inlist	= $opt_I; # for two or more libraries, you should give txt file which include the "fastq file\tprefix" with one row one library.
#my $prefix	= (defined $opt_x) ? $opt_x : "";#auspp_$start e.g. Col_r1";


if ($opt_h || ($infile eq "" && $prefix eq "")) {# && $inlist eq ""
	usage();
}
use FileHandle;
use strict;
use Cwd 'abs_path';

my $head = 100;
my $tail = 100;

my ($i,$j,$k,$m,$n,$k1,$k2,$k3,$file,$line,$count,$flag,$block,$a,$b,$end);
my (@buf,@tmp,@array,@ndir,@nlen,@rnum,@rstep);
my $ii=0;
my $jj=0;
my $kk=0;
my %snp=();
my %gnm=();
my %cro=();
my $key="";
my $dir =$infile;



#===========================================================================================================
#====================                  main
#===========================================================================================================
#print STDERR "Check the dependencies\n";
$i=readpipe("which blastn");
if ($i eq "") {
	system("cp depend/blastn $dir/.");
}
$i=readpipe("which 2bwt-builder");
if ($i eq "") {
	system("cp depend/2bwt-builder $dir/.");
}
$i=readpipe("which soap");
if ($i eq "") {
	system("cp depend/soap $dir/.");
}
$i=readpipe("which bam2wig");
if ($i eq "") {
	system("cp depend/bam2wig $dir/.");
}
$i=readpipe("which bowtie");
if ($i eq "") {
	system("cp depend/bowtie $dir/.");
}
$i=readpipe("which bowtie2");
if ($i eq "") {
	system("cp depend/bowtie2 $dir/.");
}
$i=readpipe("which bowtie-build");
if ($i eq "") {
	system("cp depend/bowtie-build $dir/.");
}
$i=readpipe("which bowtie2-build");
if ($i eq "") {
	system("cp depend/bowtie2-build $dir/.");
}
$i=readpipe("which bwa");
if ($i eq "") {
	system("cp depend/bwa $dir/.");
}
$i=readpipe("which hisat2");
if ($i eq "") {
	system("cp depend/hisat2 $dir/.");
}
$i=readpipe("which hisat2-build");
if ($i eq "") {
	system("cp depend/hisat2-build $dir/.");
}
$i=readpipe("which makeblastdb");
if ($i eq "") {
	system("cp depend/makeblastdb $dir/.");
}
$i=readpipe("which samtools");
if ($i eq "") {
	system("cp depend/samtools $dir/.");
}

#$i=readpipe("which auspp");
#if ($i eq "") {
#	system("cp bin/auspp $dir/auspp.test");
#}

system("cp bin/auspp bin/blast_m8 bin/collapseFasta bin/fastq2fasta bin/searchseq bin/searchLineACList bin/soap2sam_gl bin/trim_adaptor bin/auspp_lenDist.R bin/auspp_num.R $dir/.");

sub_end_program();


#############################################################################################################
####################################                                         ################################
####################################              "main end"                 ################################
####################################                                         ################################
#############################################################################################################
sub usage
{
	print "Program :\t$0\n";
	print "Version :\t$version\n";
	print "Author  :\tLei Gao, UC,Riverside\n";
	print "Contact :\tLei Gao <highlei\@gmail.com>\n";
	print "\nUsage:	$0 [options]\n";
	print "\t-i	<str>	input the  blast(m=8)/soap result file.";
	print " eg: blast_m8.blastn/soap.M0r2v0.soap\n";

#	print "\t-c	<str>	the parameter set for do_fold_ctFlt4miRNAInGnm.pl.";
#	print " [$fold_ctFlt]\n";

	print "\n\t-h	display this help\n";
#	print "		Note: please add quotation mark, if you input parameter in command line!\n";
	print "\nExample:\n";
	print "$0 -i blast_m8.blastn/soap.M0r2v0.soap \n";
	print ("==========================| $0  end   |=================================\n\n");

    exit(0);
}
############################################################################################################
######################                  sub_format_datetime
############################################################################################################

sub sub_format_datetime #时间子程序
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
	print STDERR ("==========================| $0  end  |==================================\n\n");
	exit(0);

}


