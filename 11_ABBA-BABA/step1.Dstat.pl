#!/usr/bin/perl -w
use strict;
use File::Basename;

die "perl $0 bindir outdir\n" unless (@ARGV==2);
my $inputdir=shift;
my $outdir=shift;

my @a=glob("$inputdir/*/*popDtest.list");

my $genotypename="/EClab/Project/Snow_leopard/Future/06.Know_new/06.Dstat_69/01.convert2eigenstrat/f4.geno";
my $snpname="/EClab/Project/Snow_leopard/Future/06.Know_new/06.Dstat_69/01.convert2eigenstrat/f4.snp";
my $indivname="/EClab/Project/Snow_leopard/Future/06.Know_new/06.Dstat_69/01.convert2eigenstrat/f4.indi.ind";
#my $poplistname="/jdfssz1/ST_EARTH/P18Z10200N0113/PMO/big_cats-16.5T/Snow_leopard-2T/Two_Pop/08.D-test/03.D_test_individuals/bin/pop.ind.list";
my $DIR="/EClab/Software/conda/Miniconda/envs/admixtools/bin";
my $qpDstat="/EClab/Software/conda/Miniconda/envs/admixtools/bin/qpDstat";
for my $list (@a){
	#my $base=basename($list);
	#my $dirname=(split/\./,$base)[0];
	my $dirname=(split/\//,$list)[-2];
	my $dir=dirname($list);
	my $poplistname=$dir."/".$dirname.".ind.list";
	`mkdir -p "$outdir/$dirname"` unless (-d "$outdir/$dirname");
	open O,">$outdir/$dirname/$dirname.parameter.file";
	print O "genotypename: $genotypename\n";
	print O "snpname: $snpname\n";
	print O "indivname: $indivname\n";
	print O "poplistname: $poplistname\n";
	print O "popfilename: $list\n";
	print O "DIR: $DIR\n";
	print O "printsd: YES\n";
	open O1,">$outdir/$dirname/step1.$dirname.calD.sh";
	print O1 "$qpDstat -p $outdir/$dirname/$dirname.parameter.file > $outdir/$dirname/$dirname.D_test.txt\n";
}

close O;
close O1;
