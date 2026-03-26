#!/usr/bin/perl -w
use strict;

die "perl $0 <bedlist> <prefix> <outdir>\n" unless (@ARGV==3);

my $bedlist=shift;
my $prefix=shift;
my $outdir=shift;

my $admix="/EClab/Software/conda/Miniconda/envs/admix/bin/admixture";


open I,"$bedlist";

while (my $line=<I>){
	chomp($line);
	my $type=(split/\//,$line)[-2];
	`mkdir -p "$outdir/$type"` unless (-d "$outdir/$type");
	for my $i (2..20){
		open O,">$outdir/$type/step1.ad$i.$prefix.sh";
		print O "cd $outdir/$type\n";
		print O "$admix -j10 --cv $line $i | tee > $outdir/$type/log$i\n";
	}
}



