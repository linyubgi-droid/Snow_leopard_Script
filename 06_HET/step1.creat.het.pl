#!/usr/bin/perl -w
use strict;

die "perl $0 <vcflist> <prefix> <length> <outdir>\n" unless (@ARGV==4);

my $vcflist=shift;
my $prefix=shift;
my $length=shift;
my $outdir=shift;

my $vcftools="/EClab/Software/conda/Miniconda/envs/normal/bin/vcftools";

open I,"$vcflist";

while (my $line=<I>){
	chomp($line);
	my ($type)=(split/\//,$line)[-3];
        `mkdir -p "$outdir/$type"` unless (-d "$outdir/$type");
        open O,">$outdir/$type/step1.$type.het.sh";
	print O "$vcftools --gzvcf $line --het --out $outdir/$type/$prefix\n";
	print O "awk \'\{print \$1\"\\t\"\$4\-\$2\"\\t\"\(\$4\-\$2\)\/$length\"\\t$prefix\"\}\' $outdir/$type/$prefix.het  > $outdir/$type/$prefix.het.txt\n";
}
