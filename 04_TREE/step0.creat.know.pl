#!/usr/bin/perl -w
use strict;

die "perl $0 <VCFlist> <outdir> <prefix>\n" unless (@ARGV==3);

my $list=shift;
my $outdir=shift;
my $prefix=shift;

my $vcftools="/EClab/Software/conda/Miniconda/envs/normal/bin/vcftools"; ## v 0.1.16
my $python="/EClab/Software/conda/Miniconda/envs/normal/bin/python"; ## 3.12.3 
my $vcf2phylip="/EClab/Software/vcf2phylip-2.7/vcf2phylip.py"; ## 2.7 
my $iqtree="/EClab/Software/conda/Miniconda/envs/normal/bin/iqtree"; ### 2.3.6 
my $phylo_tree="/EClab/Project/USER/linyu/population_script/TREE/NJ/Tree_pipeline/03_Phylogeny/bin/phylo_tree.pl"; #### BGI

open I,"$list";

while (my $line=<I>){
	chomp($line);
	my $type=(split/\//,$line)[-2];
	`mkdir -p "$outdir/$type"` unless (-d "$outdir/$type");
	open O,">$outdir/$type/step1.$type.thin.tree.sh";
	print O "cd $outdir/$type\n";
	print O "$vcftools --gzvcf $line --thin 2000 --stdout --recode --recode-INFO-all > $outdir/$type/$prefix.thin2000.vcf\n";
	print O "$python $vcf2phylip -i $outdir/$type/$prefix.thin2000.vcf --output-folder $outdir/$type --output-prefix $outdir/$type/$prefix.out\n";
	print O "$iqtree -s $outdir/$type/$prefix.out.min4.phy -nt 8 -bb 1000 -alrt 1000 -st DNA -m MFP --prefix $outdir/$type/$prefix.out\n";
	`mkdir -p "$outdir/$type/NJ"` unless (-d "$outdir/$type/NJ");
	open O1,">$outdir/$type/NJ/step2.$type.NJ.tree.sh";
	print O1 "perl $phylo_tree $outdir/$type/$prefix.out.min4.phy -type nj -d nt -b 100 -outdir $outdir/$type/NJ\n";
	`mkdir -p "$outdir/$type/Bayes"` unless (-d "$outdir/$type/Bayes");
	open O2,">$outdir/$type/Bayes/step3.$type.Bayes.tree.sh";
	print O2 "perl $phylo_tree $outdir/$type/$prefix.out.min4.phy -type bayes -data dna -outdir $outdir/$type/Bayes\n";

}

