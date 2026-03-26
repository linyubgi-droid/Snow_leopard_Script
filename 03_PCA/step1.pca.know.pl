#!/usr/bin/perl -w
use strict;

die "perl $0 <VCF.list> <Prefix> <Outdir>\n" unless (@ARGV==3);

my $vcf_list=shift;
my $prefix=shift;
my $outdir=shift;

my $bcftools="/EClab/Software/conda/Miniconda/envs/normal/bin/bcftools";
my $tabix="/EClab/Software/conda/Miniconda/envs/normal/bin/tabix";
my $vcftools="/EClab/Software/conda/Miniconda/envs/normal/bin/vcftools";
my $Rscript="/EClab/Software/conda/Miniconda/envs/ggrepel/bin/Rscript";
my $plink="/EClab/Software/conda/Miniconda/envs/normal/bin/plink"; ## PLINK v1.9.0-b.8
my $king="/EClab/Software/conda/Miniconda/envs/normal/bin/king";#KING 2.2.7
my $screen="/EClab/Project/USER/linyu/population_script/VCF/02.screen/step4.screen.pl";
#my $draw="/EClab/Project/Snow_leopard/Future/03.Re_seq/01.PCA/01.All_maf05/snow.ggplot.v4.R";
my $draw_label="/EClab/Project/USER/linyu/population_script/PCA/step1.snow.pca.label.R";
my $draw_nolabel="/EClab/Project/USER/linyu/population_script/PCA/step1.snow.pca.nolabel.R";


open I,"$vcf_list";

while (my $line=<I>){
	chomp($line);
	my ($type)=(split/\//,$line)[-2];
	`mkdir -p "$outdir/$type"` unless (-d "$outdir/$type");
	open O1,">$outdir/$type/step2.$type.pca.sh";
	print O1 "cd $outdir/$type\n";
	print O1 "$vcftools --gzvcf $line --plink --out $outdir/$type/$prefix\n";
	print O1 "$plink --allow-extra-chr --noweb --file $outdir/$type/$prefix --make-bed --out $outdir/$type/$prefix\n";
	print O1 "$plink --allow-extra-chr  --threads 6 -bfile $outdir/$type/$prefix --out $outdir/$type/$prefix --pca 10\n";
	print O1 "sed 's#SNM#UN#g' $outdir/$type/$prefix.eigenvec >$outdir/$type/$prefix.eigenvec2\n";
	open O2,">$outdir/$type/step3.$type.draw.sh";
	#	print O2 "$Rscript $draw $outdir/$type/$prefix.eigenval $outdir/$type/$prefix.eigenvec2 $outdir/$type $prefix\n";
	print O2 "$Rscript $draw_label $outdir/$type/$prefix.eigenval $outdir/$type/$prefix.eigenvec $outdir/$type $prefix\n";
	print O2 "$Rscript $draw_nolabel $outdir/$type/$prefix.eigenval $outdir/$type/$prefix.eigenvec $outdir/$type $prefix.nolabel\n";
	open O3,">$outdir/$type/step3.$type.king.sh";
	print O3 "$king -b $outdir/$type/$prefix.bed --kinship --prefix $outdir/$type/$prefix.kinship\n";
	print O3 "$king -b $outdir/$type/$prefix.bed --ibdseg --prefix $outdir/$type/$prefix.ibd\n";
}



	
	
