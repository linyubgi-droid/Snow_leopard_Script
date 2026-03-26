#!/usr/bin/perl -w
use strict;

die "perl $0 vcf_cout_list outdir \n" unless (@ARGV==2);

my $vcflist=shift;
my $outdir=shift;

#my $select_num=shift;

my $generation=100;
my $ref_count=2283352935;
my $refdir="/EClab/Project/Snow_leopard/Future/05.Know/31.slim/01.ref";
my $slim="/EClab/Project/fanyihang/software/miniconda3/bin/slim";
my $bcftools="/EClab/Software/conda/Miniconda/envs/normal/bin/bcftools";
my $java="/EClab/Software/conda/Miniconda/envs/snpEff/bin/java";
my $vcftools="/EClab/Software/conda/Miniconda/envs/normal/bin/vcftools";
my $snpeff="/EClab/Project/Snow_leopard/Future/05.Know/10.mutation_load_derived/06.ref_ancestor/snpEff.jar";
my $step1_slim="/EClab/Project/Snow_leopard/Future/06.Know_new/31.slim/00.Fanggui/SLIM5/02.slim/bin/moni.slim.simple.final.sh";
my $step1_lof="/EClab/Project/Snow_leopard/Future/06.Know_new/10.mutation_load_derived/08.SnpEff_ancestor_cds/bin/count_LoF.v3.pl";
my $step1_mis="/EClab/Project/Snow_leopard/Future/06.Know_new/10.mutation_load_derived/08.SnpEff_ancestor_cds/bin/count_mis.pl";
my $step1_lof_sum="/EClab/Project/Snow_leopard/Future/06.Know_new/31.slim/00.Fanggui/SLIM5/02.slim/bin/lof.sum.pl";
my $step1_mis_sum="/EClab/Project/Snow_leopard/Future/06.Know_new/31.slim/00.Fanggui/SLIM5/02.slim/bin/mis.sum.pl";
my $step1_het_sum="/EClab/Project/Snow_leopard/Future/06.Know_new/31.slim/00.Fanggui/SLIM5/02.slim/bin/het.sum.pl";
my $step2_damage_stat="/EClab/Project/Snow_leopard/Future/06.Know_new/31.slim/00.Fanggui/bin/slim/bin_mix/step2.count_damaging.pl";
my $step3_combine_vcf="/EClab/Project/Snow_leopard/Future/06.Know_new/31.slim/00.Fanggui/bin/slim/bin_mix/step1.com_vcf.pl";
my $step4_damage_all_stat="/EClab/Project/Snow_leopard/Future/06.Know_new/31.slim/00.Fanggui/bin/slim/bin_mix/step3.count_damaging.all.pl";
my $step5_stat_merge="/EClab/Project/Snow_leopard/Future/06.Know_new/31.slim/00.Fanggui/bin/slim/bin_mix/step4.stat_merge.pl";

open V,"$vcflist";

while (my $line=<V>){
	chomp($line);
	my @a1=(split/\s+/,$line);
	my $group=(split/\//,$a1[0])[-2];
	my $select_num=$a1[1];
	my $vcf=$a1[0];
	`mkdir -p "$outdir/$group"` unless (-d "$outdir/$group");
	open O,">$outdir/$group/step1.slim.sh";
	print O "cd $outdir/$group\n";
	print O "sh $step1_slim $vcf $group $refdir $outdir/$group $select_num $generation\n";
	open O1,">$outdir/$group/step2.$group.run.sh";
	print O1 "cd $outdir/$group\n";
        my $slim_script="slim_vcf_".$group."_gen".$generation.".slim";
        print O1 "$slim $outdir/$group/$slim_script  \> /dev/null 2\>\&1\n";
	open LOF,">$outdir/$group/lof.list\n";
	open MIS,">$outdir/$group/mis.list\n";
	open HET,">$outdir/$group/het.list\n";
	open O2,">$outdir/$group/step3.cal.sh";
	my @generations = (1,2,3,4,5,10,15,20,25,30,40,50,60,70,80,90,100);
	for my $k (@generations) {
		print O2 "$java  -jar $snpeff -no-downstream -no-upstream snow $outdir/$group/$group.gen$k.vcf | gzip -f > $outdir/$group/$group.gen$k.ann.vcf.gz\n";
		print O2 "perl $step1_lof $outdir/$group/$group.gen$k.ann.vcf.gz $outdir/$group/$group.gen$k.ann.vcf.gz.lof.count.txt  $outdir/$group/$group.gen$k.ann.vcf.gz.lof.pos.txt\n";
		print O2 "perl $step1_mis $outdir/$group/$group.gen$k.ann.vcf.gz $outdir/$group/$group.gen$k.ann.vcf.gz.mis.count.txt  $outdir/$group/$group.gen$k.ann.vcf.gz.mis.pos.txt\n";
		print O2 "$vcftools --vcf $outdir/$group/$group.gen$k.vcf --het --out $outdir/$group/$group.gen$k.vcf\n";
		print O2 "awk \'{print \$1\"\\t\"\$4\-\$2\"\\t\"\(\$4\-\$2\)\/$ref_count\"\\tSnow\"\}\' $outdir/$group/$group.gen$k.vcf.het \> $outdir/$group/$group.gen$k.vcf.het.txt\n";
		print HET "$outdir/$group/$group.gen$k.vcf.het.txt\n";
		print LOF "$outdir/$group/$group.gen$k.ann.vcf.gz.lof.count.txt\n";
		print MIS "$outdir/$group/$group.gen$k.ann.vcf.gz.mis.count.txt\n";
		#push @lof_list, "$outdir/$group/$group.gen$k.ann.vcf.gz.lof.count.txt" if $k <= $generation;
		#push @mis_list, "$outdir/$group/$group.gen$k.ann.vcf.gz.mis.count.txt" if $k <= $generation;

	}
	open SUM,">$outdir/$group/step3.sum.sh";
	print SUM "perl $step1_lof_sum $outdir/$group/lof.list $outdir/$group/lof.merge.txt\n";
	print SUM "perl $step1_mis_sum $outdir/$group/mis.list $outdir/$group/mis.merge.txt\n";
	print SUM "perl $step1_het_sum $outdir/$group/het.list $outdir/$group/het.merge.txt\n"; 
}


