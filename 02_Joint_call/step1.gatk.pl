#!usr/bin/perl -w
use strict;
use List::MoreUtils qw(uniq);
use Getopt::Long;

#######################
#####USAGE
#########################

sub USAGE{
        my $usage=<<"USAGE";

Version:1.0
Mon Mar  4 16:02:53 CST 2024    linyu\@genomics.cn
        Usage:perl $0 -conf conflist -splitwindow 10000000 -sample_count 112 -o outdir
        help          <character>   : Show Help Information
        conf          <file>        : Formulation documents
        splitwindow   <num>         : Split Window Size 10000000
        sample_count  <num>         : Sample size
	o             <character>   : outdir

USAGE
print $usage;
exit;
}


my ($conflist,$splitwindow,$sample_count,$outdir);
GetOptions("help|?" =>\&USAGE,
           "conf:s" => \$conflist,
           "splitwindow:s" => \$splitwindow,
           "sample_count:s" => \$sample_count,
           "o:s" => \$outdir,
          ) or &USAGE;

&USAGE unless ($conflist && $splitwindow && $sample_count && $outdir);


`mkdir -p "$outdir/01.split_gvcf"` unless (-d "$outdir/01.split_gvcf");
`mkdir -p "$outdir/02.joint_call"` unless (-d "$outdir/02.joint_call");
`mkdir -p "$outdir/03.personalize_filter"` unless (-d "$outdir/03.personalize_filter");

####### Step1 READ Conf ###########

my %software_paths;

open Conf,"$conflist";

while (my $line=<Conf>){
	chomp($line);
	my ($software_name,$path)=(split/\=/,$line)[0,1];
	$software_paths{$software_name} = $path;
#	print $software_paths{$software_name};
}

close Conf;

###### Step2 READ gvcfList #####

open GVCF,$software_paths{"gvcflist"};

my $gvcf_raw_list;

while (my $GVCF_line=<GVCF>){
        chomp($GVCF_line);
	$gvcf_raw_list.="\t-variant $GVCF_line";
}

###### Step3 Split REF #########

open FAI,$software_paths{"fai"};

my $chr_stat_end_vcf_list;

while (my $FAI_line=<FAI>){
        chomp($FAI_line);
	my ($chrom,$size)=(split/\t/,$FAI_line)[0,1];
	for (my $start=1;$start<$size;$start=$start+$splitwindow){
		my $end=$start + $splitwindow -1;
		$end = $size  if $end >= $size;
		open O,">$outdir/01.split_gvcf/step1.$chrom-$start-$end.sh\n";
		print O  "echo ==========start at : `date` ==========\n";
		print O join("\t",$software_paths{"gatk"},"--java-options \"-Xmx100g\" CombineGVCFs -R",$software_paths{"ref"},"-L $chrom:$start-$end $gvcf_raw_list  -O $outdir/01.split_gvcf/$chrom-$start-$end.jointcall.g.vcf.gz"),"\n";
		print O join("\t",$software_paths{"gatk"},"--java-options \"-Xmx100g  -XX:ParallelGCThreads=8 -DGATK_STACKTRACE_ON_USER_EXCEPTION=true\"  GenotypeGVCFs -R",$software_paths{"ref"},"-V $outdir/01.split_gvcf/$chrom-$start-$end.jointcall.g.vcf.gz -O $outdir/01.split_gvcf/$chrom-$start-$end.raw.vcf.gz"),"\n";
		$chr_stat_end_vcf_list.="$outdir/01.split_gvcf/$chrom-$start-$end.raw.vcf.gz\n";
		print O  "echo ==========End at : `date` ==========\n";
	}
		open O2,">$outdir/02.joint_call/chr_stat_end_vcf_list";
                print O2 "$chr_stat_end_vcf_list";
}	
	
close FAI;
close O;
close O2;

###### Step4 Jointcall & Hard Filter & SNP #########

open O,">$outdir/02.joint_call/step1.merge.hardfilter.sh";

print O  "echo ==========start at : `date` ==========\n";
print O  join("\t",$software_paths{"gatk"},"--java-options \"-Xmx50g\" MergeVcfs -R",$software_paths{"ref"},"-I $outdir/02.joint_call/chr_stat_end_vcf_list -O $outdir/02.joint_call/All.chr.vcf.gz"),"\n";
print O  join("\t",$software_paths{"gatk"},"--java-options \"-Xmx40g\" SelectVariants --select-type-to-include SNP -R",$software_paths{"ref"},"-V $outdir/02.joint_call/All.chr.vcf.gz -O $outdir/02.joint_call/All.chr.snp.vcf.gz"),"\n";
print O  join("\t",$software_paths{"gatk"},"--java-options \"-Xmx50g\"  VariantFiltration -V $outdir/02.joint_call/All.chr.snp.vcf.gz  --filter-expression  \"QUAL < 30.0 \|\| QD < 2.0 \|\| MQ < 40.0 \|\| FS > 60.0 \|\| SOR > 3.0 \|\| HaplotypeScore > 13.0 \|\| MQRankSum < -12.5 \|\| ReadPosRankSum < -8.0\" --filter-name \"snp_filter\" -O $outdir/02.joint_call/All.chr.snp.hardmarked.vcf.gz"),"\n";
print O  join("\t",$software_paths{"gatk"},"SelectVariants --exclude-filtered true -R",$software_paths{"ref"},"-V $outdir/02.joint_call/All.chr.snp.hardmarked.vcf.gz -O $outdir/02.joint_call/All.chr.snp.hardmarked.hardfil.vcf.gz"),"\n";
print O "zcat $outdir/02.joint_call/All.chr.snp.hardmarked.hardfil.vcf.gz  \| awk \'\$5\!\~\"\,\"\' \| gzip -f >$outdir/02.joint_call/All.chr.snp.hardmarked.hardfil.select_snp.vcf.gz\n";
print O  "echo ==========End at : `date` ==========\n";

close O;

###### Step5 DP Stat & Fil DP <0.25 >0.95 &  PL<20 & Missing>=80%  #########

open O,">$outdir/03.personalize_filter/step1.DP.stat.fil.sh";

print O  "echo ==========start at : `date` ==========\n";
print O "gzip -cdf $outdir/02.joint_call/All.chr.snp.hardmarked.hardfil.select_snp.vcf.gz >$outdir/03.personalize_filter/All.chr.snp.hardfil.select_snp.vcf\n";
print O join("\t",$software_paths{"bgzip"},"-c $outdir/03.personalize_filter/All.chr.snp.hardfil.select_snp.vcf > $outdir/03.personalize_filter/All.chr.snp.hardfil.select_snp.vcf.gz"),"\n";
print O join("\t",$software_paths{"tabix"},"-p vcf $outdir/03.personalize_filter/All.chr.snp.hardfil.select_snp.vcf.gz"),"\n";
print O join("\t",$software_paths{"bcftools"},"query -f \'\%CHROM\\t\%POS\\t\%DP\\n\' $outdir/03.personalize_filter/All.chr.snp.hardfil.select_snp.vcf.gz  >$outdir/03.personalize_filter/All.chr.snp.hardfil.select_snp.vcf.gz.DP.stat.txt"),"\n";
print O "less -S $outdir/03.personalize_filter/All.chr.snp.hardfil.select_snp.vcf.gz.DP.stat.txt | wc -l >$outdir/03.personalize_filter/All.chr.snp.hardfil.select_snp.vcf.gz.DP.stat.txt.line.txt\n";
print O "gzip -cf $outdir/03.personalize_filter/All.chr.snp.hardfil.select_snp.vcf.gz.DP.stat.txt > $outdir/03.personalize_filter/All.chr.snp.hardfil.select_snp.vcf.gz.DP.stat.txt.gz\n";
print O join("\t","perl",$software_paths{"step1"},"-DP $outdir/03.personalize_filter/All.chr.snp.hardfil.select_snp.vcf.gz.DP.stat.txt.gz -inputvcf $outdir/03.personalize_filter/All.chr.snp.hardfil.select_snp.vcf.gz -outstat $outdir/03.personalize_filter/All.chr.snp.hardfil.select_snp.vcf.gz.DP.0.25-0.75.stat.txt -bgzip",$software_paths{"bgzip"},"-outvcf $outdir/03.personalize_filter/All.chr.snp.hardfil.select_snp.DPfil.vcf.gz -min 0.25 -max 99.75"),"\n";
print O join("\t","perl",$software_paths{"step2"},"$outdir/03.personalize_filter/All.chr.snp.hardfil.select_snp.DPfil.vcf.gz $outdir/03.personalize_filter/All.chr.snp.hardfil.select_snp.DPfil.PLfil.vcf.gz 0.2 20 $sample_count"),"\n";
print O  "echo ==========End at : `date` ==========\n";

close O;

###### Step6 fil Missing >0.8 & indel & maf based on chromosome #########

open FAI,$software_paths{"fai"};
#print $software_paths{"vcftools"};
my $chr_fil_vcf_list;
while (my $FAI_line=<FAI>){
        chomp($FAI_line);
        my ($chrom,$size)=(split/\t/,$FAI_line)[0,1];
	open O,">$outdir/03.personalize_filter/step2.filmis.indel.$chrom.sh";
	print O  "echo ==========Start at : `date` ==========\n";
	print O join("\t",$software_paths{"vcftools"},"--gzvcf $outdir/03.personalize_filter/All.chr.snp.hardfil.select_snp.DPfil.PLfil.vcf.gz --chr $chrom --max-missing 0.8 --recode --recode-INFO-all --out $outdir/03.personalize_filter/$chrom.chr.snp.hardfil.select_snp.DPfil.PLfil.misfil"),"\n";
	print O join("\t",$software_paths{"bgzip"},"-cf $outdir/03.personalize_filter/$chrom.chr.snp.hardfil.select_snp.DPfil.PLfil.misfil.recode.vcf	>$outdir/03.personalize_filter/$chrom.chr.snp.hardfil.select_snp.DPfil.PLfil.misfil.recode.vcf.gz"),"\n";
	print O join("\t",$software_paths{"tabix"},"-p vcf $outdir/03.personalize_filter/$chrom.chr.snp.hardfil.select_snp.DPfil.PLfil.misfil.recode.vcf.gz"),"\n";
	print O join("\t",$software_paths{"vcftools"},"--gzvcf $outdir/03.personalize_filter/$chrom.chr.snp.hardfil.select_snp.DPfil.PLfil.misfil.recode.vcf.gz --chr $chrom --min-alleles 2 --max-alleles 2 --remove-indels --recode --recode-INFO-all --out $outdir/03.personalize_filter/$chrom.chr.snp.hardfil.select_snp.DPfil.PLfil.misfil.bi.indelfil"),"\n";
	print O join("\t",$software_paths{"bgzip"},"-cf $outdir/03.personalize_filter/$chrom.chr.snp.hardfil.select_snp.DPfil.PLfil.misfil.bi.indelfil.recode.vcf > $outdir/03.personalize_filter/$chrom.chr.snp.hardfil.select_snp.DPfil.PLfil.misfil.bi.indelfil.recode.vcf.gz"),"\n";
	print O join("\t",$software_paths{"tabix"},"-p vcf $outdir/03.personalize_filter/$chrom.chr.snp.hardfil.select_snp.DPfil.PLfil.misfil.bi.indelfil.recode.vcf.gz"),"\n"; 
	print O join("\t",$software_paths{"vcftools"},"--gzvcf $outdir/03.personalize_filter/$chrom.chr.snp.hardfil.select_snp.DPfil.PLfil.misfil.bi.indelfil.recode.vcf.gz --chr $chrom  --maf",$software_paths{"maf"}," --recode --out $outdir/03.personalize_filter/$chrom.chr.snp.hardfil.select_snp.DPfil.PLfil.misfil.bi.indelfil.maffil --recode-INFO-all"),"\n";	
	print O join("\t",$software_paths{"bgzip"},"-cf $outdir/03.personalize_filter/$chrom.chr.snp.hardfil.select_snp.DPfil.PLfil.misfil.bi.indelfil.maffil.recode.vcf > $outdir/03.personalize_filter/$chrom.chr.snp.hardfil.select_snp.DPfil.PLfil.misfil.bi.indelfil.maffil.recode.vcf.gz"),"\n";
	print O join("\t",$software_paths{"tabix"},"-p vcf $outdir/03.personalize_filter/$chrom.chr.snp.hardfil.select_snp.DPfil.PLfil.misfil.bi.indelfil.maffil.recode.vcf.gz"),"\n";
	print O  "echo ==========End at : `date` ==========\n";	
	$chr_fil_vcf_list.="\t$outdir/03.personalize_filter/$chrom.chr.snp.hardfil.select_snp.DPfil.PLfil.misfil.bi.indelfil.maffil.recode.vcf.gz";	
}

close FAI;
close O;

###### Step7 Merge fil vcf  #########

open O,">$outdir/03.personalize_filter/step3.merge.filvcf.sh";

print O  "echo ==========Start at : `date` ==========\n";
print O join("\t",$software_paths{"bcftools"},"concat $chr_fil_vcf_list -Oz -o $outdir/03.personalize_filter/Whole.chr.snp.hardfil.select_snp.DPfil.PLfil.misfil.bi.indelfil.maffil.recode.$sample_count.vcf.gz"),"\n";
print O join("\t",$software_paths{"tabix"},"-p vcf $outdir/03.personalize_filter/Whole.chr.snp.hardfil.select_snp.DPfil.PLfil.misfil.bi.indelfil.maffil.recode.$sample_count.vcf.gz"),"\n";
print O  "echo ==========End at : `date` ==========\n";

close O;
