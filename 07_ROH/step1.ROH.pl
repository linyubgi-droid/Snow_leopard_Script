#!/usr/bin/perl -w
use strict;
use List::MoreUtils qw(uniq);
use Getopt::Long;

#######################
#####USAGE
#########################

sub USAGE{
        my $usage=<<"USAGE";

Version:1.0
Sun Aug 31 02:30:44 PM CST 2025    linyubgi\@163.com
        Usage:perl $0 --plink plink --Rscript Rscript --step1 step1.roh_length.pl --step2 step2.roh_length_combine.pl --step3 step3.Froh_combine.pl --step4 step4.getind.overlength.pl --step5 step5.cal_froh_generation.pl --step6 step6.cal_froh_generation_Ne.pl --step7 step7.froh.generation.R --step8 step8.froh.generation.Ne.R --homozyg-window-snp 20 --homozyg-kb 10 --homozyg-density 50  -autosome_length 2246314264 -CM 1.9 -vcf test.vcf.gz -o outdir
        help               <character>   : Show Help Information
	plink	           <character>   : the path of plink (recommend v1.9.0-b.8)
	Rscript            <character>   : the path of Rscript (should install ggplot2;RColorBrewer;reshape2;dplyr;stringr;ggpubr)
	step1		   <character>   : the path of step1.roh_length.pl
	step2		   <character>   : the path of step2.roh_length_combine.pl
	step3		   <character>   : the path of step3.Froh_combine.pl
	step4	           <character>   : the path of step4.getind.overlength.pl (we set 100kb in this code)
	step5		   <character>   : the path of step5.cal_froh_generation.pl
	step6		   <character>   : the path of step6.cal_froh_generation_Ne.pl
	step7              <character>   : the path of step7.froh.generation.R
	step8              <character>   : the path of step8.froh.generation.Ne.R
	homozyg-window-snp <num>         : scanning window size 
	homozyg-kb         <num>         : min length 
	homozyg-density    <num>         : max inverse density (kb/var) 
	autosome_length    <num>         : the length of autosome
	CM		   <num>         : estimated based on a recombination rate of  cM/Mb (eg. Cat 1.9)
	vcf                <character>   : the path of vcf(.gz)
	group		   <character>   : the group name
        o                  <character>   : outdir

USAGE
print $usage;
exit;
}


my ($plink,$Rscript,$step1,$step2,$step3,$step4,$step5,$step6,$step7,$step8,$homozyg_window_snp,$homozyg_kb,$homozyg_density,$autosome_length,$CM,$vcf,$group,$outdir);
GetOptions("help|?" =>\&USAGE,
	   "plink:s" => \$plink,
	   "Rscript:s" => \$Rscript,
	   "step1:s" => \$step1,
	   "step2:s" => \$step2,
	   "step3:s" => \$step3,
	   "step4:s" => \$step4,
	   "step5:s" => \$step5,
           "step6:s" => \$step6,
	   "step7:s" => \$step7,
           "step8:s" => \$step8,
	   "homozyg-window-snp:s" => \$homozyg_window_snp,
	   "homozyg-kb:s" => \$homozyg_kb,
	   "homozyg-density:s" => \$homozyg_density,
	   "autosome_length:s" => \$autosome_length,
	   "CM:s" => \$CM,
	   "vcf:s" => \$vcf,
	   "group:s" => \$group,
           "o:s" => \$outdir,
          ) or &USAGE;

&USAGE unless ($plink && $Rscript && $step1 && $step2 && $step3 && $step4 && $step5 && $step6 && $step7 && $step8 && $homozyg_window_snp && $homozyg_kb && $homozyg_density && $autosome_length && $CM && $vcf && $group && $outdir);

`mkdir -p "$outdir/00.plink/$group"` unless (-d "$outdir/00.plink/$group");
`mkdir -p "$outdir/01.ROH_LEN/$group"` unless (-d "$outdir/01.ROH_LEN/$group");
`mkdir -p "$outdir/02.ROH_FROH/$group"` unless (-d "$outdir/02.ROH_FROH/$group");
`mkdir -p "$outdir/03.roh_ind/$group"` unless (-d "$outdir/03.roh_ind/$group");

open O1, ">$outdir/00.plink/$group/step1.$group.plink.sh";

print O1 "$plink --double-id --allow-extra-chr -vcf $vcf --make-bed --out $outdir/00.plink/$group/$group\n";
print O1 "$plink --double-id --allow-extra-chr --bfile $outdir/00.plink/$group/$group --homozyg --homozyg-window-snp $homozyg_window_snp --homozyg-kb $homozyg_kb --homozyg-density $homozyg_density --out $outdir/00.plink/$group/$group\n";
print O1 "sed \-i \'s\/\^\[\[\:blank\:\]\]\*\/\/\' $outdir/00.plink/$group/$group.hom\n";

my @LEN_arry;
my @FROH_arry;

open O2, ">$outdir/01.ROH_LEN/$group/step2.$group.LEN.sh";


print O2 "perl $step1 $outdir/00.plink/$group/$group.hom 100000 999999999999999 $group $outdir/01.ROH_LEN/$group/$group.ROH.100kb.Length.txt\n";
print O2 "perl $step1 $outdir/00.plink/$group/$group.hom 1 100000 $group $outdir/01.ROH_LEN/$group/$group.ROH.1bp-100kb.Length.txt\n";
print O2 "perl $step1 $outdir/00.plink/$group/$group.hom 100000 1000000 $group $outdir/01.ROH_LEN/$group/$group.ROH.100kb-1MB.Length.txt\n";
print O2 "perl $step1 $outdir/00.plink/$group/$group.hom 1000000 3000000 $group $outdir/01.ROH_LEN/$group/$group.ROH.1MB-3MB.Length.txt\n";
print O2 "perl $step1 $outdir/00.plink/$group/$group.hom 1000000 5000000 $group $outdir/01.ROH_LEN/$group/$group.ROH.1MB-5MB.Length.txt\n";
print O2 "perl $step1 $outdir/00.plink/$group/$group.hom 3000000 5000000 $group $outdir/01.ROH_LEN/$group/$group.ROH.3MB-5MB.Length.txt\n";
print O2 "perl $step1 $outdir/00.plink/$group/$group.hom 5000000 10000000 $group $outdir/01.ROH_LEN/$group/$group.ROH.5MB-10MB.Length.txt\n";
print O2 "perl $step1 $outdir/00.plink/$group/$group.hom 10000000 15000000 $group $outdir/01.ROH_LEN/$group/$group.ROH.10MB-15MB.Length.txt\n";
print O2 "perl $step1 $outdir/00.plink/$group/$group.hom 15000000 20000000 $group $outdir/01.ROH_LEN/$group/$group.ROH.15MB-20MB.Length.txt\n";
print O2 "perl $step1 $outdir/00.plink/$group/$group.hom 1000000 999999999999999 $group $outdir/01.ROH_LEN/$group/$group.ROH.1MB.Length.txt\n";
print O2 "perl $step1 $outdir/00.plink/$group/$group.hom 2000000 999999999999999 $group $outdir/01.ROH_LEN/$group/$group.ROH.2MB.Length.txt\n";
print O2 "perl $step1 $outdir/00.plink/$group/$group.hom 3000000 999999999999999 $group $outdir/01.ROH_LEN/$group/$group.ROH.3MB.Length.txt\n";
print O2 "perl $step1 $outdir/00.plink/$group/$group.hom 5000000 999999999999999 $group $outdir/01.ROH_LEN/$group/$group.ROH.5MB.Length.txt\n";
print O2 "perl $step1 $outdir/00.plink/$group/$group.hom 10000000 999999999999999 $group $outdir/01.ROH_LEN/$group/$group.ROH.10MB.Length.txt\n";
print O2 "perl $step1 $outdir/00.plink/$group/$group.hom 20000000 999999999999999 $group $outdir/01.ROH_LEN/$group/$group.ROH.20MB.Length.txt\n";

push @LEN_arry,"$outdir/01.ROH_LEN/$group/$group.ROH.1bp-100kb.Length.txt";
push @LEN_arry,"$outdir/01.ROH_LEN/$group/$group.ROH.100kb-1MB.Length.txt";
push @LEN_arry,"$outdir/01.ROH_LEN/$group/$group.ROH.1MB-3MB.Length.txt";
push @LEN_arry,"$outdir/01.ROH_LEN/$group/$group.ROH.1MB-5MB.Length.txt";
push @LEN_arry,"$outdir/01.ROH_LEN/$group/$group.ROH.5MB-10MB.Length.txt";
push @LEN_arry,"$outdir/01.ROH_LEN/$group/$group.ROH.10MB-15MB.Length.txt";
push @LEN_arry,"$outdir/01.ROH_LEN/$group/$group.ROH.15MB-20MB.Length.txt";
push @LEN_arry,"$outdir/01.ROH_LEN/$group/$group.ROH.1MB.Length.txt";
push @LEN_arry,"$outdir/01.ROH_LEN/$group/$group.ROH.2MB.Length.txt";
push @LEN_arry,"$outdir/01.ROH_LEN/$group/$group.ROH.3MB.Length.txt";
push @LEN_arry,"$outdir/01.ROH_LEN/$group/$group.ROH.5MB.Length.txt";
push @LEN_arry,"$outdir/01.ROH_LEN/$group/$group.ROH.10MB.Length.txt";
push @LEN_arry,"$outdir/01.ROH_LEN/$group/$group.ROH.20MB.Length.txt";

open O2_1,">$outdir/01.ROH_LEN/$group/$group.ROH.LENG.list";

print O2_1 join("\n",@LEN_arry);

open O2_2,">$outdir/01.ROH_LEN/$group/step3.$group.combineLEN.sh";

print O2_2 "perl $step2 $outdir/01.ROH_LEN/$group/$group.ROH.LENG.list $outdir/01.ROH_LEN/$group/$group.ROH.All.LENG.txt\n";




open O3,">$outdir/02.ROH_FROH/$group/step3.$group.FROH.sh";

print O3 "awk \'NR\=\=1 \{print \$1\"\\t\"\$2\"\\t\"\$NF\"\-FROH\"\} NR!=1 {print \$1\"\\t\"\$2\"\\t\"\$NF\/$autosome_length\}\' $outdir/01.ROH_LEN/$group/$group.ROH.100kb.Length.txt \> $outdir/02.ROH_FROH/$group/$group.ROH.100kb.FROH.txt\n";
print O3 "awk \'NR\=\=1 \{print \$1\"\\t\"\$2\"\\t\"\$NF\"\-FROH\"\} NR!=1 {print \$1\"\\t\"\$2\"\\t\"\$NF\/$autosome_length\}\' $outdir/01.ROH_LEN/$group/$group.ROH.100kb-1MB.Length.txt \> $outdir/02.ROH_FROH/$group/$group.ROH.100kb-1MB.FROH.txt\n";
print O3 "awk \'NR\=\=1 \{print \$1\"\\t\"\$2\"\\t\"\$NF\"\-FROH\"\} NR!=1 {print \$1\"\\t\"\$2\"\\t\"\$NF\/$autosome_length\}\' $outdir/01.ROH_LEN/$group/$group.ROH.1MB.Length.txt \> $outdir/02.ROH_FROH/$group/$group.ROH.1MB.FROH.txt\n";
print O3 "awk \'NR\=\=1 \{print \$1\"\\t\"\$2\"\\t\"\$NF\"\-FROH\"\} NR!=1 {print \$1\"\\t\"\$2\"\\t\"\$NF\/$autosome_length\}\' $outdir/01.ROH_LEN/$group/$group.ROH.2MB.Length.txt \> $outdir/02.ROH_FROH/$group/$group.ROH.2MB.FROH.txt\n";
print O3 "awk \'NR\=\=1 \{print \$1\"\\t\"\$2\"\\t\"\$NF\"\-FROH\"\} NR!=1 {print \$1\"\\t\"\$2\"\\t\"\$NF\/$autosome_length\}\' $outdir/01.ROH_LEN/$group/$group.ROH.3MB.Length.txt \> $outdir/02.ROH_FROH/$group/$group.ROH.3MB.FROH.txt\n";
print O3 "awk \'NR\=\=1 \{print \$1\"\\t\"\$2\"\\t\"\$NF\"\-FROH\"\} NR!=1 {print \$1\"\\t\"\$2\"\\t\"\$NF\/$autosome_length\}\' $outdir/01.ROH_LEN/$group/$group.ROH.5MB.Length.txt \> $outdir/02.ROH_FROH/$group/$group.ROH.5MB.FROH.txt\n";
print O3 "awk \'NR\=\=1 \{print \$1\"\\t\"\$2\"\\t\"\$NF\"\-FROH\"\} NR!=1 {print \$1\"\\t\"\$2\"\\t\"\$NF\/$autosome_length\}\' $outdir/01.ROH_LEN/$group/$group.ROH.10MB.Length.txt \> $outdir/02.ROH_FROH/$group/$group.ROH.10MB.FROH.txt\n";
print O3 "awk \'NR\=\=1 \{print \$1\"\\t\"\$2\"\\t\"\$NF\"\-FROH\"\} NR!=1 {print \$1\"\\t\"\$2\"\\t\"\$NF\/$autosome_length\}\' $outdir/01.ROH_LEN/$group/$group.ROH.20MB.Length.txt \> $outdir/02.ROH_FROH/$group/$group.ROH.20MB.FROH.txt\n";

push @FROH_arry,"$outdir/02.ROH_FROH/$group/$group.ROH.100kb.FROH.txt";
push @FROH_arry,"$outdir/02.ROH_FROH/$group/$group.ROH.1MB.FROH.txt";
push @FROH_arry,"$outdir/02.ROH_FROH/$group/$group.ROH.2MB.FROH.txt";
push @FROH_arry,"$outdir/02.ROH_FROH/$group/$group.ROH.3MB.FROH.txt";
push @FROH_arry,"$outdir/02.ROH_FROH/$group/$group.ROH.5MB.FROH.txt";
push @FROH_arry,"$outdir/02.ROH_FROH/$group/$group.ROH.10MB.FROH.txt";
push @FROH_arry,"$outdir/02.ROH_FROH/$group/$group.ROH.20MB.FROH.txt";

open O3_1,">$outdir/02.ROH_FROH/$group/$group.Froh.pct.list";

print O3_1 join("\n",@FROH_arry);

open O3_2,">$outdir/02.ROH_FROH/$group/step4.$group.combine.froh.sh";

print O3_2 "perl $step3 $outdir/02.ROH_FROH/$group/$group.Froh.pct.list $outdir/02.ROH_FROH/$group/$group.Froh.all.txt\n";

`mkdir -p "$outdir/03.roh_ind/$group"` unless (-d "$outdir/03.roh_ind/$group");

my $homozyg_roh_length=100000;

open O4_1,">$outdir/03.roh_ind/$group/step1.getind.roh.sh";
print O4_1 "perl $step4 $outdir/00.plink/$group/$group.hom $homozyg_roh_length $group $outdir/03.roh_ind/$group\n";

`mkdir -p "$outdir/04.froh_generation/$group"` unless (-d "$outdir/04.froh_generation/$group");

open O4_2,">$outdir/04.froh_generation/$group/step1.generation.froh.sh";
print O4_2 "perl $step5 $outdir/03.roh_ind/$group/$group.roh_hom.list $CM $autosome_length $group $outdir/04.froh_generation/$group/$group.froh.generation.txt $outdir/04.froh_generation/$group/$group.froh.generation.stat.txt\n";
print O4_2 "perl $step6 $outdir/03.roh_ind/$group/$group.roh_hom.list $CM $autosome_length $group $outdir/04.froh_generation/$group/$group.froh.generation.Ne.txt $outdir/04.froh_generation/$group/$group.froh.generation.Ne.stat.txt\n";

open O4_3,">$outdir/04.froh_generation/$group/step2.generation.froh.draw.sh";
print O4_3 "$Rscript $step7 $outdir/04.froh_generation/$group $group.froh.generation.txt $group.froh.generation.stat.txt $group\n";
print O4_3 "$Rscript $step8 $outdir/04.froh_generation/$group $group.froh.generation.Ne.txt $group.froh.generation.Ne.stat.txt $group\n";
