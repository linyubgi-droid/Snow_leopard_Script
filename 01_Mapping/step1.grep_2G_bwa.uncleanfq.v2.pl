#!/usr/bin/perl -w
use strict;
die "perl $0 <list> <ref.fa> <outdir>\n" unless (@ARGV==3);

my $bwa="/EClab/Software/conda/Miniconda/envs/normal/bin/bwa";
my $samtools="/EClab/Software/conda/Miniconda/envs/normal/bin/samtools";
my $fastp="/EClab/Software/conda/Miniconda/envs/normal/bin/fastp";
#my $REF="/EClab/Project/Snow_leopard/Post_editing/00.raw/Snow_leopard/00.REF/G_scaffolds.fasta";
my $java="/EClab/Software/conda/Miniconda/envs/picard/bin/java";
my $MkDup="/EClab/Software/conda/Miniconda/envs/picard/bin/../share/picard-3.3.0-0/picard.jar";
my $thread=16;
my $list=shift;
my $REF=shift;
my $outdir=shift;
my $adaptor1=shift || "AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA";
my $adaptor2=shift || "AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG";
my $bamtools="/EClab/Software/conda/Miniconda/envs/normal/bin/bamtools";

open I,"$list";

while (my $line=<I>){
	chomp($line);
	my ($Group,$Sample,$lane,$fq1,$fq2)=(split/\t/,$line)[1,2,3,4,5];
	#	my $RG_tag="\@RG\\tID\:$Sample\\tPL\:ILLUMINA\\tPU\:$Sample\\tLB\:$lane\\tSM\:$Sample\\tCN\=BGI";
	my $RG_tag="\@RG\\tID:$Sample\\tPL:BGI\\tLB:$lane\\tSM:$Sample";
	`mkdir -p "$outdir/01.raw/$Group/$Sample/$lane"` unless (-e "$outdir/01.raw/$Group/$Sample/$lane");
	`mkdir -p "$outdir/01.raw/$Group/$Sample/$lane/tmp"` unless (-e "$outdir/01.raw/$Group/$Sample/$lane/tmp");
	open O,">$outdir/01.raw/$Group/$Sample/$lane/step1.$Group-$Sample-$lane.qc.stat.sh";
	#	print O "#zcat $fq1 |head -n 80000000 |gzip ->$outdir/01.raw/$Group/$Sample/$lane/$Group-$Sample-$lane.fq1.gz\n";
	#print O "#zcat $fq2 |head -n 80000000 |gzip ->$outdir/01.raw/$Group/$Sample/$lane/$Group-$Sample-$lane.fq2.gz\n";
	print O  "echo ==========start at : `date` ==========\n";
	print O "ln -sf $fq1 $outdir/01.raw/$Group/$Sample/$lane/$Group-$Sample-$lane.fq1.gz\n";
	print O "ln -sf $fq2 $outdir/01.raw/$Group/$Sample/$lane/$Group-$Sample-$lane.fq2.gz\n";
	#	print O "$fastp -i $outdir/01.raw/$Group/$Sample/$lane/$Group-$Sample-$lane.fq1.gz -I $outdir/01.raw/$Group/$Sample/$lane/$Group-$Sample-$lane.fq2.gz -o $outdir/01.raw/$Group/$Sample/$lane/$Group-$Sample-$lane.fq1.clean.gz -O $outdir/01.raw/$Group/$Sample/$lane/$Group-$Sample-$lane.fq2.clean.gz --qualified_quality_phred=5    --unqualified_percent_limit=50 --n_base_limit=10 --adapter_sequence=\"$adaptor1\"    --adapter_sequence_r2=\"$adaptor2\"    --disable_trim_poly_g --thread=16 -j $outdir/01.raw/$Group/$Sample/$lane/$Group-$Sample-$lane.json -h $outdir/01.raw/$Group/$Sample/$lane/$Group-$Sample-$lane.html -R $Group-$Sample-$lane.Clean\n";
	print O  "echo ==========end at : `date` ==========\n";
	open O2,">$outdir/01.raw/$Group/$Sample/$lane/step2.$Group-$Sample-$lane.aln.sh";
	print O2  "echo ==========start at : `date` ==========\n";
	#	print O2 "$bwa mem -t $thread -M -R \'$RG_tag\' $REF $outdir/01.raw/$Group/$Sample/$lane/$Group-$Sample-$lane.fq1.clean.gz $outdir/01.raw/$Group/$Sample/$lane/$Group-$Sample-$lane.fq2.clean.gz \| awk \'\{if\(\!\/\^\@\/ \&\& and\(\$2\,4\) \=\= 4\)\{\$5\=0\; \$6\=\"\*\"\; gsub\(\/ \/\,\"\\t\"\,\$0\)\;\} print\}\' \| $samtools view \-Shb \-T $REF \- \> $outdir//01.raw/$Group/$Sample/$lane/$Group-$Sample-$lane.raw.bam\n";
	print O2 "$bwa mem -t $thread -M -R \'$RG_tag\' $REF $outdir/01.raw/$Group/$Sample/$lane/$Group-$Sample-$lane.fq1.gz $outdir/01.raw/$Group/$Sample/$lane/$Group-$Sample-$lane.fq2.gz \| $samtools view -bhS -o $outdir//01.raw/$Group/$Sample/$lane/$Group-$Sample-$lane.raw.bam\n";
	print O2 "$bamtools stats -in $outdir/01.raw/$Group/$Sample/$lane/$Group-$Sample-$lane.raw.bam \> $outdir/01.raw/$Group/$Sample/$lane/$Group-$Sample-$lane.raw.bam.qc\n";
	print O2 "echo raw.bam.finish\n";

	print O2 "$samtools sort -@ $thread -o $outdir/01.raw/$Group/$Sample/$lane/$Group-$Sample-$lane.raw.sorted.bam $outdir/01.raw/$Group/$Sample/$lane/$Group-$Sample-$lane.raw.bam\n";
	print O2 "$samtools index $outdir/01.raw/$Group/$Sample/$lane/$Group-$Sample-$lane.raw.sorted.bam\n";
	print O2 "$bamtools stats -in $outdir/01.raw/$Group/$Sample/$lane/$Group-$Sample-$lane.raw.sorted.bam \> $outdir/01.raw/$Group/$Sample/$lane/$Group-$Sample-$lane.raw.sorted.bam.qc\n";
	print O2 "echo raw.sorted.bam.finish\n";
	
	print O2 "$java -jar $MkDup MarkDuplicates \\\n";
	print O2 "I=$outdir/01.raw/$Group/$Sample/$lane/$Group-$Sample-$lane.raw.sorted.bam \\\n";
        print O2 "O=$outdir/01.raw/$Group/$Sample/$lane/$Group-$Sample-$lane.raw.sorted.mkdup.bam \\\n";
        print O2 "METRICS_FILE=$outdir/01.raw/$Group/$Sample/$lane/$Group-$Sample-$lane.raw.sorted.mkdup.bam.mat \\\n";
        print O2 "REMOVE_DUPLICATES=false\n";
	#MAX_FILE_HANDLES=1000 VALIDATION_STRINGENCY=SILENT\n";
	print O2 "$samtools index $outdir/01.raw/$Group/$Sample/$lane/$Group-$Sample-$lane.raw.sorted.mkdup.bam\n";
	print O2 "$samtools flagstat $outdir/01.raw/$Group/$Sample/$lane/$Group-$Sample-$lane.raw.sorted.mkdup.bam > $outdir/01.raw/$Group/$Sample/$lane/$Group-$Sample-$lane.raw.sorted.mkdup.bam.stat.txt\n";
	print O2 "$bamtools stats -in $outdir/01.raw/$Group/$Sample/$lane/$Group-$Sample-$lane.raw.sorted.mkdup.bam \> $outdir/01.raw/$Group/$Sample/$lane/$Group-$Sample-$lane.raw.sorted.mkdup.bam.qc\n";
	print O2 "rm $outdir//01.raw/$Group/$Sample/$lane/$Group-$Sample-$lane.raw.bam\n";
	print O2 "rm $outdir//01.raw/$Group/$Sample/$lane/$Group-$Sample-$lane.raw.bam.bai\n";
	print O2 "rm $outdir/01.raw/$Group/$Sample/$lane/$Group-$Sample-$lane.raw.sorted.bam\n";
	print O2 "rm $outdir/01.raw/$Group/$Sample/$lane/$Group-$Sample-$lane.raw.sorted.bam.bai\n";
	print O2  "echo ==========end at : `date` ==========\n";

	#print O "$java -Xmx10G -Djava.io.tmpdir=$outdir/01.raw/$Group/$Sample/$lane/tmp -XX:-UseGCOverheadLimit -jar $MkDup CleanSam \\\n";
	#print O "I=$outdir/01.raw/$Group/$Sample/$lane/$Group-$Sample-$lane.sort.mkdup.bam \\\n";
	#print O "O=$outdir/01.raw/$Group/$Sample/$lane/$Group-$Sample-$lane.sort.mkdup.adjust.bam\n";
	#print O "$samtools index $outdir/01.raw/$Group/$Sample/$lane/$Group-$Sample-$lane.sort.mkdup.adjust.bam\n";
	#print O "$samtools flagstat $outdir/01.raw/$Group/$Sample/$lane/$Group-$Sample-$lane.sort.mkdup.adjust.bam > $outdir/01.raw/$Group/$Sample/$lane/$Group-$Sample-$lane.sort.mkdup.adjust.bam.stat.txt\n";
}


close I;
close O;
close O2;




