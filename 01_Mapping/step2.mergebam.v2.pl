#!/usr/bin/perl -w 
use strict;
die "perl $0 bamlist ref.fa outdir genomesize\n" unless (@ARGV==4);
my $bamlist=shift;
my $ref=shift;
my $outdir=shift;
my $genome_size=shift;

my $step1_cov="/EClab/Project/Snow_leopard/Back/population_script/Lin/align/aln/bin/pipeline/step1.stat.coverage.pl";
my $java="/EClab/Software/conda/Miniconda/envs/picard/bin/java";
my $picard="/EClab/Software/conda/Miniconda/envs/picard/bin/../share/picard-3.3.0-0/picard.jar";
my $samtools="/EClab/Software/conda/Miniconda/envs/normal/bin/samtools";
my $gatk="/EClab/Software/conda/Miniconda/envs/normal/bin/gatk";
my $java_dir="/EClab/Software/conda/Miniconda/envs/picard/bin/";
my $tabix="/EClab/Software/conda/Miniconda/envs/normal/bin/tabix";
my $bgzip="/EClab/Software/conda/Miniconda/envs/normal/bin/bgzip";
my $BamDeal_Linux="/EClab/Software/BamDeal-0.27/bin/BamDeal_Linux";
my $bamtools="/EClab/Software/conda/Miniconda/envs/normal/bin/bamtools";

open B,"$bamlist";
my %sample_data;
my %sample_count;

while (my $line=<B>){
	chomp($line);
	my ($group,$sample,$lane)=(split/\//,$line)[-4,-3,-2];
	my $combine="$group\t$sample";	
	$sample_count{$combine}++;
	if (exists $sample_data{$combine}) {
                push @{$sample_data{$combine}}, $line;
        } else {
                $sample_data{$combine} = [$line];
        }
#	`mkdir -p "$outdir/$group/$sample"` unless (-e "$outdir/$group/$sample");
#	open O,">$outdir/$group/$sample/$lane/step1.$group.$sample.$lane.sh";
#	print O "perl $step1_cov $line $genome_size $outdir/$group/$sample/$lane/$group.$sample.$lane.cov.lane.stat.txt\n";
}	
	

foreach my $combine (keys %sample_data) {
	my ($group,$sample)=(split/\t/,$combine)[0,1];
	`mkdir -p "$outdir/$group/$sample"` unless (-e "$outdir/$group/$sample");
	open O,">$outdir/$group/$sample/step1.$group.$sample.merge.bam.sh";
	print O "echo ==========start at : `date` ==========\n";
	#	open O1,">$outdir/$group/$sample/step2.$group.$sample.gatk.sh";
	#print O1 "echo ==========start at : `date` ==========\n";
	#print O1 "export PATH\=$java_dir\:\$PATH\n";
	#print O1 "$gatk  \-\-java\-options \"\-Xmx100G\"  HaplotypeCaller  -R $ref -I $outdir/$group/$sample/$sample.merge.bam  -ERC GVCF -O $outdir/$group/$sample/$sample.g.vcf\n";
	#print O1 "$bgzip -f $outdir/$group/$sample/$sample.g.vcf\n";
	#print O1 "$tabix -p vcf $outdir/$group/$sample/$sample.g.vcf.gz\n";
	#print O1 "echo ==========end at : `date` ==========\n";
	my @F=@{$sample_data{$combine}};
	#print "@F\n";
	my $joined_string=join(",", @F);
	my @split_elements = split(",", $joined_string);
	#print "$split_elements[0]\n";
	if (exists $sample_count{$combine}){
		
		if ($sample_count{$combine}==1){
			print O "ln -sf $split_elements[0] $outdir/$group/$sample/$sample.merge.bam\n";
		#	print O "ln -s $split_elements[0].bai $outdir/$group/$sample/$sample.merge.bam.bai\n";
			print O "$samtools index $outdir/$group/$sample/$sample.merge.bam\n";
			print O "perl $step1_cov $outdir/$group/$sample/$sample.merge.bam $genome_size $outdir/$group/$sample/$sample.merge.bam.stat.txt\n";
			print O "$samtools flagstat $outdir/$group/$sample/$sample.merge.bam  > $outdir/$group/$sample/$sample.merge.bam.flagstat.txt\n";
			print O "$BamDeal_Linux statistics Coverage -i $outdir/$group/$sample/$sample.merge.bam -r $ref -o $outdir/$group/$sample/$sample.merge.bam.Coverage.txt\n";
			print O "$bamtools stats -in $outdir/$group/$sample/$sample.merge.bam \> $outdir/$group/$sample/$sample.merge.bam.qc\n";
		}elsif ($sample_count{$combine}==2){
			print O "$java -Xmx100G -Djava.io.tmpdir=$outdir/$group/$sample/tmp -XX:-UseGCOverheadLimit -jar $picard MergeSamFiles \\\n";
			print O "I=$split_elements[0] \\\n";
			print O "I=$split_elements[1] \\\n";
			print O "O=$outdir/$group/$sample/$sample.merge.bam\n";
			print O "$samtools index $outdir/$group/$sample/$sample.merge.bam\n";
			print O "perl $step1_cov $outdir/$group/$sample/$sample.merge.bam $genome_size $outdir/$group/$sample/$sample.merge.bam.stat.txt\n";
			print O "$samtools flagstat $outdir/$group/$sample/$sample.merge.bam  > $outdir/$group/$sample/$sample.merge.bam.flagstat.txt\n";
			print O "$BamDeal_Linux statistics Coverage -i $outdir/$group/$sample/$sample.merge.bam -r $ref -o $outdir/$group/$sample/$sample.merge.bam.Coverage.txt\n";
			print O "$bamtools stats -in $outdir/$group/$sample/$sample.merge.bam \> $outdir/$group/$sample/$sample.merge.bam.qc\n";
		}elsif ($sample_count{$combine}==3){
                        print O "$java -Xmx100G -Djava.io.tmpdir=$outdir/$group/$sample/tmp -XX:-UseGCOverheadLimit -jar $picard MergeSamFiles \\\n";
                        print O "I=$split_elements[0] \\\n";
                        print O "I=$split_elements[1] \\\n";
			print O "I=$split_elements[2] \\\n";
                        print O "O=$outdir/$group/$sample/$sample.merge.bam\n";
                        print O "$samtools index $outdir/$group/$sample/$sample.merge.bam\n";
			print O "perl $step1_cov $outdir/$group/$sample/$sample.merge.bam $genome_size $outdir/$group/$sample/$sample.merge.bam.stat.txt\n";
			print O "$samtools flagstat $outdir/$group/$sample/$sample.merge.bam  > $outdir/$group/$sample/$sample.merge.bam.flagstat.txt\n";
			print O "$BamDeal_Linux statistics Coverage -i $outdir/$group/$sample/$sample.merge.bam -r $ref -o $outdir/$group/$sample/$sample.merge.bam.Coverage.txt\n";
			print O "$bamtools stats -in $outdir/$group/$sample/$sample.merge.bam \> $outdir/$group/$sample/$sample.merge.bam.qc\n";
		}elsif ($sample_count{$combine}==4){
                        print O "$java -Xmx100G -Djava.io.tmpdir=$outdir/$group/$sample/tmp -XX:-UseGCOverheadLimit -jar $picard MergeSamFiles \\\n";
                        print O "I=$split_elements[0] \\\n";
                        print O "I=$split_elements[1] \\\n";
                        print O "I=$split_elements[2] \\\n";
			print O "I=$split_elements[3] \\\n";
                        print O "O=$outdir/$group/$sample/$sample.merge.bam\n";
                        print O "$samtools index $outdir/$group/$sample/$sample.merge.bam\n";
			print O "perl $step1_cov $outdir/$group/$sample/$sample.merge.bam $genome_size $outdir/$group/$sample/$sample.merge.bam.stat.txt\n";
			print O "$samtools flagstat $outdir/$group/$sample/$sample.merge.bam  > $outdir/$group/$sample/$sample.merge.bam.flagstat.txt\n";
			print O "$BamDeal_Linux statistics Coverage -i $outdir/$group/$sample/$sample.merge.bam -r $ref -o $outdir/$group/$sample/$sample.merge.bam.Coverage.txt\n";
			print O "$bamtools stats -in $outdir/$group/$sample/$sample.merge.bam \> $outdir/$group/$sample/$sample.merge.bam.qc\n";
		}elsif ($sample_count{$combine}==5){
                        print O "$java -Xmx100G -Djava.io.tmpdir=$outdir/$group/$sample/tmp -XX:-UseGCOverheadLimit -jar $picard MergeSamFiles \\\n";
                        print O "I=$split_elements[0] \\\n";
                        print O "I=$split_elements[1] \\\n";
                        print O "I=$split_elements[2] \\\n";
                        print O "I=$split_elements[3] \\\n";
			print O "I=$split_elements[4] \\\n";
                        print O "O=$outdir/$group/$sample/$sample.merge.bam\n";
                        print O "$samtools index $outdir/$group/$sample/$sample.merge.bam\n";
			print O "perl $step1_cov $outdir/$group/$sample/$sample.merge.bam $genome_size $outdir/$group/$sample/$sample.merge.bam.stat.txt\n";
			print O "$samtools flagstat $outdir/$group/$sample/$sample.merge.bam  > $outdir/$group/$sample/$sample.merge.bam.flagstat.txt\n";
			print O "$BamDeal_Linux statistics Coverage -i $outdir/$group/$sample/$sample.merge.bam -r $ref -o $outdir/$group/$sample/$sample.merge.bam.Coverage.txt\n";
			print O "$bamtools stats -in $outdir/$group/$sample/$sample.merge.bam \> $outdir/$group/$sample/$sample.merge.bam.qc\n";
                }elsif ($sample_count{$combine}==6){
                        print O "$java -Xmx100G -Djava.io.tmpdir=$outdir/$group/$sample/tmp -XX:-UseGCOverheadLimit -jar $picard MergeSamFiles \\\n";
                        print O "I=$split_elements[0] \\\n";
                        print O "I=$split_elements[1] \\\n";
                        print O "I=$split_elements[2] \\\n";
                        print O "I=$split_elements[3] \\\n";
                        print O "I=$split_elements[4] \\\n";
			print O "I=$split_elements[5] \\\n";
                        print O "O=$outdir/$group/$sample/$sample.merge.bam\n";
                        print O "$samtools index $outdir/$group/$sample/$sample.merge.bam\n";
			print O "perl $step1_cov $outdir/$group/$sample/$sample.merge.bam $genome_size $outdir/$group/$sample/$sample.merge.bam.stat.txt\n";
			print O "$samtools flagstat $outdir/$group/$sample/$sample.merge.bam  > $outdir/$group/$sample/$sample.merge.bam.flagstat.txt\n";
			print O "$BamDeal_Linux statistics Coverage -i $outdir/$group/$sample/$sample.merge.bam -r $ref -o $outdir/$group/$sample/$sample.merge.bam.Coverage.txt\n";
			print O "$bamtools stats -in $outdir/$group/$sample/$sample.merge.bam \> $outdir/$group/$sample/$sample.merge.bam.qc\n";
                }
		print O "echo ==========end at : `date` ==========\n";
        print O "echo ==========start at : `date` ==========\n";
        print O "export PATH\=$java_dir\:\$PATH\n";
        print O "$gatk  \-\-java\-options \"\-Xmx100G\"  HaplotypeCaller  -R $ref -I $outdir/$group/$sample/$sample.merge.bam  -ERC GVCF -O $outdir/$group/$sample/$sample.g.vcf\n";
        print O "$bgzip -f $outdir/$group/$sample/$sample.g.vcf\n";
        print O "$tabix -p vcf $outdir/$group/$sample/$sample.g.vcf.gz\n";
        print O "echo ==========end at : `date` ==========\n";
	}
}
	








