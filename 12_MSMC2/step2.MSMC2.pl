#!/usr/bin/perl -w
use strict;
use File::Basename;
die "perl $0 bindir phasedir outdir outdirne outdirsplit" unless (@ARGV==5);
my $inputdir=shift;
my $phasedir=shift;
my $outdir=shift;
my $outdir2=shift;
my $outdir3=shift;

my $bgzip="/EClab/Software/conda/Miniconda/envs/normal/bin/bgzip";
my $tabix="/EClab/Software/conda/Miniconda/envs/normal/bin/tabix";
my $java="/EClab/Software/conda/Miniconda/envs/normal/bin/java";
##my $beagle="/hwfssz4/BC_PUB/Software/03.Soft_ALL/beagle-5.0/beagle.28Sep18.793.jar";
#my $split_script="/jdfssz1/ST_EARTH/P18Z10200N0113/PMO/big_cats-16.5T/Snow_leopard-2T/Two_Pop/05.MSMC/01.phased_vcf/bin/split_individual.pl";
my $python="/EClab/Software/conda/Miniconda/envs/normal/bin/python";
my $generate_script="/EClab/Software/conda/Miniconda/envs/msmc/download/msmc-tools-master/generate_multihetsep.py";
my @list=glob("$inputdir/*re2.com.txt");
my $maskdir="/EClab/Project/Snow_leopard/Future/06.Know_new/14.MSMC_69/00.mask/00.raw";
my $msmc2="/EClab/Software/conda/Miniconda/envs/msmc/bin/msmc2_Linux";
my $boost_script="/EClab/Software/conda/Miniconda/envs/msmc/download/msmc-tools-master/multihetsep_bootstrap.py";

for my $l (@list){
	my $basename=basename($l);
	my $groupname=(split/\./,$basename)[0];
	open I,"$l";
	my $num=0;
	while (my $line=<I>){
		chomp($line);
		my ($One,$Two)=(split/\t/,$line)[0,1];
		$num++;
		my $Group1=$groupname;my $Group2=$groupname;
		`mkdir -p "$outdir/$groupname/G$num"` unless (-d "$outdir/$groupname/G$num");
		`mkdir -p "$outdir2/$groupname/G$num"` unless (-d "$outdir2/$groupname/G$num");	
		`mkdir -p "$outdir3/$groupname/G$num"` unless (-d "$outdir3/$groupname/G$num");
		`mkdir -p "$outdir2/$groupname/G$num-boost"` unless (-d "$outdir2/$groupname/G$num-boost"); 
		open O,">$outdir/$groupname/G$num/step1.$groupname.G$num.getinput.sh";
		open O2,">$outdir2/$groupname/G$num/step2.$groupname.G$num.cal.Ne.sh";
		open O3,">$outdir2/$groupname/G$num-boost/step1.$groupname.G$num.boost.sh";
#		open O3,">$outdir3/$groupname/G$num/step3.$groupname.G$num.cal.split.sh";	
		print O2 "$msmc2 --skipAmbiguous -I \'0\,1\,2\,3\' -o $outdir2/$groupname/G$num/$groupname.G$num -i 20 -t 6 -p \'1\*2\+15\*1\+1\*2\' \\\n";
		print O3 "$python $boost_script\t";
                print O3 "-n 20 -s 1000000 --chunks_per_chromosome 20 --nr_chromosomes 18 --seed 123 $outdir2/$groupname/G$num-boost/G$num \\\n";
#		print O3 "$msmc2 --skipAmbiguous -I 0-4,0-5,0-6,0-7,1-4,1-5,1-6,1-7,2-4,2-5,2-6,2-7,3-4,3-5,3-6,3-7 -o $outdir3/$groupname/G$num/$groupname.G$num -i 20 -t 6 -p \'1\*2\+15\*1\+1\*2\' \\\n";
		for my $i (1..18){
#			my $bed="$maskdir/chr".$i."_mask.bed.gz";
#			open O,">$outdir/$groupname/$line/step1.$groupname.$line.chr$i.split.sh";
			print O "$python $generate_script\t";
			#			print O "--mask $maskdir/$Group1/$One/$Group1.$One.chr$i.mask.bed.gz\t";
			#print O "--mask $maskdir/$Group2/$Two/$Group2.$Two.chr$i.mask.bed.gz\t";
			print O "--mask $maskdir/chr$i.mask.bed.gz\t";
                        print O "--mask $maskdir/chr$i.mask.bed.gz\t";
		        print O "$phasedir/$Group1/$One/$Group1.$One.chr$i.vcf.gz\t";
			print O "$phasedir/$Group2/$Two/$Group2.$Two.chr$i.vcf.gz\t";
		        print O ">$outdir/$groupname/G$num/chr$i.msmc.splitinput\n";
			#	print O3 "$python $boost_script\t";
			#print O3 "-n 20 -s 1000000 --chunks_per_chromosome 20 --nr_chromosomes 18 --seed 123 $outdir/$groupname/G$num-boost \\\n";
			if ($i==18){
				print O2 "$outdir/$groupname/G$num/chr$i.msmc.splitinput\n";
				print O3 "$outdir/$groupname/G$num/chr$i.msmc.splitinput\n";
			}else{
				print O2 "$outdir/$groupname/G$num/chr$i.msmc.splitinput \\\n";
				print O3 "$outdir/$groupname/G$num/chr$i.msmc.splitinput \\\n";
			}
		}

	     	for my $b (1..20){
			my $com_boost="G".$num."_".$b;
			`mkdir -p "$outdir2/$groupname/G$num-boost/$com_boost"` unless (-d "$outdir2/$groupname/G$num-boost/$com_boost");	
			open O4,">$outdir2/$groupname/G$num-boost/$com_boost/step2.$groupname.$com_boost.boost.Ne.sh";
			print O4 "$msmc2 --skipAmbiguous -I \'0\,1\,2\,3\' -o $outdir2/$groupname/G$num-boost/$com_boost/$groupname.$com_boost -i 20 -t 6 -p \'1\*2\+15\*1\+1\*2\' \\\n";
			for my $c (1..18){
				if ($c==18){
					print O4 "$outdir2/$groupname/G$num-boost/$com_boost/bootstrap_multihetsep.chr$c.txt\n";
				}else{
					print O4 "$outdir2/$groupname/G$num-boost/$com_boost/bootstrap_multihetsep.chr$c.txt \\\n";
				}
			}	
		}
					
	}
}
			
