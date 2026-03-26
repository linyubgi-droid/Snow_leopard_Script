inputvcf=${1}

groupname=${2}

ref=${3}

outdir=${4}

samplecount=${5}

generation=${6}

cat > slim_vcf_${groupname}_gen${generation}.slim << EOM



initialize() {
        initializeSLiMModelType("nonWF");
        initializeSLiMOptions(nucleotideBased=T);
        initializeMutationTypeNuc("m1", 0.5, "f", 0.0);
        initializeGenomicElementType("g1", m1, 1.0, mmJukesCantor(0.79e-8));
        defineConstant("K", 500);
for (id in 1:18, length in c(205058881,147847344,115010643,71602061,41135320,139809961,61312203,94926765,83399398,239984733,93858890,87077875,157992719,141163291,221327537,168335741,151709250,61800323), type in c("A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A","A"), symbols in c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18"), ref in c("${ref}/chr1.NA.fa","${ref}/chr2.NA.fa","${ref}/chr3.NA.fa","${ref}/chr4.NA.fa","${ref}/chr5.NA.fa","${ref}/chr6.NA.fa","${ref}/chr7.NA.fa","${ref}/chr8.NA.fa","${ref}/chr9.NA.fa","${ref}/chr10.NA.fa","${ref}/chr11.NA.fa","${ref}/chr12.NA.fa","${ref}/chr13.NA.fa","${ref}/chr14.NA.fa","${ref}/chr15.NA.fa","${ref}/chr16.NA.fa","${ref}/chr17.NA.fa","${ref}/chr18.NA.fa"))
        {
                initializeChromosome(id, length, type=type, symbol=symbols);
                initializeAncestralNucleotides(ref);

                initializeGenomicElement(g1, 0, length-1);
                initializeHotspotMap(c(4.94e-9), c(length-1));
                initializeRecombinationRate(1e-8);
   }

}

1 late() {
	sim.addSubpop("p1", ${samplecount});//Set input sample size to 30
	p1.individuals.readIndividualsFromVCF("${inputvcf}", m1);
}

reproduction() {
    subpop.addCrossed(individual, subpop.sampleIndividuals(1));
}
1:100 late() {
        mut = sim.mutationsOfType(m1);
}


2:100 early() {
    p1.fitnessScaling = min(1.0, K / p1.individualCount);
}

1 late() {
	sim.outputFull();
        inds = p1.sampleIndividuals(p1.individualCount);
        inds.outputIndividuals();
        inds.outputIndividualsToVCF("${outdir}/${groupname}.gen1"+".vcf");
}

2 late() {
	sim.outputFull();
        inds = p1.sampleIndividuals(p1.individualCount);
        inds.outputIndividuals();	
        inds.outputIndividualsToVCF("${outdir}/${groupname}.gen2"+".vcf");
}

3 late() {
	sim.outputFull();
        inds = p1.sampleIndividuals(p1.individualCount);
        inds.outputIndividuals();
        inds.outputIndividualsToVCF("${outdir}/${groupname}.gen3"+".vcf");
}

4 late() {
	sim.outputFull();
        inds = p1.sampleIndividuals(p1.individualCount);
        inds.outputIndividuals();
        inds.outputIndividualsToVCF("${outdir}/${groupname}.gen4"+".vcf");
}

5 late() {
	sim.outputFull();
        inds = p1.sampleIndividuals(p1.individualCount);
        inds.outputIndividuals();
        inds.outputIndividualsToVCF("${outdir}/${groupname}.gen5"+".vcf");
}

10 late() {
	sim.outputFull();
        inds = p1.sampleIndividuals(p1.individualCount);
        inds.outputIndividuals();
        inds.outputIndividualsToVCF("${outdir}/${groupname}.gen10"+".vcf");
}

15 late() {
	sim.outputFull();
        inds = p1.sampleIndividuals(p1.individualCount);
        inds.outputIndividuals();
        inds.outputIndividualsToVCF("${outdir}/${groupname}.gen15"+".vcf");
}

20 late() {
	sim.outputFull();
        inds = p1.sampleIndividuals(p1.individualCount);
        inds.outputIndividuals();
	inds.outputIndividualsToVCF("${outdir}/${groupname}.gen20"+".vcf");
}

25 late() {
	sim.outputFull();
        inds = p1.sampleIndividuals(p1.individualCount);
        inds.outputIndividuals();
        inds.outputIndividualsToVCF("${outdir}/${groupname}.gen25"+".vcf");
}

30 late() {
	sim.outputFull();
        inds = p1.sampleIndividuals(p1.individualCount);
        inds.outputIndividuals();
        inds.outputIndividualsToVCF("${outdir}/${groupname}.gen30"+".vcf");
}

40 late() {
	sim.outputFull();
        inds = p1.sampleIndividuals(p1.individualCount);
        inds.outputIndividuals();
        inds.outputIndividualsToVCF("${outdir}/${groupname}.gen40"+".vcf");
}

50 late() {
	sim.outputFull();
        inds = p1.sampleIndividuals(p1.individualCount);
        inds.outputIndividuals();
        inds.outputIndividualsToVCF("${outdir}/${groupname}.gen50"+".vcf");
}

60 late() {
	sim.outputFull();
        inds = p1.sampleIndividuals(p1.individualCount);
        inds.outputIndividuals();
        inds.outputIndividualsToVCF("${outdir}/${groupname}.gen60"+".vcf");
}

70 late() {
	sim.outputFull();
        inds = p1.sampleIndividuals(p1.individualCount);
        inds.outputIndividuals();
        inds.outputIndividualsToVCF("${outdir}/${groupname}.gen70"+".vcf");
}

80 late() {
	sim.outputFull();
        inds = p1.sampleIndividuals(p1.individualCount);
        inds.outputIndividuals();
        inds.outputIndividualsToVCF("${outdir}/${groupname}.gen80"+".vcf");
}

90 late() {
	sim.outputFull();
        inds = p1.sampleIndividuals(p1.individualCount);
        inds.outputIndividuals();
        inds.outputIndividualsToVCF("${outdir}/${groupname}.gen90"+".vcf");
}

100 late() {
	sim.outputFull();
        inds = p1.sampleIndividuals(p1.individualCount);
        inds.outputIndividuals();
        inds.outputIndividualsToVCF("${outdir}/${groupname}.gen100"+".vcf");
}




EOM


