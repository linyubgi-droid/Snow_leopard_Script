#!/bin/bash

# 检查参数数量，需要至少10个参数
if [ $# -lt 10 ]; then
    echo "ERROR：Need 10"
    echo "使用方法：$0 <prefix> <run_time> <run_time_count> <split_time> <split_time_count> <split_time_count_half> <now_time> <pop1_size> <pop2_size> <breed_count>"
    echo ""
    echo "参数说明："
    echo "  1. prefix: 前缀"
    echo "  2. run_time: 预跑时间"
    echo "  3. run_time_count: 预跑时间数目"
    echo "  4. split_time: 分化时间"
    echo "  5. split_time_count: 分化时间数目 6500"
    echo "  6. split_time_count_half: 分化时间计数的一半 3000"
    echo "  7. now_time: 当前时间"
    echo "  8. pop1_size: 种群1大小"
    echo "  9. pop2_size: 种群2大小"
    echo "  10. breed_count: 繁殖复制计数"
    exit 1
fi
prefix=${1}
run_time=${2}
run_time1=$((run_time + 1))
run_time2=$((run_time + 2))
run_time_count=${3}
split_time=${4}
split_time1=$((split_time + 1))
split_time2=$((split_time + 2))
split_time3=$((split_time + 3))
split_time_count=${5}
split_time_count_half=${6}
split_cut_run_time=$((split_time - run_time2))
now_time=${7}
now_time_1=$((now_time - 1))
now_time1=$((now_time + 1))
now_time2=$((now_time + 2))
pop1_size=${8}
pop2_size=${9}
breed_count=${10}
now_time_cut_split_time3=$((now_time - split_time3))
future_time=$((now_time1 + 100))
future_time1=$((future_time + 1))
cat > ${prefix}.slim << EOM

initialize() {
   initializeSLiMModelType("nonWF");
	//initializeSLiMOptions(keepPedigrees=T);
	defineConstant("geneLength", 1699);	
	initializeSex("A");
	defineConstant("simID", getSeed());
	//initializeTreeSeq();
	initializeMutationRate(0.49e-8);
	
	//set DFE parameters with additive mutations
	//enforce dominance parameters with mutationEffect function below
	initializeMutationType("m1", 0.5, "g",-0.2, 0.33);
	initializeMutationType("m2", 0.5, "f",0.0);
	initializeGenomicElementType("g1", c(m1,m2), c(2.31,1.0));
	m1.mutationStackPolicy = "f";
	m2.mutationStackPolicy = "f";
	
	//number of genes on each autosome from saola annotations
	gene_vec=c(1082,1286,1552,715,703,1193,1156,730,472,1373,869,673,944,1303,1868,1667,1184,1129);
	defineConstant("seqLength", sum(gene_vec)*geneLength);
	
	gene_num=sum(gene_vec);
	
	for (i in 1:gene_num){
		initializeGenomicElement(g1, ((i-1)*geneLength)+(i-1), (i*geneLength)+(i-2) );
	}
	
	rates=NULL;
	
	//assume no recombination within genes, a rate of 1e-3 between genes, and free recombination between chroms
	for (i in 1:(size(gene_vec)-1)){
		rates=c(rates, 0, rep(c(1e-3, 0), asInteger(gene_vec[i-1]-1)), 0.5);
	}
	rates=c(rates, 0, rep(c(1e-3, 0), asInteger(gene_vec[size(gene_vec)-1]-1)));
	
	ends=NULL;
	for (i in 1:gene_num){
		ends=c(ends, (i*geneLength)+(i-2), (i*geneLength)+(i-1));
	}
	ends=ends[0:(size(ends)-2)];
	initializeRecombinationRate(rates, ends);

}

1:${run_time} reproduction() {
	N = sim.getValue("N");
	for (i in seqLen(N))
	{
		firstParent = p1.sampleIndividuals(1, sex='F');
		secondParent = p1.sampleIndividuals(1, exclude=firstParent, sex='M');
		p1.addCrossed(firstParent, secondParent);
	}
	self.active = 0;
}

early()
{
        // parents die; offspring survive proportional to fitness
        inds = sim.subpopulations.individuals;
        inds[inds.age > 0].fitnessScaling = 0.0;
}

// Saola Demography
// can make this demography arbitrarily complex as described in the SLiM manual
1 first() {
	sim.addSubpop("p1", ${run_time_count});
	sim.setValue("N", ${run_time_count});
// convert fixed mutations to substitutions for computational efficiency
	m1.convertToSubstitution = F;
	m2.convertToSubstitution = T;
	m1.mutationStackPolicy = "f";
	m2.mutationStackPolicy = "f";

	writeFile("saola_stair_perez_out.csv", "POP,gen,popSize,heterozygosity,genLoadSlim,inbLoad,genLoad,realizedLoad,lethal_seg,sublet_seg,strDel_seg,modDel_seg,wkDel_seg,lethal_fix,sublet_fix,strDel_fix,modDel_fix,wkDel_fix,males,females", append=F);
	writeFile("load_decomposed_byh.txt", "POP,gen,popSize,heterozygosity,genLoadSlim,inbLoad,genLoad,realizedLoad,lethal_seg,sublet_seg,strDel_seg,modDel_seg,wkDel_seg,lethal_fix,sublet_fix,strDel_fix,modDel_fix,wkDel_fix,males,females", append=F);

	writeFile("saola_breeding_first.csv", "POP,gen,popSize,heterozygosity,genLoadSlim,inbLoad,genLoad,realizedLoad,lethal_seg,sublet_seg,strDel_seg,modDel_seg,wkDel_seg,lethal_fix,sublet_fix,strDel_fix,modDel_fix,wkDel_fix,males,females", append=F);
	writeFile("saola_breeding_early.csv", "POP,gen,popSize,heterozygosity,genLoadSlim,inbLoad,genLoad,realizedLoad,lethal_seg,sublet_seg,strDel_seg,modDel_seg,wkDel_seg,lethal_fix,sublet_fix,strDel_fix,modDel_fix,wkDel_fix,males,females", append=F);
	writeFile("saola_breeding_late.csv", "POP,gen,popSize,heterozygosity,genLoadSlim,inbLoad,genLoad,realizedLoad,lethal_seg,sublet_seg,strDel_seg,modDel_seg,wkDel_seg,lethal_fix,sublet_fix,strDel_fix,modDel_fix,wkDel_fix,males,females", append=F);
	writeFile("saola_hets_counter.csv", paste("generation", "population", "individual", "neutral", "weakly_deleterious", "moderately_deleterious", "strongly_deleterious", "sublethal", "lethal", sep=","), append=F);
}

// draw h from a uniform based on s
mutation(m1) {

	mut.setValue("dom", runif(n=1,min=0,max=exp(7.6*mut.selectionCoeff)));
	
	return T;
}

// use mutationEffect to enforce h-s relationship - note that this is very slow
mutationEffect(m1) {
	if (homozygous)
		return 1.0 + mut.selectionCoeff;
	else
		return 1.0 + mut.getValue("dom") * mut.selectionCoeff;
}



//  first bottleneck
${run_time1}:${split_time} reproduction() {
        N = sim.getValue("N");
        for (i in seqLen(N))
        {
                firstParent = p1.sampleIndividuals(1, sex='F');
                secondParent = p1.sampleIndividuals(1, exclude=firstParent, sex='M');
                p1.addCrossed(firstParent, secondParent);
        }
        self.active = 0;
}

${run_time2}:${split_time} first() {
        t = sim.cycle - ${run_time2};
        //exponential decay of log(p2 minimum last size / p2 initial size) divided by 167 generations until the present.
	exp_rate = log(${split_time_count}/${run_time_count})/${split_cut_run_time};
        p1_size = asInteger(round(exp(exp_rate * t) * ${run_time_count}));
        sim.setValue("N", p1_size);
        if (sim.cycle == ${split_time}){
                sim.setValue("p1_size",asInteger(p1_size));
        }
}
// Split
${split_time} late() {
        // p1 is North and p2 is Central
        sim.addSubpop("p2",0).setValue("C",asInteger(${split_time_count_half}));
        p2.takeMigrants(p1.sampleIndividuals(asInteger((${split_time_count_half})/2), sex="M"));
        p2.takeMigrants(p1.sampleIndividuals(asInteger((${split_time_count_half})/2), sex="F"));
        sim.setValue("N", ${split_time_count_half});
   sim.setValue("C", asInteger(${split_time_count_half}));
}

${split_time1} reproduction() {
        N = sim.getValue("N");
        for (i in seqLen(N))
        {
                firstParent = p1.sampleIndividuals(1, sex='F');
                secondParent = p1.sampleIndividuals(1, exclude=firstParent, sex='M');
                p1.addCrossed(firstParent, secondParent);
        }
        C = sim.getValue("C");
        for (j in seqLen(C))
        {
                firstParent2 = p2.sampleIndividuals(1, sex='F');
                secondParent2 = p2.sampleIndividuals(1, exclude=firstParent2, sex='M');
                p2.addCrossed(firstParent2, secondParent2);
        }
        self.active = 0;
}


${split_time2}:${now_time_1} reproduction() {
        N = sim.getValue("N");
        for (i in seqLen(N))
        {
                firstParent = p1.sampleIndividuals(1, sex='F');
                secondParent = p1.sampleIndividuals(1, exclude=firstParent, sex='M');
                p1.addCrossed(firstParent, secondParent);
        }
        C = sim.getValue("C");
        for (j in seqLen(C))
        {
                firstParent2 = p2.sampleIndividuals(1, sex='F');
                secondParent2 = p2.sampleIndividuals(1, exclude=firstParent2, sex='M');
                p2.addCrossed(firstParent2, secondParent2);
        }
        self.active = 0;
}

//Second bottleneck
${split_time3}:${now_time} first() {
	t = sim.cycle - ${split_time3};
	//exponential decay of log(p2 minimum last size / p2 initial size) divided by 167 generations until the present.
	exp_rate1 = log(${pop1_size}/${split_time_count_half})/${now_time_cut_split_time3};
	exp_rate2 = log(${pop2_size}/${split_time_count_half})/${now_time_cut_split_time3};
	p1_size = asInteger(round(exp(exp_rate1 * t) * ${split_time_count_half}));
	p2_size = asInteger(round(exp(exp_rate2 * t) * ${split_time_count_half}));
	sim.setValue("N", p1_size);
   	sim.setValue("C", p2_size);
	if (sim.cycle == ${now_time}){
		sim.setValue("p1_size",asInteger(p1_size));
		sim.setValue("p2_size",asInteger(p2_size));
	}
}
// for last generation create MANY (1000) inds per eac of the 2 population population
// else we had not enough individauls to simulate all 9 breeding programs simultaneously
// it is not realistic but it does not matter, each of the 9 breeding programs are alternative universes :)
${now_time} reproduction() {
	N = ${breed_count};
	for (i in seqLen(N))
	{
		firstParent = p1.sampleIndividuals(1, sex='F');
		secondParent = p1.sampleIndividuals(1, exclude=firstParent, sex='M');
		p1.addCrossed(firstParent, secondParent);
	}
	C = ${breed_count};
	for (j in seqLen(C))
	{
		firstParent2 = p2.sampleIndividuals(1, sex='F');
		secondParent2 = p2.sampleIndividuals(1, exclude=firstParent2, sex='M');
		p2.addCrossed(firstParent2, secondParent2);
	}
	self.active = 0;
}


//Admixture breeding program
// create 9 populations (3 founding proportions x 3 founders number) 
${now_time} late() {


	// ADMIXED AT EQUAL PROPORTIONS
 	sim.addSubpop("p4",0);
        p4.name="admix4";
        p4.takeMigrants(p1.sampleIndividuals(1, sex="M"));
        p4.takeMigrants(p1.sampleIndividuals(1, sex="F"));
        p4.takeMigrants(p2.sampleIndividuals(1, sex="M"));
        p4.takeMigrants(p2.sampleIndividuals(1, sex="F"));
        // sim.setValue("K",p3.individualCount);
        catn("Migrants in p4: " + p4.individualCount);
       
	sim.addSubpop("p5",0);
	p5.name="admix8";
	p5.takeMigrants(p1.sampleIndividuals(2, sex="M"));
	p5.takeMigrants(p1.sampleIndividuals(2, sex="F"));
	p5.takeMigrants(p2.sampleIndividuals(2, sex="M"));
	p5.takeMigrants(p2.sampleIndividuals(2, sex="F"));
	// sim.setValue("K",p3.individualCount);
	catn("Migrants in p5: " + p5.individualCount);

	sim.addSubpop("p6",0);
	p6.name = "admix12";
	p6.takeMigrants(p1.sampleIndividuals(3, sex="M"));
	p6.takeMigrants(p1.sampleIndividuals(3, sex="F"));
	p6.takeMigrants(p2.sampleIndividuals(3, sex="M"));
	p6.takeMigrants(p2.sampleIndividuals(3, sex="F"));
	// value to keep count on which subpopulations are alive
	catn("Migrants in p6: " + p6.individualCount);
 
	sim.addSubpop("p7",0);
        p7.name="admix16";
        p7.takeMigrants(p1.sampleIndividuals(4, sex="M"));
        p7.takeMigrants(p1.sampleIndividuals(4, sex="F"));
        p7.takeMigrants(p2.sampleIndividuals(4, sex="M"));
        p7.takeMigrants(p2.sampleIndividuals(4, sex="F"));
        // sim.setValue("K",p3.individualCount);
        catn("Migrants in p7: " + p7.individualCount);
	
	sim.addSubpop("p8",0);
        p8.name="admix20";
        p8.takeMigrants(p1.sampleIndividuals(5, sex="M"));
        p8.takeMigrants(p1.sampleIndividuals(5, sex="F"));
        p8.takeMigrants(p2.sampleIndividuals(5, sex="M"));
        p8.takeMigrants(p2.sampleIndividuals(5, sex="F"));
        // sim.setValue("K",p3.individualCount);
        catn("Migrants in p8: " + p8.individualCount);

	sim.addSubpop("p9",0);
	p9.name = "admix24";
	p9.takeMigrants(p1.sampleIndividuals(6, sex="M"));
	p9.takeMigrants(p1.sampleIndividuals(6, sex="F"));
	p9.takeMigrants(p2.sampleIndividuals(6, sex="M"));
	p9.takeMigrants(p2.sampleIndividuals(6, sex="F"));
	// sim.setValue("K",p3.individualCount);
	catn("Migrants in p9: " + p9.individualCount);
	
	sim.addSubpop("p10",0);
        p10.name="admix28";
        p10.takeMigrants(p1.sampleIndividuals(7, sex="M"));
        p10.takeMigrants(p1.sampleIndividuals(7, sex="F"));
        p10.takeMigrants(p2.sampleIndividuals(7, sex="M"));
        p10.takeMigrants(p2.sampleIndividuals(7, sex="F"));
        // sim.setValue("K",p3.individualCount);
        catn("Migrants in p10: " + p10.individualCount);

	sim.addSubpop("p11",0);
	p11.name="admix32";
        p11.takeMigrants(p1.sampleIndividuals(8, sex="M"));
        p11.takeMigrants(p1.sampleIndividuals(8, sex="F"));
        p11.takeMigrants(p2.sampleIndividuals(8, sex="M"));
        p11.takeMigrants(p2.sampleIndividuals(8, sex="F"));
        // sim.setValue("K",p3.individualCount);
        catn("Migrants in p11: " + p11.individualCount);

	// AM
	sim.addSubpop("p12",0);
	p12.name = "AM4";
	p12.takeMigrants(p2.sampleIndividuals(2, sex="M"));
	p12.takeMigrants(p2.sampleIndividuals(2, sex="F"));
	// sim.setValue("K",p3.individualCount);
	catn("Migrants in p12: " + p12.individualCount);

	sim.addSubpop("p13",0);
	p13.name = "AM8";
	p13.takeMigrants(p2.sampleIndividuals(4, sex="M"));
	p13.takeMigrants(p2.sampleIndividuals(4, sex="F"));
	// sim.setValue("K",p3.individualCount);
	catn("Migrants in p13: " + p13.individualCount);

	sim.addSubpop("p14",0);
	p14.name = "AM12";
	p14.takeMigrants(p2.sampleIndividuals(6, sex="M"));
	p14.takeMigrants(p2.sampleIndividuals(6, sex="F"));
	// sim.setValue("K",p3.individualCount);
	catn("Migrants in p14: " + p14.individualCount);

	sim.addSubpop("p15",0);
        p15.name = "AM16";
        p15.takeMigrants(p2.sampleIndividuals(8, sex="M"));
        p15.takeMigrants(p2.sampleIndividuals(8, sex="F"));
        // sim.setValue("K",p3.individualCount);
        catn("Migrants in p15: " + p15.individualCount);

	sim.addSubpop("p16",0);
        p16.name = "AM20";
        p16.takeMigrants(p2.sampleIndividuals(10, sex="M"));
        p16.takeMigrants(p2.sampleIndividuals(10, sex="F"));
        // sim.setValue("K",p3.individualCount);
        catn("Migrants in p16: " + p16.individualCount);

	sim.addSubpop("p17",0);
        p17.name = "AM24";
        p17.takeMigrants(p2.sampleIndividuals(12, sex="M"));
        p17.takeMigrants(p2.sampleIndividuals(12, sex="F"));
        // sim.setValue("K",p3.individualCount);
        catn("Migrants in p17: " + p17.individualCount);

	sim.addSubpop("p18",0);
        p18.name = "AM28";
        p18.takeMigrants(p2.sampleIndividuals(14, sex="M"));
        p18.takeMigrants(p2.sampleIndividuals(14, sex="F"));
        // sim.setValue("K",p3.individualCount);
        catn("Migrants in p18: " + p18.individualCount);

	sim.addSubpop("p19",0);
        p19.name = "AM32";
        p19.takeMigrants(p2.sampleIndividuals(16, sex="M"));
        p19.takeMigrants(p2.sampleIndividuals(16, sex="F"));
        // sim.setValue("K",p3.individualCount);
        catn("Migrants in p19: " + p19.individualCount);

	// NORTHERN
	sim.addSubpop("p20",0);
	p20.name = "AH4";
	p20.takeMigrants(p1.sampleIndividuals(2, sex="M"));
	p20.takeMigrants(p1.sampleIndividuals(2, sex="F"));
	// sim.setValue("K",p10.individualCount);
	catn("Migrants in p20: " + p20.individualCount);

	sim.addSubpop("p21",0);
        p21.name = "AH8";
        p21.takeMigrants(p1.sampleIndividuals(4, sex="M"));
        p21.takeMigrants(p1.sampleIndividuals(4, sex="F"));
        // sim.setValue("K",p10.individualCount);
        catn("Migrants in p21: " + p21.individualCount);


	sim.addSubpop("p22",0);
	p22.name = "AH12";
	p22.takeMigrants(p1.sampleIndividuals(6, sex="M"));
	p22.takeMigrants(p1.sampleIndividuals(6, sex="F"));
	// sim.setValue("K",p3.individualCount);
	catn("Migrants in p22: " + p22.individualCount);

	sim.addSubpop("p23",0);
        p23.name = "AH16";
        p23.takeMigrants(p1.sampleIndividuals(8, sex="M"));
        p23.takeMigrants(p1.sampleIndividuals(8, sex="F"));
        // sim.setValue("K",p10.individualCount);
        catn("Migrants in p23: " + p23.individualCount);

	sim.addSubpop("p24",0);
        p24.name = "AH20";
        p24.takeMigrants(p1.sampleIndividuals(10, sex="M"));
        p24.takeMigrants(p1.sampleIndividuals(10, sex="F"));
        // sim.setValue("K",p10.individualCount);
        catn("Migrants in p24: " + p24.individualCount);

	sim.addSubpop("p25",0);
        p25.name = "AH24";
        p25.takeMigrants(p1.sampleIndividuals(12, sex="M"));
        p25.takeMigrants(p1.sampleIndividuals(12, sex="F"));
        // sim.setValue("K",p10.individualCount);
        catn("Migrants in p25: " + p25.individualCount);


	sim.addSubpop("p26",0);
        p26.name = "AH28";
        p26.takeMigrants(p1.sampleIndividuals(14, sex="M"));
        p26.takeMigrants(p1.sampleIndividuals(14, sex="F"));
        // sim.setValue("K",p10.individualCount);
        catn("Migrants in p26: " + p26.individualCount);
	
	sim.addSubpop("p27",0);
        p27.name = "AH32";
        p27.takeMigrants(p1.sampleIndividuals(16, sex="M"));
        p27.takeMigrants(p1.sampleIndividuals(16, sex="F"));
        // sim.setValue("K",p10.individualCount);
        catn("Migrants in p27: " + p27.individualCount);


	p4.setValue("survives", "Yes");
	p5.setValue("survives", "Yes");
	p6.setValue("survives", "Yes");
	p7.setValue("survives", "Yes");
	p8.setValue("survives", "Yes");
	p9.setValue("survives", "Yes");
	p10.setValue("survives", "Yes");
	p11.setValue("survives", "Yes");
	p12.setValue("survives", "Yes");
	p13.setValue("survives", "Yes");
        p14.setValue("survives", "Yes");
        p15.setValue("survives", "Yes");
        p16.setValue("survives", "Yes");
        p17.setValue("survives", "Yes");
	p18.setValue("survives", "Yes");
        p19.setValue("survives", "Yes");
        p20.setValue("survives", "Yes");
        p21.setValue("survives", "Yes");
        p22.setValue("survives", "Yes");
	p23.setValue("survives", "Yes");
        p24.setValue("survives", "Yes");
        p25.setValue("survives", "Yes");
        p26.setValue("survives", "Yes");
        p27.setValue("survives", "Yes");


	// value to keep track if subpopulations are active
	p1.setValue("active", F);
	p2.setValue("active", F);
	p4.setValue("active", T);
	p5.setValue("active", T);
	p6.setValue("active", T);
	p7.setValue("active", T);
	p8.setValue("active", T);
	p9.setValue("active", T);
	p10.setValue("active", T);
	p11.setValue("active", T);
	p12.setValue("active", T);
	p13.setValue("active", T);
        p14.setValue("active", T);
        p15.setValue("active", T);
        p16.setValue("active", T);
        p17.setValue("active", T);
        p18.setValue("active", T);
        p19.setValue("active", T);
        p20.setValue("active", T);
        p21.setValue("active", T);
	p22.setValue("active", T);
        p23.setValue("active", T);
        p24.setValue("active", T);
        p25.setValue("active", T);
        p26.setValue("active", T);
        p27.setValue("active", T);


}

${now_time1}:${future_time} reproduction() {
	// reproduce all subpopulations with popsize > 0, producing 2 offspring per individual
	for (spop in sim.subpopulations){

		if (!spop.getValue("active"))
			next;

		popsize = spop.individualCount * 2;

		if(popsize == 0)
			next;

		for (i in seqLen(popsize)){
			firstParent = spop.sampleIndividuals(1, sex='F');
			secondParent = spop.sampleIndividuals(1, sex='M');
			spop.addCrossed(firstParent, secondParent);
		}	
	}
	self.active=0;

}

${now_time1} first() {
	sim.killIndividuals(p1.individuals);
	sim.killIndividuals(p2.individuals);
}


1:${now_time} early(){
	//muts=sim.mutations;
	muts=sim.mutationsOfType(m1);
//	currentids = muts.id;
	
	// output statistics every 1000 generations and other conditions
	if ((sim.cycle % 1000 == 0 |  sim.cycle == 1 | sim.cycle == ${run_time1} | sim.cycle == ${run_time2} | (sim.cycle >= ${run_time2})) & (size(p1.individuals) > 0)) {
		// North: p1 estimations
		

		// print heterozygosity in different categories for sample of 5 individuals
		for(ind in p1.sampleIndividuals(5)) {
	
			ind_muts = ind.genomes.mutations;
			counts = ind.genomes.mutationCountsInGenomes(ind_muts);
			hets = counts[counts==1];
			s = ind_muts.selectionCoeff;
			s_hets = s[counts==1];

			neutral_het = hets[s_hets== 0];
			lethal_het = hets[s_hets == -1];
			sublet_het = hets[s_hets > -1 & s_hets <= -0.1 ];
			strDel_het = hets[s_hets <= -0.01 & s_hets > -0.1 ];
			modDel_het = hets[s_hets <= -0.001 & s_hets > -0.01 ];
			wkDel_het = hets[s_hets < 0 & s_hets > -0.001];
			writeFile("saola_hets_counter.csv", paste(sim.cycle, "AH", ind.index, sum(neutral_het), sum(wkDel_het), sum(modDel_het), sum(strDel_het), sum(sublet_het), sum(lethal_het), sep=","), append=T);
	}

		// Calculate gender ratio, if any are extinct, stop sim.
		females = size(p1.individuals.sex[p1.individuals.sex=='F']);
		males = size(p1.individuals.sex[p1.individuals.sex=='M']);
		
		//get segregating muts
		freq = p1.genomes.mutationFrequenciesInGenomes(muts);
		seg = muts[freq != 1.0 & freq > 0];

		lethal_seg = seg[seg.selectionCoeff==-1];
		sublet_seg = seg[seg.selectionCoeff > -1 & seg.selectionCoeff <= -0.1 ];
		strDel_seg = seg[seg.selectionCoeff <= -0.01 & seg.selectionCoeff > -0.1 ];
		modDel_seg = seg[seg.selectionCoeff <= -0.001 & seg.selectionCoeff > -0.01 ];
		wkDel_seg = seg[seg.selectionCoeff > -0.001];
		
		//get fixed muts
		fixed = muts[freq == 1.0];
		
		lethal_fix = fixed[fixed.selectionCoeff==-1];
		sublet_fix = fixed[fixed.selectionCoeff > -1 & fixed.selectionCoeff <= -0.1 ];
		strDel_fix = fixed[fixed.selectionCoeff <= -0.01 & fixed.selectionCoeff > -0.1 ];
		modDel_fix = fixed[fixed.selectionCoeff <= -0.001 & fixed.selectionCoeff > -0.01 ];
		wkDel_fix = fixed[fixed.selectionCoeff > -0.001];
		
		// get load as 1 - fitness
		// SLiM calculates fitness multiplicatively each generation and caches it
		genetic_load = 1-mean(p1.cachedFitness(NULL));
		
		// calculate inbreeding load
		q = p1.genomes.mutationFrequenciesInGenomes(muts);
		s = -muts.selectionCoeff;

		// replace mutations with s>1.0 with 1.0 (can happen when drawing from gamma distribution)
		s[s>1.0]=1.0;
	

		// get h for each mutation
		// note that this will not work if changing h using fitness callbacks
		//h=muts.mutationType.dominanceCoeff;
		//
		if(muts.size() == 0){
		h=float(0);
		} else {
		h=sapply(muts, "applyValue.getValue('dom');");
		}

		// calculate and print to file loads decomposed by H bin
		bin1 = h >= 0 & h <= 0.01;
		bin2 = h > 0.01 & h <= 0.05;
		bin3 = h > 0.05 & h <= 0.1;
		bin4 = h > 0.1 & h <= 0.2;
		bin5 = h > 0.2 & h <= 0.3;
		bin6 = h > 0.3 & h <= 0.4;
		bin7 = h > 0.4 & h <= 0.5;
		bin8 = h > 0.5 & h <= 0.6;
		bin9 = h > 0.6 & h <= 0.7;
		bin10 = h > 0.7 & h <= 0.8;
		bin11 = h > 0.8 & h <= 0.9;
		bin12 = h > 0.9 & h <= 1;

		genLoad1 = sum(q[bin1] * s[bin1]);
		genLoad2 = sum(q[bin2] * s[bin2]);
		genLoad3 = sum(q[bin3] * s[bin3]);
		genLoad4 = sum(q[bin4] * s[bin4]);
		genLoad5 = sum(q[bin5] * s[bin5]);
		genLoad6 = sum(q[bin6] * s[bin6]);
		genLoad7 = sum(q[bin7] * s[bin7]);
		genLoad8 = sum(q[bin8] * s[bin8]);
		genLoad9 = sum(q[bin9] * s[bin9]);
		genLoad10 = sum(q[bin10] * s[bin10]);
		genLoad11 = sum(q[bin11] * s[bin11]);
		genLoad12 = sum(q[bin12] * s[bin12]);


		relative_load1 = sum(q[bin1]^2*s[bin1])+2*sum(q[bin1]*(1-q[bin1])*s[bin1]*h[bin1]);
		relative_load2 = sum(q[bin2]^2*s[bin2])+2*sum(q[bin2]*(1-q[bin2])*s[bin2]*h[bin2]);
		relative_load3 = sum(q[bin3]^2*s[bin3])+2*sum(q[bin3]*(1-q[bin3])*s[bin3]*h[bin3]);
		relative_load4 = sum(q[bin4]^2*s[bin4])+2*sum(q[bin4]*(1-q[bin4])*s[bin4]*h[bin4]);
		relative_load5 = sum(q[bin5]^2*s[bin5])+2*sum(q[bin5]*(1-q[bin5])*s[bin5]*h[bin5]);
		relative_load6 = sum(q[bin6]^2*s[bin6])+2*sum(q[bin6]*(1-q[bin6])*s[bin6]*h[bin6]);
		relative_load7 = sum(q[bin7]^2*s[bin7])+2*sum(q[bin7]*(1-q[bin7])*s[bin7]*h[bin7]);
		relative_load8 = sum(q[bin8]^2*s[bin8])+2*sum(q[bin8]*(1-q[bin8])*s[bin8]*h[bin8]);
		relative_load9 = sum(q[bin9]^2*s[bin9])+2*sum(q[bin9]*(1-q[bin9])*s[bin9]*h[bin9]);
		relative_load10 = sum(q[bin10]^2*s[bin10])+2*sum(q[bin10]*(1-q[bin10])*s[bin10]*h[bin10]);
		relative_load11 = sum(q[bin11]^2*s[bin11])+2*sum(q[bin11]*(1-q[bin11])*s[bin11]*h[bin11]);
		relative_load12 = sum(q[bin12]^2*s[bin12])+2*sum(q[bin12]*(1-q[bin12])*s[bin12]*h[bin12]);

		writeFile("load_decomposed_byh.txt", "AH," + sim.cycle + "," + p1.individuals.size() + "," + genLoad1 + "," + genLoad2 + "," + genLoad3 + "," + genLoad4 + "," + genLoad5 + "," + genLoad6 + "," + genLoad7 + "," + genLoad8 + "," + genLoad9 + "," + genLoad10 + "," + genLoad11 + "," + genLoad12 + "," + relative_load1 + "," + relative_load2 + "," + relative_load3 + "," + relative_load4 + "," + relative_load5 + "," + relative_load6 + "," + relative_load7 + "," + relative_load8 + "," + relative_load9 + "," + relative_load10 + "," + relative_load11 + "," + relative_load12, append=T);

		//print(h);
		// calculate number of diploid lethal equivalents (B or inbreeding load)
		// equation from Morton et al 1956
		// note that this equation assumes random mating
		inbreeding_load = sum(q*s)-sum(q^2*s)-2*sum(q*(1-q)*s*h);
		genetic_pop = sum(q*s);
		relative_load = sum(q^2*s)+2*sum(q*(1-q)*s*h);
		
		heterozygosity = calcHeterozygosity(p1.genomes);
		
		writeFile("saola_stair_perez_out.csv", "AH" + "," + sim.cycle + "," + p1.individuals.size() + "," + heterozygosity + "," + genetic_load + "," + inbreeding_load + "," + genetic_pop + "," + relative_load + "," + lethal_seg.length()  + "," + sublet_seg.length() + "," + strDel_seg.length() + "," + modDel_seg.length() + "," + wkDel_seg.length() + "," + lethal_fix.length() + "," + sublet_fix.length() + "," + strDel_fix.length() + "," + modDel_fix.length() + "," + wkDel_fix.length() + "," + males + "," + females, append=T);
		
		//Central: p2 estimations 
		if (sim.cycle >= ${split_time1}) {
			if(size(p2.individuals) > 0){
			// Calculate gender ratio, if any are extinct, stop sim.
			females2 = size(p2.individuals.sex[p2.individuals.sex=='F']);
			males2 = size(p2.individuals.sex[p2.individuals.sex=='M']);
			

		// print heterozygosity in different categories for sample of 5 individuals
		for(ind in p2.sampleIndividuals(5)) {
	
			ind_muts = ind.genomes.mutations;
			counts = ind.genomes.mutationCountsInGenomes(ind_muts);
			hets = counts[counts==1];
			s = ind_muts.selectionCoeff;
			s_hets = s[counts==1];

			neutral_het = hets[s_hets== 0];
			lethal_het = hets[s_hets == -1];
			sublet_het = hets[s_hets > -1 & s_hets <= -0.1 ];
			strDel_het = hets[s_hets <= -0.01 & s_hets > -0.1 ];
			modDel_het = hets[s_hets <= -0.001 & s_hets > -0.01 ];
			wkDel_het = hets[s_hets < 0 & s_hets > -0.001];
			writeFile("saola_hets_counter.csv", paste(sim.cycle, "AM", ind.index, sum(neutral_het), sum(wkDel_het), sum(modDel_het), sum(strDel_het), sum(sublet_het), sum(lethal_het), sep=","), append=T);
	}

			//p2 popSize
			p2_size = sim.subpopulations.individuals.size()-p1.individuals.size();
			
			//get segregating muts
			freq2 = p2.genomes.mutationFrequenciesInGenomes(muts);
			seg2 = muts[freq2 != 1.0 & freq2 > 0];
			

			lethal_seg2 = seg2[seg2.selectionCoeff==-1];
			sublet_seg2 = seg2[seg2.selectionCoeff > -1 & seg2.selectionCoeff <= -0.1 ];
			strDel_seg2 = seg2[seg2.selectionCoeff <= -0.01 & seg2.selectionCoeff > -0.1 ];
			modDel_seg2 = seg2[seg2.selectionCoeff <= -0.001 & seg2.selectionCoeff > -0.01 ];
			wkDel_seg2 = seg2[seg2.selectionCoeff > -0.001];
			
			//get fixed muts
			fixed2 = muts[freq2 == 1.0];
			
			lethal_fix2 = fixed2[fixed2.selectionCoeff==-1];
			sublet_fix2 = fixed2[fixed2.selectionCoeff > -1 & fixed2.selectionCoeff <= -0.1 ];
			strDel_fix2 = fixed2[fixed2.selectionCoeff <= -0.01 & fixed2.selectionCoeff > -0.1 ];
			modDel_fix2 = fixed2[fixed2.selectionCoeff <= -0.001 & fixed2.selectionCoeff > -0.01 ];
			wkDel_fix2 = fixed2[fixed2.selectionCoeff > -0.001];
			
			// get load as 1 - fitness
			// SLiM calculates fitness multiplicatively each generation and caches it
			genetic_load2 = 1-mean(p2.cachedFitness(NULL));
			
			// calculate inbreeding load
			q2 = p2.genomes.mutationFrequenciesInGenomes(muts);
			s2 = -muts.selectionCoeff;
			
			// replace mutations with s>1.0 with 1.0 (can happen when drawing from gamma distribution)
			s2[s2>1.0]=1.0;
			
			// get h for each mutation
			// note that this will not work if changing h using fitness callbacks
			//h2 = muts.mutationType.dominanceCoeff;
			//h2 = sapply(muts, "applyValue.getValue('dom');");

			if(muts.size() == 0){
			h2=float(0);
			} else {
			  h2=sapply(muts, "applyValue.getValue('dom');");
				 }

			
		// calculate and print to file loads decomposed by H bin
		bin1 = h2 >= 0 & h2 <= 0.01;
		bin2 = h2 > 0.01 & h2 <= 0.05;
		bin3 = h2 > 0.05 & h2 <= 0.1;
		bin4 = h2 > 0.1 & h2 <= 0.2;
		bin5 = h2 > 0.2 & h2 <= 0.3;
		bin6 = h2 > 0.3 & h2 <= 0.4;
		bin7 = h2 > 0.4 & h2 <= 0.5;
		bin8 = h2 > 0.5 & h2 <= 0.6;
		bin9 = h2 > 0.6 & h2 <= 0.7;
		bin10 = h2 > 0.7 & h2 <= 0.8;
		bin11 = h2 > 0.8 & h2 <= 0.9;
		bin12 = h2 > 0.9 & h2 <= 1;

		genLoad1 = sum(q2[bin1] * s2[bin1]);
		genLoad2 = sum(q2[bin2] * s2[bin2]);
		genLoad3 = sum(q2[bin3] * s2[bin3]);
		genLoad4 = sum(q2[bin4] * s2[bin4]);
		genLoad5 = sum(q2[bin5] * s2[bin5]);
		genLoad6 = sum(q2[bin6] * s2[bin6]);
		genLoad7 = sum(q2[bin7] * s2[bin7]);
		genLoad8 = sum(q2[bin8] * s2[bin8]);
		genLoad9 = sum(q2[bin9] * s2[bin9]);
		genLoad10 = sum(q2[bin10] * s2[bin10]);
		genLoad11 = sum(q2[bin11] * s2[bin11]);
		genLoad12 = sum(q2[bin12] * s2[bin12]);

		relative_load1 = sum(q2[bin1]^2*s2[bin1])+2*sum(q2[bin1]*(1-q2[bin1])*s2[bin1]*h2[bin1]);
		relative_load2 = sum(q2[bin2]^2*s2[bin2])+2*sum(q2[bin2]*(1-q2[bin2])*s2[bin2]*h2[bin2]);
		relative_load3 = sum(q2[bin3]^2*s2[bin3])+2*sum(q2[bin3]*(1-q2[bin3])*s2[bin3]*h2[bin3]);
		relative_load4 = sum(q2[bin4]^2*s2[bin4])+2*sum(q2[bin4]*(1-q2[bin4])*s2[bin4]*h2[bin4]);
		relative_load5 = sum(q2[bin5]^2*s2[bin5])+2*sum(q2[bin5]*(1-q2[bin5])*s2[bin5]*h2[bin5]);
		relative_load6 = sum(q2[bin6]^2*s2[bin6])+2*sum(q2[bin6]*(1-q2[bin6])*s2[bin6]*h2[bin6]);
		relative_load7 = sum(q2[bin7]^2*s2[bin7])+2*sum(q2[bin7]*(1-q2[bin7])*s2[bin7]*h2[bin7]);
		relative_load8 = sum(q2[bin8]^2*s2[bin8])+2*sum(q2[bin8]*(1-q2[bin8])*s2[bin8]*h2[bin8]);
		relative_load9 = sum(q2[bin9]^2*s2[bin9])+2*sum(q2[bin9]*(1-q2[bin9])*s2[bin9]*h2[bin9]);
		relative_load10 = sum(q2[bin10]^2*s2[bin10])+2*sum(q2[bin10]*(1-q2[bin10])*s2[bin10]*h2[bin10]);
		relative_load11 = sum(q2[bin11]^2*s2[bin11])+2*sum(q2[bin11]*(1-q2[bin11])*s2[bin11]*h2[bin11]);
		relative_load12 = sum(q2[bin12]^2*s2[bin12])+2*sum(q2[bin12]*(1-q2[bin12])*s2[bin12]*h2[bin12]);

		writeFile("load_decomposed_byh.txt", "AM," + sim.cycle + "," + p2.individuals.size() + "," + genLoad1 + "," + genLoad2 + "," + genLoad3 + "," + genLoad4 + "," + genLoad5 + "," + genLoad6 + "," + genLoad7 + "," + genLoad8 + "," + genLoad9 + "," + genLoad10 + "," + genLoad11 + "," + genLoad12 + "," + relative_load1 + "," + relative_load2 + "," + relative_load3 + "," + relative_load4 + "," + relative_load5 + "," + relative_load6 + "," + relative_load7 + "," + relative_load8 + "," + relative_load9 + "," + relative_load10 + "," + relative_load11 + "," + relative_load12, append=T);


		// calculate number of diploid lethal equivalents (B or inbreeding load)
		// equation from Morton et al 1956
		// note that this equation assumes random mating
		inbreeding_load2 = sum(q2*s2)-sum(q2^2*s2)-2*sum(q2*(1-q2)*s2*h2);
		genetic_pop2 = sum(q2*s2);
		relative_load2 = sum(q2^2*s2)+2*sum(q2*(1-q2)*s2*h2);

		heterozygosity = calcHeterozygosity(p2.genomes);

		writeFile("saola_stair_perez_out.csv", "AM" + "," + sim.cycle + "," + p2_size + "," + heterozygosity + "," + genetic_load2 + "," + inbreeding_load2 + "," + genetic_pop2 + "," + relative_load2 + "," + lethal_seg2.length()  + "," + sublet_seg2.length() + "," + strDel_seg2.length() + "," + modDel_seg2.length() + "," + wkDel_seg2.length() + "," + lethal_fix2.length() + "," + sublet_fix2.length() + "," + strDel_fix2.length() + "," + modDel_fix2.length() + "," + wkDel_fix2.length() + "," + males2 + "," + females2, append=T);
			}
		}
	}
}

${now_time1}:${future_time} first(){

	muts = sim.mutationsOfType(m1);
	for (subpop in sim.subpopulations){

		if (!subpop.getValue("active"))
			next;

		if (subpop.individuals.size() == 0)
			next;


		females = size(subpop.individuals.sex[subpop.individuals.sex=='F']);
		males = size(subpop.individuals.sex[subpop.individuals.sex=='M']);

		//get segregating muts
		freq = subpop.genomes.mutationFrequenciesInGenomes(muts);
		seg = muts[freq != 1.0 & freq > 0];

		lethal_seg = seg[seg.selectionCoeff==-1];
		sublet_seg = seg[seg.selectionCoeff > -1 & seg.selectionCoeff <= -0.1 ];
		strDel_seg = seg[seg.selectionCoeff <= -0.01 & seg.selectionCoeff > -0.1 ];
		modDel_seg = seg[seg.selectionCoeff <= -0.001 & seg.selectionCoeff > -0.01 ];
		wkDel_seg = seg[seg.selectionCoeff > -0.001];
		
		//get fixed muts
		fixed = muts[freq == 1.0];
		
		lethal_fix = fixed[fixed.selectionCoeff==-1];
		sublet_fix = fixed[fixed.selectionCoeff > -1 & fixed.selectionCoeff <= -0.1 ];
		strDel_fix = fixed[fixed.selectionCoeff <= -0.01 & fixed.selectionCoeff > -0.1 ];
		modDel_fix = fixed[fixed.selectionCoeff <= -0.001 & fixed.selectionCoeff > -0.01 ];
		wkDel_fix = fixed[fixed.selectionCoeff > -0.001];
		
		// get load as 1 - fitness
		// SLiM calculates fitness multiplicatively each generation and caches it
		genetic_load = 1-mean(subpop.cachedFitness(NULL));
		
		// calculate inbreeding load
		q = subpop.genomes.mutationFrequenciesInGenomes(muts);
		s = -muts.selectionCoeff;

		// replace mutations with s>1.0 with 1.0 (can happen when drawing from gamma distribution)
		s[s>1.0]=1.0;
	

		// get h for each mutation
		// note that this will not work if changing h using fitness callbacks
		//h=muts.mutationType.dominanceCoeff;
		//
		if(muts.size() == 0){
		h=float(0);
		} else {
		h=sapply(muts, "applyValue.getValue('dom');");
		}
		inbreeding_load = sum(q*s)-sum(q^2*s)-2*sum(q*(1-q)*s*h);
		genetic_pop = sum(q*s);
		relative_load = sum(q^2*s)+2*sum(q*(1-q)*s*h);
		
		heterozygosity = calcHeterozygosity(subpop.genomes);

		
		writeFile("saola_breeding_first.csv", subpop.name + "," + sim.cycle + "," + subpop.individuals.size() + "," + heterozygosity + "," + genetic_load + "," + inbreeding_load + "," + genetic_pop + "," + relative_load + "," + lethal_seg.length()  + "," + sublet_seg.length() + "," + strDel_seg.length() + "," + modDel_seg.length() + "," + wkDel_seg.length() + "," + lethal_fix.length() + "," + sublet_fix.length() + "," + strDel_fix.length() + "," + modDel_fix.length() + "," + wkDel_fix.length() + "," + males + "," + females, append=T);
		
	}

} 

${now_time1}:${future_time} early(){

	muts = sim.mutationsOfType(m1);
	for (subpop in sim.subpopulations){

		if (!subpop.getValue("active"))
			next;

		if (subpop.individuals.size() == 0)
			next;
		

		females = size(subpop.individuals.sex[subpop.individuals.sex=='F']);
		males = size(subpop.individuals.sex[subpop.individuals.sex=='M']);

		//get segregating muts
		freq = subpop.genomes.mutationFrequenciesInGenomes(muts);
		seg = muts[freq != 1.0 & freq > 0];

		lethal_seg = seg[seg.selectionCoeff==-1];
		sublet_seg = seg[seg.selectionCoeff > -1 & seg.selectionCoeff <= -0.1 ];
		strDel_seg = seg[seg.selectionCoeff <= -0.01 & seg.selectionCoeff > -0.1 ];
		modDel_seg = seg[seg.selectionCoeff <= -0.001 & seg.selectionCoeff > -0.01 ];
		wkDel_seg = seg[seg.selectionCoeff > -0.001];

		//get fixed muts
		fixed = muts[freq == 1.0];
		
		lethal_fix = fixed[fixed.selectionCoeff==-1];
		sublet_fix = fixed[fixed.selectionCoeff > -1 & fixed.selectionCoeff <= -0.1 ];
		strDel_fix = fixed[fixed.selectionCoeff <= -0.01 & fixed.selectionCoeff > -0.1 ];
		modDel_fix = fixed[fixed.selectionCoeff <= -0.001 & fixed.selectionCoeff > -0.01 ];
		wkDel_fix = fixed[fixed.selectionCoeff > -0.001];
		
		// get load as 1 - fitness
		// SLiM calculates fitness multiplicatively each generation and caches it
		genetic_load = 1-mean(subpop.cachedFitness(NULL));
		
		// calculate inbreeding load
		q = subpop.genomes.mutationFrequenciesInGenomes(muts);
		s = -muts.selectionCoeff;

		// replace mutations with s>1.0 with 1.0 (can happen when drawing from gamma distribution)
		s[s>1.0]=1.0;
	

		// get h for each mutation
		// note that this will not work if changing h using fitness callbacks
		//h=muts.mutationType.dominanceCoeff;
		//
		if(muts.size() == 0){
		h=float(0);
		} else {
		h=sapply(muts, "applyValue.getValue('dom');");
		}
		inbreeding_load = sum(q*s)-sum(q^2*s)-2*sum(q*(1-q)*s*h);
		genetic_pop = sum(q*s);
		relative_load = sum(q^2*s)+2*sum(q*(1-q)*s*h);

		heterozygosity = calcHeterozygosity(subpop.genomes);

		writeFile("saola_breeding_early.csv", subpop.name + "," + sim.cycle + "," + subpop.individuals.size() + "," + heterozygosity + "," + genetic_load + "," + inbreeding_load + "," + genetic_pop + "," + relative_load + "," + lethal_seg.length()  + "," + sublet_seg.length() + "," + strDel_seg.length() + "," + modDel_seg.length() + "," + wkDel_seg.length() + "," + lethal_fix.length() + "," + sublet_fix.length() + "," + strDel_fix.length() + "," + modDel_fix.length() + "," + wkDel_fix.length() + "," + males + "," + females, append=T);
		
	}

} 



${now_time1}:${future_time} late(){

	muts = sim.mutationsOfType(m1);
	for (subpop in sim.subpopulations){

		if (!subpop.getValue("active"))
			next;

		if (subpop.individuals.size() == 0){
			// for cases where all individuals died during fitness seleciton
			subpop.setValue("active", F);
			subpop.setValue("survives", "No");
			sim.killIndividuals(subpop.individuals);
			next;
		}

		females = size(subpop.individuals.sex[subpop.individuals.sex=='F']);
		males = size(subpop.individuals.sex[subpop.individuals.sex=='M']);

		//get segregating muts
		freq = subpop.genomes.mutationFrequenciesInGenomes(muts);
		seg = muts[freq != 1.0 & freq > 0];
		
		lethal_seg = seg[seg.selectionCoeff==-1];
		sublet_seg = seg[seg.selectionCoeff > -1 & seg.selectionCoeff <= -0.1 ];
		strDel_seg = seg[seg.selectionCoeff <= -0.01 & seg.selectionCoeff > -0.1 ];
		modDel_seg = seg[seg.selectionCoeff <= -0.001 & seg.selectionCoeff > -0.01 ];
		wkDel_seg = seg[seg.selectionCoeff > -0.001];

		//get fixed muts
		fixed = muts[freq == 1.0];
		
		lethal_fix = fixed[fixed.selectionCoeff==-1];
		sublet_fix = fixed[fixed.selectionCoeff > -1 & fixed.selectionCoeff <= -0.1 ];
		strDel_fix = fixed[fixed.selectionCoeff <= -0.01 & fixed.selectionCoeff > -0.1 ];
		modDel_fix = fixed[fixed.selectionCoeff <= -0.001 & fixed.selectionCoeff > -0.01 ];
		wkDel_fix = fixed[fixed.selectionCoeff > -0.001];
		
		// get load as 1 - fitness
		// SLiM calculates fitness multiplicatively each generation and caches it
		genetic_load = 1-mean(subpop.cachedFitness(NULL));
		
		// calculate inbreeding load
		q = subpop.genomes.mutationFrequenciesInGenomes(muts);
		s = -muts.selectionCoeff;

		// replace mutations with s>1.0 with 1.0 (can happen when drawing from gamma distribution)
		s[s>1.0]=1.0;
	

		// get h for each mutation
		// note that this will not work if changing h using fitness callbacks
		//h=muts.mutationType.dominanceCoeff;
		//
		if(muts.size() == 0){
		h=float(0);
		} else {
		h=sapply(muts, "applyValue.getValue('dom');");
		}
		inbreeding_load = sum(q*s)-sum(q^2*s)-2*sum(q*(1-q)*s*h);
		genetic_pop = sum(q*s);
		relative_load = sum(q^2*s)+2*sum(q*(1-q)*s*h);
		
		heterozygosity = calcHeterozygosity(subpop.genomes);

		writeFile("saola_breeding_late.csv", subpop.name + "," + sim.cycle + "," + subpop.individuals.size() +  "," + heterozygosity + "," + genetic_load + "," + inbreeding_load + "," + genetic_pop + "," + relative_load + "," + lethal_seg.length()  + "," + sublet_seg.length() + "," + strDel_seg.length() + "," + modDel_seg.length() + "," + wkDel_seg.length() + "," + lethal_fix.length() + "," + sublet_fix.length() + "," + strDel_fix.length() + "," + modDel_fix.length() + "," + wkDel_fix.length() + "," + males + "," + females, append=T);
		
		if(subpop.individuals.size() > 1000){
			// inactivate populations that group too much, these will survive
			subpop.setValue("active", F);
			sim.killIndividuals(subpop.individuals);
		}

		if(females == 0 | males == 0){
			// inactivate populations that die and set survival to NO
			subpop.setValue("active", F);
			subpop.setValue("survives", "No");
			sim.killIndividuals(subpop.individuals);
		}
	}

} 



${future_time1} first(){

	// last step
	writeFile("survival_log.txt", "admix4 admix8 admix12 admix16 admix20 admix24 admix28 admix32 AM4 AM8 AM12 AM16 AM20 AM24 AM28 AM32 AH4 AH8 AH12 AH16 AH20 AH24 AH28 AH32", append=F);
	writeFile("survival_log.txt", p4.getValue("survives") + " " +  p5.getValue("survives") + " " + p6.getValue("survives") + " " + p7.getValue("survives") + " " +  p8.getValue("survives") + " " + p9.getValue("survives") + " " + p10.getValue("survives") + " " +  p11.getValue("survives") + " " + p12.getValue("survives") + " " + p13.getValue("survives") + " " +  p14.getValue("survives") + " " + p15.getValue("survives") + " " + p16.getValue("survives") + " " +  p17.getValue("survives") + " " + p18.getValue("survives") + " " + p19.getValue("survives") + " " +  p20.getValue("survives") + " " + p21.getValue("survives") + " " + p22.getValue("survives") + " " +  p23.getValue("survives") + " " + p24.getValue("survives")  + " " + p25.getValue("survives") + " " +  p26.getValue("survives") + " " + p27.getValue("survives"), append=T);


}

EOM
