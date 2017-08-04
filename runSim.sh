#/bin/bash
sim=$1
#Extract 1K random SNPs
plink19 --allow-no-sex --bfile geno/GreenlandSiblings.autosome.geno.0.01.withAdmix --thin-count 1000 --out plink/1K.sim.$sim.SNPs --make-bed

#Add betas and extract 95% tails
Rscript trunc.sim/addBetas.R plink/1K.sim.$sim.SNPs.bim trunc.sim/sim.$sim.truncSNPs

#Simulate phenotypes for each h
for h in 0.4 0.6 0.8
do
    gcta64 --bfile geno/GreenlandSiblings.autosome.geno.0.01.withAdmix --simu-qt --simu-causal-loci trunc.sim/sim.$sim.truncSNPs.withBetas --simu-hsq $h --out trunc.sim/sim.$sim.h.$h.gctaSim
done

#Calculate local h estimates from geno only
##All genos first
plink19 --allow-no-sex --bfile geno/GreenlandSiblings.autosome.geno.0.01.withAdmix --extract trunc.sim/sim.$sim.truncSNPs --recode A --out plink/sim.$sim.all
#Extra important step for getting CAF for causal SNPs:
plink19 --allow-no-sex --bfile geno/GreenlandSiblings.autosome.geno.0.01.withAdmix --extract trunc.sim/sim.$sim.truncSNPs --freq --out plink/sim.$sim.all
Rscript ./creating.own.pheno.R  plink/sim.$sim.all.raw trunc.sim/sim.$sim.truncSNPs.withBetas plink/sim.$sim.all.frq trunc.sim/sim.$sim.all.genoSum

##Proportions for both danes and inuit
for admix in 0.5.to.0.6 0.6.to.0.7 0.7.to.0.8
do
    for pop in inuit danes
    do
	plink19 --allow-no-sex --keep-allele-order --bfile geno/GreenlandSiblings.autosome.geno.0.01.withAdmix --extract trunc.sim/sim.$sim.truncSNPs --keep plink/${admix}.$pop --recode A --out plink/sim.$sim.${admix}.$pop
	Rscript ./creating.own.pheno.R plink/sim.$sim.${admix}.$pop.raw trunc.sim/sim.$sim.truncSNPs.withBetas plink/sim.$sim.all.frq trunc.sim/sim.$sim.${admix}.$pop.genoSum
    done
done
##Proportions only fitting for inuit
for admix in 0.8.to.0.9 0.9.to.0.99 0.99.to.1
do
    pop=inuit
    plink19 --allow-no-sex --keep-allele-order --bfile geno/GreenlandSiblings.autosome.geno.0.01.withAdmix --extract trunc.sim/sim.$sim.truncSNPs --keep plink/${admix}.$pop --recode A --out plink/sim.$sim.${admix}.$pop
    Rscript ./creating.own.pheno.R plink/sim.$sim.${admix}.$pop.raw trunc.sim/sim.$sim.truncSNPs.withBetas plink/sim.$sim.all.frq trunc.sim/sim.$sim.${admix}.$pop.genoSum	    
done

#All genosums and phenos are now created for each h - time to collect and plot  
