#/bin/bash
sim=$1

#Using same sim as for intervals - only calculate other genoSums

##No all genos first - they are calculated in the intervals

##Proportions for both danes and inuit
for admix in 0.5 0.6 0.7
do
    for pop in inuit danes
    do
	plink19 --keep-allele-order --allow-no-sex --bfile geno/GreenlandSiblings.autosome.geno.0.01.withAdmix --extract trunc.sim/sim.$sim.truncSNPs --keep plink/${admix}.$pop --recode A --out plink/sim.$sim.${admix}.$pop
	Rscript ./creating.own.pheno.R plink/sim.$sim.${admix}.$pop.raw trunc.sim/sim.$sim.truncSNPs.withBetas plink/sim.$sim.all.frq trunc.sim/sim.$sim.${admix}.$pop.genoSum
    done
done
##Proportions only fitting for inuit
for admix in 0.8 0.9 0.99
do
    pop=inuit
    plink19 --keep-allele-order --allow-no-sex --bfile geno/GreenlandSiblings.autosome.geno.0.01.withAdmix --extract trunc.sim/sim.$sim.truncSNPs --keep plink/${admix}.$pop --recode A --out plink/sim.$sim.${admix}.$pop
    Rscript ./creating.own.pheno.R plink/sim.$sim.${admix}.$pop.raw trunc.sim/sim.$sim.truncSNPs.withBetas plink/sim.$sim.all.frq trunc.sim/sim.$sim.${admix}.$pop.genoSum	    
done

#All genosums and phenos are now created for each h - time to collect and plot  

