#!/bin/sh
#SBATCH --job-name=rfmix_${admix}  # The job name.
#SBATCH --cpus-per-task=1       # The number of cpu cores to use per task
#SBATCH --partition=Batch          # partition to use for the job

echo " create par files, convert them, and run RFMix on it to output an accuracy"
basepath="/my/home/gchu/realdata/india/extracted_pops"

chrm=$1
gen=$2
admix=$3
prop=$4
numsnps=$5

echo ${admix} > admixed_${admix}${chrm}.txt

echo "WestEur" > ancestrylist_${admix}${chrm}.txt
echo "ASI" >> ancestrylist_${admix}${chrm}.txt

cat admixed_${admix}${chrm}.txt > poplist_${admix}${chrm}.txt
cat ancestrylist_${admix}${chrm}.txt >> poplist_${admix}${chrm}.txt

extract() {
	
	local SSS=$1
	local poplist=$2
	local chrm=$3
	local isref=$4
	local numsnps=$5

	FILE=${SSS}.phind
	if [ ! -f "${FILE}" ]; then
		scp ${basepath}/${SSS}.phind .
		sed -i 's/GBR/WestEur/g' ${SSS}.phind
		sed -i 's/Jarwa/ASI/g' ${SSS}.phind
		sed -i 's/Onge/ASI/g' ${SSS}.phind
                sed -i 's/CHB/ASI/g' ${SSS}.phind
	fi 
	
	badsnpname="badsnpname${chrm}_leave300k"
	if [ ! -f "${badsnpname}" ]; then 
                echo "filtering out the position 0 snps"
                awk ' { if($3=="0") print $1 } ' ${basepath}/${SSS}.phsnp > ${badsnpname} 
		scp ${basepath}/${SSS}.phsnp ${SSS}.phsnp
                wc -l "${badsnpname}" > t${chrm}
                num_zero=$( awk '{ print $1 }' t$chrm ) 
                rm t${chrm}

                echo "generating badsnpname_leave300k"
		wc -l "${SSS}.phsnp" > t$chrm 
		totalsnps=$( awk '{ print $1 }' t$chrm ) 
		rm t${chrm}
		echo "${SSS}.phsnp has ${totalsnps} snps"
		numBadSnps=$(( ${totalsnps} - ${num_zero} - ${numsnps} ))
		echo "${SSS}.phsnp has ${totalsnps} snps, with ${num_zero} snps at position 0, so I will leave out ${numBadSnps} in order to get the desired ${numsnps}" 
        	shuf -n ${numBadSnps} ${SSS}.phsnp >> ${badsnpname}
	fi 

	echo "indivname: ${SSS}.phind" > par_extract_${isref}${admix}${chrm}
        echo "snpname: ${SSS}.phsnp" >> par_extract_${isref}${admix}${chrm}
        echo "genotypename: ${basepath}/${SSS}.phgeno" >> par_extract_${isref}${admix}${chrm}
        echo "outputformat: EIGENSTRAT" >> par_extract_${isref}${admix}${chrm}
        echo "genotypeoutname: admix_ref_data/ga_chr${chrm}_${isref}_300k.phgeno" >> par_extract_${isref}${admix}${chrm}
        echo "indivoutname: admix_ref_data/ga_chr${chrm}_${isref}_300k.phind" >> par_extract_${isref}${admix}${chrm}
        echo "snpoutname: admix_ref_data/ga_chr${chrm}_${isref}_300k.phsnp" >> par_extract_${isref}${admix}${chrm}
        echo "phasedmode: haploid" >> par_extract_${isref}${admix}${chrm}
        echo "poplistname: ${poplist}" >> par_extract_${isref}${admix}${chrm}
	echo "badsnpname: ${badsnpname}" >> par_extract_${isref}${admix}${chrm}
	echo "chrom: ${chrm}" >> par_extract_${isref}${admix}${chrm}
	/my/home/pmoorjani/bin/convertf -p par_extract_${isref}${admix}${chrm}
}

SSS="ga_chr${chrm}" # orig chrm file
# extract reference pops
FILE="admix_ref_data/ga_chr${chrm}_ref_300k.phgeno"
mkdir ${outdir}
if [ ! -f "${FILE}" ]; then
	echo "extracting reference pops"
	extract ${SSS} ancestrylist_${admix}${chrm}.txt ${chrm} "ref" ${numsnps}
fi

# extract admixed pops
outdir="admix_ref_data"
if [ ! -f "admix_ref_data/ga_chr${chrm}_${admix}_300k.phgeno" ]; then
        echo "extracting admixed pops"
	extract ${SSS} admixed_${admix}${chrm}.txt ${chrm} ${admix} ${outdir}
fi 
# run rfmix conversion
echo "running rfmix conversion on ${admix} set" 
basepath="/my/home/gchu/realdata/india/rfmix_300k_again/admix_ref_data/"
output_tag="_${admix}_${chrm}_300k"

echo "basepath:${basepath}" > par_${admix}_${chrm}
echo "chrm:${chrm}" >> par_${admix}_${chrm}
echo "query_root: ga_chr${chrm}_${admix}_300k" >> par_${admix}_${chrm}
echo "ref_root: ga_chr${chrm}_ref_300k" >> par_${admix}_${chrm}
echo "poplist:poplist_${admix}${chrm}.txt" >> par_${admix}_${chrm}
echo "ancestrylist:ancestrylist_${admix}${chrm}.txt" >> par_${admix}_${chrm}
echo "output_tag:${output_tag}" >> par_${admix}_${chrm}
echo "admixedanc: False" >> par_${admix}_${chrm}
echo "windows: 0.01"

if [ ! -f "outputPopPhased${chrm}_${admix}_${chrm}_300k.0.ForwardBackward.txt" ]; then
	echo "Running conversion..."
	python2 /my/home/gchu/lai/RFMix_v1.5.4/conversion_real.py -p par_${admix}_${chrm}
	echo "Running RFMix..."
	python2 RunRFMix.py PopPhased data/alleles${chrm}${output_tag}.txt data/classes${chrm}${output_tag}.txt data/markerLocations${chrm}${output_tag}.txt  -o outputPopPhased${chrm}${output_tag} -w 0.01 --forward-backward -G ${gen} 
fi

echo "DONE."


