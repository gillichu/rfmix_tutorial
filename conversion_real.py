import argparse
import numpy as np

def parse_par(args):
     filename = args.par
     params = {}
     with open(filename) as f:
          for line in f:
               parts = line.split(":")
               parts = [part.strip() for part in parts]
               # print(parts)
               if True not in [bool(part) for part in parts]:
                    continue
               params[parts[0].strip()] = parts[1].strip()
     for param in ['output_tag', 'basepath', 'chrm', 'query_root', 'ref_root', 'poplist', 'ancestrylist', 'admixedanc']:
          assert param in params.keys(), "Must give required parameter: " + param
     params['basepath'] = str(params['basepath'])
     params['chrm'] = str(params['chrm'])
     params['query_root'] = str(params['query_root'])
     params['ref_root'] = str(params['ref_root'])
     params['poplist'] = str(params['poplist'])
     params['ancestrylist'] = str(params['ancestrylist'])
     params['output_tag'] = str(params['output_tag'])
     params['admixedanc'] = True if params["admixedanc"] == "True" else False
# for param in ['basepath', 'chrm', 'query_root', 'ref_root', 'output_tag']:
          # print(param, ":", params[param])
     return params

parser = argparse.ArgumentParser(description="Run a simulation given the parameters in the par file.")
parser.add_argument("-p", help=".par file to attach", dest="par", type=str, required=True)
parser.set_defaults(func=parse_par)
args = parser.parse_args()
params = args.func(args)

query_ind = params['basepath'] + params['query_root'] + ".phind"
query_geno = params['basepath'] + params['query_root'] + ".phgeno"
query_snp = params['basepath'] + params['query_root'] + ".phsnp" 

ref_ind = params['basepath'] + params['ref_root'] + ".phind"
ref_geno = params['basepath'] + params['ref_root'] + ".phgeno"
poplist = np.loadtxt(params['poplist'], dtype=str)
popdict = dict()

# map to index number in dictionary
for i, pop in enumerate(poplist):
     popdict[pop] = i

all_snps = dict()
with open(query_snp, "r") as f:
    z = f.readlines()
    for i, snpline in enumerate(z):
        snpline = snpline.split()
        snp_id, snp_chrm = snpline[0], snpline[1]
        if (snp_id in all_snps): 
            print(snp_id)
        if snp_chrm == params['chrm']:
            all_snps[snp_id] = i

print("Saved %s snps" % len(all_snps))

## row per snp # col per haplotype
## vals are 0 or 1
## no space between alleles
allelename = "data/" + "alleles" +  params['chrm'] + params['output_tag'] + ".txt"
with open(allelename, "w") as f:
    with open(query_geno, "r") as f2:
      q = f2.readlines()
      with open(query_snp, "r") as f3:
          z = f3.readlines()
          with open(ref_geno, "r") as f4:
              y = f4.readlines()
              for i, snpline in enumerate(z):
                  snpline = snpline.split()
                  if snpline[1] == params['chrm']:
                      snp_id = snpline[0]
                      snp_idx = all_snps[snp_id]
                      line = q[snp_idx].replace('\n', '') + y[snp_idx].replace('\n', '') + "\n"
                      f.write(line)

## one row, one col per haplotype
classesname = "data/" +'classes' +  params['chrm'] + params['output_tag'] + ".txt"
with open(classesname, "w") as f:
    strline = ""
    with open(ref_ind, "r") as f2:
        q = f2.readlines()
        with open(query_ind, "r") as f3:
            z = f3.readlines()
            for ind in z:
                ind = ind.split()
                pop = ind[2]
                strline += str(popdict[pop]) + " "
            for ind in q:
                ind = ind.split()
                pop = ind[2]
                strline += str(popdict[pop]) + " "
            f.write(strline + "\n") #print(strline)


## one row per SNP, one col
## expects them to be in order, RFMix expects centimorgans
markersname = "data/" +'markerLocations'  +  params['chrm'] + params['output_tag'] + ".txt"
with open(markersname, "w") as f:
    genpos_arr = []
    with open(query_snp, "r") as f2:
        z = f2.readlines()
        for snpline in z:
            snpline = snpline.split()
            if snpline[1] == params['chrm']:
                genpos = float(snpline[2]) * 100
                strline = str(genpos)
                genpos_arr.append(float(genpos))
                f.write(strline + "\n") # print(strline)

print("Done!")



