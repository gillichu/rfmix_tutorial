import subprocess
import os

def to_txt_file(poplist, filename):
    with open(filename, "w+") as f:
        for pop in poplist:
            f.write(pop + "\n")

def run_conversion(parfile):
    if not os.path.exists("data"):
        os.mkdir("data")
    subprocess.call(["python2", "conversion_real.py", "-p", parfile])

def make_conversion_par(parfile, chrm, query_root, ref_root, poplist, ancestrylist, output_tag):
    with open(parfile, "w+") as f:
        f.write("basepath:./\n")
        f.write("chrm:%s\n" % str(chrm))
        f.write("query_root:%s\n" % str(query_root))
        f.write("ref_root:%s\n" % str(ref_root))
        f.write("poplist:%s\n" % str(poplist))
        f.write("ancestrylist:%s\n" % str(ancestrylist))
        f.write("output_tag:%s\n" % str(output_tag))
        f.write("admixedanc: False\n")
        f.write("windows: 0.2\n") # old parameter, unused, ignore here

def run_rfmix(parfile, chrm, output_tag, w_size, logfile):
    with open(logfile, "w+") as f:
        subprocess.call(["python2", "RunRFMix.py", "PopPhased", "data/alleles%s%s.txt" % (chrm, output_tag), "data/classes%s%s.txt" % (chrm, output_tag), "data/markerLocations%s%s.txt" % (chrm, output_tag), "-o", "outputPopPhased%s%s" % (chrm, output_tag), "-w", str(w_size), "--forward-backward"], stdout=f, stderr=f)
    print("Done running RFMix.")

