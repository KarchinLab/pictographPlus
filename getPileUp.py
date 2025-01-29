import yaml
from yaml import load
try:
    from yaml import CLoader as Loader
except ImportError:
    from yaml import Loader

import argparse
import os

def generateHetPos(haplotype_vcf, output_dir):
    output_dir = output_dir.rstrip("/") + "/pileup"
    makeDir(output_dir)
    
    outputFile = open(f"{output_dir}/germline_het.txt", "w")
    outFile = open(f"{output_dir}/germline_het_pos.txt", "w")
    outputFile.write("chrom\tpos\tref\talt\trefCount\taltCount\n")
    
    with open(haplotype_vcf, "r") as f:
        for line in f:
            line = line.rstrip()
            if line.startswith("#"):
                continue
            line = line.split("\t")
            if len(line[3]) == 1 and len(line[4]) == 1:
                if "0" in line[9].split(":")[0] and "1" in line[9].split(":")[0]:
                    count = line[9].split(":")[1].split(",")
                    outputFile.write(f"{line[0]}\t{line[1]}\t{line[3]}\t{line[4]}\t{count[0]}\t{count[1]}\n")
                    outFile.write(f"{line[0]} {line[1]}\n")
    
    outputFile.close()
    outFile.close()

def makeDir(directory):
    try:
        os.mkdir(directory)
    except FileExistsError:
        pass

def getPileup(tumor_bam_files, fasta, pos_file, output_dir, submit):
    output_dir = output_dir.rstrip("/") + "/pileup"
    makeDir(output_dir)
    bashDir = output_dir + "/bash"
    makeDir(bashDir)

    bashFilePath = f"{bashDir}/pileup.sh"
        
    with open(bashFilePath, "w") as bashFile:
        bashFile.write("#!/bin/bash\n\n")
        bashFile.write(f"#SBATCH --job-name=pu\n")
        bashFile.write(f"#SBATCH --output={bashDir}/pu.out\n")
        bashFile.write("#SBATCH --mem=16g\n")
        bashFile.write("#SBATCH --time=1980:00:00\n\n")
        bashFile.write("\n\ndate\n\n")

    for tumor_bam in tumor_bam_files:
        tumor_name = os.path.basename(tumor_bam).replace(".bam", "")
        bashFile.write(f"\nsamtools mpileup -f {fasta} -l {pos_file} {tumor_bam} > {output_dir}/{tumor_name}_pileup.txt\n")
        bashFile.write(f"\nawk '{{ref=$3;ref_count=gsub(/\\./, \"\", $5) + gsub(/,/, \"\", $5); a=gsub(/[aA]/, \"\", $5); t=gsub(/[tT]/, \"\", $5); c=gsub(/[cC]/, \"\", $5); g=gsub(/[gG]/, \"\", $5); print $1, $2, \"Ref:\", ref, ref_count, \"A:\", a, \"T:\", t, \"C:\", c, \"G:\", g }}' {output_dir}/{tumor_name}_pileup.txt > {output_dir}/{tumor_name}_pileup_summary.txt\n")
    
    bashFile.write("\n\ndate\n\n")
        
    if submit:
        os.system(f"sbatch {bashFilePath}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generate pileup information from HaplotypeCaller and tumor BAM files")
    parser.add_argument("-v", "--haplotype_vcf", help="Path to haplotype VCF file", required=True)
    parser.add_argument("-b", "--tumor_bam_files", help="Paths to tumor BAM files", nargs='+', required=True)
    parser.add_argument("-o", "--output_dir", help="Output directory", required=True)
    parser.add_argument("-f", "--fasta", help="Reference genome fasta file", required=True)
    parser.add_argument("--submit", help="Submit job to scheduler", default=False, action=argparse.BooleanOptionalAction)
    args = parser.parse_args()
    
    generateHetPos(args.haplotype_vcf, args.output_dir)
    pos_file = f"{args.output_dir.rstrip('/')}/pileup/germline_het_pos.txt"
    getPileup(args.tumor_bam_files, args.fasta, pos_file, args.output_dir, args.submit)
