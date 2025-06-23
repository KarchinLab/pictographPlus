import os
import argparse

def make_dir(directory):
    """
    Safely create a directory if it does not exist.
    """
    try:
        os.mkdir(directory)
    except FileExistsError:
        pass

def generate_het_positions(haplotype_vcf, output_dir, normal_sample):
    """
    Parse a haplotype VCF and extract germline heterozygous sites.
    Write two files:
      1) <normal_sample>_germline_het.txt  (detailed info: chrom, pos, ref, alt, refCount, altCount)
      2) germline_het_pos.txt              (just chrom pos for mpileup)
    """
    pileup_dir = os.path.join(output_dir.rstrip("/"), "pileup")
    make_dir(pileup_dir)

    # e.g. Normal sample specific file name
    germline_het_file = os.path.join(pileup_dir, f"{normal_sample}_germline_het.txt")
    germline_het_pos_file = os.path.join(pileup_dir, "germline_het_pos.txt")

    with open(germline_het_file, "w") as out1, open(germline_het_pos_file, "w") as out2:
        # Header in the normal-sample file
        out1.write("chrom\tpos\tref\talt\trefCount\taltCount\n")

        with open(haplotype_vcf, "r") as f:
            for line in f:
                line = line.strip()
                if line.startswith("#"):
                    continue
                parts = line.split("\t")
                # VCF columns: CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, SAMPLE
                chrom, pos, ref, alt = parts[0], parts[1], parts[3], parts[4]

                # We'll look in parts[9], the genotype data for the normal sample (assuming 1 sample VCF)
                # For multi-sample VCF, you'd have to adjust the index accordingly.
                sample_info = parts[9].split(":")  # e.g. GT:AD:...
                genotype = sample_info[0]         # e.g. "0/1"
                
                # Make sure we only handle single-nucleotide REF/ALT
                if len(ref) == 1 and len(alt) == 1:
                    # Check that we have a heterozygote genotype (0/1 or 0|1, etc.)
                    if "0" in genotype and "1" in genotype:
                        # 'AD' field is typically index 1 if FORMAT is GT:AD:DP:...
                        # This might vary depending on your VCF's FORMAT order.
                        if len(sample_info) > 1:
                            ad_counts = sample_info[1].split(",")  # e.g. ["10", "12"]
                            if len(ad_counts) == 2:
                                ref_count, alt_count = ad_counts
                                out1.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{ref_count}\t{alt_count}\n")
                                out2.write(f"{chrom} {pos}\n")

def run_pileup_and_awk(bam_file, fasta, pos_file, output_dir):
    """
    For a given BAM, run samtools mpileup and then an awk command to summarize
    read counts for each base at each position.
    """
    tumor_name = os.path.basename(bam_file).replace(".bam", "")
    pileup_txt = os.path.join(output_dir, f"{tumor_name}_pileup.txt")
    pileup_summary = os.path.join(output_dir, f"{tumor_name}_pileup_summary.txt")

    # 1) Run samtools mpileup
    mpileup_cmd = (
        f"samtools mpileup -f {fasta} -l {pos_file} {bam_file} > {pileup_txt}"
    )
    print("Running:", mpileup_cmd)
    os.system(mpileup_cmd)

    # 2) Run awk to parse read counts
    # Here, we interpret each line from mpileup:
    #   col1=chrom, col2=pos, col3=refBase, col4=depth, col5=bases, col6=quality
    # The AWK logic you used:
    #   ref = $3
    #   ref_count = gsub(/\./, "", $5) + gsub(/,/, "", $5)
    #   a = gsub(/[aA]/, "", $5)
    #   t = gsub(/[tT]/, "", $5)
    #   c = gsub(/[cC]/, "", $5)
    #   g = gsub(/[gG]/, "", $5)
    #   print chrom, pos, "Ref:", ref, ref_count, "A:", a, "T:", t, "C:", c, "G:", g
    awk_cmd = (
        f"awk '{{ref=$3;ref_count=gsub(/\\./, \"\", $5) + gsub(/,/, \"\", $5); "
        f"a=gsub(/[aA]/, \"\", $5); t=gsub(/[tT]/, \"\", $5); c=gsub(/[cC]/, \"\", $5); g=gsub(/[gG]/, \"\", $5); "
        f"print $1, $2, \"Ref:\", ref, ref_count, \"A:\", a, \"T:\", t, \"C:\", c, \"G:\", g }}' "
        f"{pileup_txt} > {pileup_summary}"
    )
    print("Running:", awk_cmd)
    os.system(awk_cmd)

def process_and_combine(normal_sample, tumor_samples, output_dir, minreads, vaf_range):
    """
    Process the normal-sample germline het file, filter based on minreads & VAF,
    then parse each tumor's pileup_summary to match those sites, and finally
    write out a combined germline_SNV.csv.
    """
    pileup_dir = os.path.join(output_dir.rstrip("/"), "pileup")
    
    # 1) Make sure the normal file from generate_het_positions exists
    normal_het_file_name = f"{normal_sample}_germline_het.txt"
    if normal_het_file_name not in os.listdir(pileup_dir):
        print(f"[ERROR] Expected normal file {normal_het_file_name} not found in {pileup_dir}. Exiting.")
        exit(1)
    
    normal_het_file = os.path.join(pileup_dir, normal_het_file_name)
    
    # We'll store every site that passes minreads + VAF in `hetList`
    # We'll also store a dictionary: hetDict[sampleName][pos] = [refCount, altCount, refBase, altBase] (for normal)
    # For tumor, it will be  [refCount, altCount].
    hetList = []
    hetDict = {normal_sample: {}}

    # 2) Read & filter normal germline_het.txt
    with open(normal_het_file, "r") as f:
        header = f.readline()  # skip header line "chrom pos ref alt refCount altCount"
        for line in f:
            line = line.strip().split("\t")
            chrom, pos, ref_base, alt_base = line[0], line[1], line[2], line[3]
            ref_count, alt_count = int(line[4]), int(line[5])

            # Filter by minreads
            if ref_count >= minreads and alt_count >= minreads:
                tot = ref_count + alt_count
                vaf = alt_count / tot
                # Filter by VAF range
                if vaf_range[0] <= vaf <= vaf_range[1]:
                    pos_key = f"{chrom}-{pos}"
                    hetList.append(pos_key)
                    # Store relevant info for normal
                    hetDict[normal_sample][pos_key] = [ref_count, alt_count, ref_base, alt_base]

    # 3) Now process each tumor
    #    We'll parse the <tumor_name>_pileup_summary.txt that was generated by run_pileup_and_awk()
    for tumor in tumor_samples:
        tumor_name = os.path.basename(tumor).replace(".bam", "")
        tumor_pileup_summary = os.path.join(pileup_dir, f"{tumor_name}_pileup_summary.txt")
        
        # Initialize the dictionary for this tumor
        hetDict[tumor_name] = {}

        # The columns in <tumor_name>_pileup_summary.txt look like:
        # [chr, pos, "Ref:", REF, ref_count, "A:", a_count, "T:", t_count, "C:", c_count, "G:", g_count]
        # indexes: 0=chr, 1=pos, 2="Ref:", 3=REF, 4=ref_count, 5="A:", 6=a_count, 7="T:", 8=t_count, ...
        idx_dict = {"A": 6, "T": 8, "C": 10, "G": 12}
        
        if not os.path.exists(tumor_pileup_summary):
            print(f"[WARNING] Tumor pileup summary not found: {tumor_pileup_summary}")
            continue
        
        with open(tumor_pileup_summary, "r") as tf:
            for line in tf:
                parts = line.strip().split()
                if len(parts) < 13:
                    continue  # skip incomplete lines

                chrom = parts[0]
                position = parts[1]
                ref_base = parts[3]
                ref_count = int(parts[4])
                
                pos_key = f"{chrom}-{position}"
                # Only consider the site if it was a previously identified normal het site
                if pos_key in hetList:
                    # Check that the reference base from tumor matches reference base from normal
                    normal_ref_base = hetDict[normal_sample][pos_key][2].upper()
                    if ref_base.upper() != normal_ref_base:
                        raise ValueError(
                            f"Reference mismatch for position {pos_key} "
                            f"(normal={normal_ref_base}, tumor={ref_base})"
                        )
                    alt_base = hetDict[normal_sample][pos_key][3]  # from normal
                    alt_index = idx_dict.get(alt_base.upper())
                    if alt_index is not None:
                        alt_count = int(parts[alt_index])
                        # If the tumor also meets minreads on ref & alt
                        if ref_count >= minreads and alt_count >= minreads:
                            # Store in dictionary
                            hetDict[tumor_name][pos_key] = [ref_count, alt_count]

    # 4) Determine the set of positions present in all samples
    #    We only want positions that appear in normal & *all* tumors
    all_pos_temp = []
    for sample_name in hetDict.keys():
        if len(all_pos_temp) == 0:
            all_pos_temp = list(hetDict[sample_name].keys())
        else:
            all_pos_temp = list(set(all_pos_temp).intersection(list(hetDict[sample_name].keys())))

    # among those, only keep the ones that originally passed normal filters (hetList)
    final_positions = [p for p in hetList if p in all_pos_temp]

    # 5) Write final germline_SNV.csv
    # pictograph_dir = os.path.join(output_dir.rstrip("/"), "pictograph")
    # make_dir(pictograph_dir)
    output_file = os.path.join(pileup_dir, "germline_SNV.csv")

    with open(output_file, "w") as out:
        # Write header
        out.write("chroms,position,ref,alt,germline_ref,germline_alt")
        for sample_name in hetDict.keys():
            if sample_name == normal_sample:
                continue
            out.write(f",{sample_name}_ref,{sample_name}_alt")
        out.write("\n")

        # Write rows
        for pos in final_positions:
            chrom, position = pos.split("-")
            ref_base = hetDict[normal_sample][pos][2]
            alt_base = hetDict[normal_sample][pos][3]
            germline_ref_count = hetDict[normal_sample][pos][0]
            germline_alt_count = hetDict[normal_sample][pos][1]

            row = [
                chrom,
                position,
                ref_base,
                alt_base,
                str(germline_ref_count),
                str(germline_alt_count)
            ]
            for sample_name in hetDict.keys():
                if sample_name == normal_sample:
                    continue
                ref_alt_counts = hetDict[sample_name][pos]
                row.append(str(ref_alt_counts[0]))
                row.append(str(ref_alt_counts[1]))

            out.write(",".join(row) + "\n")

    print(f"[INFO] Final germline_SNV.csv written to: {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Generate pileups and combine SNV calls.")
    parser.add_argument("-v", "--haplotype_vcf", required=True,
                        help="Path to haplotype VCF file (normal sample).")
    parser.add_argument("-b", "--tumor_bam_files", nargs='+', required=True,
                        help="Paths to one or more tumor BAM files.")
    parser.add_argument("-o", "--output_dir", required=True,
                        help="Output directory.")
    parser.add_argument("-f", "--fasta", required=True,
                        help="Reference genome fasta file.")
    parser.add_argument("--minreads", type=int, default=3,
                        help="Minimum read count for both ref and alt to keep a site.")
    parser.add_argument("--vaf", type=float, nargs=2, default=[0.1, 0.9],
                        help="Minimum and maximum VAF for normal heterozygous sites.")
    args = parser.parse_args()

    # 1) Generate the heterozygous positions from the haplotype VCF
    normal_sample_name = "normal"
    generate_het_positions(args.haplotype_vcf, args.output_dir, normal_sample_name)

    # 2) Run samtools mpileup + AWK for each tumor
    pileup_dir = os.path.join(args.output_dir.rstrip("/"), "pileup")
    pos_file = os.path.join(pileup_dir, "germline_het_pos.txt")

    # For each tumor, run mpileup/awk
    for bam in args.tumor_bam_files:
        run_pileup_and_awk(bam, args.fasta, pos_file, pileup_dir)

    # 3) Process the normal germline_het.txt and tumor pileups, then combine into final germline_SNV.csv
    process_and_combine(
        normal_sample=normal_sample_name,
        tumor_samples=args.tumor_bam_files,
        output_dir=args.output_dir,
        minreads=args.minreads,
        vaf_range=args.vaf
    )

if __name__ == "__main__":
    main()
