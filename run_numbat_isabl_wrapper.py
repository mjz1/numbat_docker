from datetime import datetime
import pandas as pd
import os
import isabl_utils
import subprocess
import argparse
from pathlib import Path

parser = argparse.ArgumentParser(description='Wrapper for running numbat on isabl projects')

# Required Arguments
parser.add_argument("--outdir", help="Output directory", type=str)
parser.add_argument("--isabl_project", help="Isabl project name")
parser.add_argument("--isabl_assay", help="Isabl assay")
parser.add_argument("--mem", help="Amount of memory to use per core", default=8)
parser.add_argument("--cores", help="Number of cores", default=4)
parser.add_argument("--combine_patients", help="Flag to combine patient samples. Allowable values: single, combine, both", default="both")
parser.add_argument("--combine_mem", default=30, help="For combination runs how much memory to use per core")

# Parse the arguments
args = parser.parse_args()

# Process args
outdir = args.outdir
isabl_project = args.isabl_project
isabl_assay = args.isabl_assay
mem = args.mem
cores = args.cores
combine_patients = args.combine_patients
combine_mem = args.combine_mem

print(isabl_assay)
print(isabl_project)


if combine_patients not in {'single', 'combine', 'both'}:
    print("combine_patients argument incorrect:")
    parser.print_help()
    quit()
    
if isabl_project==None or isabl_assay==None:
    parser.print_help()
    quit("ERROR: Please provide isabl project and assay")


# Get needed filepaths to run the pipeline
numbat_py = "/home/zatzmanm/work/repos/run_numbat/run_numbat.py"
numbat_rscript = f"/home/zatzmanm/work/repos/run_numbat/numbat.R"
numbat_refdir = "/home/zatzmanm/work/references/numbat/"
numbat_img = "/home/zatzmanm/work/images/numbat/numbat_0.1.1-2022-04-25-a4e75ca9f783.sif"
pileup_script = f"{numbat_refdir}/numbat/pileup_and_phase.R"

genome_ver = "hg38"
gmap = f"/src/Eagle_v2.4.1/tables/genetic_map_{genome_ver}_withX.txt.gz"  # Within the dockerimage
panel_dir = f"{numbat_refdir}/1000G_{genome_ver}/"
snp_vcf = f"{numbat_refdir}/genome1K.phase3.SNP_AF5e2.chr1toX.{genome_ver}.vcf.gz"

# 0. Get the cohort data from isabl and build a dict with the inputs
file_loc = isabl_utils.get_paths(isabl_assay, isabl_project, details=True)

file_loc['sample_id'] = file_loc['target_sample'].str[0]

# For each patient want names of each sample
pt_dict = {}
# Loop over each unique samples
for i, samp in enumerate(set(file_loc['sample_id'])):
    print(i, samp)
    patient = list(set(file_loc.loc[file_loc['sample_id'] == samp, 'individual']))[0]
    bam = file_loc.loc[(file_loc['sample_id'] == samp) & (file_loc['file_type'] == 'bam'), 'path'].iloc[0]
   
    # Initialize the patient and sample dicts if it doesn't already exist
    if patient not in pt_dict:
        pt_dict[patient] = {}
    if samp not in pt_dict[patient]:
        pt_dict[patient][samp] = {}
    
    # BAM file    
    pt_dict[patient][samp]['bam'] = bam
    
    # Cellranger 10x files
    filt_mtx_dir = f"{os.path.dirname(bam)}/filtered_feature_bc_matrix/"
    pt_dict[patient][samp]['filt_mtx_dir'] = filt_mtx_dir
    pt_dict[patient][samp]['barcodes'] = f"{filt_mtx_dir}/barcodes.tsv.gz"
    
# Now loop over the dict and submit numbat for each patient
for i, pt in enumerate(pt_dict.keys()):
    samples = []
    bams = []
    barcodes = []
    mtx_dirs = []
    outdir_pt = f"{outdir}/{pt}"
    
    Path(outdir_pt).mkdir(parents=True, exist_ok=True)
    
    # Loop over each sample and construct the inputs
    for j, samp in enumerate(pt_dict[pt].keys()):
        # Create lists of information for combined run
        samples.append(samp)
        bams.append(pt_dict[pt][samp]['bam'])
        barcodes.append(pt_dict[pt][samp]['barcodes'])
        mtx_dirs.append(pt_dict[pt][samp]['filt_mtx_dir'])
        
        outdir_samp = f"{outdir_pt}/{samp}/"
        
        if combine_patients in {'single', 'both'}:
            # Run each sample individually
            cmd=f"""python {numbat_py} \
                    --numbat_img {numbat_img} \
                    --numbat_refdir {numbat_refdir} \
                    --pileup_script {pileup_script} \
                    --numbat_rscript {numbat_rscript} \
                    --gmap {gmap} \
                    --panel_dir {panel_dir} \
                    --snp_vcf {snp_vcf} \
                    --mem {mem} \
                    --walltime 72:00 \
                    --patient {samp} \
                    --samples {samp} \
                    --bams {pt_dict[pt][samp]['bam']} \
                    --barcodes {pt_dict[pt][samp]['barcodes']} \
                    --mtxdirs {pt_dict[pt][samp]['filt_mtx_dir']} \
                    --outdir {outdir_samp} \
                    --genome_ver {genome_ver} \
                    --cores {cores}"""
                    
            subprocess.run(cmd, shell=True, check=True)
        
    if combine_patients in {'combine', 'both'}:
        # Check if there is more than one sample
        if len(samples) > 1:
            print(f"Running samples: {samples} for patient {pt} in combined mode")
    
             # Run samples merged across patient
            cmd=f"""python {numbat_py} \
                        --numbat_img {numbat_img} \
                        --numbat_refdir {numbat_refdir} \
                        --pileup_script {pileup_script} \
                        --numbat_rscript {numbat_rscript} \
                        --gmap {gmap} \
                        --panel_dir {panel_dir} \
                        --snp_vcf {snp_vcf} \
                        --mem {combine_mem} \
                        --walltime 144:00 \
                        --patient {pt} \
                        --samples {','.join(samples)} \
                        --bams {','.join(bams)} \
                        --barcodes {','.join(barcodes)} \
                        --mtxdirs {','.join(mtx_dirs)} \
                        --outdir {outdir_pt}/combined \
                        --genome_ver {genome_ver} \
                        --cores {cores}"""
                    
        
            subprocess.run(cmd, shell=True, check=True)
        