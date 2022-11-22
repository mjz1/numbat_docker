from datetime import datetime
import pandas as pd
import os
import isabl_utils
import subprocess
import argparse
from pathlib import Path

parser = argparse.ArgumentParser(description='Wrapper for running numbat on isabl projects', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Required Arguments
parser.add_argument("--outdir", help="Output directory", type=str)
parser.add_argument("--isabl_project", help="Isabl project name")
parser.add_argument("--isabl_assay", help="Isabl assay")
parser.add_argument("--numbat_img", help="The numbat image file", default="/home/zatzmanm/work/images/numbat/numbat-rbase_latest.sif")
parser.add_argument("--pileup_script", help="The numbat preprocessing pileup and phasing Rscript", default="/numbat/inst/bin/pileup_and_phase.R")
parser.add_argument("--numbat_rscript", help="The Rscript to run numbat", default="/home/zatzmanm/work/repos/run_numbat/numbat.R")
parser.add_argument("--numbat_py", help="The python script that wraps numbat", default="/home/zatzmanm/work/repos/run_numbat/run_numbat.py")
parser.add_argument("--mem", help="Amount of memory to use per core", default=8)
parser.add_argument("--cores", help="Number of cores", default=4)
parser.add_argument("--combine_patients", help="Flag to combine patient samples. Allowable values: single, combine, both", choices=["single","combine", "both"], default="both")
parser.add_argument("--combine_mem", default=30, help="For combination runs how much memory to use per core")
parser.add_argument("--walltime", help="amount of walltime to use", default="48:00")
parser.add_argument("--combine_walltime", help="amount of walltime to use in combined mode", default="72:00")
parser.add_argument("--high_purity", help="Flag to detect and exclude regions of clonal deletions/LOH before running Numbat. Recommended for cell line data or high-purity tumors", action='store_true', default=False)
parser.add_argument("--dry_run", help="Flag to simply print all samples that will run, but not run them", action='store_true', default=False)
parser.add_argument("--target_column", help="Column in the isabl df to use as the sample id. Usually 'target_sample' or 'target_aliquot' when samples have multiple associated samples.", default='target_sample')
parser.add_argument("--trans", help="Numbat HMM transmission probability", default=1e-5, type=float)
parser.add_argument("--gamma", help="Numbat overdispersion parameter in allele counts", default=20)
parser.add_argument("--min_cells", help="Numbat minimum number of cells for which an pseudobulk HMM will be run", default=20)
parser.add_argument("--min_LLR", help="Numbat minimum log-likelihood ratio threshold to filter CNVs by. ", default=50)
parser.add_argument("--init_k", help="Number of clusters in the initial clustering", default=3, type = "integer")
parser.add_argument("--max_iter", help="Maximum number of iterations to run the phyologeny optimization", default=2, type = "integer")
parser.add_argument("--max_entropy", help="Entropy threshold to filter CNVs", default=0.5)


# Parse the arguments
args = parser.parse_args()

# Process args
outdir = args.outdir
isabl_project = args.isabl_project
numbat_img = args.numbat_img
pileup_script = args.pileup_script
numbat_py = args.numbat_py
isabl_assay = args.isabl_assay
mem = args.mem
cores = args.cores
combine_patients = args.combine_patients
combine_mem = args.combine_mem
high_purity = args.high_purity
numbat_rscript = args.numbat_rscript
numbat_py = args.numbat_py
dry_run = args.dry_run
target_column = args.target_column
walltime = args.walltime
combine_walltime = args.combine_walltime
trans = args.trans
gamma = args.gamma
min_cells = args.min_cells
min_LLR = args.min_LLR
init_k = args.init_k
max_iter = args.max_iter
max_entropy = args.max_entropy

# Set high purity flag
if high_purity:
    high_purity="--high_purity"
else:
    high_purity=""

print(f"Looking for assay: {isabl_assay} in project: {isabl_project}")


if combine_patients not in {'single', 'combine', 'both'}:
    print("combine_patients argument incorrect:")
    parser.print_help()
    quit()
    
if isabl_project==None or isabl_assay==None:
    parser.print_help()
    quit("ERROR: Please provide isabl project and assay")

# 0. Get the cohort data from isabl and build a dict with the inputs
file_loc = isabl_utils.get_paths(isabl_assay, isabl_project, details=True)

file_loc['sample_id'] = file_loc[target_column].str[0]

# For each patient want names of each sample
pt_dict = {}
# Loop over each unique samples
for i, samp in enumerate(set(file_loc['sample_id'])):
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
    
    if dry_run==True:
        for j, samp in enumerate(pt_dict[pt].keys()):
            # Create lists of information for combined run
            print(f"Donor: {pt} -- Sample: {samp}")
            
        continue

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
                    --pileup_script {pileup_script} \
                    --numbat_rscript {numbat_rscript} \
                    --mem {mem} \
                    --walltime {walltime} \
                    --patient {samp} \
                    --samples {samp} \
                    --bams {pt_dict[pt][samp]['bam']} \
                    --barcodes {pt_dict[pt][samp]['barcodes']} \
                    --mtxdirs {pt_dict[pt][samp]['filt_mtx_dir']} \
                    --outdir {outdir_samp} \
                    --trans {trans} \
                    --gamma {gamma} \
                    --min_cells {min_cells} \
                    --min_LLR {min_LLR} \
                    --max_iter {max_iter} \
                    --init_k {init_k} \
                    --max_entropy {max_entropy} \  
                    --cores {cores} \
                    {high_purity}"""
                    
            subprocess.run(cmd, shell=True, check=False)
            
        
    if combine_patients in {'combine', 'both'}:
        # Check if there is more than one sample
        if len(samples) > 1:
            print(f"Running samples: {samples} for patient {pt} in combined mode")
    
             # Run samples merged across patient
            cmd=f"""python {numbat_py} \
                        --numbat_img {numbat_img} \
                        --pileup_script {pileup_script} \
                        --numbat_rscript {numbat_rscript} \
                        --mem {combine_mem} \
                        --walltime {combine_walltime} \
                        --patient {pt} \
                        --samples {','.join(samples)} \
                        --bams {','.join(bams)} \
                        --barcodes {','.join(barcodes)} \
                        --mtxdirs {','.join(mtx_dirs)} \
                        --outdir {outdir_pt}/combined \
                        --trans {trans} \
                        --gamma {gamma} \
                        --min_cells {min_cells} \
                        --min_LLR {min_LLR} \
                        --max_iter {max_iter} \
                        --init_k {init_k} \
                        --max_entropy {max_entropy} \  
                        --cores {cores} \
                        {high_purity}"""
                            
            subprocess.run(cmd, shell=True, check=False)
        