import argparse
import subprocess
import os

# Improvements
## Split apart the genotyping and R into seperate jobs: 
## genotyping is light on memory but faster with more cores
## Numbat.R is more memory hungry, especially for multi-sample patients
## Add logic to skip genotyping if already completed

parser = argparse.ArgumentParser(description='Run the Numbat allele specific scRNA copy number pipeline')

# Required Arguments
parser.add_argument("--numbat_img", help="The numbat image file", type=str)
parser.add_argument("--numbat_refdir", help="Directory of numbat reference files")
parser.add_argument("--pileup_script", help="The numbat preprocessing pileup and phasing Rscript")
parser.add_argument("--numbat_rscript", help="The Rscript to run numbat")
parser.add_argument("--gmap", help="Location of the Eagle genetic map for phasing")
parser.add_argument("--panel_dir", help="Directory of the 1000g reference panel")
parser.add_argument("--snp_vcf", help="SNP vcf")
parser.add_argument("--patient", help="Patient identifier")
parser.add_argument("--samples", help="Patient sample names")
parser.add_argument("--bams", help="comma seperated list of sample bams")
parser.add_argument("--barcodes", help="comma seperated list of sample barcode files")
parser.add_argument("--mtxdirs", help="comma seperated list of matrix folders")
parser.add_argument("--outdir", help="output directory")
parser.add_argument("--UMItag", help="UMItag option for pileup script", default = "Auto")
parser.add_argument("--cellTAG", help="cellTAG option for pileup script", default = "CB")
parser.add_argument("--genome_ver", help="Genome version (hg19, hg38)", default = "hg38")
parser.add_argument("--cores", help="number of cores to use", default=4)
parser.add_argument("--mem", help="amount of memory to use", default=8)
parser.add_argument("--walltime", help="amount of walltime to use", default="48:00")
parser.add_argument("--trans", help="Numbat HMM transmission probability", default=1e-5)
parser.add_argument("--gamma", help="Numbat overdispersion parameter in allele counts", default=20)
parser.add_argument("--min_cells", help="Numbat minimum number of cells for which an pseudobulk HMM will be run", default=20)
parser.add_argument("--min_LLR", help="Numbat minimum log-likelihood ratio threshold to filter CNVs by. ", default=50)


# Parse the arguments
args = parser.parse_args()

# Process args
numbat_img = args.numbat_img
numbat_refdir = args.numbat_refdir
pileup_script = args.pileup_script
numbat_rscript = args.numbat_rscript
gmap = args.gmap
panel_dir = args.panel_dir
snp_vcf = args.snp_vcf
patient = args.patient
samples = args.samples
bams = args.bams
barcodes = args.barcodes
mtxdirs = args.mtxdirs
outdir = args.outdir
UMItag = args.UMItag
cellTAG = args.cellTAG
genome_ver = args.genome_ver
cores = args.cores
mem = args.mem
walltime = args.walltime
trans = args.trans
gamma = args.gamma
min_cells = args.min_cells
min_LLR = args.min_LLR


# Create the preprocessing command
preprocess_cmd = f"""Rscript {pileup_script} \
    --label {patient} \
    --samples {samples} \
    --bams {bams} \
    --barcodes {barcodes} \
    --gmap {gmap} \
    --snpvcf {snp_vcf} \
    --paneldir {panel_dir} \
    --outdir {outdir} \
    --ncores {cores} \
    --UMItag {UMItag} \
    --cellTAG {cellTAG}"""

    
numbat_r_cmd = f"""Rscript {numbat_rscript} \
    --patient {patient} \
    --samples {samples} \
    --mtxdirs {mtxdirs} \
    --outdir {outdir} \
    --genome_ver {genome_ver} \
    --cores {cores} \
    --trans {trans} \
    --gamma {gamma} \
    --min_cells {min_cells} \
    --min_LLR {min_LLR}"""
    
# Open a run_numbat.sh script and print the commands to it
sh_numbat = f"{outdir}/scripts/run_numbat.sh"
os.makedirs(os.path.dirname(sh_numbat), exist_ok=True)
with open(sh_numbat, "w") as f:
    print(preprocess_cmd, file=f)
    print(numbat_r_cmd, file=f)
    f.close()

# Now submit the bash script in the container
bsub_cmd=f"module purge; module load singularity/3.7.1; bsub \
        -J {patient}_numbat \
        -R \"rusage[mem={mem}]\" \
        -R \"select[type==CentOS7]\" \
        -n {cores} \
        -W {walltime} \
        -o {outdir}/out \
        -e {outdir}/err \
        singularity exec \
        --bind /juno:/juno \
        {numbat_img} sh {sh_numbat}"

numbat_bsub = f"{outdir}/scripts/bsub_numbat.sh"

with open(numbat_bsub, "w") as f:
    print(f"{bsub_cmd}", file=f)
    f.close()       

subprocess.run(f"sh {numbat_bsub}", shell=True, check = True)
