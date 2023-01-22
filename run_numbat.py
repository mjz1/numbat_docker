import argparse
import subprocess
import os

parser = argparse.ArgumentParser(description='Run the Numbat allele specific scRNA copy number pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Required Arguments
parser.add_argument("--numbat_img", help="The numbat image file", default=argparse.SUPPRESS)
parser.add_argument("--pileup_script", help="The numbat preprocessing pileup and phasing Rscript", default="/numbat/inst/bin/pileup_and_phase.R")
parser.add_argument("--numbat_rscript", help="The Rscript to run numbat", default=argparse.SUPPRESS)
parser.add_argument("--gmap", help="Location of the Eagle genetic map for phasing", default="/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz")
parser.add_argument("--panel_dir", help="Directory of the 1000g reference panel", default="/data/1000G_hg38/")
parser.add_argument("--snp_vcf", help="SNP vcf", default="/data/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf")
parser.add_argument("--patient", help="Patient name. ", default=argparse.SUPPRESS)
parser.add_argument("--samples", help="Comma seperated list of samples", default=argparse.SUPPRESS)
parser.add_argument("--bams", help="Comma seperated list of sample bams.", default=argparse.SUPPRESS)
parser.add_argument("--barcodes", help="Comma seperated list of sample barcode files", default=argparse.SUPPRESS)
parser.add_argument("--mtxdirs", help="Comma seperated list of matrix folders (filtered_feature_bc_matrix folder from 10X cellranger)", default=argparse.SUPPRESS)
parser.add_argument("--outdir", help="output directory", default=argparse.SUPPRESS)
parser.add_argument("--UMItag", help="UMItag option for pileup script", default = "Auto")
parser.add_argument("--cellTAG", help="cellTAG option for pileup script", default = "CB")
parser.add_argument("--genome_ver", help="Genome version (hg19, hg38)", default = "hg38", choices=["hg19","hg38"])
parser.add_argument("--cores", help="number of cores to use", default=4, type=int)
parser.add_argument("--mem", help="amount of memory to use", default=16, type=int)
parser.add_argument("--walltime", help="amount of walltime to use", default="48:00")
parser.add_argument("--trans", help="Numbat HMM transmission probability", default=1e-5, type=float)
parser.add_argument("--gamma", help="Numbat overdispersion parameter in allele counts", default=20)
parser.add_argument("--min_cells", help="Numbat minimum number of cells for which an pseudobulk HMM will be run", default=50)
parser.add_argument("--min_LLR", help="Numbat minimum log-likelihood ratio threshold to filter CNVs by. ", default=5)
parser.add_argument("--init_k", help="Number of clusters in the initial clustering", default=10, type=int)
parser.add_argument("--max_iter", help="Maximum number of iterations to run the phyologeny optimization", default=2, type=int)
parser.add_argument("--max_entropy", help="Entropy threshold to filter CNVs", default=0.5)
parser.add_argument("--multi_allelic", help="Flag to enable multi-allelic calling", action='store_true', default=True)
parser.add_argument("--high_purity", help="Flag to detect and exclude regions of clonal deletions/LOH before running Numbat. Recommended for cell line data or high-purity tumors", action='store_true', default=False)


# Parse the arguments
args = parser.parse_args()

# Process args
numbat_img = args.numbat_img
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
init_k = args.init_k
max_iter = args.max_iter
max_entropy = args.max_entropy
multi_allelic = args.multi_allelic
high_purity = args.high_purity


pileup_ncores = 8

# Set multi_allelic flag
if multi_allelic:
    multi_allelic="--multi_allelic"
else:
    multi_allelic=""
    
# Set high purity flag
if high_purity:
    high_purity="--high_purity"
else:
    high_purity=""

# Create the preprocessing command
pileup_cmd = f"""Rscript {pileup_script} \
    --label {patient} \
    --samples {samples} \
    --bams {bams} \
    --barcodes {barcodes} \
    --gmap {gmap} \
    --snpvcf {snp_vcf} \
    --paneldir {panel_dir} \
    --outdir {outdir} \
    --ncores {pileup_ncores} \
    --UMItag {UMItag} \
    --cellTAG {cellTAG}"""
    

# Run numbat command    
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
    --min_LLR {min_LLR} \
    --max_iter {max_iter} \
    --init_k {init_k} \
    --max_entropy {max_entropy} \
    {multi_allelic} \
    {high_purity}"""
    

# Create bash runscripts script and print the commands to them
sh_pileup = f"{outdir}/scripts/run_pileup.sh"
sh_numbat = f"{outdir}/scripts/run_numbat.sh"

# Create directory for scripts
os.makedirs(os.path.dirname(sh_pileup), exist_ok=True)

# Print pileup command
with open(sh_pileup, "w") as f:
    print(pileup_cmd, file=f)
    f.close()

# Print run numbat command
with open(sh_numbat, "w") as f:
    print(numbat_r_cmd, file=f)
    f.close()

# Pileup bsub commands
bsub_pileup=f"module purge; module load singularity/3.7.1; bsub \
        -J {patient}_numbat_pileups \
        -R \"rusage[mem=4]\" \
        -R \"select[type==CentOS7]\" \
        -n {pileup_ncores} \
        -W 24:00 \
        -o {outdir}/pileups.out \
        -e {outdir}/pileups.err \
        singularity exec \
        --no-home \
        --bind /juno:/juno \
        --bind /work:/work \
        {numbat_img} sh {sh_pileup}"

bsub_pileup_sh = f"{outdir}/scripts/bsub_pileups.sh"
with open(bsub_pileup_sh, "w") as f:
    print(f"{bsub_pileup}", file=f)
    f.close()      
        
bsub_numbat=f"module purge; module load singularity/3.7.1; bsub \
        -J {patient}_numbat_r \
        -w 'done({patient}_numbat_pileups)' \
        -R \"rusage[mem={mem}]\" \
        -R \"select[type==CentOS7]\" \
        -n {cores} \
        -W {walltime} \
        -o {outdir}/numbat.out \
        -e {outdir}/numbat.err \
        singularity exec \
        --no-home \
        --bind /juno:/juno \
        --bind /work:/work \
        {numbat_img} sh {sh_numbat}"

bsub_numbat_sh = f"{outdir}/scripts/bsub_numbat.sh"
with open(bsub_numbat_sh, "w") as f:
    print(f"{bsub_numbat}", file=f)
    f.close() 

# TODO: Add output paths
final_sample = samples.split(",")[-1]
pileup_out = f"{outdir}/{final_sample}_allele_counts.tsv.gz"

numbat_out = f"{outdir}/numbat_out.rda"


# Submit commands depending on existing outputs
if os.path.exists(numbat_out):
    print(f"Numbat already completed for {patient}. Skipping...")
elif os.path.exists(pileup_out):
    # Need to remove the job dependency
    # This can be done better
    #####
    bsub_numbat=f"module purge; module load singularity/3.7.1; bsub \
        -J {patient}_numbat_r \
        -R \"rusage[mem={mem}]\" \
        -R \"select[type==CentOS7]\" \
        -n {cores} \
        -W {walltime} \
        -o {outdir}/numbat.out \
        -e {outdir}/numbat.err \
        singularity exec \
        --no-home \
        --bind /juno:/juno \
        --bind /work:/work \
        {numbat_img} sh {sh_numbat}"

    bsub_numbat_sh = f"{outdir}/scripts/bsub_numbat.sh"
    with open(bsub_numbat_sh, "w") as f:
        print(f"{bsub_numbat}", file=f)
        f.close() 
    #####
    
    print(f"Running numbat for {patient}")
    subprocess.run(f"sh {bsub_numbat_sh}", shell=True, check = True)
else:
    # Neither pilups or numbat completed.
    print(f"Running pileups and numbat for {patient}")
    subprocess.run(f"sh {bsub_pileup_sh}", shell=True, check = True)
    subprocess.run(f"sh {bsub_numbat_sh}", shell=True, check = True)
