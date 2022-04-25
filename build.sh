# Build docker image
docker build -t numbat:0.1.1 --platform linux/x86_64 -f numbat.Dockerfile .

mkdir singularity
output_path=$(pwd)/singularity

# Convert to singularity image
docker run --platform linux/x86_64 \
    --volume=/var/run/docker.sock:/var/run/docker.sock \
    --volume=/Users/mzatzman/repos/numbat_docker/singularity:/output \
    --privileged --tty --rm quay.io/singularity/docker2singularity \
    numbat:0.1.1


rsync -ar --progress ~/repos/numbat_docker/singularity/numbat_*.sif zatzmanm@juno.mskcc.org:/home/zatzmanm/work/images/numbat