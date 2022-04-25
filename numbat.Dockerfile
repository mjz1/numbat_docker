FROM rocker/r-ver:4.1.2

# Install base utilities
RUN apt-get update && \
    apt-get install -y build-essential && \
    apt-get install -y wget && \
    apt-get install -y libz-dev && \
    apt-get install -y libbz2-dev && \
    apt-get install -y liblzma-dev && \
    apt-get install -y libcurl4-openssl-dev && \
    apt-get install -y libxml2-dev && \
    apt-get install -y libglpk-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# install conda
ENV CONDA_DIR /opt/conda

RUN wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda

ENV PATH=$CONDA_DIR/bin:$PATH

RUN conda config --add channels conda-forge
RUN conda config --add channels bioconda

# Install cellsnp-lite
RUN conda install -y cellsnp-lite

# install eagle
RUN mkdir src && \
    cd src && \
    wget -q https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/Eagle_v2.4.1.tar.gz && \
    tar -xvzf Eagle_v2.4.1.tar.gz

ENV PATH=/src/Eagle_v2.4.1/:$PATH

# Get 1000g SNP reference sets
# We opt not to include these in the image to keep size down
# RUN mkdir ref && cd ref
# RUN wget -q https://sourceforge.net/projects/cellsnp/files/SNPlist/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz
# RUN wget -q https://sourceforge.net/projects/cellsnp/files/SNPlist/genome1K.phase3.SNP_AF5e2.chr1toX.hg19.vcf.gz
# Get 1000g reference panel (phasing?)
# RUN wget -q http://pklab.med.harvard.edu/teng/data/1000G_hg38.zip
# RUN wget -q http://pklab.med.harvard.edu/teng/data/1000G_hg19.zip

# Install necessary R packages
RUN Rscript -e 'install.packages("devtools", quiet = TRUE)'
RUN Rscript -e 'install.packages("BiocManager", quiet = TRUE)'
RUN Rscript -e 'BiocManager::install()'
RUN Rscript -e 'install.packages("Seurat")'

# install numbat
RUN Rscript -e 'devtools::install_github("https://github.com/kharchenkolab/numbat", quiet = TRUE)'