FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="25846baaf1384172ae96e327c06051c49a227f9b2dceba7181c58c2f33a8e121"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: workflow/envs/iqtree_env.yaml
#   prefix: /conda-envs/275e768c6f61d7d02adfbcf9c09fef1d
#   channels:
#     - bioconda
#   dependencies:
#     - bioconda::iqtree=2.2.0_beta
RUN mkdir -p /conda-envs/275e768c6f61d7d02adfbcf9c09fef1d
COPY workflow/envs/iqtree_env.yaml /conda-envs/275e768c6f61d7d02adfbcf9c09fef1d/environment.yaml

# Conda environment:
#   source: workflow/envs/muscle_env.yaml
#   prefix: /conda-envs/aa39c1a9de06616bc46a5eccb739bf0e
#   channels:
#     - bioconda
#   dependencies:
#     - bioconda::muscle=3.8.1551
RUN mkdir -p /conda-envs/aa39c1a9de06616bc46a5eccb739bf0e
COPY workflow/envs/muscle_env.yaml /conda-envs/aa39c1a9de06616bc46a5eccb739bf0e/environment.yaml

# Conda environment:
#   source: workflow/envs/openjdk_env.yaml
#   prefix: /conda-envs/02373219334870eb0ad3f91206043108
#   channels:
#     - conda-forge
#   dependencies:
#     - conda-forge::openjdk=10.0.2
RUN mkdir -p /conda-envs/02373219334870eb0ad3f91206043108
COPY workflow/envs/openjdk_env.yaml /conda-envs/02373219334870eb0ad3f91206043108/environment.yaml

# Conda environment:
#   source: workflow/envs/raxml_env.yaml
#   prefix: /conda-envs/202009c5e81d3324de81c25721b22ede
#   channels:
#     - bioconda
#   dependencies:
#     - bioconda::raxml=8.2.12
RUN mkdir -p /conda-envs/202009c5e81d3324de81c25721b22ede
COPY workflow/envs/raxml_env.yaml /conda-envs/202009c5e81d3324de81c25721b22ede/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/275e768c6f61d7d02adfbcf9c09fef1d --file /conda-envs/275e768c6f61d7d02adfbcf9c09fef1d/environment.yaml && \
    mamba env create --prefix /conda-envs/aa39c1a9de06616bc46a5eccb739bf0e --file /conda-envs/aa39c1a9de06616bc46a5eccb739bf0e/environment.yaml && \
    mamba env create --prefix /conda-envs/02373219334870eb0ad3f91206043108 --file /conda-envs/02373219334870eb0ad3f91206043108/environment.yaml && \
    mamba env create --prefix /conda-envs/202009c5e81d3324de81c25721b22ede --file /conda-envs/202009c5e81d3324de81c25721b22ede/environment.yaml && \
    mamba clean --all -y
