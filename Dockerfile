FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="04053bac5affd345f6138e6c34a9f6dd2d39abc016a736ac97cd51a3577a074f"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: workflow/envs/iqtree_env.yaml
#   prefix: /conda-envs/79d5f95a98f35db6ea804c88387f5fdb
#   channels:
#     - bioconda
#     - conda-forge
#   dependencies:
#     - conda-forge::libgcc-ng>=9.4.0
#     - bioconda::iqtree=2.2.0_beta
RUN mkdir -p /conda-envs/79d5f95a98f35db6ea804c88387f5fdb
COPY workflow/envs/iqtree_env.yaml /conda-envs/79d5f95a98f35db6ea804c88387f5fdb/environment.yaml

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

RUN mamba env create --prefix /conda-envs/79d5f95a98f35db6ea804c88387f5fdb --file /conda-envs/79d5f95a98f35db6ea804c88387f5fdb/environment.yaml && \
    mamba env create --prefix /conda-envs/aa39c1a9de06616bc46a5eccb739bf0e --file /conda-envs/aa39c1a9de06616bc46a5eccb739bf0e/environment.yaml && \
    mamba env create --prefix /conda-envs/02373219334870eb0ad3f91206043108 --file /conda-envs/02373219334870eb0ad3f91206043108/environment.yaml && \
    mamba env create --prefix /conda-envs/202009c5e81d3324de81c25721b22ede --file /conda-envs/202009c5e81d3324de81c25721b22ede/environment.yaml && \
    mamba clean --all -y
