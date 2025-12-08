# Use Python 3.10 slim
FROM python:3.10-slim

# 1. Install System Dependencies
# Added: pkg-config (CRITICAL FIX for rpy2 compilation)
RUN apt-get update && apt-get install -y --no-install-recommends \
    r-base \
    build-essential \
    cmake \
    git \
    pkg-config \
    gfortran \
    libblas-dev \
    liblapack-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    libglpk-dev \
    libgmp-dev \
    libmpfr-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfontconfig1-dev \
    libjpeg-dev \
    libpng-dev \
    libtiff-dev \
    libdeflate-dev \
    liblzma-dev \
    libzstd-dev \
    libbz2-dev \
    libtirpc-dev \
    && rm -rf /var/lib/apt/lists/*

# 2. Set working directory
WORKDIR /app

# 3. Install Python Dependencies
# We upgrade pip/setuptools first to ensure wheel building works for rpy2
COPY requirements.txt .
RUN pip install --no-cache-dir --upgrade pip setuptools wheel && \
    pip install --no-cache-dir -r requirements.txt

# 4. Install R Dependencies (Sequential Order)

# Step A: Bioconductor (ComplexHeatmap)
RUN R -e "install.packages('BiocManager', repos='http://cran.rstudio.com/'); \
          BiocManager::install('ComplexHeatmap', ask=FALSE, update=FALSE)"

# Step B: SeuratObject (Must be installed BEFORE Seurat)
RUN R -e "install.packages('SeuratObject', repos='http://cran.rstudio.com/')"

# Step C: Seurat (The heavy lifter)
RUN R -e "install.packages('Seurat', repos='http://cran.rstudio.com/')"

# Step D: Arrow & Data Wrangling
RUN R -e "Sys.setenv(NOT_CRAN='true'); install.packages('arrow', repos='http://cran.rstudio.com/')"

# Step E: Visualization & Utilities
RUN R -e "install.packages(c( \
    'dplyr', 'tibble', 'jsonlite', 'stringr', 'Matrix', 'tidyr', \
    'ggplot2', 'forcats', 'patchwork', 'reshape2', 'ggh4x', \
    'duckdb', 'glue', 'RColorBrewer', 'circlize', 'ragg', \
    'scales', 'data.table', 'svglite', 'sp', 'cowplot' \
    ), repos='http://cran.rstudio.com/')"

# 5. Copy Application Code
COPY . .

# 6. Run Configuration
EXPOSE 8050
CMD ["gunicorn", "--bind", "0.0.0.0:8050", "app:server", "--workers", "4", "--timeout", "300"]