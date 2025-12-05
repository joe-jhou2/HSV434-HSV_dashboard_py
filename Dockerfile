# Use an official Python runtime as a parent image
FROM python:3.10-slim

# 1. Install System Dependencies
# We added:
# - libhdf5-dev: Required for Seurat (hdf5r)
# - libglpk-dev: Required for Seurat (igraph)
# - libharfbuzz-dev, libfribidi-dev: Required for text shaping (ragg/systemfonts)
# - libfontconfig1-dev: Required for svglite
# - libjpeg-dev, libpng-dev: Image processing
# - cmake: Required for building arrow
RUN apt-get update && apt-get install -y --no-install-recommends \
    r-base \
    build-essential \
    libxml2-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    libcairo2-dev \
    libpango1.0-dev \
    libxt-dev \
    libhdf5-dev \
    libglpk-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfontconfig1-dev \
    libjpeg-dev \
    libpng-dev \
    cmake \
    git \
    && rm -rf /var/lib/apt/lists/*

# 2. Set working directory
WORKDIR /app

# 3. Install Python Dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# 4. Install R Dependencies (Split into steps to prevent memory crashes)

# Step 4A: Install BiocManager and Bioconductor packages (ComplexHeatmap)
RUN R -e "install.packages('BiocManager', repos='http://cran.rstudio.com/'); \
          BiocManager::install('ComplexHeatmap', ask=FALSE, update=FALSE)"

# Step 4B: Install Heavy CRAN Packages (Seurat, Arrow) separately
# We define NOT_CRAN=true for Arrow to allow it to download pre-compiled binaries (much faster)
RUN R -e "Sys.setenv(NOT_CRAN='true'); install.packages(c('arrow', 'Seurat'), repos='http://cran.rstudio.com/')"

# Step 4C: Install Remaining Visualization & Utility Packages
RUN R -e "install.packages(c( \
    'dplyr', 'tibble', 'jsonlite', 'stringr', 'Matrix', 'tidyr', \
    'ggplot2', 'forcats', 'patchwork', 'reshape2', 'ggh4x', \
    'duckdb', 'glue', 'RColorBrewer', 'circlize', 'ragg', \
    'scales', 'data.table', 'svglite', 'SeuratObject', 'sp', 'cowplot' \
    ), repos='http://cran.rstudio.com/')"

# 5. Copy the rest of the application code
COPY . .

# 6. Expose the port
EXPOSE 8050

# 7. Define the command to run the app
CMD ["gunicorn", "--bind", "0.0.0.0:8050", "app:server", "--workers", "4", "--timeout", "300"]