# HSV Genital Skin Biopsy scRNAseq Explorer

**By: Joe Hou**

![Workflow](www/images/Work_Flow_ver2.png)

## Overview

This Shiny application provides an interactive interface for exploring 10X scRNAseq data derived from HSV+ skin biopsy samples.
[**https://hsvdashboard.joehou.net/**](https://hsvdashboard.joehou.net/).\
As the manuscript is currently under review for publication, the raw and processed data remain undisclosed.

### About the cohort:
This study involved healthy adults with confirmed herpes simplex virus type 2 (HSV-2) seropositivity, verified via Western blot analysis. All participants tested negative for HIV. In this longitudinal HSV study, we enrolled 17 participants and collected genital skin biopsies at three distinct time points to capture HSV shedding and healing phases.

### About the time points:
Samples were collected at three key stages: the first, labeled "Prior," was collected before any visible lesions appeared. When lesions became visible, samples were labeled as "Lesion." Finally, samples taken 8 weeks after lesion appearance are referred to as "Post."

### About the dataset:
This release includes datasets with refined mapping for T cells, Myeloid cells and Visium data. Users can freely select from the available datasets for exploration. We recommend using the "CellType_Level3" annotation for visualization, as it provides the most detailed and granular clustering. While "CellType_Level1" offers a broader categorization, "CellType_Level3" represents the final, detailed clustering, with cell types identified at a higher level of specificity. Full definitions of each cell cluster and cell type will be provided in the manuscript upon publication.

## Features and Functionalities

Two main tabs are provided for **scRNAseq** and **spatial Visium** data.

For **scRNAseq** data:

1. **Cell Cluster Composition**:
   - **Goal:** Analyze the composition of cell clusters across time points and individuals, with some clusters showing specificity to time points and variability across participants.
   
2. **Gene Expression Patterns**:
   - **Goal:** Investigate gene expression patterns across different cell clusters, including gene distribution and expression variability across time points.

For **Visium** data:

1. **Gene Expression Patterns**:
   - **Goal:** Investigate gene expression distribution and intensity across tissue.
   
2. **Cell Type De-convolution**:
   - **Goal:** To further understand cell type composition per spatial location.

### Interactive Exploration on Cell Type (Cluster Discovery Tab)
- **Cell Type and Status Identification**: Users can view all available cell types and subjects, with the flexibility to select specific ones for detailed exploration.
- **Dynamic UMAP Visualization**: Initial displays show holistic UMAP results, which can be refined by time point. Selecting different subjects and cell types dynamically updates the UMAP to show specific clusters.
- **Cell Type Statistics**: Displays percentages and counts of cell types for selected or all subjects.
- **Cell Type by Subject**: Visualizes the percentage of each cell type within a given sample using stacked bar plots to understand the heterogeneity of cell type across samples.

### Interactive Exploration on Gene Expression (Gene Discovery Tab)
- **Feature Exploration**: Allows for the selection of cell types, subjects, and specific genes of interest.
- **At Single Cell Level**:
  - **Feature Gene Highlight**: The expression of selected gene projected to UMAP overlay with cell type annotation, vividly display gene expression pattern across cell types.
  - **Heatmap**: Shows selected gene expression intensity per cell, grouped by cell types, statuses, and subject, aiding in providing gene details per sample.
- **At Cell Type Level**:
  - **Violin Plots**: The selected gene expression was summarized into cell type level and shows average expression levels across cell types and statuses, aiding in understanding temporal dynamics.
  - **Feature Percentage Plots**: Summarizes gene expression at the sample/subject level, indicating the proportion of cells expressing a specific gene within a cell type.
  - **Dot Plots**: Dot size indicates the percentage of cells expressing each gene, while color intensity represents gene expression levels.

## Gallery
Here are some visual showcases of the application features:

|  |  |
|------------------------------------|------------------------------------|
| ![Demo 1](assets/images/demo1.png) | ![Demo 3](assets/images/demo3.png) |
| ![Demo 5](assets/images/demo5.png) | ![Demo 7](assets/images/demo7.png) |

### Version Release
Dec-5-2025: release ver 1.0  
