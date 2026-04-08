# EchoPipe
Tutorial on how to use the EchoPipe script.
The scripts are command-line based, written in python and uses conda.
Conda environment, and installations are found below.

For guidance of how to use the script, follow the Example, or read the manual.

## Flowchart of the workflow of EchoPipe:

<img width="982" height="1024" alt="Figure_1_workflow" src="https://github.com/user-attachments/assets/f944a13a-1524-4534-bfc7-1915753dec80" />  
  
Flowchart of the EchoPipe's workflow. Each color represents a module of the workflow and the corresponding script. Figure 1. EchoPipe modular workflow for iterative database creation and curation. The pipeline architecture is divided into four primary stages: Reference Template Generation (light yellow), Sequence Retrieval and Database Creation (light orange), Diagnostic Curation (orange), and Database Completion/Evaluation (red). The iterative feedback loop (bottom center) enables users to continuously update curated databases with novel accessions from public repositories without duplicating previous computational efforts.  



# Dependencies and installation
The script requires miniconda, and works on both Ubuntu (WSL) and OS.

Miniconda
https://docs.anaconda.com/free/miniconda/index.html

```
conda install conda=24.5.0
```

## Ubuntu (WSL)
Download the miniconda environment EchoPipe, environment_240930.yml, following instruction underneath.

```
wget https://raw.githubusercontent.com/EivindStensrud/EchoPipe/refs/heads/main/environment_240930.yml
conda env create -f environment_240930.yml

conda activate EchoPipe

```


## OS and Windows
Download the miniconda environment EchoPipe, environment_OS.yml, following instruction underneath for Ubuntu.

```
wget https://raw.githubusercontent.com/EivindStensrud/EchoPipe/refs/heads/main/environment_240930_OS.yml
conda env create -f environment_240930_OS.yml

conda activate EchoPipe

```

# Using EchoPipe
Follow the tutorial where the amphibian batra primers and Scandinavian amphibian species group was used.  
https://github.com/EivindStensrud/EchoPipe/blob/main/Tutorial_amphibians.md  
