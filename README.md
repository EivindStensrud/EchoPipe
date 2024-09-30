# EchoPipe
Tutorial on how to use the EchoPipe script.
The scripts are command-line based, written in python and uses conda.
Conda environment, and installations are found below.

For guidance of how to use the script, follow the Example, or read the manual.

## Flowchart of the workflow of EchoPipe:


![Workflow_echopipe](https://github.com/EivindStensrud/EchoPipe/assets/83813403/e4a0ca27-31a7-48c9-95ab-b05cfbc6d086)




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
wget https://raw.githubusercontent.com/EivindStensrud/EchoPipe/main/environment.yml
conda env create -f environment_240930.yml

conda activate EchoPipe

```


## OS and Windows
Download the miniconda environment EchoPipe, environment_OS.yml, following instruction underneath for Ubuntu.

```
wget https://raw.githubusercontent.com/EivindStensrud/EchoPipe/main/environment_OS.yml
conda env create -f environment_OS.yml

conda activate EchoPipe

```
