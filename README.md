# AiCE: High-fitness mutation prediction tool
<img width="1433" alt="image" src="https://github.com/user-attachments/assets/abd98244-7e87-4ada-87e5-6e9d74387930" />
The repository is an official implementation of Harnessing structural and evolutionary constraints to enhance protein evolution using inverse folding models.

AiCE is an approach that optimizes protein function by incorporating structural and evolutionary constraints into the process of AI-assisted mutation nomination. It is compatible with widely used protein inverse folding models such as [ProteinMPNN](https://github.com/dauparas/ProteinMPNN), [LigandMPNN](https://github.com/dauparas/LigandMPNN), [ESM-IF1](https://github.com/facebookresearch/esm/tree/main/examples/inverse_folding), [SaProt](https://github.com/westlake-repl/SaProt), and others. A demo for nominating high-fitness mutations using AiCE-ProteinMPNN is provided in this repository.

<details open><summary><b>Table of contents</b></summary>
  
- [Overview](#overview)
- [Requirements](#requirements)
- [Installation](#installation)
- [Optional Dependencies](#optional-dependencies)
- [Usage](#usage)
  - [1. Single Mutation Prediction](#1-single-mutation-prediction)
  - [2. LD Matrix Construction](#2-ld-matrix-construction)
  - [3. SCA Matrix Construction](#3-sca-matrix-construction)
  - [4. Multi-mutation Prediction](#4-multi-mutation-prediction)
- [Citing This Work](#citing-this-work)
- [Credits](#credits)
</details>

## Overview
This method nominates mutations based on sampling inverse folding protein sequences. It works with common protein inverse folding models including ProteinMPNN, LigandMPNN, ESM-IF1, and SaProt. A demo based on AiCE-ProteinMPNN for nominating high-fitness mutations is provided in this repository.

## Requirements
To run AiCE, you will need:
- **Python:** Version ≥ 3.8
- **Libraries:** PyTorch, NumPy, SciPy, and Pandas
- **Biological Sequence Analysis:** Biopython
- **PDB File Handling:** ProDy

Dependencies can be installed directly using the provided `requirements.txt` file.

## Installation
Clone the repository and set up your environment:
```bash
# Clone the AiCE repository
git clone https://github.com/ScorpioLea/AiCE
cd AiCE

# Setup your conda environment
conda create -n AiCE python=3.11
conda activate AiCE
pip install -r requirements.txt

# Install plink
wget -c https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20241022.zip
unzip -d scripts/plink/ plink_linux_x86_64_20241022.zip
```

## Optional Dependencies

- **Inverse Folding Model:** 
An inverse folding model is required to output structure-compatible sequences from a given protein structure. For demonstration purposes, we use ProteinMPNN—a lightweight model based on graph neural networks. A pre-deployed version of ProteinMPNN is provided in the scripts folder.
- **Secondary Structure Prediction:** 
The DSSP algorithm is used to predict the protein secondary structure. The repository includes the mkdssp module (version 4.4.7).
- **Linkage Disequilibrium (LD) Calculation:** 
Plink is used to calculate the LD score. We provide a deployment workflow for plink version v1.9.0-b.7.7. Note that plink v2.0 is not compatible with our workflow by default; you may need to modify scripts/02.caculated_ld.py to use plink v2.0 or later.
- **Evolutionary Coupling Analysis:**
The repository contains a modified version of the pySCA module (originally from pySCA) to calculate amino acid evolutionary coupling effects.

## Usage
A demo notebook (AiCE_demo.ipynb) is provided for a simple demonstration. Change to the example directory to get started:
```
cd example/
```
The scripts in this repository use relative paths; you may modify them according to your specific requirements.

###1-single-mutation-prediction
Run the following script to nominate single mutations using a protein inverse folding model:
```
bash ../scripts/01.single_mut_prediction.sh <scripts_dir> <input_folder> <beta> <gamma> [output_folder]
```
- ** <scripts_dir>**: Directory containing the necessary sub-scripts (by default, the scripts folder).
	•	<input_folder>: Folder containing input structure files (PDB/mmCIF). The script automatically searches for these files and outputs the nominated single mutations to [output_folder] using the same file prefix as the structure file.
	•	<beta> and <gamma>: Screening thresholds for global occurrence and flexible region occurrence, respectively. We recommend 0.8 and 0.5 as general thresholds (AiCE filtering).
	•	[output_folder]: (Optional) Folder for storing output results; the default is ../output.
