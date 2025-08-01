# AiCE: High-fitness mutation prediction tool
<img width="1433" alt="image" src="https://github.com/user-attachments/assets/abd98244-7e87-4ada-87e5-6e9d74387930" />

The repository is an official implementation of [Advancing protein evolution with inverse folding models integrating structural and evolutionary constraints](https://www.cell.com/cell/abstract/S0092-8674(25)00680-4?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867425006804%3Fshowall%3Dtrue).

AiCE is an approach that optimizes protein function by incorporating structural and evolutionary constraints into the process of AI-assisted mutation nomination. It is compatible with widely used protein inverse folding models such as [ProteinMPNN](https://github.com/dauparas/ProteinMPNN), [LigandMPNN](https://github.com/dauparas/LigandMPNN), [ESM-IF1](https://github.com/facebookresearch/esm/tree/main/examples/inverse_folding), [SaProt](https://github.com/westlake-repl/SaProt), and others. A demo for nominating high-fitness (HF) mutations using AiCE-ProteinMPNN is provided in this repository.

<details open><summary><b>Table of contents</b></summary>
  
- [Overview](#overview)
- [Requirements](#requirements)
- [Installation](#installation)
- [Optional Dependencies](#optional-dependencies)
- [Usage](#usage)
  - [1. Single mutation nomination](#1-single-mutation-nomination)
  - [2. LD matrix construction](#2-ld-matrix-construction)
  - [3. SCA matrix construction](#3-sca-matrix-construction)
  - [4. Multi-mutation nomination](#4-multi-mutation-nomination)
- [Citing this work](#citing-this-work)
- [Credits](#credits)
</details>

## Overview
This method nominates mutations based on sampling inverse folding protein sequences. It works with common protein inverse folding models. A demonstration using AiCE-ProteinMPNN for **monomeric** proteins is provided in this repository. For mutation prediction in protein complexes, the script can be modified or customized accordingly, as the underlying principle is similar. Alternatively, the complex structure can be decomposed into multiple substructures, and mutations can be predicted for each separately. This strategy effectively reduces computational demands. For detailed methodology, please refer to our paper.

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
pip3 install -r requirements.txt

# Installing mkdssp

The `mkdssp` program is required for secondary structure annotation in AiCE.

## ✅ Option 1: Use the pre-installed script (if available)

This repository includes a pre-downloaded `mkdssp` binary in the `scripts/` directory. To make it executable, run:

```bash
sudo chmod 755 scripts/mkdssp
# Note: The mkdssp executable provided in the `AiCE/scripts/` directory may fail to run properly on certain systems due to compiler or library version mismatches.
## ✅ Option 2: Compile DSSP manually from source
git clone https://github.com/PDB-REDO/dssp.git
cd dssp
cmake -S . -B build
cmake --build build
cmake --install build
rm scripts/mkdssp
ln -sf $(pwd)/src/mkdssp scripts/mkdssp
## ✅ Option 3: Install via Conda
conda install conda-forge::dssp
After installation, create a symbolic link from the conda-installed binary to the scripts/ directory:
ln -sf ~/anaconda3/envs/AiCE/bin/mkdssp /path/to/AiCE/scripts/mkdssp

# Install plink
wget -c https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20241022.zip
unzip -d scripts/plink/ plink_linux_x86_64_20241022.zip
rm plink_linux_x86_64_20241022.zip
```

## Optional Dependencies

- **Inverse folding models:** 
An inverse folding model is required to output structure-compatible sequences from a given protein structure. For demonstration purposes, we use ProteinMPNN—a lightweight model based on graph neural networks. A pre-deployed version of ProteinMPNN is provided in the `scripts` folder.
- **Secondary structure prediction:** 
The DSSP algorithm is used to predict the protein secondary structure. The repository includes the mkdssp module (version **4.4.7**).
- **Linkage Disequilibrium (LD) calculation:** 
[Plink](https://www.cog-genomics.org/plink/) is used to calculate the LD score. We provide a deployment workflow for plink version **v1.9.0-b.7.7**. Note that plink **v2.0** is not compatible with our workflow by default; you may need to modify `scripts/02.caculated_ld.py` to use plink **v2.0** or later.
- **Evolutionary coupling analysis (SCA) :**
The repository contains a modified version of the pySCA module (originally from [pySCA](https://github.com/reynoldsk/pySCA)) to calculate amino acid evolutionary coupling effects.

## Usage
A demo notebook (`AiCE_demo.ipynb`) is provided for a simple demonstration. Change to the example directory to get started:
```
cd example/
```
The scripts in this repository use relative paths; you may modify them according to your specific requirements.

### 1. Single mutation nomination
Run the following script to nominate single mutations using a protein inverse folding model:
```
bash ../scripts/01.single_mut_prediction.sh <scripts_dir> <input_folder> <beta> <gamma> [output_folder]
```
- `<scripts_dir>`: Directory containing the necessary sub-scripts (by default, the `scripts` folder).
- `<input_folder>`: Folder containing input structure files (`PDB/mmCIF file`). The script automatically searches for these files and outputs the nominated single mutations to `[output_folder]` using the same file prefix as the structure file.
- `<beta>` and `<gamma>`: Screening thresholds for global occurrence and flexible region occurrence, respectively. We recommend **0.8** and **0.5** as general thresholds ("AiCE filtering").
- `[output_folder]`: (**Optional**) Folder for storing output results; the default is `../output`.

Example:
```
bash ../scripts/01.single_mut_prediction.sh ../scripts ./ 0.8 0.5
```
**Note**: Use bash (not sh) to execute the script to avoid unnecessary errors.

An alternative script automatically recommends `<beta>` and `<gamma>` values based on the input structure:
```
bash scripts/01.single_mut_Auto_prediction.sh <input_folder> [output_folder]
```
Additionally, the `scripts/inverse_MPNN.sh` provides a ProteinMPNN-based inverse folding workflow. You can adjust parameters such as **num_seq_per_target** and **sampling_temp** to specify the number of output sequences and the sampling temperature.

### 2. LD matrix construction
Construct the LD matrix based on the inverse folding output sequences:
```
python ../scripts/02.caculated_ld.py <seq_dir> <output_ld_dir>
```
- `<seq_dir>`: Directory containing the inverse folding sequences with a **.fa** extension.
- `<output_ld_dir>`: Directory where the LD matrix files will be saved.

The script automatically searches for **.fa** files in `<seq_dir>`, predicts the LD matrix, and outputs files with the same prefix as the input. Output files include:
- `.ld`: Linkage disequilibrium matrix (derived from pseudo-reverse translated sequences)
- `.vcf`: File recording mutation information

Example:
```
python ../scripts/02.caculated_ld.py ../output/ ../output
```
For more details, please refer to the accompanying manuscript.

### 3. SCA matrix construction
Generate the Statistical Coupling Analysis (SCA) matrix:
```
bash ../scripts/03.caculated_sca.sh <script_dir> <input_dir> <output_dir>
```
- `<input_dir>`: Directory containing the inverse folding sequences with a **.fa** extension.
- `<script_dir>`: Directory containing the sub-scripts for generating the evolutionary coupling matrix.
- `<output_ld_dir>`: Directory where the output files will be stored.

The script automatically searches for **.fa** files in `<input_dir>`, calls the necessary sub-scripts in `<script_dir>`, and outputs the results with the same file prefix as the input. Output files include:
- `.sca_matrix.tsv`: Amino acid evolutionary coupling matrix.
- `.db`: Binary file.

Example:
```
bash ../scripts/03.caculated_sca.sh ../scripts/pySCA/ ../output ../output
```
The folder `../scripts/pySCA/` contains the modified pySCA scripts.

### 4. Multi mutation nomination
Nominate multi-mutations using the LD and SCA matrices:
```
bash ../scripts/04.com_mut_prediction.sh <script_dir> <input_dir> <number-or-list> <output_dir>
```
•	The script automatically searches `<input_dir>` for the **.fa**, **.sca_matrix.tsv**, **.ld**, **.comb**, and **.vcf** files produced in steps 1, 2, and 3.

•	It outputs the SCA and LD scores for multi-mutations to files ending in **.sca.result** and **.ld.result**, respectively. The output file prefix will match the corresponding input file.

- `<number-or-list>`: If a number is provided, the script will iterate through all mutation types of that order；If a list is provided (e.g., "1 3 5"), only the scores for the specified mutation combinations will be output.
  	
	•	Each line in the list represents a mutation combination.
	
 	•	Positions are space-separated using 1-based indexing.
 
 	•	Lines that are empty or contain unparseable non-integer characters will be skipped or marked as invalid (output as `None	NaN`).

The output file format is as follows:
```
Mutation Type	Mean Pairwise score: score	Log Mean Pairwise score: score	Logical Flag (0/1)
```
The fourth column indicates whether the mutation combination is recommended based on the screening thresholds (1 for recommended, 0 for not recommended). By default, the screening thresholds are set to 0.5 for LD scores and 0.9 for SCA score percentiles. You can customize these thresholds by providing an additional `-t` parameter (see **lines 106–124** in `scripts/04.com_mut_prediction.sh`) directly.

Example:
```
bash ../scripts/04.com_mut_prediction.sh ../scripts ../output/ 2 ../output
```

This command will iterate over all double mutation combinations and nominate HF mutations.

## Citing this work
If you use this code, please cite:
```
@article{Li2025AiCE,
  title={Harnessing structural and evolutionary constraints to enhance protein evolution using inverse folding models},
  author={Fei, Hongyuan and Li, Yunjia and Liu, Yijing and Wei, Jingjing and Chen, Aojie and Gao, Caixia},
  journal={Cell},
  year={2025},
  publisher={Cell Press}
```
## Credits
This repository incorporates code from:

- **[ProteinMPNN](https://github.com/dauparas/ProteinMPNN)**
- **[Plink](https://www.cog-genomics.org/plink/)**
- **[pySCA](https://github.com/reynoldsk/pySCA)**

## Patent Notice

This project is licensed under the MIT License. Please note that portions of this project are covered by one or more patents. 

Using this software for research or educational purposes is permitted under the MIT License.  
However, **commercial use of the patented invention requires explicit permission from the patent holder**.  
For licensing inquiries, please contact [cxgao@genetics.ac.cn].

By using this software, you acknowledge and agree to comply with all applicable patent laws and regulations.
