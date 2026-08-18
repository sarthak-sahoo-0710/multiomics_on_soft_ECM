# Multi-omics investigation reveals molecular determinants of cancer cell evolution on soft extracellular matrix

## Overview
This repository contains the in-silico mechanism-based model associated with the manuscript **"Multi-omics investigation reveals molecular determinants of cancer cell evolution on soft extracellular matrix."** 

Cancer cell adaptation to physical tumor microenvironments, specifically soft extracellular matrices (ECM), imposes a selection pressure leading to the enrichment of high-fitness genetic variants. While ancestral populations exhibit a stressed cell state (cell cycle arrest, distinct metabolic shifts) upon short-duration exposure to soft ECM, sustained culture selects for a robust proliferative phenotype characterized by a silenced stress response, reduced chromatin accessibility, and *de novo* DNA methylation (including *CDH1*/E-cadherin promoter hypermethylation).

The mathematical model provided here demonstrates how these molecular differences—specifically elevated MYBL2 and FAK levels, alongside high YAP1 nuclear localization—serve as salient features of genetic clones capable of FAK upregulation and subsequent survival in mechanically soft microenvironments.

## Repository Contents
* **Mathematical Model**: Scripts and source files for the in-silico mechanism-based model evaluating FAK upregulation, YAP1 nuclear localization, and MYBL2 dynamics. *(Note: Omics datasets such as RNA-seq, ATAC-seq, and RRBS-seq are hosted in their respective public repositories as detailed in the manuscript).*

## Usage
*(Provide instructions on how to run the mathematical model, required dependencies, and expected outputs here.)*

```bash
# Example standard execution
python run_model.py
```

## Authors and Affiliations
* **Sarthak Sahoo**<sup>1,#</sup>
* **Sejal Khanna**<sup>1,#</sup>
* **Ting-Ching Wang**<sup>2</sup>
* **Tanmay P. Lele**<sup>2,3,4,*</sup>
* **Mohit Kumar Jolly**<sup>1,*</sup>

<sup>1</sup> Department of Bioengineering, Indian Institute of Science, Bengaluru, Karnataka, India – 560012<br>
<sup>2</sup> Artie McFerrin Department of Chemical Engineering, Texas A&M University, College Station, TX 77843<br>
<sup>3</sup> Department of Biomedical Engineering, Texas A&M University, College Station, TX 77843<br>
<sup>4</sup> Department of Translational Medical Sciences, Texas A&M University, Houston, TX 77030<br>
<sup>#</sup> Equal contribution<br>
<sup>*</sup> Co-corresponding authors: mkjolly@iisc.ac.in (M.K.J); tanmay.lele@tamu.edu (T.P.L)

## Citation
If you use this model in your research, please cite the associated publication:

```bibtex
@article{sahoo202Xmultiomics,
  title={Multi-omics investigation reveals molecular determinants of cancer cell evolution on soft extracellular matrix},
  author={Sahoo, Sarthak and Khanna, Sejal and Wang, Ting-Ching and Lele, Tanmay P. and Jolly, Mohit Kumar},
  journal={TBD},
  year={202X}
}
```
