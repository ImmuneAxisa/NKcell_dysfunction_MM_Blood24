This set of scripts extract relevant embeddings from R Seurat objects, saves it to csv files to run PHATE and MELD. PHATE and MELD were run on an HPC cluster with slurm scheduler.

# Prerequisites

* install [MELD](https://github.com/KrishnaswamyLab/MELD) (and [phate](https://github.com/KrishnaswamyLab/PHATE)) in a conda environment
* Adapt the path in the `phate_meld_sbatcher_w_args_v2.sh` script to your computing environment, notably if working on an HPC cluster

# Scripts

* `phate_meld_wrapper_v2_NK_analysis.R` loads the Seurat objects, saves embeddings as csv and calls `phate_meld_sbatcher_w_args_v2.sh`.
* `phate_meld_sbatcher_w_args_v2.sh` generates the scripts that will submitted to the slurms scheduler. Those will call the python scripts to run MELD and PHATE.
* `PHATE_v2.py` and `MELD_v2.1.py` will run the analyses and output csv files with PHATE embeddings and MELD likelihood and density scores.