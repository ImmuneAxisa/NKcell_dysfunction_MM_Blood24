#!/bin/bash

###############################################################
######################### set up ##############################
###############################################################


# command line arguments
# column name for the metadata table to be used in MELD
sample_labels=$1       # eg "treatment"
cluster=$2             # eg "cluster_label"
out_suffix=$3          # basename for tha output tables
input_csv_file=$4      # input files basenames (has to finish in _pca.csv and _meta.csv)


# Raw data and working directories:
base_dir="PATH/TO/YOUR/PROJECT/DIR/"
data_dir="${base_dir}/data"
input_dir="${base_dir}/results/NK_analysis/phate_meld_input_csv"

# Output directories:
slurm_dir="${base_dir}/slurm"
logs_dir="${base_dir}/logs"
doc_dir="${base_dir}/doc"
results_dir="${base_dir}/results/NK_analysis/phate_meld_output"


# Create any missing directories
directories=( $data_dir $logs_dir $doc_dir $slurm_dir \
  $input_dir $results_dir)
for directory in ${directories[@]}; do
  if [ ! -d $directory ]; then
    echo "Creating directory ${directory}"
    mkdir -p $directory
  fi
done

###############################################################
######################### sbatch ##############################
###############################################################

# Sample table prefixes: 
# structure must be <prefix>_pca.csv and <prefix>_meta.csv
# for both files the first column must be the cell IDs eg barcodes

# Array of samples with prefixes to be processed
#samples=(CD4 CD8 CD8NK B M)
samples=$input_csv_file
#samples=(B)




for sample in ${samples[@]} #array is deprecated but still functional with one argument
do
    # format inputs and check they exists
    pca=${input_dir}/${sample}_pca.csv
    meta=${input_dir}/${sample}_meta.csv
    
    if [ ! -f "$pca" ]; then
        echo "$pca not found"
        exit
    fi
    
    if [ ! -f "$meta" ]; then
        echo "$meta not found"
        exit
    fi
    # modify name for output if wanted
    sample=${sample}_${out_suffix}
    
    echo "processing $sample"
    
    #slurm scripts
    phate_script=${slurm_dir}/${sample}_phate.slurm
    meld_vfc_script=${slurm_dir}/${sample}_meld_vfc.slurm
  
  
############ PHATE ##############

    echo \
"#!/bin/bash
#SBATCH --partition=general
#SBATCH --job-name=${sample}.phate
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=119gb
#SBATCH --time=48:00:00
#SBATCH -o ${logs_dir}/${sample}.phate.out
#SBATCH -e ${logs_dir}/${sample}.phate.err

module purge
module load miniconda

source activate MELD_env

python -u ~/project/MM/src/phate_meld/PHATE_v2.py \\
    -i ${pca} \\
    -o ${results_dir} \\
    -p ${sample}


source deactivate

sstat \${SLURM_JOBID}.batch


" > $phate_script

#sbatch $phate_script


############ MELD ##############

    echo \
"#!/bin/bash
#SBATCH --partition=general
#SBATCH --job-name=${sample}.meld_vfc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=119gb
#SBATCH --time=48:00:00
#SBATCH -o ${logs_dir}/${sample}.meld_vfc.out
#SBATCH -e ${logs_dir}/${sample}.meld_vfc.err

module purge
module load miniconda

#conda init

source activate MELD_env

python -u ~/project/MM/src/phate_meld/MELD_v2.1.py \\
    -i ${pca} \\
    -m ${meta} \\
    -l ${sample_labels} \\
    -c ${cluster} \\
    -o ${results_dir} \\
    -p ${sample}

source deactivate

sstat \${SLURM_JOBID}.batch

" > $meld_vfc_script

#sbatch $meld_vfc_script


done