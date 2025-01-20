import loompy
import os

data_dir = '/mnt/SERVER-CRCT-STORAGE/CRCT13/Ruchan/01_NK_project/99_back_up/cellranger_counts_human/'

files = [os.path.join(data_dir, f'{sample}/velocyto/{sample}.loom')
         for sample in ['H21b', 'H21m', 'H22b', 'H22m', 'H23b', 'H23m', 'H24b', 'H24m',
                        'H25b', 'H25m', 'H26b', 'H26m', 'H27b', 'H27m', 'H28b', 'H28m',
                        'H29b', 'H29m', 'H30b', 'H30m', 'T22467b', 'T22467m',
                        'T23491b_NK', 'T23491b_T', 'T23491m_NK', 'T23491m_T']]

loompy.combine(files,"data/merged.loom",key="Accession")

# From the Velocyto run, required quantification of intronic reads and exonic reads are already abundant in the project's folder.  

import scvelo as sc
# Import merged loom file, just created above
adata = sc.read_loom("data/merged.loom")
adata.var_names_make_unique()
adata.var.index.is_unique
# Calculate the nuclear RNA fraction using the spliced and unspliced reads' matrices
exon_sum = adata.layers['spliced'].sum(axis=1)
intron_sum = adata.layers['unspliced'].sum(axis=1)
nuclear_fraction = intron_sum/(exon_sum + intron_sum)
### alternative (poor precision) method ### 
# sc_obj$nuc_frac <- (sc_obj$nCount_unspliced) / (big_sc$nCount_spliced + big_sc$nCount_unspliced)
