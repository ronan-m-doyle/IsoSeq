# synnovis_isolate_sequencing_pipeline
Bacterial isolate sequencing bioinformatics pipeline.

## Install

```
git clone https://github.com/ronan-m-doyle/synnovis_isolate_sequencing_pipeline.git
cd synnovis_isolate_sequencing_pipeline
mamba env create -f isolate_pipeline_env.yml -n isoseq_pipeline
mamba create -n abricate -c conda-forge -c bioconda abricate
mamba create -n abritamr -c bioconda abritamr
mamba create -n mlst_contigs -c conda-forge -c bioconda mlst
wget "https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_08_GB_20251015.tar.gz"
```
