# IsoSeq
Bacterial isolate sequencing bioinformatics pipeline.

## Install

```
git clone https://github.com/ronan-m-doyle/IsoSeq.git
cd IsoSeq/

mamba env create -f isoseq_pipeline_env.yml -n isoseq_pipeline
mamba create -n abricate -c conda-forge -c bioconda abricate
mamba create -n abritamr -c bioconda abritamr
mamba create -n mlst_contigs -c conda-forge -c bioconda mlst
mamba create -n gubbins -c bioconda gubbins
mamba create -n clair3 -c bioconda clair3

wget "https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_08_GB_20251015.tar.gz"
mkdir kraken_reference/
tar -xzvf k2_pluspf_08_GB_20251015.tar.gz -C kraken_reference/
rm k2_pluspf_08_GB_20251015.tar.gz
```

## Run

```
mamba activate isoseq_pipeline
bash isoseq_pipeline.sh --samplesheet samplesheet.csv --threads 1 --rundir example_run_directory/
```
