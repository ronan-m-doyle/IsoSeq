# IsoSeq
Bacterial and Candida isolate sequencing bioinformatics pipeline.

## Install

```
git clone https://github.com/ronan-m-doyle/IsoSeq.git
cd IsoSeq/

mamba env create -f isoseq_pipeline_env.yml -n isoseq_pipeline
mamba create -n mlst_contigs -c conda-forge -c bioconda mlst=2.32.2
mamba create -n gubbins -c bioconda gubbins=3.4.3
mamba create -n clair3 -c bioconda clair3=1.2.0
mamba create -n chroquetas -c bioconda -c nmquijada chroquetas=1.0.1

git clone https://github.com/rrwick/Porechop.git
cd Porechop
make
./porechop-runner.py -h
cd ..
cp scripts/adapters.py Porechop/porechop/

wget "https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_16_GB_20251015.tar.gz"
mkdir databases/kraken_reference_16/
tar -xzvf k2_pluspf_16_GB_20251015.tar.gz -C databases/kraken_reference_16/
rm k2_pluspf_16_GB_20251015.tar.gz
```

## Run

```
mamba activate isoseq_pipeline
bash isoseq_pipeline.sh --samplesheet samplesheet.csv --threads 1 --rundir example_run_directory/
```
