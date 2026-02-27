#!/usr/bin/env bash
set -uo pipefail
IFS=$'\n\t'

# -------------------------------------------------------------
# Refactored Isolate Sequencing Pipeline
# Required arguments: --run_directory, --output_directory, --samplesheet
# -------------------------------------------------------------

# Parse command-line options
run_directory=""
output_directory=""
SAMPLESHEET=""

while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --run_directory)
            run_directory="$2"
            shift; shift ;;
        --output_directory)
            output_directory="$2"
            shift; shift ;;
        --samplesheet)
            SAMPLESHEET="$2"
            shift; shift ;;
        *)
            echo "Unknown option: $1"; exit 1 ;;
    esac
done

# Validate required arguments
if [[ -z "$run_directory" || -z "$output_directory" || -z "$SAMPLESHEET" ]]; then
    echo "❌ ERROR: You must specify --run_directory, --output_directory, and --samplesheet"
    exit 1
fi

# Resolve absolute paths to avoid issues after changing directories
run_directory=$(realpath "$run_directory")
output_directory=$(realpath "$output_directory")
SAMPLESHEET=$(realpath "$SAMPLESHEET")

# Check if samplesheet exists
if [[ ! -f "$SAMPLESHEET" ]]; then
    echo "❌ ERROR: Samplesheet not found at $SAMPLESHEET"
    exit 1
fi

# Move into the output directory
mkdir -p "$output_directory"
cd "$output_directory"

# -------------------------------------------------------------
# Conda environments
# -------------------------------------------------------------
ENV_MLST="mlst_contigs"
ENV_PROKKA="prokka"
ENV_AMRFINDER="amrfinder"

# -------------------------------------------------------------
# Species dictionary: maps organism name → parameters
# -------------------------------------------------------------
declare -A species_dict_read_mlst=(
    ["Acinetobacter baumannii"]="Acinetobacter_baumannii"
    ["Escherichia coli"]="Escherichia_coli"
    ["Klebsiella pneumoniae"]="Klebsiella_pneumoniae"
    ["Staphylococcus aureus"]="Staphylococcus_aureus"
    ["Enterococcus faecium"]="Enterococcus_faecium"
    ["Enterococcus faecalis"]="Enterococcus_faecalis"
)

declare -A species_dict_mlst_contigs=(
    ["Acinetobacter baumannii"]="abaumannii_2"
    ["Escherichia coli"]="ecoli"
    ["Klebsiella pneumoniae"]="klebsiella"
    ["Staphylococcus aureus"]="saureus"
    ["Enterococcus faecium"]="efaecium"
    ["Enterococcus faecalis"]="efaecalis"
)

declare -A species_dict_prokka=(
    ["Acinetobacter baumannii"]="Acinetobacter_baumannii"
    ["Escherichia coli"]="Escherichia"
    ["Klebsiella pneumoniae"]="Klebsiella_pneumoniae"
    ["Staphylococcus aureus"]="Staphylococcus_aureus"
    ["Enterococcus faecium"]="Enterococcus_faecium"
    ["Enterococcus faecalis"]="Enterococcus_faecalis"
)

declare -A species_dict_busco=(
    ["Acinetobacter baumannii"]="acinetobacter_odb12"
    ["Escherichia coli"]="enterobacteriaceae_odb12"
    ["Klebsiella pneumoniae"]="enterobacteriaceae_odb12"
    ["Staphylococcus aureus"]="staphylococcus_odb12"
    ["Enterococcus faecium"]="enterococcus_odb12"
    ["Enterococcus faecalis"]="enterococcus_odb12"
)

declare -A species_dict_size=(
    ["Acinetobacter baumannii"]="3.9m"
    ["Escherichia coli"]="4.6m"
    ["Klebsiella pneumoniae"]="5.2m"
    ["Staphylococcus aureus"]="2.8m"
    ["Enterococcus faecium"]="2.5m"
    ["Enterococcus faecalis"]="2.8m"
)

# -------------------------------------------------------------
# Read samplesheet and populate arrays
# -------------------------------------------------------------
SAMPLES=()
BARCODES=()
ORGANISMS=()

while IFS=, read -r barcode sample_name organism; do
    # Skip header or empty lines
    if [[ "$barcode" == "barcode" || -z "$barcode" ]]; then
        continue
    fi

    if [[ -z "${species_dict_read_mlst[$organism]:-}" ]]; then
        echo "❌ ERROR: Unknown organism '$organism' in samplesheet. Please add to dictionary."
        exit 1
    fi

    BARCODES+=("$barcode")
    SAMPLES+=("$sample_name")
    ORGANISMS+=("$organism")

done < "$SAMPLESHEET"

# -------------------------------------------------------------
# Step 1: Concatenate barcode fastq files to sample.fastq.gz
# -------------------------------------------------------------
for idx in "${!BARCODES[@]}"; do
    barcode="${BARCODES[$idx]}"
    sample="${SAMPLES[$idx]}"
    src_pattern="$run_directory/*/*/fastq_pass/barcode${barcode}/*.fastq.gz"
    out_file="${sample}.fastq.gz"

    shopt -s nullglob
    files=( $src_pattern )
    if [ ${#files[@]} -eq 0 ]; then
        echo "WARNING: No fastq.gz files found for barcode${barcode}, skipping."
        continue
    fi
    echo "Concatenating ${#files[@]} files for barcode ${barcode} -> $out_file"
    cat ${files[@]} > "$out_file"
done

# -------------------------------------------------------------
# Step 2: Trim adapters with Porechop and NanoStat QC
# -------------------------------------------------------------
mkdir -p trimmed qc
for sample in "${SAMPLES[@]}"; do
    fq="${sample}.fastq.gz"
    out_trim="trimmed/${sample}.fastq.gz"
    echo "Trimming $fq -> $out_trim"
    porechop -i "$fq" -o "$out_trim" -t 20 --no_split
    echo "Running NanoStat for $out_trim"
    NanoStat --fastq "$out_trim" -n "qc/${sample}_nanostat.txt" -t 20
done

# Clean up non-trimmed fastq files
rm *.fastq.gz

# -------------------------------------------------------------
# Step 3: Kraken2 classification
# -------------------------------------------------------------
mkdir -p kraken2
KRAKEN2_DB="/home/ronan/reference_genomes/kraken2_db/"
KRKN_THREADS=12
for fq in trimmed/*.fastq.gz; do
    base="$(basename "$fq" .fastq.gz)"
    kraken_out="kraken2/${base}.k2output"
    kraken_report="kraken2/${base}.k2report"
    echo "Running Kraken2 for $fq"
    kraken2 --db "$KRAKEN2_DB" --threads $KRKN_THREADS --output "$kraken_out" --report "$kraken_report" --use-names "$fq"
done

PY_KRAKEN_SUMMARY="/home/ronan/useful_scripts/kraken2_top5_species.py"
[ -f "$PY_KRAKEN_SUMMARY" ] && python3 "$PY_KRAKEN_SUMMARY"

# -------------------------------------------------------------
# Step 4: Read-based MLST with Krocus
# -------------------------------------------------------------
mkdir -p mlst_reads
for idx in "${!SAMPLES[@]}"; do
    sample="${SAMPLES[$idx]}"
    organism="${ORGANISMS[$idx]}"
    read_mlst="${species_dict_read_mlst[$organism]}"
    trimmed_fq="trimmed/${sample}.fastq.gz"
    out_tsv="mlst_reads/${sample}_krocusmlst.tsv"
    echo "Running Krocus for $sample (${read_mlst})"
    krocus -o "$out_tsv" "/home/ronan/isolate_sequencing/mlst_krocus/${read_mlst}/" "$trimmed_fq"
done

# -------------------------------------------------------------
# Step 5: Assembly with autoautocycler
# -------------------------------------------------------------
mkdir -p assemblies
for idx in "${!SAMPLES[@]}"; do
    sample="${SAMPLES[$idx]}"
    size="${species_dict_size[${ORGANISMS[$idx]}]}"
    trimmed_fq="trimmed/${sample}.fastq.gz"
    outdir="assemblies/${sample}"
    mkdir -p "$outdir"
    echo "Assembling $sample -> $outdir"
    autoautocycler.sh -o "$outdir/" -t 20 -s "$size" -a "flye raven" "$trimmed_fq"
    cp "$outdir/${sample}.fasta" "assemblies/${sample}.fasta" 2>/dev/null || true
    cp "$outdir/${sample}.gfa" "assemblies/${sample}.gfa" 2>/dev/null || true
done

# -------------------------------------------------------------
# Step 6: Map reads to assemblies, BAM, indexing, Qualimap
# -------------------------------------------------------------
for asm in assemblies/*.fasta; do
    sample="$(basename "$asm" .fasta)"
    trimmed_fq="trimmed/${sample}.fastq.gz"
    out_bam="assemblies/${sample}.bam"
    echo "Mapping $trimmed_fq to $asm -> $out_bam"
    minimap2 -x map-ont --secondary no -a -t 20 "$asm" "$trimmed_fq" | samtools sort -O BAM -o "$out_bam"
    samtools index -M "$out_bam"
    mkdir -p "qc/${sample}"
    qualimap bamqc -bam "$out_bam" -outdir "qc/${sample}" --java-mem-size=16G
done

# -------------------------------------------------------------
# Step 7: QUAST
# -------------------------------------------------------------
for asm in assemblies/*.fasta; do
    sample="$(basename "$asm" .fasta)"
    quast -o "qc/${sample}_quast" -t 20 -m 100 -l "$sample" "$asm"
done

# -------------------------------------------------------------
# Step 8: Bandage
# -------------------------------------------------------------
for gfa in assemblies/*.gfa; do
    sample="$(basename "$gfa" .gfa)"
    Bandage image "$gfa" "qc/${sample}.svg"
done

# -------------------------------------------------------------
# Step 9: MLST (contigs)
# -------------------------------------------------------------
for idx in "${!SAMPLES[@]}"; do
    sample="${SAMPLES[$idx]}"
    organism="${ORGANISMS[$idx]}"
    scheme="${species_dict_mlst_contigs[$organism]}"
    echo "Running MLST (contigs) for $sample ($scheme)"
    conda run -n "$ENV_MLST" mlst --scheme "$scheme" "assemblies/${sample}.fasta" > "${sample}_contigmlst.tsv"
done

# -------------------------------------------------------------
# Step 10: Prokka annotation
# -------------------------------------------------------------
for idx in "${!SAMPLES[@]}"; do
    sample="${SAMPLES[$idx]}"
    organism="${ORGANISMS[$idx]}"
    prokka_species="${species_dict_prokka[$organism]}"
    echo "Running Prokka for $sample ($prokka_species)"
    conda run -n "$ENV_PROKKA" prokka --outdir "assemblies/${sample}/" --force --species "$prokka_species" \
        --prefix "$sample" --strain "$sample" --kingdom Bacteria --cpus 20 --gcode 11 "assemblies/${sample}.fasta"
done

# -------------------------------------------------------------
# Step 11: AMRFinder
# -------------------------------------------------------------
for idx in "${!SAMPLES[@]}"; do
    sample="${SAMPLES[$idx]}"
    organism="${ORGANISMS[$idx]}"
    amr_species="${species_dict_prokka[$organism]}"
    echo "Running AMRFinder for $sample ($amr_species)"
    conda run -n "$ENV_AMRFINDER" amrfinder -n "assemblies/${sample}.fasta" -O "$amr_species" --threads 20 \
        -o "${sample}_resistance.tsv" -g "assemblies/${sample}/${sample}.gff" -p "assemblies/${sample}/${sample}.faa" -a prokka
done

# -------------------------------------------------------------
# Step 12: BUSCO assessment
# -------------------------------------------------------------
mkdir -p busco_results
for idx in "${!SAMPLES[@]}"; do
    sample="${SAMPLES[$idx]}"
    organism="${ORGANISMS[$idx]}"
    busco_lineage="${species_dict_busco[$organism]}"
    assembly_fasta="assemblies/${sample}.fasta"
    outdir="qc/busco_results/"

    echo "Running BUSCO for $sample ($busco_lineage)"
    busco -i "$assembly_fasta" -f -m genome -l "$busco_lineage" -c 20 --out_path "$outdir" -o "$sample"
done

rm -rf busco_downloads/

# -------------------------------------------------------------
# Step 13: MultiQC report
# -------------------------------------------------------------
currentPWD=$(basename "$output_directory")
mkdir -p multiqc_report
multiqc -i "$currentPWD" -f -v -o multiqc_report/ --no-ai -c ~/isolate_sequencing/multiqc_config.yaml .

# -------------------------------------------------------------
# Step 14: Copy respective genomes to specifc organism tree folders
# -------------------------------------------------------------

for idx in "${!SAMPLES[@]}"; do
    sample="${SAMPLES[$idx]}"
    organism="${ORGANISMS[$idx]}"
    read_mlst="${species_dict_read_mlst[$organism]}"
    trimmed_fq="trimmed/${sample}.fastq.gz"
    assembly_fasta="assemblies/${sample}.fasta"
    organism_tree_folder="/home/ronan/isolate_sequencing/${read_mlst}_tree/"
    echo "copying $sample assembly fasta and read fastq files to $read_mlst tree folder"
    mkdir "$organism_tree_folder"
    cp "$assembly_fasta" "$organism_tree_folder/${sample}.fasta"
    cp "$trimmed_fq" "$organism_tree_folder/${sample}.fastq.gz"
done
