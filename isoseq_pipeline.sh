#!/usr/bin/env bash
set -euo pipefail
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
# Species dictionary: maps organism name → parameters
# -------------------------------------------------------------
declare -A species_dict_read_mlst=(
["Acinetobacter baumannii"]="Acinetobacter_baumannii"
["Escherichia coli"]="Escherichia_coli"
)


declare -A species_dict_mlst_contigs=(
["Acinetobacter baumannii"]="abaumannii_2"
["Escherichia coli"]="ecoli"
)


declare -A species_dict_prokka=(
["Acinetobacter baumannii"]="Acinetobacter_baumannii"
["Escherichia coli"]="Escherichia"
)


declare -A species_dict_busco=(
["Acinetobacter baumannii"]="acinetobacter_odb12"
["Escherichia coli"]="enterobacteriaceae_odb12"
)


declare -A species_dict_size=(
["Acinetobacter baumannii"]="3.9m"
["Escherichia coli"]="4.6m"
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

# ---------------------------
# Other configuration (static)
# ---------------------------
KRAKEN2_DB="/home/ronan/reference_genomes/kraken2_db/"
KRKN_THREADS=12
COMMON_THREADS=20
MLST_KROCUS_DIR_BASE="/home/ronan/isolate_sequencing/mlst_krocus"
MULTIQC_CONFIG="~/isolate_sequencing/multiqc_config.yaml"
PY_KRAKEN_SUMMARY="/home/ronan/useful_scripts/kraken2_top5_species.py"

# Output directories
DIR_TRIMMED="trimmed"
DIR_KRAKEN="kraken2"
DIR_MLST_READS="mlst_reads"
DIR_ASSEMBLIES="assemblies"
DIR_QC="qc"
DIR_PROKKA="prokka"
DIR_AMRFINDER="amrfinder"
DIR_BUSCO="busco"
DIR_MULTIQC="multiqc_report"

# Tools
CMD_PORECHOP="porechop"
CMD_NANOSTAT="NanoStat"
CMD_KRAKEN2="kraken2"
CMD_KROCUS="krocus"
CMD_AUTOAUTO="autoautocycler.sh"
CMD_MINIMAP2="minimap2"
CMD_SAMTOOLS="samtools"
CMD_QUALIMAP="qualimap"
CMD_QUAST="quast"
CMD_BANDAGE="Bandage"
CMD_BUSCO="busco"
CMD_MULTIQC="multiqc"

# ---------------------------
# Helper functions and checks
# ---------------------------
log(){ echo "[\$(date -Iseconds)] \$*"; }
check_cmd(){
  local c
  for c in "$@"; do
    if ! command -v "$c" >/dev/null 2>&1; then
      echo "ERROR: required command '$c' not found in PATH. Please install or update PATH." >&2
      exit 1
    fi
  done
}

safe_mkdir(){ mkdir -p "$@"; }
trap 'echo "ERROR at line $LINENO. Last command: $BASH_COMMAND"' ERR

log "Starting pre-flight checks"
check_cmd bash gzip zcat cat sed head pigz "$CMD_PORECHOP" "$CMD_NANOSTAT" "$CMD_KRAKEN2" python3 "$CMD_MINIMAP2" "$CMD_SAMTOOLS" "$CMD_QUALIMAP" "$CMD_QUAST" "$CMD_BANDAGE" "$CMD_BUSCO" "$CMD_MULTIQC"
log "All required commands present"

safe_mkdir "$DIR_TRIMMED" "$DIR_KRAKEN" "$DIR_MLST_READS" "$DIR_ASSEMBLIES" "$DIR_QC" "$DIR_PROKKA" "$DIR_AMRFINDER" "$DIR_BUSCO" "$DIR_MULTIQC"

# ---------------------------
# Step 1: Concatenate barcodes -> named sample.fastq.gz
# ---------------------------
log "Concatenating barcodes from $run_directory"
for barcode in "${!barcode_to_sample[@]}"; do
  sample=${barcode_to_sample[$barcode]}
  src_pattern="${run_directory}/barcode${barcode}/*.fastq.gz"
  out_file="${sample}.fastq.gz"
  log "Processing barcode $barcode -> sample $sample"
  shopt -s nullglob
  files=( $src_pattern )
  if [ ${#files[@]} -eq 0 ]; then
    log "WARNING: no fastq.gz files found for barcode${barcode}. Skipping."
    continue
  fi
  cat "${run_directory}/barcode${barcode}"/*.fastq.gz > "$out_file"
  log "Wrote $out_file"
done

# ---------------------------
# Step 2: Trim adapters with Porechop and compute NanoStat
# ---------------------------
log "Trimming and NanoStat generation"
for sample in "${samples_order[@]}"; do
  fq="${sample}.fastq.gz"
  if [ ! -f "$fq" ]; then
    log "Raw FASTQ for \$sample not found (\$fq). Skipping trimming."
    continue
  fi
  base="$sample"
  out_trim="${DIR_TRIMMED}/${base}.trim.fastq.gz"
  log "Trimming: $fq -> $out_trim"
  "$CMD_PORECHOP" -i "$fq" -o "$out_trim" -t 20 --no_split
  log "Running NanoStat on trimmed reads"
  "$CMD_NANOSTAT" --fastq "$out_trim" -n "${DIR_QC}/${base}_nanostat" -t 20
done

# ---------------------------
# Step 3: Kraken2 classification
# ---------------------------
log "Running Kraken2 classification (sequential)"
for fq in ${DIR_TRIMMED}/*.trim.fastq.gz; do
  [ -e "$fq" ] || continue
  base="$(basename "$fq" .trim.fastq.gz)"
  out_report="${DIR_KRAKEN}/${base}.k2report"
  out_output="${DIR_KRAKEN}/${base}.k2output"
  log "Kraken2: $fq -> report $out_report"
  "$CMD_KRAKEN2" --db "$KRAKEN2_DB" --threads "$KRKN_THREADS" --output "$out_output" --report "$out_report" --use-names "$fq"
done

log "Summarising Kraken2 results with python script (if present)"
if [ -x "$PY_KRAKEN_SUMMARY" ] || [ -f "$PY_KRAKEN_SUMMARY" ]; then
  python3 "$PY_KRAKEN_SUMMARY"
else
  log "Warning: Kraken summary script not found at $PY_KRAKEN_SUMMARY. Skipping."
fi

# ---------------------------
# Step 4: Read-based MLST with Krocus
# ---------------------------
log "Running Krocus (MLST from reads)"
for sample in "${samples_order[@]}"; do
  trimmed_fq="${DIR_TRIMMED}/${sample}.trim.fastq.gz"
  if [ ! -f "$trimmed_fq" ]; then
    log "Trimmed FASTQ not found for $sample (expected $trimmed_fq). Skipping Krocus."
    continue
  fi
  organism="${sample_to_read_mlst[$sample]}"
  out_tsv="${DIR_MLST_READS}/${sample}_krocusmlst.tsv"
  log "Krocus: sample=$sample organism=$organism"
  pigz -d -c "$trimmed_fq" | head -n 40000 | "$CMD_KROCUS" -o "$out_tsv" "${MLST_KROCUS_DIR_BASE}/${organism}/" -
  log "Wrote $out_tsv"
done

# ---------------------------
# Step 5: Assembly (autoautocycler)
# ---------------------------
log "Running autoautocycler.sh for each sample"
for sample in "${samples_order[@]}"; do
  trimmed_fq="${DIR_TRIMMED}/${sample}.trim.fastq.gz"
  if [ ! -f "$trimmed_fq" ]; then
    log "Trimmed FASTQ not found for $sample. Skipping assembly."
    continue
  fi
  size_param="${sample_to_size[$sample]}"
  outdir="${DIR_ASSEMBLIES}/${sample}_autocycler"
  log "Assembly: sample=$sample size=$size_param -> $outdir"
  safe_mkdir "$outdir"
  "$CMD_AUTOAUTO" -o "$outdir/" -t "$COMMON_THREADS" -s "$size_param" -a "flye raven" "$trimmed_fq"
  if [ -f "$outdir/${sample}.fasta" ]; then
    cp "$outdir/${sample}.fasta" "${DIR_ASSEMBLIES}/${sample}.fasta"
  fi
  if [ -f "$outdir/${sample}.gfa" ]; then
    cp "$outdir/${sample}.gfa" "${DIR_ASSEMBLIES}/${sample}.gfa"
  fi
done

# ---------------------------
# Step 6: Map reads to assemblies -> BAM, index, Qualimap
# ---------------------------
log "Mapping trimmed reads to assemblies and running Qualimap"
for asm in ${DIR_ASSEMBLIES}/*.fasta; do
  [ -e "$asm" ] || continue
  sample_name="$(basename "$asm" .fasta)"
  trimmed_fq="${DIR_TRIMMED}/${sample_name}.trim.fastq.gz"
  if [ ! -f "$trimmed_fq" ]; then
    log "Trimmed reads for $sample_name not found. Skipping mapping."
    continue
  fi
  out_bam="${DIR_QC}/${sample_name}.bam"
  log "minimap2 mapping: $trimmed_fq -> $asm"
  "$CMD_MINIMAP2" -x map-ont --secondary no -a -t "$COMMON_THREADS" "$asm" "$trimmed_fq" | "$CMD_SAMTOOLS" sort -O BAM -o "$out_bam"
  log "Indexing BAM: $out_bam"
  "$CMD_SAMTOOLS" index -M "$out_bam"
  log "Running Qualimap on $out_bam"
  safe_mkdir "${DIR_QC}/${sample_name}"
  "$CMD_QUALIMAP" bamqc -bam "$out_bam" -outdir "${DIR_QC}/${sample_name}" --java-mem-size=16G
done

# ---------------------------
# Step 7: QUAST on assemblies
# ---------------------------
log "Running QUAST for each assembly"
for asm in ${DIR_ASSEMBLIES}/*.fasta; do
  [ -e "$asm" ] || continue
  sample_name="$(basename "$asm" .fasta)"
  outdir="${DIR_QC}/${sample_name}_quast"
  log "QUAST: $asm -> $outdir"
  safe_mkdir "$outdir"
  "$CMD_QUAST" -o "$outdir" -t "$COMMON_THREADS" -m 100 -L "$asm"
done

# ---------------------------
# Step 8: Bandage images for GFAs
# ---------------------------
log "Generating Bandage images for any .gfa files"
for gfa in ${DIR_ASSEMBLIES}/*.gfa; do
  [ -e "$gfa" ] || continue
  sample_name="$(basename "$gfa" .gfa)"
  out_svg="${DIR_QC}/${sample_name}.svg"
  log "Bandage image: $gfa -> $out_svg"
  "$CMD_BANDAGE" image "$gfa" "$out_svg"
done

# ---------------------------
# Step 9: MLST (contigs) with mlst
# ---------------------------
log "Running MLST on contigs"
for sample in "${samples_order[@]}"; do
  fasta="${DIR_ASSEMBLIES}/${sample}_autocycler.fasta"
  if [ ! -f "$fasta" ]; then
    fasta="${DIR_ASSEMBLIES}/${sample}.fasta"
  fi
  if [ ! -f "$fasta" ]; then
    log "No assembly fasta found for $sample. Skipping mlst."
    continue
  fi
  scheme="${sample_to_mlst_contigs[$sample]}"
  out_tsv="${DIR_QC}/${sample}_contigmlst.tsv"
  log "MLST(contigs): sample=$sample scheme=$scheme"
  conda run -n mlst_contigs mlst --scheme "$scheme" "$fasta" > "$out_tsv" || log "mlst returned non-zero for $sample"
done

# ---------------------------
# Step 10: Prokka annotation
# ---------------------------
log "Running Prokka annotations"
for sample in "${samples_order[@]}"; do
  fasta="${DIR_ASSEMBLIES}/${sample}_autocycler.fasta"
  if [ ! -f "$fasta" ]; then
    fasta="${DIR_ASSEMBLIES}/${sample}.fasta"
  fi
  if [ ! -f "$fasta" ]; then
    log "No assembly fasta for $sample. Skipping Prokka."
    continue
  fi
  org="${sample_to_prokka[$sample]}"
  outdir="${DIR_PROKKA}/${sample}_prokka"
  safe_mkdir "$outdir"
  log "Prokka: sample=$sample species=$org outdir=$outdir"
  conda run -n prokka prokka --outdir "$outdir" --force --species "$org" --prefix "${sample}" --kingdom Bacteria --cpus "$COMMON_THREADS" --gcode 11 "$fasta"
  if [ -f "$outdir/${sample}.gff" ]; then
    cp "$outdir/${sample}.gff" "${DIR_PROKKA}/${sample}.gff"
  fi
  if [ -f "$outdir/${sample}.faa" ]; then
    cp "$outdir/${sample}.faa" "${DIR_PROKKA}/${sample}.faa"
  fi
done

# ---------------------------
# Step 11: AMRFinder
# ---------------------------
log "Running AMRFinder"
for sample in "${samples_order[@]}"; do
  fasta="${DIR_ASSEMBLIES}/${sample}_autocycler.fasta"
  if [ ! -f "$fasta" ]; then
    fasta="${DIR_ASSEMBLIES}/${sample}.fasta"
  fi
  if [ ! -f "$fasta" ]; then
    log "No assembly fasta for $sample. Skipping AMRFinder."
    continue
  fi
  org="${sample_to_prokka[$sample]}"
  prokka_dir="${DIR_PROKKA}/${sample}_prokka"
  gff="${prokka_dir}/${sample}.gff"
  faa="${prokka_dir}/${sample}.faa"
  out_tsv="${DIR_AMRFINDER}/${sample}_amrfinder.tsv"
  if [ ! -f "$gff" ] || [ ! -f "$faa" ]; then
    log "Prokka outputs missing for $sample (gff/faa). Skipping AMRFinder."
    continue
  fi
  log "AMRFinder: sample=$sample org=$org"
  conda run -n amrfinder amrfinder -n "$fasta" -O "$org" --threads "$COMMON_THREADS" -o "$out_tsv" -g "$gff" -p "$faa" -a prokka || log "AMRFinder non-zero for $sample"
done

# ---------------------------
# Step 12: BUSCO
# ---------------------------
log "Running BUSCO"
for sample in "${samples_order[@]}"; do
  fasta="${DIR_ASSEMBLIES}/${sample}_autocycler.fasta"
  if [ ! -f "$fasta" ]; then
    fasta="${DIR_ASSEMBLIES}/${sample}.fasta"
  fi
  if [ ! -f "$fasta" ]; then
    log "No assembly fasta for $sample. Skipping BUSCO."
    continue
  fi
  lineage="${sample_to_busco[$sample]}"
  outpath="${DIR_BUSCO}/${sample}"
  safe_mkdir "$outpath"
  log "BUSCO: sample=$sample lineage=$lineage"
  "$CMD_BUSCO" -i "$fasta" -f -m genome -l "$lineage" -c "$COMMON_THREADS" --out_path "$outpath" -o "${sample}"
done

if [ -d "busco_downloads" ]; then
  log "Removing busco_downloads/ directory"
  rm -rf busco_downloads/
fi

# ---------------------------
# Step 13: MultiQC
# ---------------------------
log "Running MultiQC to aggregate results"
"$CMD_MULTIQC" -i "$PWD" -f -v -o "$DIR_MULTIQC" --no-ai -c "$MULTIQC_CONFIG" . || log "MultiQC returned non-zero"

log "Pipeline complete. Outputs summarized in: $PWD/$DIR_MULTIQC"

exit 0
