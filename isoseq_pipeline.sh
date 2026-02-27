#!/usr/bin/env bash
# ----------------------------------------------------------------------
# Isolate Sequencing Pipeline
# ----------------------------------------------------------------------
IFS=$'\n\t'

# -------------------------- CONFIG ----------------------------
KRAKEN2_DB="${HOME}/reference_genomes/kraken2_db/"
PY_KRAKEN_SUMMARY="${HOME}/useful_scripts/kraken2_top5_species.py"
MULTIQC_CONFIG="${HOME}/isolate_sequencing/multiqc_config.yaml"
TREE_BASE="${HOME}/isolate_sequencing"
KROCUS_DB="${HOME}/isolate_sequencing/mlst_krocus/"

# -------------------------- ARG PARSE --------------------------
THREADS=1
SAMPLESHEET=""
RUN_DIR=""

while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --samplesheet) SAMPLESHEET="$2"; shift; shift ;;
        --threads) THREADS="$2"; shift; shift ;;
        --rundir) RUN_DIR="$2"; shift; shift ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

if [[ -z "$SAMPLESHEET" || -z "$THREADS" || -z "$RUN_DIR" ]]; then
    echo "❌ ERROR: Must specify --samplesheet, --threads, --rundir"
    exit 1
fi

RUN_DIR=$(realpath "$RUN_DIR")
[[ -d "$RUN_DIR" ]] || { echo "Run directory not found"; exit 1; }

SAMPLESHEET=$(realpath "$SAMPLESHEET")
[[ -f "$SAMPLESHEET" ]] || { echo "Samplesheet not found"; exit 1; }

# -------------------------- ENVS -------------------------------
ENV_MLST="mlst_contigs"
ENV_PROKKA="prokka"
ENV_ABRITAMR="abritamr"

# ----------------------- DICTIONARIES --------------------------
declare -A species_dict_read_mlst=(
    ["Acinetobacter baumannii"]="Acinetobacter_baumannii"
    ["Escherichia coli"]="Escherichia_coli"
    ["Klebsiella pneumoniae"]="Klebsiella_pneumoniae"
    ["Staphylococcus aureus"]="Staphylococcus_aureus"
    ["Enterococcus faecium"]="Enterococcus_faecium"
    ["Enterococcus faecalis"]="Enterococcus_faecalis"
    ["Citrobacter koseri"]="Citrobacter_freundii"
    ["Enterobacter hormaechei"]="Enterobacter_cloacae"
    ["Pseudomonas aeruginosa"]="Pseudomonas_aeruginosa"
)

declare -A species_dict_mlst_contigs=(
    ["Acinetobacter baumannii"]="abaumannii_2"
    ["Escherichia coli"]="ecoli"
    ["Klebsiella pneumoniae"]="klebsiella"
    ["Staphylococcus aureus"]="saureus"
    ["Enterococcus faecium"]="efaecium"
    ["Enterococcus faecalis"]="efaecalis"
    ["Citrobacter koseri"]="cfreundii"
    ["Enterobacter hormaechei"]="ecloacae"
    ["Pseudomonas aeruginosa"]="paeruginosa"
)

declare -A species_dict_amrfinder=(
    ["Acinetobacter baumannii"]="Acinetobacter_baumannii"
    ["Escherichia coli"]="Escherichia"
    ["Klebsiella pneumoniae"]="Klebsiella_pneumoniae"
    ["Staphylococcus aureus"]="Staphylococcus_aureus"
    ["Enterococcus faecium"]="Enterococcus_faecium"
    ["Enterococcus faecalis"]="Enterococcus_faecalis"
    ["Citrobacter koseri"]="Citrobacter_freundii"
    ["Enterobacter hormaechei"]="Enterobacter_cloacae"
    ["Pseudomonas aeruginosa"]="Pseudomonas_aeruginosa"
)

declare -A species_dict_busco=(
    ["Acinetobacter baumannii"]="acinetobacter_odb12"
    ["Escherichia coli"]="enterobacteriaceae_odb12"
    ["Klebsiella pneumoniae"]="enterobacteriaceae_odb12"
    ["Staphylococcus aureus"]="staphylococcus_odb12"
    ["Enterococcus faecium"]="enterococcus_odb12"
    ["Enterococcus faecalis"]="enterococcus_odb12"
    ["Citrobacter koseri"]="enterobacteriaceae_odb12"
    ["Enterobacter hormaechei"]="enterobacter_odb12"
    ["Pseudomonas aeruginosa"]="pseudomonas_odb12"
)

declare -A species_dict_size=(
    ["Acinetobacter baumannii"]="3.9m"
    ["Escherichia coli"]="4.6m"
    ["Klebsiella pneumoniae"]="5.2m"
    ["Staphylococcus aureus"]="2.8m"
    ["Enterococcus faecium"]="2.5m"
    ["Enterococcus faecalis"]="2.8m"
    ["Citrobacter koseri"]="4.9m"
    ["Enterobacter hormaechei"]="5m"
    ["Pseudomonas aeruginosa"]="6.3m"
)

# ---------------------- READ SAMPLESHEET -----------------------
SAMPLES=(); BARCODES=(); ORGANISMS=()
while IFS=, read -r barcode sample organism; do
    [[ "$barcode" == "barcode" || -z "$barcode" ]] && continue
    if [[ -z "${species_dict_read_mlst[$organism]:-}" ]]; then
        echo "❌ Unknown organism: $organism"; exit 1
    fi
    SAMPLES+=("$sample")
    BARCODES+=("$barcode")
    ORGANISMS+=("$organism")
done < "$SAMPLESHEET"

# ---------------------- PIPELINE LOG DIR -----------------------
mkdir -p pipeline_logs
RUN_TIMESTAMP=$(date '+%Y%m%d_%H%M%S')
MASTER_SUMMARY="pipeline_logs/pipeline_summary_${RUN_TIMESTAMP}.tsv"
echo -e "Sample\tStatus\tFailed_Step\tTotal_Time(s)" > "$MASTER_SUMMARY"

# ============================================================== 
# PER-SAMPLE PIPELINE 
# ============================================================== 

for idx in "${!SAMPLES[@]}"; do
    sample="${SAMPLES[$idx]}"
    barcode="${BARCODES[$idx]}"
    organism="${ORGANISMS[$idx]}"

    sample_dir="${sample}"
    mkdir -p "$sample_dir/logs"
    LOGFILE="${sample_dir}/logs/pipeline.log"

    log() { echo -e "[$(date '+%F %T')] $*" | tee -a "$LOGFILE"; }

    log "========== Processing sample: $sample ($organism) =========="
    start_total=$(date +%s)
    status="OK"; failed_step="-"

    # ---------------- STEP 1: Concatenate ----------------
    step="concat"; step_start=$(date +%s)
    log "▶ Step $step"
    src_pattern="${RUN_DIR}/**/fastq_pass/barcode${barcode}/*.fastq.gz"
    shopt -s globstar nullglob
    files=( $src_pattern )
    if (( ${#files[@]} == 0 )); then
        log "⚠️ No FASTQ found for barcode${barcode}"
        status="FAILED"; failed_step="$step"
        echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"
        continue
    fi
    log "Concatenating ${#files[@]} fastq files → ${sample_dir}/${sample}.fastq.gz"
    if !  cat "${files[@]}" > "${sample_dir}/${sample}.fastq.gz"; then
       log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; continue
    fi
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- STEP 2: Porechop & NanoStat ----------------
    step="porechop_nanostat"; step_start=$(date +%s)
    log "▶ Step $step"
    mkdir -p ${sample_dir}/trimmed ${sample_dir}/qc
    fq="${sample_dir}/${sample}.fastq.gz"; out_trim="${sample_dir}/trimmed/${sample}.fastq.gz"
    if ! porechop -i "$fq" -o "$out_trim" -t "$THREADS" --no_split 2>>"$LOGFILE"; then
        log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; continue
    fi
    NanoStat --fastq "$out_trim" -n "${sample_dir}/qc/${sample}_nanostat.txt" -t "$THREADS" >>"$LOGFILE" 2>&1
    rm -f "$fq"
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- STEP 3: Kraken2 ----------------
    step="kraken2"; step_start=$(date +%s)
    log "▶ Step $step"
    mkdir -p ${sample_dir}/kraken2
    kraken_out="${sample_dir}/kraken2/kraken2_output.tsv"
    kraken_report="${sample_dir}/kraken2/kraken2_report.tsv"
    if ! kraken2 --db "$KRAKEN2_DB" --threads "$THREADS" --output "$kraken_out" --report "$kraken_report" --use-names "${sample_dir}/trimmed/${sample}.fastq.gz" >> "$LOGFILE" 2>&1; then
        log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; continue
    fi
    cd ${sample_dir}/
    [[ -f "$PY_KRAKEN_SUMMARY" ]] && python3 "$PY_KRAKEN_SUMMARY"
    cd ..
    kraken_summary="${sample_dir}/kraken2_top5_taxa.csv"
    if [[ -f "$kraken_summary" ]]; then
        top_hit=$(awk -F, 'NR>1 && $1 {print $2; exit}' "$kraken_summary")
        if [[ "$top_hit" != "$organism" ]]; then
            log "❌ Species mismatch: expected $organism, got $top_hit"; status="FAILED"; failed_step="species_mismatch"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; continue
        else
                log "✅ Species match confirmed for $sample: $top_hit (expected: $organism)"
            fi
        else
        log "⚠️ Kraken summary file not found for $sample"; continue
    fi
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- STEP 4: Read-based MLST ----------------
    step="mlst_reads"; step_start=$(date +%s)
    log "▶ Step $step"
    mkdir -p ${sample_dir}/mlst_reads
    rm -f ${sample_dir}/"mlst_reads/krocusmlst.tsv" # Remove any previous runs, baisically an overwrite
    pigz -dc "${sample_dir}/trimmed/${sample}.fastq.gz" > "${sample_dir}/mlst_reads/${sample}.fastq" # decompress fastq
    head -n200000 "${sample_dir}/mlst_reads/${sample}.fastq" > "${sample_dir}/mlst_reads/${sample}_krocus.fastq" # Take first 50,000 reads only
    rm -f "${sample_dir}/mlst_reads/${sample}.fastq" # Remove intermediary file
    read_mlst="${species_dict_read_mlst[$organism]}"
    if ! krocus -o "${sample_dir}/mlst_reads/krocusmlst.tsv" "${KROCUS_DB}/${read_mlst}/" "${sample_dir}/mlst_reads/${sample}_krocus.fastq" 2>>"$LOGFILE"; then
        log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; continue
    fi
    rm -f "${sample_dir}/mlst_reads/${sample}_krocus.fastq" # Remove other intermediary file
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- STEP 5: Assembly ----------------
    step="assembly"; step_start=$(date +%s)
    log "▶ Step $step"
    outdir="${sample_dir}/assemblies"; mkdir -p "$outdir"
    size="${species_dict_size[$organism]}"
    if ! autoautocycler.sh -o "$outdir/" -t "$THREADS" -s "$size" -a "flye raven" "${sample_dir}/trimmed/${sample}.fastq.gz" 2>>"$LOGFILE"; then
        log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; continue
    fi
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- STEP 6: Coverage depth  ----------------
    step="mapping_depth"; step_start=$(date +%s)
    log "▶ Step $step"
    asm="${sample_dir}/assemblies/${sample}.fasta"; depth="${sample_dir}/qc/genome_coverage_depth.txt"
    if ! minimap2 -x map-ont --secondary no -a -t "$THREADS" "$asm" "${sample_dir}/trimmed/${sample}.fastq.gz" | samtools sort -O BAM | samtools coverage - | cut -f 7 | head -n 2 > "$depth" 2>>"$LOGFILE"; then
        log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; continue
    fi
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- STEP 7: QUAST ----------------
    step="quast"; step_start=$(date +%s)
    log "▶ Step $step"
    if ! quast -o "${sample_dir}/qc/quast" -t "$THREADS" -m 100 -l "$sample" "$asm" >>"$LOGFILE" 2>&1; then
        log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; continue
    fi
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- STEP 8: Bandage ----------------
    step="bandage"; step_start=$(date +%s)
    log "▶ Step $step"
    gfa="${sample_dir}/assemblies/${sample}.gfa"
    if ! [[ -f "$gfa" ]] && Bandage image "$gfa" "${sample_dir}/qc/assembly_image.svg" 2>>"$LOGFILE"; then
        log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; continue
    fi
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- STEP 9: MLST contigs ----------------
    step="mlst_contigs"; step_start=$(date +%s)
    log "▶ Step $step"
    scheme="${species_dict_mlst_contigs[$organism]}"
    if ! conda run -n "$ENV_MLST" mlst --scheme "$scheme" "$asm" > "${sample_dir}/contig_mlst.tsv" 2>>"$LOGFILE"; then
        log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; continue
    fi
    sequence_type=$(awk '{print $3; exit}' "${sample_dir}/contig_mlst.tsv")
    log "Multi-locus Sequence Type found for $top_hit: ST${sequence_type}" 
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- STEP 10: abritamr (AMRFINDER) ----------------
    step="amrfinder"; step_start=$(date +%s)
    log "▶ Step $step"
    amr_species="${species_dict_amrfinder[$organism]}"
    if ! conda run -n "$ENV_ABRITAMR" abritamr run -c "${sample_dir}/assemblies/${sample}.fasta" -px "$sample_dir" -j "$THREADS" -sp "$amr_species" 2>>"$LOGFILE"; then
        log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; continue
    fi
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- STEP 11: BUSCO ----------------
    step="busco"; step_start=$(date +%s)
    log "▶ Step $step"
    mkdir -p "${sample_dir}/qc/busco_results"
    lineage="${species_dict_busco[$organism]}"
    if ! busco -i "$asm" -f -m genome -l "$lineage" -c "$THREADS" --out_path "${sample_dir}/qc/busco_results" -o "$sample" 2>>"$LOGFILE"; then
        log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; continue
    fi
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- STEP 12: MultiQC ----------------
    step="multiqc"; step_start=$(date +%s)
    log "▶ Step $step"
    if ! multiqc -i "$sample" -f -v -o ${sample_dir}/multiqc_report/ --no-ai -c "$MULTIQC_CONFIG" "$sample_dir"/ >>"$LOGFILE" 2>&1; then
        log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; continue
    fi
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- STEP 13: Copy genomes for ST specifc tree ----------------
    step="copy_genomes"; step_start=$(date +%s)
    log "▶ Step $step"
    species="${species_dict_read_mlst[$organism]}"
    species_tree_dir="${TREE_BASE}/${species}_trees"
    mlst_tree_dir="$species_tree_dir/ST${sequence_type}"
    mkdir -p "$species_tree_dir"
    mkdir -p "$mlst_tree_dir"
    cp "$asm" "$mlst_tree_dir/${sample}.fasta" 2>>"$LOGFILE" || log "⚠️ copy fasta failed"
    cp "${sample_dir}/trimmed/${sample}.fastq.gz" "$mlst_tree_dir/${sample}.fastq.gz" 2>>"$LOGFILE" || log "⚠️ copy fastq failed"
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- Sample completion ----------------
    end_total=$(date +%s)
    total_time=$((end_total-start_total))
    echo -e "${sample}\t${status}\t${failed_step}\t${total_time}" >> "${MASTER_SUMMARY}"
    log "✅ Sample $sample completed in ${total_time}s"

done

log "Pipeline complete. Master summary saved in $MASTER_SUMMARY"
