#!/usr/bin/env bash
# ----------------------------------------------------------------------
# Isolate Sequencing Pipeline - 8h stage - (Initial read-based MLST and AMR)
# ----------------------------------------------------------------------
IFS=$'\n\t'

source "$(dirname "$0")/isoseq_pipeline_config.sh"

SAMPLESHEET="$1"
THREADS="$2"
RUN_DIR="$3"
MASTER_SUMMARY="$4"
RUN_NAME="$5"

# ---------------------- READ SAMPLESHEET -----------------------
SAMPLES=(); BARCODES=(); ORGANISMS=(); COLLECTIONDATE=()
while IFS=, read -r barcode sample organism collection; do
    [[ "$barcode" == "barcode" || -z "$barcode" ]] && continue
    if [[ -z "${species_dict_size[$organism]:-}" ]]; then
        echo "❌ Unknown organism: $organism"; exit 1
    fi
    SAMPLES+=("$sample")
    BARCODES+=("$barcode")
    ORGANISMS+=("$organism")
    COLLECTIONS+=("$collection")
done < "$SAMPLESHEET"

# ============================================================== 
# PER-SAMPLE PIPELINE 
# ============================================================== 

for idx in "${!SAMPLES[@]}"; do 
{
    sample="${SAMPLES[$idx]}"
    barcode="${BARCODES[$idx]}"
    organism="${ORGANISMS[$idx]}"
    collection="${COLLECTIONS[$idx]}"

    sample_dir="results/${sample}/8h"
    mkdir -p "${sample_dir}/logs"
    LOGFILE="${sample_dir}/logs/pipeline.log"

    mkdir -p "${sample_dir}/samplesheet"
    cp -v "$SAMPLESHEET" "${sample_dir}/samplesheet/samplesheet.csv"
    touch "${sample_dir}/samplesheet/${RUN_NAME}"

    log() { echo -e "[$(date '+%F %T')] $*"; }
    exec 3>&1 4>&2
    exec > >(tee -a "$LOGFILE") 2>&1

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
    if ! cat "${files[@]}" > "${sample_dir}/${sample}_raw.fastq.gz"; then
       log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; continue
    fi
    ontime --to 8h -o "${sample_dir}/${sample}.fastq.gz" "${sample_dir}/${sample}_raw.fastq.gz"  # Use ontime to filter for only reads for this step
    rm -vf "${sample_dir}/${sample}_raw.fastq.gz"
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- STEP 2: Porechop & NanoStat ----------------
    step="porechop_nanostat"; step_start=$(date +%s)
    log "▶ Step $step"
    mkdir -p ${sample_dir}/trimmed ${sample_dir}/qc
    fq="${sample_dir}/${sample}.fastq.gz"; out_trim="${sample_dir}/trimmed/${sample}.fastq.gz"
    if ! porechop -i "$fq" -o "$out_trim" -t "$THREADS" --no_split > /dev/tty; then
        log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; continue
    fi
    NanoStat --fastq "$out_trim" -n "${sample_dir}/qc/${sample}_nanostat.txt" -t "$THREADS"
    rm -vf "$fq"
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- STEP 3: Kraken2 ----------------
    step="kraken2"; step_start=$(date +%s)
    log "▶ Step $step"
    mkdir -p ${sample_dir}/kraken2
    kraken_out="${sample_dir}/kraken2/kraken2_output.tsv"
    kraken_report="${sample_dir}/kraken2/kraken2_report.tsv"
    if ! kraken2 --db "$KRAKEN2_DB" --threads "$THREADS" --output "$kraken_out" --report "$kraken_report" --use-names "${sample_dir}/trimmed/${sample}.fastq.gz"; then
        log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; continue
    fi
    python3 "$PY_KRAKEN_SUMMARY" --input "${sample_dir}/kraken2/kraken2_report.tsv" --output "${sample_dir}/kraken2_top5_taxa.csv"
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
    rm -vf ${sample_dir}/"mlst_reads/krocusmlst.tsv" # Remove any previous runs, baisically an overwrite
    read_mlst="${species_dict_read_mlst[$organism]}"
    if ! krocus -k "19" -o "${sample_dir}/mlst_reads/krocusmlst.tsv" "${KROCUS_DB}/${read_mlst}/" "${sample_dir}/trimmed/${sample}.fastq.gz"; then
        log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; continue
    fi
    tail -n1 "${sample_dir}/mlst_reads/krocusmlst.tsv" > "${sample_dir}/mlst_result.tsv"
    sequence_type=$(awk '{print $1; exit}' "${sample_dir}/mlst_result.tsv")
    log "Multi-locus Sequence Type found for $top_hit: ST${sequence_type}" 
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- STEP 5: abricate for AMR gene prediction ----------------
    step="abricate"; step_start=$(date +%s)
    log "▶ Step $step"
    if ! conda run -n "$ENV_ABRICATE" abricate --threads "$THREADS" "${sample_dir}/trimmed/${sample}.fastq.gz" > "$sample_dir"/amr.tsv; then
        log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; continue
    fi
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- Stage complete ----------------
    end_total=$(date +%s)
    total_time=$((end_total-start_total))
    status="8h COMPLETE"
    echo -e "${sample}\t${status}\t${failed_step}\t${total_time}" >> "${MASTER_SUMMARY}"
    log "✅ Sample $sample 8h analysis completed in ${total_time}s"
    exec 1>&3 2>&4
    exec 3>&- 4>&-
}; done
