#!/usr/bin/env bash
# ----------------------------------------------------------------------
# Isolate Candida Sequencing Pipeline v1.0 - 8h stage - (Initial read-based Typing and AMR)
# ----------------------------------------------------------------------
IFS=$'\n\t'

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/isoseq_pipeline_config.sh"

SAMPLESHEET="$1"
THREADS="$2"
RUN_DIR="$3"
MASTER_SUMMARY="$4"
RUN_NAME="$5"

# ---------------------- READ SAMPLESHEET -----------------------
SAMPLES=(); BARCODES=(); ORGANISMS=(); COLLECTIONDATE=()
while IFS=, read -r barcode sample organism collection || [[ -n "$barcode" ]]; do
    [[ "$barcode" == "barcode" || -z "$barcode" ]] && continue
    # "Unknown" is an explicitly allowed value and bypasses the dict lookup below.
    # Any other unrecognized organism still hard-fails the run.
    if [[ "$organism" != "Unknown" && -z "${species_dict_size[$organism]:-}" ]]; then
        echo "❌ Unknown organism: $organism"; exit 1
    fi
    SAMPLES+=("$sample")
    BARCODES+=("$barcode")
    ORGANISMS+=("$organism")
    COLLECTIONDATE+=("$collection")
done < "$SAMPLESHEET"

# ============================================================== 
# PER-SAMPLE PIPELINE 
# ============================================================== 

process_sample() {
    local idx="$1"
    sample="${SAMPLES[$idx]}"
    barcode="${BARCODES[$idx]}"
    organism="${ORGANISMS[$idx]}"
    collection="${COLLECTIONDATE[$idx]}"

    sample_dir="/data/IsoSeq_results/${sample}/8h"
    mkdir -p "${sample_dir}/logs"
    LOGFILE="${sample_dir}/logs/pipeline.log"

    mkdir -p "${sample_dir}/samplesheet"
    cp -v "$SAMPLESHEET" "${sample_dir}/samplesheet/samplesheet.csv"
    touch "${sample_dir}/samplesheet/${RUN_NAME}"

    log() { echo -e "[$(date '+%F %T')] $*"; }
    exec 3>&1 4>&2
    exec > >(tee -a "$LOGFILE") 2>&1
    trap 'exec 1>&3 2>&4; exec 3>&- 4>&-' RETURN

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
        return
    fi
    log "Concatenating ${#files[@]} fastq files → ${sample_dir}/${sample}.fastq.gz"
    if ! printf '%s\0' "${files[@]}" | xargs -0 cat > "${sample_dir}/${sample}_raw.fastq.gz"; then
        log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; return
    fi
    ontime --to 8h -o "${sample_dir}/${sample}.fastq.gz" "${sample_dir}/${sample}_raw.fastq.gz"  # Use ontime to filter for only reads for this step
    rm -vf "${sample_dir}/${sample}_raw.fastq.gz"
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- STEP 2: Porechop & NanoStat ----------------
    step="porechop_nanostat"; step_start=$(date +%s)
    log "▶ Step $step"
    mkdir -p ${sample_dir}/trimmed ${sample_dir}/qc
    fq="${sample_dir}/${sample}.fastq.gz"; out_trim="${sample_dir}/trimmed/${sample}.fastq.gz"
    if ! python ${SCRIPT_DIR}/Porechop/porechop-runner.py -i "$fq" -b "${sample_dir}" -t "$THREADS" --no_split --discard_unassigned; then
        log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; return
    fi
    cp -v "${sample_dir}/BC${barcode}.fastq.gz" "$out_trim" # Copy porechop trimmed and demultiplexed fastq to trimmed dir
    if ! NanoStat --fastq "$out_trim" -n "${sample_dir}/qc/${sample}_nanostat.txt" -t "$THREADS"; then
        log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; return
    fi
    rm -vf "${sample_dir}/BC${barcode}.fastq.gz" # Clean up leftover fastq
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- STEP 3: Kraken2 ----------------
    step="kraken2"; step_start=$(date +%s)
    log "▶ Step $step"
    mkdir -p ${sample_dir}/kraken2
    kraken_out="${sample_dir}/kraken2/kraken2_output.tsv"
    kraken_report="${sample_dir}/kraken2/kraken2_report.tsv"
    if ! kraken2 --db "$KRAKEN2_DB" --threads "$THREADS" --output "$kraken_out" --report "$kraken_report" --use-names "$out_trim"; then
        log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; return
    fi
    python3 "$PY_KRAKEN_SUMMARY" --input "${sample_dir}/kraken2/kraken2_report.tsv" --output "${sample_dir}/kraken2_top5_taxa.csv"
    kraken_summary="${sample_dir}/kraken2_top5_taxa.csv"
    if [[ -f "$kraken_summary" ]]; then
        top_hit=$(awk -F, 'NR>1 && $1 {print $2; exit}' "$kraken_summary")
        if [[ "$organism" == "Unknown" ]]; then
            log "ℹ️ Organism marked Unknown — skipping species-match check (Kraken2 top hit: $top_hit)"
        elif [[ "$top_hit" != "$organism" ]]; then
            log "❌ Species mismatch: expected $organism, got $top_hit"; status="FAILED"; failed_step="species_mismatch"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; return
        else
            log "✅ Species match confirmed for $sample: $top_hit (expected: $organism)"

            # ---- Extra check: require top-hit read fraction > 15% ----
            top_percent=$(awk -F, -v hit="$top_hit" 'NR>1 && $2==hit {print $3; exit}' "$kraken_summary")
            if [[ -z "$top_percent" ]]; then
                log "⚠️ Could not determine read fraction for $top_hit"
                status="FAILED"; failed_step="kraken2_fraction_missing"
                echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"
                return
            fi
            if ! awk -v p="$top_percent" 'BEGIN{exit !(p > 15)}'; then
                log "❌ Read fraction for $top_hit too low: ${top_percent}% (threshold: >15%)"
                status="FAILED"; failed_step="kraken2_low_fraction"
                echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"
                return
            fi
            log "✅ Read fraction check passed for $top_hit: ${top_percent}% (>15%)"
        fi
    else
        log "⚠️ Kraken summary file not found for $sample"; status="FAILED"; failed_step="kraken2_summary_missing"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; return
    fi
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- STEP 4: Auriclass Candidozyma auris clade prediction ----------------
    step="auriclass"; step_start=$(date +%s)
    if [[ "$organism" == "Unknown" ]]; then
        log "⏭️  Skipping Step $step (organism marked $organism — reference to use)"
    else
        log "▶ Step $step"
        if ! auriclass -n "$sample" -o "${sample_dir}/auriclass_report.tsv" --log_file_path "${sample_dir}/logs/auriclass_log.txt" --verbose --fastq "$out_trim"; then
            log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; return
        fi
        clade=$(awk -F'\t' 'NR>1 {print $2}' "${sample_dir}/auriclass_report.tsv")
        log "C. auris clade: ${clade}"
        step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"
    fi

    # ---------------- STEP 5: Chroquetas Fungal AMR prediction ----------------
    step="chroquetas"; step_start=$(date +%s)
    if [[ "$organism" == "Unknown" ]]; then
        log "⏭️  Skipping Step $step (organism marked $organism — no reference to use)"
    else
        log "▶ Step $step"
        amr_species="${species_dict_amrfinder[$organism]}"
        if ! conda run -n "$ENV_CHROQUETAS" ChroQueTas.sh -g "$out_trim" -s "$amr_species" -o "${sample_dir}/amr_output/" -t "$THREADS"; then
            log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; continue
        fi
    cp -v "${sample_dir}/amr_output/${sample}.ChroQueTaS.AMR_summary.txt" "${sample_dir}/amr_summary.txt"
    fi
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- Stage complete ----------------
    end_total=$(date +%s)
    total_time=$((end_total-start_total))
    status="8h COMPLETE"
    echo -e "${sample}\t${status}\t${failed_step}\t${total_time}" >> "${MASTER_SUMMARY}"
    log "✅ Sample $sample 8h analysis completed in ${total_time}s"
}

for idx in "${!SAMPLES[@]}"; do
    process_sample "$idx"
done