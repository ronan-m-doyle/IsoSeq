#!/usr/bin/env bash
# ----------------------------------------------------------------------
# Isolate Sequencing Pipeline v1.0 - 72h stage - (Assembly, QC, MLST, AMR)
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

    sample_dir="/data/IsoSeq_results/${sample}/72h"
    mkdir -p "${sample_dir}/logs"
    LOGFILE="${sample_dir}/logs/pipeline.log"

    mkdir -p "${sample_dir}/samplesheet"
    cp -v "$SAMPLESHEET" "${sample_dir}/samplesheet/samplesheet.csv"
    touch "${sample_dir}/samplesheet/${RUN_NAME}"

    log() { echo -e "[$(date '+%F %T')] $*"; }
    exec 3>&1 4>&2
    exec > >(tee -a "$LOGFILE") 2>&1
    trap 'exec 1>&3 2>&4; exec 3>&- 4>&-' RETURN

    if [[ -f "/data/IsoSeq_results/${sample}/8h/kraken2_top5_taxa.csv" ]]; then
        top_hit=$(awk -F, 'NR>1 && $1 {print $2; exit}' "/data/IsoSeq_results/${sample}/8h/kraken2_top5_taxa.csv")
        if [[ "$organism" == "Unknown" ]]; then
            log "ℹ️ Organism marked Unknown — skipping 8h species consistency check (8h Kraken2 top hit: $top_hit)"
        elif [[ "$top_hit" != "$organism" ]]; then
            log "8h analysis failed for $sample - sample skipped"; status="FAILED"; failed_step="8h_species_mismatch"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; return
            fi
    fi


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
    ontime --from "48h1ns" --to "72h" -o "${sample_dir}/${sample}.fastq.gz" "${sample_dir}/${sample}_raw.fastq.gz"  # Use ontime to filter for only reads for this step
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
    cat "${sample_dir}/BC${barcode}.fastq.gz" "/data/IsoSeq_results/${sample}/24h/trimmed/${sample}.fastq.gz" > "$out_trim" # concat prev timepoint trimmed with latest
    if ! NanoStat --fastq "$out_trim" -n "${sample_dir}/qc/${sample}_nanostat.txt" -t "$THREADS"; then
        log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; return
    fi
    rm -vf "${sample_dir}/BC${barcode}.fastq.gz" # Clean up leftover fastq
    rm -vf "${sample_dir}/${sample}.fastq.gz" # Clean up leftover fastq
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"
    # ---------------- STEP 3: Assembly ----------------
    step="assembly"; step_start=$(date +%s)
    log "▶ Step $step"
    outdir="${sample_dir}/assemblies"; mkdir -p "$outdir"
    if [[ "$organism" == "Unknown" ]]; then
        log "ℹ️ Organism marked Unknown — running assembly without genome size estimate (-s)"
        if ! bash ${SCRIPT_DIR}/scripts/autoautocycler.sh -o "$outdir/" -t "$THREADS" -c "2" -a "metamdbg myloasm" "${sample_dir}/trimmed/${sample}.fastq.gz"; then
            log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; return
        fi
    else
        size="${species_dict_size[$organism]}"
        if ! bash ${SCRIPT_DIR}/scripts/autoautocycler.sh -o "$outdir/" -t "$THREADS" -c "2" -s "$size" -a "metamdbg myloasm" "${sample_dir}/trimmed/${sample}.fastq.gz"; then
            log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; return
        fi
    fi
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- STEP 4: Coverage depth  ----------------
    step="mapping_depth"; step_start=$(date +%s)
    log "▶ Step $step"
    asm="${sample_dir}/assemblies/${sample}.fasta"; depth="${sample_dir}/qc/genome_coverage_depth.txt"
    if ! minimap2 -x lr:hq --secondary no -a -t "$THREADS" "$asm" "${sample_dir}/trimmed/${sample}.fastq.gz" | samtools sort -O BAM | samtools coverage - | cut -f 7 | head -n 2 > "$depth"; then
        log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; return
    fi
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- STEP 5: QUAST ----------------
    step="quast"; step_start=$(date +%s)
    log "▶ Step $step"
    if ! quast -o "${sample_dir}/qc/quast" -t "$THREADS" -m 100 -l "$sample" "$asm"; then
        log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; return
    fi
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- STEP 6: Bandage ----------------
    step="bandage"; step_start=$(date +%s)
    log "▶ Step $step"
    gfa="${sample_dir}/assemblies/${sample}/autocycler_out/consensus_assembly.gfa"
    if ! Bandage image "$gfa" "${sample_dir}/qc/assembly_image.svg"; then
        log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; return
    fi
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- STEP 7: BUSCO ----------------
    step="busco"; step_start=$(date +%s)
    if [[ "$organism" == "Unknown" ]]; then
        log "⏭️  Skipping Step $step (organism marked Unknown — no BUSCO lineage to use)"
    else
        log "▶ Step $step"
        mkdir -p "${sample_dir}/qc/busco_results"
        lineage="${species_dict_busco[$organism]}"
        if ! busco -i "$asm" -f -m genome -l "$lineage" -c "$THREADS" --out_path "${sample_dir}/qc/busco_results" -o "$sample"; then
            log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; return
        fi
    fi
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- STEP 8: Chroquetas Fungal AMR prediction ----------------
    step="chroquetas"; step_start=$(date +%s)
    log "▶ Step $step"
    if [[ "$organism" == "Unknown" ]]; then
        log "⏭️  Skipping Step $step (organism marked $organism — no reference to use)"
    else
        amr_species="${species_dict_amrfinder[$organism]}"
        if ! conda run -n "$ENV_CHROQUETAS" ChroQueTas.sh -g "${sample_dir}/assemblies/${sample}.fasta" -s "$amr_species" -o "${sample_dir}/amr_output/" -t "$THREADS"; then
            log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; continue
        fi
    cp -v "${sample_dir}/amr_output/${sample}.ChroQueTaS.AMR_summary.txt" "${sample_dir}/amr_summary.txt"
    fi
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- STEP 9: Auriclass Candidozyma auris clade prediction ----------------
    step="auriclass"; step_start=$(date +%s)
    if [[ "$organism" == "Unknown" ]]; then
        log "⏭️  Skipping Step $step (organism marked $organism — reference to use)"
    else
        log "▶ Step $step"
        if ! auriclass -n "$sample" -o "${sample_dir}/auriclass_report.tsv" --log_file_path "${sample_dir}/logs/auriclass_log.txt" --verbose --fasta "${sample_dir}/assemblies/${sample}.fasta"; then
            log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; return
        fi
        clade=$(awk -F'\t' 'NR>1 {print $2}' "${sample_dir}/auriclass_report.tsv")
        log "C. auris clade: ${clade}"
        step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"
    fi

    # ---------------- Stage complete ----------------
    end_total=$(date +%s)
    total_time=$((end_total-start_total))
    status="72h COMPLETE"
    echo -e "${sample}\t${status}\t${failed_step}\t${total_time}" >> "${MASTER_SUMMARY}"
    log "✅ Sample $sample 72h analysis completed in ${total_time}s"
}

for idx in "${!SAMPLES[@]}"; do
    process_sample "$idx"
done
