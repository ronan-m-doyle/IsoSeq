#!/usr/bin/env bash
# -------------------------------------------------------------
# Refactored Isolate Sequencing Pipeline (stable logging version)
# -------------------------------------------------------------
IFS=$'\n\t'

# -------------------------- CONFIG ----------------------------
KRAKEN2_DB="/home/ronan/reference_genomes/kraken2_db/"
PY_KRAKEN_SUMMARY="/home/ronan/useful_scripts/kraken2_top5_species.py"
MULTIQC_CONFIG="$HOME/isolate_sequencing/multiqc_config.yaml"
TREE_BASE="/home/ronan/isolate_sequencing"

# -------------------------- LOGGING ----------------------------
LOGFILE="pipeline.log"
SUMMARY_FILE="pipeline_summary.tsv"
log() { echo -e "[$(date '+%F %T')] $*" | tee -a "$LOGFILE"; }

# -------------------------- ARG PARSE --------------------------
run_directory=""
output_directory=""
SAMPLESHEET=""
THREADS=20

while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --run_directory) run_directory="$2"; shift; shift ;;
        --output_directory) output_directory="$2"; shift; shift ;;
        --samplesheet) SAMPLESHEET="$2"; shift; shift ;;
        --threads) THREADS="$2"; shift; shift ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

if [[ -z "$run_directory" || -z "$output_directory" || -z "$SAMPLESHEET" ]]; then
    log "❌ ERROR: Must specify --run_directory, --output_directory, --samplesheet"
    exit 1
fi

run_directory=$(realpath "$run_directory")
output_directory=$(realpath "$output_directory")
SAMPLESHEET=$(realpath "$SAMPLESHEET")

[[ -f "$SAMPLESHEET" ]] || { log "❌ ERROR: Samplesheet not found"; exit 1; }

mkdir -p "$output_directory"
cd "$output_directory" || { log "❌ Cannot cd to $output_directory"; exit 1; }

log "Pipeline started:"
log " Run dir:  $run_directory"
log " Output:   $output_directory"
log " Samplesheet: $SAMPLESHEET"
log " Threads:  $THREADS"

# -------------------------- ENVS -------------------------------
ENV_MLST="mlst_contigs"
ENV_PROKKA="prokka"
ENV_AMRFINDER="amrfinder"

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

declare -A species_dict_prokka=(
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
        log "❌ Unknown organism: $organism"; exit 1
    fi
    BARCODES+=("$barcode"); SAMPLES+=("$sample"); ORGANISMS+=("$organism")
done < "$SAMPLESHEET"

declare -A sample_status sample_step sample_time_start sample_time_end

# ============================================================== 
# PER-SAMPLE PIPELINE 
# ============================================================== 
for idx in "${!SAMPLES[@]}"; do
    sample="${SAMPLES[$idx]}"
    barcode="${BARCODES[$idx]}"
    organism="${ORGANISMS[$idx]}"
    log "-------------------------------------------------------------"
    log "▶ Processing sample: $sample ($organism)"
    start_total=$(date +%s)
    sample_status["$sample"]="OK"
    sample_step["$sample"]="started"

    # ---------- Step 1: Concatenate ----------
    step_start=$(date +%s)
    src_pattern="$run_directory/*/*/fastq_pass/barcode${barcode}/*.fastq.gz"
    out_file="${sample}.fastq.gz"
    shopt -s nullglob; files=( $src_pattern )
    if (( ${#files[@]} == 0 )); then
        log "⚠️ No FASTQ for barcode${barcode}"
        sample_status["$sample"]="FAILED"; sample_step["$sample"]="concat"; continue
    fi
    log "Concatenating ${#files[@]} → $out_file"
    if ! cat "${files[@]}" > "$out_file" 2>>"$LOGFILE"; then
        log "❌ concat failed"; sample_status["$sample"]="FAILED"; sample_step["$sample"]="concat"; continue
    fi
    step_end=$(date +%s); log "✅ Step concat done in $((step_end-step_start))s"

    # ---------- Step 2: Porechop and Nanostat ----------
    mkdir -p trimmed qc
    step_start=$(date +%s)
    fq="${sample}.fastq.gz"; out_trim="trimmed/${sample}.fastq.gz"
    if ! porechop -i "$fq" -o "$out_trim" -t "$THREADS" --no_split 2>>"$LOGFILE"; then
        log "❌ Porechop failed"; sample_status["$sample"]="FAILED"; sample_step["$sample"]="porechop"; continue
    fi
    NanoStat --fastq "$out_trim" -n "qc/${sample}_nanostat.txt" -t "$THREADS" >>"$LOGFILE" 2>&1
    rm -f "$fq"
    step_end=$(date +%s); log "✅ Step trim done in $((step_end-step_start))s"

    # ---------- Step 3: Kraken2 ----------
    mkdir -p kraken2
    step_start=$(date +%s)
    kraken_out="kraken2/${sample}.k2output"
    kraken_report="kraken2/${sample}.k2report"
    if ! kraken2 --db "$KRAKEN2_DB" --threads "$THREADS" \
        --output "$kraken_out" --report "$kraken_report" --use-names "trimmed/${sample}.fastq.gz" >>"$LOGFILE" 2>&1; then
        log "❌ Kraken2 failed"
        sample_status["$sample"]="FAILED"
        sample_step["$sample"]="kraken2"
        continue
    fi

    # ---------- Step 3.1: Kraken species check ----------
    [[ -f "$PY_KRAKEN_SUMMARY" ]] && python3 "$PY_KRAKEN_SUMMARY" 2>>"$LOGFILE"

    kraken_summary="kraken2_top5_per_sample.csv"
    if [[ -f "$kraken_summary" ]]; then
        # Get top hit (first line for that sample)
        top_hit=$(awk -F, -v s="$sample" 'NR>1 && $1==s {print $2; exit}' "$kraken_summary")
        if [[ -z "$top_hit" ]]; then
            log "⚠️ No top hit found in Kraken summary for $sample"
        else
            log "🔍 Kraken top hit for $sample: $top_hit (expected: $organism)"
            if [[ "$top_hit" != "$organism" ]]; then
                log "❌ Species mismatch for $sample — expected '$organism' but found '$top_hit'"
                sample_status["$sample"]="FAILED"
                sample_step["$sample"]="species_mismatch"
                # Skip remaining steps for this sample
                continue
            else
                log "✅ Species match confirmed for $sample"
            fi
        fi
    else
        log "⚠️ Kraken summary file not found for $sample"
    fi

    step_end=$(date +%s)
    log "✅ Step kraken2 done in $((step_end-step_start))s"

    # ---------- Step 4: Read-based MLST ----------
    mkdir -p mlst_reads
    rm -f "mlst_reads/${sample}_krocusmlst.tsv"
    step_start=$(date +%s)
    read_mlst="${species_dict_read_mlst[$organism]}"
    if ! krocus -o "mlst_reads/${sample}_krocusmlst.tsv" "/home/ronan/isolate_sequencing/mlst_krocus/${read_mlst}/" "trimmed/${sample}.fastq.gz" 2>>"$LOGFILE"; then
        log "❌ Krocus failed"; sample_status["$sample"]="FAILED"; sample_step["$sample"]="krocus"; continue
    fi
    step_end=$(date +%s); log "✅ Step krocus done in $((step_end-step_start))s"

    # ---------- Step 5: Assembly ----------
    mkdir -p assemblies
    step_start=$(date +%s)
    outdir="assemblies/${sample}"; mkdir -p "$outdir"
    size="${species_dict_size[$organism]}"
    if ! autoautocycler.sh -o "$outdir/" -t "$THREADS" -s "$size" -a "flye raven" "trimmed/${sample}.fastq.gz" 2>>"$LOGFILE"; then
        log "❌ Assembly failed"; sample_status["$sample"]="FAILED"; sample_step["$sample"]="assembly"; continue
    fi
    cp "$outdir/${sample}.fasta" "assemblies/${sample}.fasta" 2>/dev/null || true
    cp "$outdir/${sample}.gfa" "assemblies/${sample}.gfa" 2>/dev/null || true
    step_end=$(date +%s); log "✅ Step assembly done in $((step_end-step_start))s"

    # ---------- Step 6: Mapping ----------
    step_start=$(date +%s)
    asm="assemblies/${sample}.fasta"; bam="assemblies/${sample}.bam"
    if ! minimap2 -x map-ont --secondary no -a -t "$THREADS" "$asm" "trimmed/${sample}.fastq.gz" 2>>"$LOGFILE" \
        | samtools sort -O BAM -o "$bam" 2>>"$LOGFILE"; then
        log "❌ Mapping failed"; sample_status["$sample"]="FAILED"; sample_step["$sample"]="mapping"; continue
    fi
    samtools index -M "$bam" 2>>"$LOGFILE"
    mkdir -p "qc/${sample}"
    qualimap bamqc -bam "$bam" -outdir "qc/${sample}" --java-mem-size=16G 2>>"$LOGFILE"
    step_end=$(date +%s); log "✅ Step mapping done in $((step_end-step_start))s"

    # ---------- Step 7: QUAST ----------
    step_start=$(date +%s)
    quast -o "qc/${sample}_quast" -t "$THREADS" -m 100 -l "$sample" "$asm" >>"$LOGFILE" 2>&1
    step_end=$(date +%s); log "✅ Step quast done in $((step_end-step_start))s"

    # ---------- Step 8: Bandage ----------
    step_start=$(date +%s)
    gfa="assemblies/${sample}.gfa"
    [[ -f "$gfa" ]] && Bandage image "$gfa" "qc/${sample}.svg" 2>>"$LOGFILE"
    step_end=$(date +%s); log "✅ Step bandage done in $((step_end-step_start))s"

    # ---------- Step 9: MLST contigs ----------
    step_start=$(date +%s)
    scheme="${species_dict_mlst_contigs[$organism]}"
    conda run -n "$ENV_MLST" mlst --scheme "$scheme" "$asm" > "${sample}_contigmlst.tsv" 2>>"$LOGFILE"
    step_end=$(date +%s); log "✅ Step mlst_contigs done in $((step_end-step_start))s"

    # ---------- Step 10: Prokka ----------
    step_start=$(date +%s)
    prokka_species="${species_dict_prokka[$organism]}"
    if ! conda run -n "$ENV_PROKKA" prokka --outdir "assemblies/${sample}/" --force \
        --species "$prokka_species" --prefix "$sample" --strain "$sample" \
        --kingdom Bacteria --cpus "$THREADS" --gcode 11 "$asm" 2>>"$LOGFILE"; then
        log "❌ Prokka failed"; sample_status["$sample"]="FAILED"; sample_step["$sample"]="prokka"; continue
    fi
    step_end=$(date +%s); log "✅ Step prokka done in $((step_end-step_start))s"

    # ---------- Step 11: AMRFinder ----------
    step_start=$(date +%s)
    amr_species="${species_dict_prokka[$organism]}"
    conda run -n "$ENV_AMRFINDER" amrfinder -n "$asm" -O "$amr_species" \
        --threads "$THREADS" -o "${sample}_resistance.tsv" \
        -g "assemblies/${sample}/${sample}.gff" -p "assemblies/${sample}/${sample}.faa" -a prokka 2>>"$LOGFILE"
    step_end=$(date +%s); log "✅ Step amrfinder done in $((step_end-step_start))s"

    # ---------- Step 12: BUSCO ----------
    step_start=$(date +%s)
    mkdir -p qc/busco_results
    lineage="${species_dict_busco[$organism]}"
    if ! busco -i "$asm" -f -m genome -l "$lineage" -c "$THREADS" \
        --out_path "qc/busco_results" -o "$sample" 2>>"$LOGFILE"; then
        log "❌ BUSCO failed"; sample_status["$sample"]="FAILED"; sample_step["$sample"]="busco"; continue
    fi
    rm -rf busco_downloads/ 2>>"$LOGFILE"
    step_end=$(date +%s); log "✅ Step busco done in $((step_end-step_start))s"

    # ---------- Step 14: Copy genomes ----------
    step_start=$(date +%s)
    species="${species_dict_read_mlst[$organism]}"
    tree_dir="${TREE_BASE}/${species}_tree"
    mkdir -p "$tree_dir"
    cp "$asm" "$tree_dir/${sample}.fasta" 2>>"$LOGFILE" || log "⚠️ copy fasta failed"
    cp "trimmed/${sample}.fastq.gz" "$tree_dir/${sample}.fastq.gz" 2>>"$LOGFILE" || log "⚠️ copy fastq failed"
    step_end=$(date +%s); log "✅ Step copy done in $((step_end-step_start))s"

    # ---------- Sample complete ----------
    end_total=$(date +%s)
    sample_time_start["$sample"]=$start_total
    sample_time_end["$sample"]=$end_total
    log "✅ Sample $sample completed in $((end_total-start_total))s"
done

# ==============================================================
# Step 13: MultiQC (run once at end)
# ==============================================================
log "STEP 13: MultiQC report"
mkdir -p multiqc_report
currentPWD=$(basename "$output_directory")
multiqc -i "$currentPWD" -f -v -o multiqc_report/ --no-ai -c "$MULTIQC_CONFIG" . >>"$LOGFILE" 2>&1

# ==============================================================
# FINAL SUMMARY
# ==============================================================
log "-------------------------------------------------------------"
log "Pipeline summary:"
echo -e "Sample\tStatus\tFailed_Step\tTotal_Time(s)" > "$SUMMARY_FILE"
for sample in "${SAMPLES[@]}"; do
    status="${sample_status[$sample]:-OK}"
    step="${sample_step[$sample]:--}"
    if [[ "$status" == "OK" ]]; then
        total=$(( sample_time_end[$sample]-sample_time_start[$sample] ))
        log "✅ $sample — OK — ${total}s"
        echo -e "${sample}\t${status}\t${step}\t${total}" >> "$SUMMARY_FILE"
    else
        log "❌ $sample — failed at step ${step}"
        echo -e "${sample}\tFAILED\t${step}\t-" >> "$SUMMARY_FILE"
    fi
done
log "Summary written to $SUMMARY_FILE"
log "Pipeline finished at $(date '+%F %T')"
