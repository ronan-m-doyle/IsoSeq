#!/usr/bin/env bash
# ----------------------------------------------------------------------
# Isolate Sequencing Pipeline
# ----------------------------------------------------------------------
IFS=$'\n\t'

# -------------------------- CONFIG ----------------------------
KRAKEN2_DB="${HOME}/reference_genomes/kraken2_pluspf_db/"
PY_KRAKEN_SUMMARY="${HOME}/useful_scripts/kraken2_top5_species.py"
MULTIQC_CONFIG="${HOME}/isolate_sequencing/multiqc_config.yaml"
TREE_BASE="${HOME}/isolate_sequencing"
KROCUS_DB="${HOME}/isolate_sequencing/mlst_krocus/"
CLAIR3_MODEL="${HOME}/miniforge3/envs/clair3/bin/models/r1041_e82_400bps_sup_v430_bacteria_finetuned/"

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
    ["Klebsiella oxytoca"]="Klebsiella_oxytoca"
    ["Legionella pneumophila"]="Legionella_pneumophila"
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
    ["Klebsiella oxytoca"]="koxytoca"
    ["Legionella pneumophila"]="lpneumophila"
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
    ["Klebsiella oxytoca"]="Klebsiella_oxytoca"
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
    ["Klebsiella oxytoca"]="enterobacteriaceae_odb12"
    ["Legionella pneumophila"]="legionellaceae_odb12"
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
    ["Klebsiella oxytoca"]="5.9m"
    ["Legionella pneumophila"]="3.4m"
)

# ---------------------- READ SAMPLESHEET -----------------------
SAMPLES=(); BARCODES=(); ORGANISMS=()
while IFS=, read -r barcode sample organism; do
    [[ "$barcode" == "barcode" || -z "$barcode" ]] && continue
    if [[ -z "${species_dict_size[$organism]:-}" ]]; then
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
{
    sample="${SAMPLES[$idx]}"
    barcode="${BARCODES[$idx]}"
    organism="${ORGANISMS[$idx]}"

    sample_dir="${sample}"
    mkdir -p "$sample_dir/logs"
    LOGFILE="${sample_dir}/logs/pipeline.log"

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
    if !  cat "${files[@]}" > "${sample_dir}/${sample}.fastq.gz"; then
       log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; continue
    fi
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
    rm -vf ${sample_dir}/"mlst_reads/krocusmlst.tsv" # Remove any previous runs, baisically an overwrite
    pigz -dc "${sample_dir}/trimmed/${sample}.fastq.gz" > "${sample_dir}/mlst_reads/${sample}.fastq" # decompress fastq
    head -n200000 "${sample_dir}/mlst_reads/${sample}.fastq" > "${sample_dir}/mlst_reads/${sample}_krocus.fastq" # Take first 50,000 reads only
    rm -vf "${sample_dir}/mlst_reads/${sample}.fastq" # Remove intermediary file
    read_mlst="${species_dict_read_mlst[$organism]}"
    if ! krocus -k "19" -o "${sample_dir}/mlst_reads/krocusmlst.tsv" "${KROCUS_DB}/${read_mlst}/" "${sample_dir}/mlst_reads/${sample}_krocus.fastq"; then
        log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; continue
    fi
    rm -vf "${sample_dir}/mlst_reads/${sample}_krocus.fastq" # Remove other intermediary file
    tail -n1 "${sample_dir}/mlst_reads/krocusmlst.tsv" > "${sample_dir}/mlst_result.tsv"
    sequence_type=$(awk '{print $1; exit}' "${sample_dir}/mlst_result.tsv")
    log "Multi-locus Sequence Type found for $top_hit: ST${sequence_type}" 
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- STEP 5: Assembly ----------------
    step="assembly"; step_start=$(date +%s)
    log "▶ Step $step"
    outdir="${sample_dir}/assemblies"; mkdir -p "$outdir"
    size="${species_dict_size[$organism]}"
    if ! autoautocycler.sh -o "$outdir/" -t "$THREADS" -s "$size" -a "flye raven" "${sample_dir}/trimmed/${sample}.fastq.gz"; then
        log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; continue
    fi
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- STEP 6: Coverage depth  ----------------
    step="mapping_depth"; step_start=$(date +%s)
    log "▶ Step $step"
    asm="${sample_dir}/assemblies/${sample}.fasta"; depth="${sample_dir}/qc/genome_coverage_depth.txt"
    if ! minimap2 -x map-ont --secondary no -a -t "$THREADS" "$asm" "${sample_dir}/trimmed/${sample}.fastq.gz" | samtools sort -O BAM | samtools coverage - | cut -f 7 | head -n 2 > "$depth"; then
        log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; continue
    fi
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- STEP 7: QUAST ----------------
    step="quast"; step_start=$(date +%s)
    log "▶ Step $step"
    if ! quast -o "${sample_dir}/qc/quast" -t "$THREADS" -m 100 -l "$sample" "$asm"; then
        log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; continue
    fi
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- STEP 8: Bandage ----------------
    step="bandage"; step_start=$(date +%s)
    log "▶ Step $step"
    gfa="${sample_dir}/assemblies/${sample}.gfa"
    if ! Bandage image "$gfa" "${sample_dir}/qc/assembly_image.svg"; then
        log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; continue
    fi
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- STEP 9: BUSCO ----------------
    step="busco"; step_start=$(date +%s)
    log "▶ Step $step"
    mkdir -p "${sample_dir}/qc/busco_results"
    lineage="${species_dict_busco[$organism]}"
    if ! busco -i "$asm" -f -m genome -l "$lineage" -c "$THREADS" --out_path "${sample_dir}/qc/busco_results" -o "$sample"; then
        log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; continue
    fi
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- STEP 10: abritamr (AMRFINDER) ----------------
    step="amrfinder"; step_start=$(date +%s)
    log "▶ Step $step"
    amr_species="${species_dict_amrfinder[$organism]}"
    if ! conda run -n "$ENV_ABRITAMR" abritamr run -c "${sample_dir}/assemblies/${sample}.fasta" -px "$sample_dir" -j "$THREADS" -sp "$amr_species"; then
        log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; continue
    fi
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- STEP 11: MLST contigs ----------------
    step="mlst_contigs"; step_start=$(date +%s)
    log "▶ Step $step"
    mkdir -p ${sample_dir}/mlst_contigs
    scheme="${species_dict_mlst_contigs[$organism]}"
    if ! conda run -n "$ENV_MLST" mlst --scheme "$scheme" "${sample_dir}/assemblies/${sample}.fasta" > "${sample_dir}/mlst_contigs/mlst_result_contigs.tsv"; then
        log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; continue
    fi
    # Add a step that if krocus MLST result is ND, then use contig mlst result instead
    if [[ -f "${sample_dir}/mlst_result.tsv" ]]; then
        krocus_result=$(awk '{print $1; exit}' "${sample_dir}/mlst_result.tsv")
        if ! [[ "$krocus_result" =~ ^[0-9]+$ ]]; then
            log "❌ Krocus result is ST \"ND\" so using contig mlst result instead"
            cp -v "${sample_dir}/mlst_contigs/mlst_result_contigs.tsv" "${sample_dir}/mlst_result.tsv"
            sequence_type=$(awk '{print $3; exit}' "${sample_dir}/mlst_contigs/mlst_result_contigs.tsv")
            log "✅ MLST result has succesfully changed from ST${krocus_result} to ST${sequence_type}"
        else
                log "✅ Krocus result is present so do not need to use contig mlst result"
            fi
        else
        log "⚠️ Krocus mlst result not found for sample $sample"; continue
    fi
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- STEP 12: Copy genomes for ST specifc tree ----------------
    step="copy_genomes"; step_start=$(date +%s)
    log "▶ Step $step"
    species="${species_dict_read_mlst[$organism]}"
    species_tree_dir="${TREE_BASE}/${species}_trees"
    # Adds a step in MLST result exists
    if [[ -f "${sample_dir}/mlst_result.tsv" ]]; then
        log "✅ MLST result exists."
        else
        log "⚠️ MLST result does not exist, ST set to \"none\"."
        sequence_type=`echo "_NONE"`
    fi
    mlst_tree_dir="${species_tree_dir}/ST${sequence_type}"
    mkdir -p "$species_tree_dir"
    mkdir -p "$mlst_tree_dir"
    cp -v "${sample_dir}/assemblies/${sample}.fasta" "${mlst_tree_dir}/${sample}.fasta" || log "⚠️ copy fasta failed"
    cp -v "${sample_dir}/trimmed/${sample}.fastq.gz" "${mlst_tree_dir}/${sample}.fastq.gz" || log "⚠️ copy fastq failed"
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"
    
    # ---------------- STEP 13: Picking Reference for ST specifc tree ----------------
    step="tree_reference_pick"; step_start=$(date +%s)
    # Sort by sample date and choose the oldest sample as reference
    REF_SAMPLE=`ls -1 "$mlst_tree_dir" | egrep "^[0-9]{1,6}.*.fasta" | grep -v "consensus" | sort -n | head -n 1 | sed 's/.fasta//'`
    log "$REF_SAMPLE chosen as reference sample for mapping as root of tree."
    # Take the first sequence in the fasta. It has to be the full genome 
    if [[ `grep ">1" "${mlst_tree_dir}/${REF_SAMPLE}.fasta" | grep -o "circular=true"` != "circular=true" ]]; then
    log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"
    continue
    else
        log "First sequence of $REF_SAMPLE is a completed bacterial assembly, continue with mapping"
    fi
    grep -a1 ">1" "${mlst_tree_dir}/${REF_SAMPLE}.fasta" > "${mlst_tree_dir}/ref_full_genome.fasta"
    REF_FASTA="${mlst_tree_dir}/ref_full_genome.fasta"
    samtools faidx $REF_FASTA
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- STEP 14: Mapping against tree reference ----------------
    step="reference_mapping"; step_start=$(date +%s)
    bam="${sample}.bam"
    log "Mapping $sample to $REF_SAMPLE and outputting $bam";
    if ! minimap2 -x map-ont -a -t "$THREADS" "$REF_FASTA" "$mlst_tree_dir/${sample}.fastq.gz" \
        | samtools view -h -F 4 \
        | samtools sort -O BAM -o "${mlst_tree_dir}/${bam}"; then
        log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; continue
    fi
    samtools index "${mlst_tree_dir}/${bam}"
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- STEP 15: Variant calling with Clair3 ----------------
    step="variant_calling"; step_start=$(date +%s)
    log "Calling variants for $sample against $REF_SAMPLE and outputting VCF";
    if ! conda run -n clair3 run_clair3.sh \
            --bam_fn="${mlst_tree_dir}/${bam}" \
            --ref_fn="$REF_FASTA" \
            --threads="$THREADS" \
            --platform="ont" \
            --model_path="$CLAIR3_MODEL" \
            --output="${mlst_tree_dir}/clair3/${sample}" \
            --sample_name="$sample" \
            --min_coverage="10" \
            --include_all_ctgs \
            --haploid_precise \
            --no_phasing_for_fa \
            --enable_long_indel; then
        log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; continue
    fi
    mv -v "${mlst_tree_dir}/clair3/${sample}/merge_output.vcf.gz" "${mlst_tree_dir}/${sample}_full.vcf.gz"
    # Filter variants to split MNVs to SNVs and keep only varints that pass filter
    bcftools norm -a -m - "${mlst_tree_dir}/${sample}_full.vcf.gz" |
    bcftools norm -a -d "none" |
    bcftools view -v "snps" -f "PASS" -O "z" > "${mlst_tree_dir}/${sample}_filtered.vcf.gz"
    bcftools index -f "${mlst_tree_dir}/${sample}_filtered.vcf.gz"
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- STEP 16: Generate consensus and multi-alignment FASTA ----------------
    step="consensus"; step_start=$(date +%s)
    log "Generating consensus sequence for $sample";
    bedtools genomecov -ibam "${mlst_tree_dir}/${bam}" -bga | awk '$4 < 10' > "${mlst_tree_dir}/${sample}.bed" # Generate coverage less than 10 bed file
    if ! bcftools consensus -f "$REF_FASTA" -m "${mlst_tree_dir}/${sample}.bed" "${mlst_tree_dir}/${sample}_filtered.vcf.gz" > "${mlst_tree_dir}/${sample}_consensus.fasta"; then
        log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; continue
    fi # Generate consensus, masking low cov regions
    sed -i "s/>.*$/>${sample}/" "${mlst_tree_dir}/${sample}_consensus.fasta" # Change name of sequences in consensus fasta file to sample name
    cat ${mlst_tree_dir}/*_consensus.fasta > ${mlst_tree_dir}/alignment.fasta
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- STEP 17: Filter recombination and generate tree ----------------
    step="recombination"; step_start=$(date +%s)
    # Remove intermediate folders that are made by a failed gubbins run
    rm -vrf ${PWD}/tmp*
    if ! conda run -n gubbins run_gubbins.py --threads "$THREADS" \
            -p "${mlst_tree_dir}/ST${sequence_type}" \
            --first-tree-builder "iqtree-fast"  \
            --tree-builder "iqtree" \
            -o "$REF_SAMPLE" \
            --first-model "JC" \
            --model "GTRGAMMA" \
            -v \
            -f "90" \
            ${mlst_tree_dir}/alignment.fasta; then
        log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; continue
    fi
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- STEP 18: Generate tree image ----------------
    step="tree"; step_start=$(date +%s)
    conda run -n gubbins plot_gubbins.R -t "${mlst_tree_dir}/ST${sequence_type}.node_labelled.final_tree.tre" \
        -r "${mlst_tree_dir}/ST${sequence_type}.recombination_predictions.gff" \
        -o "${sample_dir}/ST${sequence_type}.node_labelled.final_tree.png" \
        --taxon-label-size "2" \
        --tree-width "10" \
        --show-taxa \
        --tree-axis-expansion "100"
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"

    # ---------------- STEP 20: MultiQC ----------------
        step="multiqc"; step_start=$(date +%s)
    log "▶ Step $step"
    if ! multiqc -i "$sample" -f -v -o ${sample_dir}/multiqc_report/ --no-ai -c "$MULTIQC_CONFIG" "$sample_dir"/; then
        log "❌ Step $step failed"; status="FAILED"; failed_step="$step"; echo -e "${sample}\t${status}\t${failed_step}\t-" >> "${MASTER_SUMMARY}"; continue
    fi
    step_end=$(date +%s); log "✅ Step $step done in $((step_end-step_start))s"
    
    # ---------------- Sample completion ----------------
    end_total=$(date +%s)
    total_time=$((end_total-start_total))
    echo -e "${sample}\t${status}\t${failed_step}\t${total_time}" >> "${MASTER_SUMMARY}"
    log "✅ Sample $sample completed in ${total_time}s"
    exec 1>&3 2>&4
    exec 3>&- 4>&-
}; done

echo "Pipeline complete. Master summary saved in $MASTER_SUMMARY"
