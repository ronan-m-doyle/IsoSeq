#!/usr/bin/env bash
# ----------------------------------------------------------------------
# Isolate Sequencing Pipeline Config v0.9
# ----------------------------------------------------------------------
IFS=$'\n\t'

# -------------------------- CONFIG ----------------------------
KRAKEN2_DB="/data/IsoSeq/databases/kraken_reference/"
PY_KRAKEN_SUMMARY="/data/IsoSeq/scripts/kraken2_top5_species.py"
TREE_BASE="/data/IsoSeq_results/trees/"
KROCUS_DB="/data/IsoSeq/databases/mlst_krocus/"
CLAIR3_MODEL="${HOME}/miniforge3/envs/clair3/bin/models/r1041_e82_400bps_sup_v500/"
MLST_BLAST_DB="/data/IsoSeq/databases/mlst/blast/mlst.fa"
MLST_DB="/data/IsoSeq/databases/mlst/pubmlst/"

# -------------------------- ENVS -------------------------------
ENV_MLST="mlst_contigs"
ENV_GUBBINS="gubbins"
ENV_CLAIR3="clair3"

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
    ["Streptococcus pneumoniae"]="Streptococcus_pneumoniae"
    ["Klebsiella aerogenes"]="Klebsiella_aerogenes"
    ["Corynebacterium diphtheriae"]="Corynebacterium_diphtheriae"
    ["Listeria monocytogenes"]="Listeria_monocytogenes"
)

declare -A species_dict_mlst_contigs=(
    ["Acinetobacter baumannii"]="abaumannii"
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
    ["Streptococcus pneumoniae"]="spneumoniae"
    ["Klebsiella aerogenes"]="kaerogenes"
    ["Corynebacterium diphtheriae"]="diphtheria_3"
    ["Listeria monocytogenes"]="listeria_2"
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
    ["Streptococcus pneumoniae"]="Streptococcus_pneumoniae"
    ["Klebsiella aerogenes"]="Klebsiella_pneumoniae"
    ["Corynebacterium diphtheriae"]="Corynebacterium_diphtheriae"
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
    ["Streptococcus pneumoniae"]="streptococcaceae_odb12.2"
    ["Klebsiella aerogenes"]="enterobacteriaceae_odb12"
    ["Corynebacterium diphtheriae"]="corynebacterium_odb12.2"
    ["Listeria monocytogenes"]="listeria_odb12.2"
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
    ["Streptococcus pneumoniae"]="2.0m"
    ["Klebsiella aerogenes"]="5.3m"
    ["Corynebacterium diphtheriae"]="2.5m"
    ["Listeria monocytogenes"]="2.9m"
)
