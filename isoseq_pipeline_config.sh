#!/usr/bin/env bash
# ----------------------------------------------------------------------
# Isolate Sequencing Pipeline Config - Version 0.6
# ----------------------------------------------------------------------
IFS=$'\n\t'

# -------------------------- CONFIG ----------------------------
KRAKEN2_DB="kraken_reference/"
PY_KRAKEN_SUMMARY="scripts/kraken2_top5_species.py"
TREE_BASE="trees/"
KROCUS_DB="mlst_krocus/"
CLAIR3_MODEL="${HOME}/miniforge3/envs/clair3/bin/models/r1041_e82_400bps_sup_v500/"

# -------------------------- ENVS -------------------------------
ENV_MLST="mlst_contigs"
ENV_ABRITAMR="abritamr"
ENV_ABRICATE="abricate"
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
