

PATH=VTD_SCR_TestD_D3_F1_pathogenic_variants.csv
NEUTRAL=VTD_SCR_TestD_D3_F2_neutral_variants.csv
ID_MAP=uniprot_id_to_name.tsv
OUT=merged.fa
./format_variants.py --pathogenic $PATH --neutral $NEUTRAL --id_map $ID_MAP --outfile $OUT
