set -e
python3 final_tables.py > log_final.txt 2>&1
for K in 3 4; do python3 analyse_longrun.py $K > log_long$K.txt 2>&1; done
python3 analyse_tuned.py 3 > log_tuned3.txt 2>&1
python3 compile_tables.py > log_compile.txt 2>&1
