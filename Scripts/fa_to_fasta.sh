mkdir -p MAGs/A
mkdir -p MAGs/B

for file in 5_BIN-REASSEMBLY/A/reassembled_bins/*.fa;do base=$(basename $file .fa);cp "$file" "MAGs/A/$base.fasta";done
for file in 5_BIN-REASSEMBLY/B/reassembled_bins/*.fa;do base=$(basename $file .fa);cp "$file" "MAGs/B/$base.fasta";done
