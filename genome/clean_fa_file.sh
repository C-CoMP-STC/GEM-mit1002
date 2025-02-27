# Add the string "MTI1002_anvio_prot_seq_" before each sequence ID in the fastaa file amd save it as a new file
# This is necessary for CarveMe to work
sed 's/^>/>MTI1002_anvio_prot_seq_/' genome/MIT1002_anvio_prot_seqs.fa > genome/MIT1002_anvio_prot_seqs_clean.fa