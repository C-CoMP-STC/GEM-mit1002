# divide_batches.py
'''
Split up the FASTA file of all protein sequences into multiple smaller files
with equal numbers of sequences
'''

import sys

batch_size = int(sys.argv[1])

# accumulate batch_size seqs, then write to a new file
with open('2738541267_aagenesequences.fa', 'r') as f_in:
    seqs = 0
    batch = 1
    batch_lines = list()
    for in_line in f_in:
        if in_line.startswith('>'):
            # if we're at the batch size we want to write the lines we have
            # in the list before adding the header for the next sequence
            if seqs == batch_size:
                with open(f'batched_input_seqs/batch_{batch}.fa', 'w') as f_out:
                    for out_line in batch_lines:
                        f_out.write(out_line)
                # reset list of lines and sequence counter, increment batch
                batch_lines = list()
                seqs = 0
                batch += 1
            # now we can add the header of the next sequence to lines because
            # we're either in the middle of a batch or the start of a new one
            batch_lines.append(in_line)
            seqs += 1
        else:
            batch_lines.append(in_line)
