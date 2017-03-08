# fc_consensus is in FALCON/

N?=0
F?=./foo

new:
	../n.exe >| cnew.fasta < $F
go:
	../go-consensus --output_multi --min_idt 0.70 --min_cov 1 --max_n_read 20000 --n_core $N >| cnew.fasta < $F
old:
	fc_consensus --output_multi --min_idt 0.70 --min_cov 1 --max_n_read 20000 --n_core $N >| cns_00001.fasta < $F
comp:
	diff cnew.fasta cns_00001.fasta
