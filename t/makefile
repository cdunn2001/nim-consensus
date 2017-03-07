# fc_consensus is in FALCON/

N?=0
F?=./foo

new:
	../n.exe >| cnew.fasta < foo
go:
	../go-consensus --output_multi --min_idt 0.70 --min_cov 1 --max_n_read 20000 --n_core $N >| cnew.fasta < foo
old:
	fc_consensus --output_multi --min_idt 0.70 --min_cov 1 --max_n_read 20000 --n_core $N >| cns_00001.fasta < foo
comp:
	diff cnew.fasta cns_00001.fasta
