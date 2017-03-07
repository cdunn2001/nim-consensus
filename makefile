do: n.exe
	make -C t # LA4Falcon | ./n.exe

foo.nim: foo.h
	c2nim foo.h --out:foo.nim
common.nim:
	c2nim --header --cdecl foo.h --out:foo.nim
run-%: %.nim
	nim -d:debug --out:$*.exe c $<
	./$*.exe
%.exe: %.nim
	nim -d:debug --out:$*.exe c $<

# Not needed, but possible.
LIBS:=DW_banded.o falcon.o kmer_lookup.o
lib: ${LIBS}
	${CC} -shared -o libcommon.so ${LIBS}
