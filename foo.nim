#{.deadCodeElim: on.}
proc poo*() {.cdecl, importc, header:"../foo.h".}

poo()
