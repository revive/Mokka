name := tbeam
G4TARGET := $(name)
G4EXLIB := true

.PHONY: all
all: lib 

CPPFLAGS += -I../../Kernel/include -I../CGA/include 

ifdef GEAR
CPPFLAGS += -I../MokkaGear/include
endif


include ../../Kernel/GNUmakefile.common
include $(G4INSTALL)/config/binmake.gmk
