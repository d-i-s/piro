# Makefile to build class 'pan~' for Pure Data.
# Needs Makefile.pdlibbuilder as helper makefile for platform-dependent build
# settings and rules.

# library name
lib.name = piro

# common sources

common.sources := \
./lib/m_memory.c \
./lib/z_array.c \
./lib/z_fft.c

hirmeasure~ := \
./lib/z_utils.c \
./lib/z_core.c 

hmulticonvolve~ := \
./lib/u_partition_convolve.c \
./lib/u_time_domain_convolve.c\
./lib/u_zero_latency_convolve.c \
./lib/u_output_channel_convolve.c \
./lib/u_multi_channel_convolve.c 

hirmanip := \
./lib/z_matrix_math.c \
./lib/z_utils.c \
./lib/z_core.c 

# input source file (class name == source file basename)
irmeasure~.class.sources := irmeasure~.c $(hirmeasure~)
multiconvolve~.class.sources := multiconvolve~.c $(hmulticonvolve~)
irmanip.class.sources := irmanip.c $(hirmanip)

# all extra files to be included in binary distribution of the library
datafiles = irmeasure~-help.pd multiconvolve~-help.pd irmanip-help.pd

# extra dirs
datadirs = snd \


cflags += -fno-finite-math-only

ldflags +=

suppress-wunused += -Wunused-variable

# include Makefile.pdlibbuilder from submodule directory 'pd-lib-builder'
PDLIBBUILDER_DIR=./pd-lib-builder/
include $(PDLIBBUILDER_DIR)/Makefile.pdlibbuilder
