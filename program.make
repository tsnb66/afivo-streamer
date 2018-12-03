# This is a template to build an afivo-streamer program

ifndef NDIM
$(error NDIM is not set)
endif

ifndef MAIN_DIR
$(error MAIN_DIR is not set)
endif

AFIVO_DIR := $(MAIN_DIR)/afivo
LIBDIRS := $(MAIN_DIR)/lib_$(NDIM)d $(AFIVO_DIR)/lib_$(NDIM)d	\
$(AFIVO_DIR)/external_libraries/silo/lib
INCDIRS := $(TARGET_DIR) $(MAIN_DIR)/lib_$(NDIM)d $(AFIVO_DIR)/lib_$(NDIM)d
LIBS := streamer afivo silo
PROG := streamer_$(NDIM)d

.PHONY: all clean

all: $(PROG)

clean:
	$(RM) *.o *.mod $(PROG)

vpath %.f90 $(MAIN_DIR)/src

# Include compilation rules
include $(AFIVO_DIR)/src/makerules.make

# Optionally include a local makefile
-include local.make

FFLAGS += -DNDIM=$(NDIM)

$(PROG): streamer.f90
	$(FC) -o $@ $(filter %.f90 %.o, $^) $(FFLAGS) $(addprefix -I,$(INCDIRS)) \
	$(addprefix -L,$(LIBDIRS)) $(addprefix -l,$(LIBS))

# Dependencies
$(PROG): $(MAIN_DIR)/lib_$(NDIM)d/libstreamer.a
$(PROG): $(AFIVO_DIR)/lib_$(NDIM)d/libafivo.a
$(PROG): m_user.o
