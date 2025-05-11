EXTENSION = pg_dphyp

MODULE_big = pg_dphyp
OBJS = $(WIN32RES) \
	src/pg_dphyp.o \
	src/unionset.o \
	src/dphyp_generic.o \
	src/dphyp_simple.o \
	src/simplebms.o

PG_CPPFLAGS += -Isrc/include

ifdef USE_PGXS
PG_CONFIG := pg_config
PGXS := $(shell $(PG_CONFIG) --pgxs)
include $(PGXS)
else
subdir = contrib/pg_dphyp
top_builddir = ../..
include $(top_builddir)/src/Makefile.global
include $(top_srcdir)/contrib/contrib-global.mk
endif
