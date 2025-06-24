EXTENSION = pg_dphyp

MODULE_big = pg_dphyp
OBJS = $(WIN32RES) pg_dphyp.o

DATA = pg_dphyp--0.1.sql

REGRESS = join
REGRESS_OPTS = --temp-config $(top_srcdir)/contrib/pg_dphyp/regress.conf

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
