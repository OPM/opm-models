SUBDIRS = . material nonlinear twophase

# the standard debug streams are put into the libdune
noinst_LTLIBRARIES = libdumux.la
SUBLIBS = shape /libshapefunct ns.la /libonedinndgrid.la

libdumux_la_SOURCES =
libdumux_la_LIBADD = $(SUBLIBS)

dumuxdir = $(includedir)/dune/dumux

include $(top_srcdir)/am/global-rules
