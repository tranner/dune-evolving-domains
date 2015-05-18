GRIDTYPE = ALUGRID_SIMPLEX
GRIDDIM = 2
WORLDDIM = 3
POLORDER = 1
WANT_ISTL = 1

CLEANFILES = *.log

TESTS = main
check_PROGRAMS = main

AM_CPPFLAGS += $(ALL_PKG_CPPFLAGS) $(DUNEMPICPPFLAGS) -DPOLORDER=$(POLORDER) -DWANT_ISTL=$(WANT_ISTL)
AM_LDFLAGS = $(ALL_PKG_LDFLAGS) $(DUNEMPILDFLAGS)
LDADD = $(ALL_PKG_LIBS) $(DUNEMPILIBS)

main_SOURCES = main.cc


LESSON = 12
SOLUTION = 1

EXTRA_DIST =  main.cc deformation.hh dirichletconstraints.hh elliptic.hh femscheme.hh heat.hh heatmodel.hh heatscheme.hh model.hh probleminterface.hh rhs.hh temporalprobleminterface.hh heat2.dox
include $(top_srcdir)/am/global-rules
