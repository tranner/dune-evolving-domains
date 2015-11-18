
set(SUBDIRS  data output source)
set(noinst_PROGRAMS  dune_evolving_domains)
add_executable(dune_evolving_domains dune_evolving_domains.cc)
add_dune_mpi_flags(dune_evolving_domains )
target_link_libraries(dune_evolving_domains ${DUNE_LIBS})
add_dune_alugrid_flags(dune_evolving_domains )
# The libraries have to be given in reverse order (most basic libraries
# last).  Also, due to some misunderstanding, a lot of libraries include the
# -L option in LDFLAGS instead of LIBS -- so we have to include the LDFLAGS
# here as well.
# The following lines are comment out as the translation is not clear. TODO review
#
# set(dune_evolving_domains_ LDADD
#        ${DUNE_LDFLAGS} ${DUNE_LIBS}
#        ${ALUGRID_LDFLAGS} ${ALUGRID_LIBS}
#        ${ALBERTA_LDFLAGS} ${ALBERTA_LIBS}
#        ${AMIRAMESH_LDFLAGS} ${AMIRAMESH_LIBS}
#        ${UG_LDFLAGS} ${UG_LIBS}
#        ${DUNEMPILIBS}
#        ${LDADD})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(dune_evolving_domains_ LDFLAGS ${AM_LDFLAGS}
#        ${DUNEMPILDFLAGS}
#        ${UG_LDFLAGS}
#        ${AMIRAMESH_LDFLAGS}
#        ${ALBERTA_LDFLAGS}
#        ${ALUGRID_LDFLAGS}
#        ${DUNE_LDFLAGS})
# don't follow the full GNU-standard
# we need automake 1.9
set(AUTOMAKE_OPTIONS  foreign 1.9)
# pass most important options when "make distcheck" is used
set(DISTCHECK_CONFIGURE_FLAGS  --with-dune-common=${DUNE_COMMON_ROOT} --with-dune-geometry=${DUNE_GEOMETRY_ROOT} --with-dune-grid=${DUNE_GRID_ROOT} --with-dune-alugrid=${DUNE_ALUGRID_ROOT} --with-dune-istl=${DUNE_ISTL_ROOT} --with-dune-localfunctions=${DUNE_LOCALFUNCTIONS_ROOT} --with-dune-fem=${DUNE_FEM_ROOT}  CXX="${CXX}" CC="${CC}")
set(EXTRA_DIST  CMakeLists.txt)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
foreach(i ${SUBDIRS})
  if(${i} STREQUAL "test")
    set(opt EXCLUDE_FROM_ALL)
  endif(${i} STREQUAL "test")
  add_subdirectory(${i} ${opt})
  unset(opt)
endforeach(i ${SUBDIRS})
