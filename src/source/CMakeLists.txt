set(GRIDTYPE  ALUGRID_SIMPLEX)
set(GRIDDIM  3)
set(WORLDDIM  3)
set(POLORDER  1)
set(WANT_ISTL  1)

set(CLEANFILES  *.log)

set(check_PROGRAMS  main-surface main-coupled)

SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -DPOLORDER=${POLORDER} -DWANT_ISTL=${WANT_ISTL}")

# list(APPEND AM_CPPFLAGS ${ALL_PKG_CPPFLAGS} ${DUNEMPICPPFLAGS} -DPOLORDER=${POLORDER} -DWANT_ISTL=${WANT_ISTL})
# set(AM_LDFLAGS  ${ALL_PKG_LDFLAGS} ${DUNEMPILDFLAGS})
# set(LDADD  ${ALL_PKG_LIBS} ${DUNEMPILIBS})

add_executable(main-surface main-surface.cc)
target_link_dune_default_libraries( main-surface )
add_dune_alberta_flags(GRIDDIM ${GRIDDIM} main-surface)

add_executable(main-coupled main-coupled.cc)
target_link_dune_default_libraries( main-coupled )
add_dune_alberta_flags(GRIDDIM ${GRIDDIM} main-surface)

set(EXTRA_DIST   main.cc deformation.hh dirichletconstraints.hh elliptic.hh femscheme.hh heat.hh heatmodel.hh heatscheme.hh model.hh probleminterface.hh rhs.hh temporalprobleminterface.hh heat2.dox)

#add_directory_test_target(_test_target)
#add_dependencies(${_test_target} ${TESTS} ${COMPILE_XFAIL_TESTS})
