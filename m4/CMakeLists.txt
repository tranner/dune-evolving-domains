set(M4FILES  dune-evolving-domains.m4)
set(aclocaldir  ${CMAKE_INSTALL_DATADIR}/dune/aclocal)
set(aclocal_DATA  ${M4FILES})
set(EXTRA_DIST  ${M4FILES} CMakeLists.txt)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${aclocal_DATA} DESTINATION ${aclocaldir})