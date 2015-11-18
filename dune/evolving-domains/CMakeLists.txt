set(evolving_domainsincludedir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/evolving-domains)
set(evolving_domainsinclude_HEADERS  evolving-domains.hh)
set(EXTRA_DIST  CMakeLists.txt)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${evolving_domainsinclude_HEADERS} DESTINATION ${evolving_domainsincludedir})