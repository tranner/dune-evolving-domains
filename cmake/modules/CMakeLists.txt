set(MODULES  DuneEvolvingDomainsMacros.cmake)
set(modulesdir  ${CMAKE_INSTALL_DATADIR}/dune/cmake/modules)
install(FILES ${dist_modules_DATA} DESTINATION ${dist_modulesdir})