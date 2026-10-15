# htslib: bundled 1.24 release tarball (default) or the system library.
#
# Either way the result is one target, htslib::hts, that carries the include directories
# (the sources include "sam.h" etc. directly, as the htslib tree lays them out, and also
# <htslib/...>) and the libraries libhts needs. Targets link htslib::hts and nothing else.

# htslib's thread pool and samstat's reader/worker threads both need pthreads.
find_package(Threads REQUIRED)

if(SAMSTAT_USE_SYSTEM_HTSLIB)
  find_package(PkgConfig REQUIRED)
  # 1.22 is the first release that no longer fetches CRAM reference sequences from EBI's server
  # when REF_PATH is unset: an older library would send a CRAM's reference MD5s over the network.
  pkg_check_modules(HTSLIB REQUIRED IMPORTED_TARGET GLOBAL htslib>=1.22)
  set_property(TARGET PkgConfig::HTSLIB APPEND PROPERTY
    INTERFACE_INCLUDE_DIRECTORIES "${HTSLIB_INCLUDEDIR}/htslib")
  add_library(htslib::hts ALIAS PkgConfig::HTSLIB)
else()
  include(ExternalProject)

  # libhts.a needs zlib, the CRAM codec libraries, pthreads and libm; htslib::hts carries them,
  # so they follow libhts.a on every link line that uses it.
  find_package(ZLIB REQUIRED)
  find_package(BZip2 REQUIRED)
  find_package(LibLZMA REQUIRED)

  # htslib's Makefile needs GNU make, whatever generator builds samstat (with Ninja,
  # CMAKE_MAKE_PROGRAM is ninja; on the BSDs "make" is BSD make).
  find_program(SAMSTAT_GNU_MAKE NAMES gmake make DOC "GNU make, used to build the bundled htslib")
  if(NOT SAMSTAT_GNU_MAKE)
    message(FATAL_ERROR "GNU make (gmake or make) is required to build the bundled htslib")
  endif()
  execute_process(COMMAND ${SAMSTAT_GNU_MAKE} --version
    RESULT_VARIABLE gnu_make_rc OUTPUT_VARIABLE gnu_make_version ERROR_QUIET)
  if(NOT gnu_make_rc EQUAL 0 OR NOT gnu_make_version MATCHES "GNU Make")
    message(FATAL_ERROR "${SAMSTAT_GNU_MAKE} is not GNU make, which the bundled htslib needs; "
      "install GNU make or set SAMSTAT_GNU_MAKE to it")
  endif()

  set(htslib_PREFIX ${CMAKE_BINARY_DIR}/htslib)
  set(htslib_SOURCE_DIR ${htslib_PREFIX}/src/htslib)
  set(htslib_LIBRARY ${htslib_SOURCE_DIR}/libhts.a)
  set(htslib_URL ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/htslib-1.24.tar.bz2)
  set(htslib_URL_HASH SHA256=28a8de191381c7a97a35675ceac76fa1ea95e7b678d6a2e9d600a7874e4077de)

  # htslib is built with the compiler, target architecture and SDK of samstat, so that
  # libhts.a links into it (a universal or cross build would otherwise get a host libhts.a).
  set(htslib_ARCH_FLAGS "")
  if(APPLE)
    foreach(arch IN LISTS CMAKE_OSX_ARCHITECTURES)
      string(APPEND htslib_ARCH_FLAGS " -arch ${arch}")
    endforeach()
    if(CMAKE_OSX_SYSROOT AND IS_DIRECTORY "${CMAKE_OSX_SYSROOT}")
      string(APPEND htslib_ARCH_FLAGS " -isysroot ${CMAKE_OSX_SYSROOT}")
    endif()
    if(CMAKE_OSX_DEPLOYMENT_TARGET)
      string(APPEND htslib_ARCH_FLAGS " -mmacosx-version-min=${CMAKE_OSX_DEPLOYMENT_TARGET}")
    endif()
  endif()

  # bz2 and lzma are enabled (CRAM codecs); libcurl/s3 stay off, so there is no network access.
  # htslib is built independently of the samstat target flags (e.g. sanitizers); the global
  # CMAKE_C_FLAGS / CMAKE_EXE_LINKER_FLAGS (e.g. a ThreadSanitizer build) do reach it.
  set(htslib_CFLAGS "-O2 -fPIC ${CMAKE_C_FLAGS}${htslib_ARCH_FLAGS}")
  set(htslib_LDFLAGS "${CMAKE_EXE_LINKER_FLAGS}${htslib_ARCH_FLAGS}")
  set(htslib_CONFIGURE_ARGS
    "CC=${CMAKE_C_COMPILER}"
    "CFLAGS=${htslib_CFLAGS}"
    "LDFLAGS=${htslib_LDFLAGS}"
    --prefix=${htslib_PREFIX}
    --disable-libcurl --disable-s3 --without-libdeflate)
  # A cross build tells configure the target triple (it cannot run test programs otherwise).
  set(SAMSTAT_HTSLIB_HOST "" CACHE STRING
    "Target triple passed as --host to the bundled htslib's configure in a cross build (default: CMAKE_C_COMPILER_TARGET)")
  set(htslib_HOST "${SAMSTAT_HTSLIB_HOST}")
  if(NOT htslib_HOST)
    set(htslib_HOST "${CMAKE_C_COMPILER_TARGET}")
  endif()
  if(CMAKE_CROSSCOMPILING AND htslib_HOST)
    list(APPEND htslib_CONFIGURE_ARGS --host=${htslib_HOST})
  endif()

  ExternalProject_Add(htslib
        BUILD_IN_SOURCE 1
        URL ${htslib_URL}
        URL_HASH ${htslib_URL_HASH}
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
        PREFIX ${htslib_PREFIX}
        SOURCE_DIR ${htslib_SOURCE_DIR}
        UPDATE_COMMAND ""
        CONFIGURE_COMMAND ./configure ${htslib_CONFIGURE_ARGS}
        BUILD_COMMAND ${SAMSTAT_GNU_MAKE}
        # Ninja must know which step writes libhts.a, or it has no rule for the link input
        BUILD_BYPRODUCTS ${htslib_LIBRARY}
        INSTALL_COMMAND "")

  # An imported target may not name include directories that do not exist yet; the tarball
  # is only unpacked at build time (the extract step replaces this empty directory).
  file(MAKE_DIRECTORY ${htslib_SOURCE_DIR}/htslib)

  add_library(htslib::hts STATIC IMPORTED GLOBAL)
  set_target_properties(htslib::hts PROPERTIES
    IMPORTED_LOCATION ${htslib_LIBRARY}
    INTERFACE_INCLUDE_DIRECTORIES "${htslib_SOURCE_DIR};${htslib_SOURCE_DIR}/htslib"
    INTERFACE_LINK_LIBRARIES "ZLIB::ZLIB;BZip2::BZip2;LibLZMA::LibLZMA;Threads::Threads;$<$<BOOL:${UNIX}>:m>")
  # every target that links htslib::hts waits for the htslib build
  add_dependencies(htslib::hts htslib)
endif()
