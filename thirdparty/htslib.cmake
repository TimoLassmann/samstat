# htslib: bundled 1.24 release tarball (default) or the system library.
if(SAMSTAT_USE_SYSTEM_HTSLIB)
  find_package(PkgConfig REQUIRED)
  pkg_check_modules(HTSLIB REQUIRED IMPORTED_TARGET htslib>=1.16)
  # sources include "sam.h" etc. directly, as the bundled tree lays them out
  include_directories("${HTSLIB_INCLUDEDIR}/htslib")
else()
  include(ExternalProject)

  set(htslib_PREFIX ${CMAKE_BINARY_DIR}/htslib)
  set(htslib_URL ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/htslib-1.24.tar.bz2)
  set(htslib_URL_HASH SHA256=28a8de191381c7a97a35675ceac76fa1ea95e7b678d6a2e9d600a7874e4077de)

  # htslib is built independently of the samstat target flags (e.g. sanitizers).
  set(htslib_CFLAGS "-O2 -fPIC ${CMAKE_C_FLAGS}")

  ExternalProject_Add(htslib
        BUILD_IN_SOURCE 1
        URL ${htslib_URL}
        URL_HASH ${htslib_URL_HASH}
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
        PREFIX ${htslib_PREFIX}
        UPDATE_COMMAND ""
        CONFIGURE_COMMAND ./configure "CFLAGS=${htslib_CFLAGS}" --prefix=${CMAKE_BINARY_DIR}/htslib --disable-bz2 --disable-lzma --disable-libcurl --disable-s3 --without-libdeflate
        BUILD_COMMAND make
        INSTALL_COMMAND "")
  ExternalProject_Get_Property(htslib SOURCE_DIR)
  set(HTSLIB_SRC_DIR ${SOURCE_DIR})

  include_directories("${HTSLIB_SRC_DIR}/htslib")
endif()
