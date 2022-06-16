# build htslib
set(htslib_PREFIX ${CMAKE_BINARY_DIR}/htslib)
set (FLAGS "-fPIC")
# Enable ExternalProject CMake module
include(ExternalProject)

ExternalProject_Add(htslib
        BUILD_IN_SOURCE 1        
        URL  ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/htslib-1.15.1.tar.bz2 #[https://github.com/samtools/htslib/releases/download/1.15.1/htslib-1.15.1.tar.bz2]
        PREFIX ${htslib_PREFIX}
        UPDATE_COMMAND ""
        CONFIGURE_COMMAND autoreconf -i && ./configure --prefix=${CMAKE_BINARY_DIR}/htslib
        BUILD_COMMAND make CFLAGS=${CMAKE_C_FLAGS}
        INSTALL_COMMAND "")
ExternalProject_Get_Property(htslib SOURCE_DIR)
set(HTSLIB_SRC_DIR ${SOURCE_DIR})

include_directories("${HTSLIB_SRC_DIR}/htslib")
