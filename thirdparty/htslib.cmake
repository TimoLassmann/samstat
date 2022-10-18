# build htslib
set(htslib_PREFIX ${CMAKE_BINARY_DIR}/htslib)
set (FLAGS "-fPIC")
# Enable ExternalProject CMake module
include(ExternalProject)

set(htslib_URL ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/htslib-1.16.tar.bz2)

# calculate the MD5 sum of the file downloaded and set it in a variable
set(htslib_URL_MD5 d31777ef90d1369a52049ba0ac3c0375)

ExternalProject_Add(htslib
        BUILD_IN_SOURCE 1        
        URL ${htslib_URL}
	    URL_MD5 ${htslib_URL_MD5}
        PREFIX ${htslib_PREFIX}
        UPDATE_COMMAND ""
        CONFIGURE_COMMAND autoreconf -i && ./configure --prefix=${CMAKE_BINARY_DIR}/htslib --disable-bz2 --disable-lzma --disable-libcurl --disable-s3 
        BUILD_COMMAND make CFLAGS=${CMAKE_C_FLAGS}
        INSTALL_COMMAND "")
ExternalProject_Get_Property(htslib SOURCE_DIR)
set(HTSLIB_SRC_DIR ${SOURCE_DIR})

include_directories("${HTSLIB_SRC_DIR}/htslib")
