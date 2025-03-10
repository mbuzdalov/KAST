# Minimum cmake version
cmake_minimum_required (VERSION 3.20)

# Name of project and that it is C++ only.
project (kast CXX)

# We want to be at least C++23 because of seqan
set(CMAKE_CXX_STANDARD 23)

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search for zlib as a dependency for SeqAn.
find_package (ZLIB REQUIRED)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package (Threads REQUIRED)

SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")

# Load the SeqAn module and fail if not found.
find_package (seqan3 REQUIRED)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN3_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN3_DEFINITIONS})

# Add CXX flags found by find_package (SeqAn).
#set (CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${SEQAN3_CXX_FLAGS} -Wl,--whole-archive -lrt -pthread -lgomp -Wl,--no-as-needed -ldl")

# -Wno-interference-size: the interference size warning is triggered by seqan3/contrib/parallel/buffer_queue.hpp,
#                         and it is not something I can deal with.
set (CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${SEQAN3_CXX_FLAGS} -Wno-interference-size")

# Add executable and link against SeqAn dependencies.
add_executable (kast utils.cpp distance.cpp kast.cpp)
target_link_libraries (kast ${SEQAN3_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
add_executable (kast_tests utils.cpp distance.cpp kast_tests.cpp)
target_link_libraries (kast_tests ${SEQAN3_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
