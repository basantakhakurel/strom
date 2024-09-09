cmake_minimum_required(VERSION 3.10)
project(strom VERSION 1.0 LANGUAGES CXX)

# Set C++ standard to C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set the prefix (installation directory)
set(CMAKE_INSTALL_PREFIX "/home/basanta/projects/strom")

# Static linking option
set(CMAKE_EXE_LINKER_FLAGS "-static")

# Include directory for the NCL library
include_directories("/home/basanta/include")

# Find and link the NCL static library
find_library(NCL_LIB ncl PATHS "/home/basanta/lib/static" REQUIRED)

# Define the executable
add_executable(strom main.cpp)

# Link the NCL library to the executable
target_link_libraries(strom PRIVATE ${NCL_LIB})

# Set the installation directory for the executable
install(TARGETS strom DESTINATION ${CMAKE_INSTALL_PREFIX})

# Install the additional data file
install(FILES test.tre DESTINATION ${CMAKE_INSTALL_PREFIX})
