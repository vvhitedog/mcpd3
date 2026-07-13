function(mcpd_configure_capacity target_name include_root)
  set(MCPD_CAPACITY_MODE "32" CACHE STRING
    "Capacity precision: 32, 64, 128, or gmp")
  set_property(CACHE MCPD_CAPACITY_MODE PROPERTY STRINGS 32 64 128 gmp)

  add_library(${target_name} INTERFACE)
  target_include_directories(${target_name} INTERFACE ${include_root})

  find_path(MCPD_BOOST_INCLUDE_DIR NAMES boost/multiprecision/cpp_int.hpp
    HINTS ENV BOOST_ROOT PATH_SUFFIXES include)
  if (NOT MCPD_BOOST_INCLUDE_DIR)
    message(FATAL_ERROR
      "MCPD capacity support requires Boost.Multiprecision headers")
  endif ()
  target_include_directories(${target_name} INTERFACE
    ${MCPD_BOOST_INCLUDE_DIR})

  if (MCPD_CAPACITY_MODE STREQUAL "32")
    target_compile_definitions(${target_name} INTERFACE
      MCPD_CAPACITY_MODE_32=1)
  elseif (MCPD_CAPACITY_MODE STREQUAL "64")
    target_compile_definitions(${target_name} INTERFACE
      MCPD_CAPACITY_MODE_64=1)
  elseif (MCPD_CAPACITY_MODE STREQUAL "128")
    target_compile_definitions(${target_name} INTERFACE
      MCPD_CAPACITY_MODE_128=1)
  elseif (MCPD_CAPACITY_MODE STREQUAL "gmp")
    find_path(MCPD_GMPXX_INCLUDE_DIR NAMES gmpxx.h
      HINTS ENV GMP_ROOT PATH_SUFFIXES include)
    find_path(MCPD_GMP_INCLUDE_DIR NAMES gmp.h
      HINTS ENV GMP_ROOT PATH_SUFFIXES include include/x86_64-linux-gnu)
    find_library(MCPD_GMPXX_LIBRARY NAMES gmpxx
      HINTS ENV GMP_ROOT PATH_SUFFIXES lib lib64 lib/x86_64-linux-gnu)
    find_library(MCPD_GMP_LIBRARY NAMES gmp
      HINTS ENV GMP_ROOT PATH_SUFFIXES lib lib64 lib/x86_64-linux-gnu)
    if (NOT MCPD_GMPXX_INCLUDE_DIR OR NOT MCPD_GMP_INCLUDE_DIR OR
        NOT MCPD_GMPXX_LIBRARY OR NOT MCPD_GMP_LIBRARY)
      message(FATAL_ERROR
        "MCPD_CAPACITY_MODE=gmp requires the GMP C and C++ development files")
    endif ()
    target_compile_definitions(${target_name} INTERFACE
      MCPD_CAPACITY_MODE_GMP=1)
    target_include_directories(${target_name} INTERFACE
      ${MCPD_GMPXX_INCLUDE_DIR} ${MCPD_GMP_INCLUDE_DIR})
    target_link_libraries(${target_name} INTERFACE
      ${MCPD_GMPXX_LIBRARY} ${MCPD_GMP_LIBRARY})
  else ()
    message(FATAL_ERROR
      "MCPD_CAPACITY_MODE must be one of: 32, 64, 128, gmp")
  endif ()
endfunction()
