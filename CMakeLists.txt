cmake_minimum_required(VERSION 3.14)
project(cf-fitter VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CF_FITTER_ROOT ${CMAKE_SOURCE_DIR})
set(CF_FITTER_CONFIG_DIR ${CMAKE_SOURCE_DIR}/config CACHE PATH "Scenario/channel YAML root")

add_library(CfFitterCore STATIC
  src/YamlReader.cpp
  src/PotentialUtils.cpp
  src/PotentialBuilder.cpp
  src/ChannelSpec.cpp
  src/ChannelRegistry.cpp
  src/Scenario.cpp
  src/ModelFactory.cpp
)
set_target_properties(CfFitterCore PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(CfFitterCore PUBLIC
  ${CMAKE_SOURCE_DIR}/include
)
target_compile_definitions(CfFitterCore PUBLIC
  CF_FITTER_CONFIG_DIR="${CF_FITTER_CONFIG_DIR}"
)

add_executable(cf-calc
  calc/main.cpp
  calc/validate_cmd.cpp
  calc/make_cf_cmd.cpp
)
target_link_libraries(cf-calc PRIVATE CfFitterCore)
target_include_directories(cf-calc PRIVATE ${CMAKE_SOURCE_DIR}/calc)

add_executable(test_square_well_solver tests/test_square_well_solver.cpp)
target_link_libraries(test_square_well_solver PRIVATE CfFitterCore)
target_include_directories(test_square_well_solver PRIVATE ${CMAKE_SOURCE_DIR}/calc)

option(CF_FITTER_BUILD_ROOT "Build ROOT fit library and tests" ON)
if(CF_FITTER_BUILD_ROOT)
  find_program(ROOT_CONFIG root-config)
  if(ROOT_CONFIG)
    execute_process(COMMAND ${ROOT_CONFIG} --cflags OUTPUT_VARIABLE ROOT_CFLAGS_RAW OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(COMMAND ${ROOT_CONFIG} --libs OUTPUT_VARIABLE ROOT_LIBS_RAW OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(COMMAND ${ROOT_CONFIG} --incdir OUTPUT_VARIABLE ROOT_INCDIR OUTPUT_STRIP_TRAILING_WHITESPACE)
    separate_arguments(ROOT_CFLAGS NATIVE_COMMAND ${ROOT_CFLAGS_RAW})
    separate_arguments(ROOT_LIBS NATIVE_COMMAND ${ROOT_LIBS_RAW})
    set(ROOT_CXX_FLAGS "")
    foreach(flag ${ROOT_CFLAGS})
      if(NOT flag MATCHES "^-std=")
        list(APPEND ROOT_CXX_FLAGS ${flag})
      endif()
    endforeach()
    add_library(CfFitter SHARED src/CfFitterRoot.cpp)
    set_target_properties(CfFitter PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED ON)
    target_include_directories(CfFitter PUBLIC ${CMAKE_SOURCE_DIR}/include ${ROOT_INCDIR})
    target_compile_options(CfFitter PRIVATE ${ROOT_CXX_FLAGS})
    target_link_libraries(CfFitter PUBLIC CfFitterCore)
    target_link_libraries(CfFitter PRIVATE ${ROOT_LIBS})
    set_target_properties(CfFitter PROPERTIES OUTPUT_NAME CfFitter PREFIX lib)

    add_executable(cf-fit calc/fit_cmd.cpp)
    target_include_directories(cf-fit PRIVATE ${ROOT_INCDIR})
    target_compile_options(cf-fit PRIVATE ${ROOT_CXX_FLAGS})
    set_target_properties(cf-fit PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED ON)
    target_link_libraries(cf-fit PRIVATE CfFitter ${ROOT_LIBS})

    add_executable(test_ll_kp_consistency tests/test_ll_kp_consistency.cpp)
    target_link_libraries(test_ll_kp_consistency PRIVATE CfFitterCore)
    add_test(NAME ll_kp_consistency COMMAND test_ll_kp_consistency)
  else()
    message(WARNING "root-config not found: skipping libCfFitter.so and fit tests")
  endif()
endif()

enable_testing()
add_test(NAME square_well_solver COMMAND test_square_well_solver)
add_test(NAME validate COMMAND cf-calc validate)

install(TARGETS CfFitterCore cf-calc
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
)
install(DIRECTORY include/femto DESTINATION include)
install(DIRECTORY config DESTINATION share/cf-fitter)
install(DIRECTORY macro DESTINATION share/cf-fitter)
if(TARGET CfFitter)
  install(TARGETS CfFitter LIBRARY DESTINATION lib)
  install(TARGETS cf-fit RUNTIME DESTINATION bin)
endif()
