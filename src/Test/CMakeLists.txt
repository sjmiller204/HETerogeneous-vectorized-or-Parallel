# HETP/src/CMakeLists.txt

add_library(HETP_Test STATIC EXCLUDE_FROM_ALL
	hetp_main.F90
)
target_link_libraries(HETP_Test
	PUBLIC HETP_core
)
target_include_directories(HETP_Test
        INTERFACE ${HETP_BINARY_DIR}/mod
)
target_compile_options(HETP_Test
        PRIVATE
                ""
                $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:-ffree-line-length-none>
)

#if("hetp_test" IN_LIST HETP_EXE_TARGETS)
        add_executable(hetp_test
                hetp_main.F90
        )
        target_link_libraries(hetp_test
                PUBLIC HETP_Test
        )
        set_target_properties(hetp_test PROPERTIES
                RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
        install(TARGETS hetp_test
                RUNTIME DESTINATION ${RUNDIR}
        )
#endif()
