# HETP/src/Core/CMakeLists.txt

add_library(HETP_core STATIC EXCLUDE_FROM_ALL
        hetp_mod.F90
        mach_hetp_mod.F90
)
target_link_libraries(HETP_core
	PUBLIC HETPBuildProperties
)
target_include_directories(HETP_core
	INTERFACE ${HETP_BINARY_DIR}/mod
)
target_compile_options(HETP_core
	PRIVATE
		""
		$<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:-ffree-line-length-none>
)
