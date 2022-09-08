set(WRAP_FILE CustomAnisotropicNonbondedPluginWrapper.cpp)
set(Module_Name canplugin)

add_custom_command(
	OUTPUT "${WRAP_FILE}"
	COMMAND "${SWIG_EXECUTABLE}"
		-python -c++
		-o "${WRAP_FILE}"
		"-I${OPENMM_DIR}/include"
		"${CMAKE_CURRENT_SOURCE_DIR}/canplugin.i"
	DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/canplugin.i"
	WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

add_custom_target(PythonInstall DEPENDS "${WRAP_FILE}")
set(CUSTOMANISOTROPICNONBONDEDPLUGIN_HEADER_DIR "${CMAKE_SOURCE_DIR}/openmmapi/include")
set(CUSTOMANISOTROPICNONBONDEDPLUGIN_LIBRARY_DIR "${CMAKE_BINARY_DIR}")
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/setup.py ${CMAKE_CURRENT_BINARY_DIR}/setup.py)
add_custom_command(TARGET PythonInstall
	COMMAND "${PYTHON_EXECUTABLE}" setup.py build
	COMMAND "${PYTHON_EXECUTABLE}" setup.py install
	WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)

