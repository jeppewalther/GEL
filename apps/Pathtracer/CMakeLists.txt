ADD_EXECUTABLE(Pathtrace
	camera.cpp
	glass.cpp
	luminaire.cpp
	material.cpp
	matte.cpp
	mesh.cpp
	metal.cpp
	mt19937ar.cpp
	omni.cpp
	pathtracer.cpp
	plastic.cpp
	scene.cpp
	camera.h
	core.h
	glass.h
	luminaire.h
	material.h
	matte.h
	mesh.h
	metal.h
	mt19937ar.h
	omni.h
	pathtracer.h
	plastic.h
	scene.h
)

set(PATHTRACE_RELEASE ${LIBS_GLUT})
set(PATHTRACE_DEBUG ${LIBS_GLUTD})
  
foreach(RELEASE_LIB ${PATHTRACE_RELEASE})
	target_link_libraries(Pathtrace optimized ${RELEASE_LIB})
endforeach(RELEASE_LIB)
foreach(DEBUG_LIB ${PATHTRACE_DEBUG})
	target_link_libraries(Pathtrace debug ${DEBUG_LIB})
endforeach(DEBUG_LIB)

SET_TARGET_PROPERTIES(Pathtrace PROPERTIES XCODE_ATTRIBUTE_CLANG_CXX_LANGUAGE_STANDARD c++0x)
SET_TARGET_PROPERTIES(Pathtrace PROPERTIES XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY libc++)


TARGET_LINK_LIBRARIES(Pathtrace
	GEL
)
