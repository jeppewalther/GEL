ADD_EXECUTABLE(Raytracer
	Camera.h
	raytrace.cpp
)

set(RAYTRACER_RELEASE ${LIBS_GLUT})
set(RAYTRACER_DEBUG ${LIBS_GLUTD})
  
foreach(RELEASE_LIB ${RAYTRACER_RELEASE})
	target_link_libraries(Raytracer optimized ${RELEASE_LIB})
endforeach(RELEASE_LIB)
foreach(DEBUG_LIB ${RAYTRACER_DEBUG})
	target_link_libraries(Raytracer debug ${DEBUG_LIB})
endforeach(DEBUG_LIB)

SET_TARGET_PROPERTIES(Raytracer PROPERTIES XCODE_ATTRIBUTE_CLANG_CXX_LANGUAGE_STANDARD c++0x)
SET_TARGET_PROPERTIES(Raytracer PROPERTIES XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY libc++)

TARGET_LINK_LIBRARIES(Raytracer
	GEL
)
