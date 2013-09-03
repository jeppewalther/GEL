add_executable(CGLA-ogl
	ogl_test.cpp
)

set(CGLA-OGL_RELEASE ${LIBS_GLUT})
set(CGLA-OGL_DEBUG ${LIBS_GLUTD})
  
foreach(RELEASE_LIB ${CGLA-OGL_RELEASE})
	target_link_libraries(CGLA-ogl optimized ${RELEASE_LIB})
endforeach(RELEASE_LIB)
foreach(DEBUG_LIB ${CGLA-OGL_DEBUG})
	target_link_libraries(CGLA-ogl debug ${DEBUG_LIB})
endforeach(DEBUG_LIB)

SET_TARGET_PROPERTIES(CGLA-ogl PROPERTIES XCODE_ATTRIBUTE_CLANG_CXX_LANGUAGE_STANDARD c++0x)
SET_TARGET_PROPERTIES(CGLA-ogl PROPERTIES XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY libc++)

target_link_libraries(CGLA-ogl
	GEL
)
