add_executable(MeshEdit
	harmonics.cpp
	harmonics.h
	polarize.h
	polarize.cpp
	meshedit.cpp
	VisObj.cpp
	VisObj.h
)

set(MESHEDIT_RELEASE ${LIBS_GLUT} ${LIBS_GLEW} ${LIBS_GLCONSOLE})
set(MESHEDIT_DEBUG ${LIBS_GLUTD} ${LIBS_GLEWD} ${LIBS_GLCONSOLED})
SET_TARGET_PROPERTIES(MeshEdit PROPERTIES XCODE_ATTRIBUTE_CLANG_CXX_LANGUAGE_STANDARD c++0x)
SET_TARGET_PROPERTIES(MeshEdit PROPERTIES XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY libc++)

foreach(RELEASE_LIB ${MESHEDIT_RELEASE})
	target_link_libraries(MeshEdit optimized ${RELEASE_LIB})
endforeach(RELEASE_LIB)
foreach(DEBUG_LIB ${MESHEDIT_DEBUG})
	target_link_libraries(MeshEdit debug ${DEBUG_LIB})
endforeach(DEBUG_LIB)

target_link_libraries(MeshEdit GEL)
