ADD_EXECUTABLE(CGLA-mat
  mat_test.cpp
  )

SET_TARGET_PROPERTIES(CGLA-mat PROPERTIES XCODE_ATTRIBUTE_CLANG_CXX_LANGUAGE_STANDARD c++0x)
SET_TARGET_PROPERTIES(CGLA-mat PROPERTIES XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY libc++)


TARGET_LINK_LIBRARIES(CGLA-mat
  GEL
  )
