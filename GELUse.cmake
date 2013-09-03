#import GEL's build settings
INCLUDE(CMakeImportBuildSettings)

CMAKE_IMPORT_BUILD_SETTINGS(${GEL_BUILD_SETTINGS_FILE})

# Tell the compiler where to find GEL's header files
INCLUDE_DIRECTORIES(${GEL_INCLUDE_DIRS})

# Tell the linker where to find GEL's librarires
LINK_DIRECTORIES(${GEL_LIBRARY_DIRS})
