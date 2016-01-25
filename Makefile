DIR_RELEASE = build/Release
DIR_DEBUG   = build/Debug
DIR_PROFILE = build/Profile

DIR_CURRENT = $(PWD)
DIR_INSTALL = $(DIR_CURRENT)/target

# Choose build system for CMake
ifeq ($(OS),Windows_NT)
	CMAKE_BUILD_SYSTEM := "MSYS Makefiles"
else
	CMAKE_BUILD_SYSTEM := "Unix Makefiles"
endif

# Set common flags
CMAKE_FLAGS := -G $(CMAKE_BUILD_SYSTEM) -DCMAKE_INSTALL_PREFIX=$(DIR_INSTALL)

cmake-release:
	(mkdir -p $(DIR_RELEASE) && cd $(DIR_RELEASE) ; cmake $(CMAKE_FLAGS) -DCMAKE_BUILD_TYPE=RELEASE $(DIR_CURRENT))

cmake-debug:
	(mkdir -p $(DIR_DEBUG) && cd $(DIR_DEBUG) ; cmake $(CMAKE_FLAGS) -DCMAKE_BUILD_TYPE=DEBUG $(DIR_CURRENT))

cmake-profile:
	(mkdir -p $(DIR_PROFILE) && cd $(DIR_PROFILE) ; cmake $(CMAKE_FLAGS) -DCMAKE_BUILD_TYPE=PROFILE $(DIR_CURRENT))
