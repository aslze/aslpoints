cmake_minimum_required(VERSION 3.4)

project(aslpoints)

set(CMAKE_CXX_STANDARD 11)

if(NOT TARGET asls)
	find_package(ASL)
	if(NOT TARGET asls)
		cmake_minimum_required(VERSION 3.14...3.24)
		include(FetchContent)
		FetchContent_Declare(asl URL https://github.com/aslze/asl/archive/1.11.5.zip)
		FetchContent_MakeAvailable(asl)
	endif()
endif()

add_library(aslpoints INTERFACE)
target_include_directories(aslpoints INTERFACE .)

if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)

	set(TARGET demo)

	set(SRC
		asl/points.h
		demo.cpp
	)

	add_executable(${TARGET} ${SRC})

	target_link_libraries(${TARGET} asls aslpoints)

endif()
