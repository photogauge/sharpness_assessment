# CMakeLists.txt

# Older versions of CMake are likely to work just fine but, since
# I don't know where to cut off I just use the version I'm using
cmake_minimum_required(VERSION "3.15")

# set OpenCV_DIR variable equal to the path to the cmake
# files within the previously installed opencv program
#set(OpenCV_DIR /Users/adammcquistan/code/c++/opencv/install/lib/cmake/opencv4)

# Tell compiler to use C++ 14 features which is needed because
# Clang version is often behind in the XCode installation
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# configure the necessary common CMake environment variables
# needed to include and link the OpenCV program into this
# demo project, namely OpenCV_INCLUDE_DIRS and OpenCV_LIBS
find_package( OpenCV REQUIRED )

# tell the build to include the headers from OpenCV
include_directories( ${OpenCV_INCLUDE_DIRS} )

# specify the executable target to be built
add_library(${PROJECT_NAME} STATIC sharpness.cpp s3tools.cpp s3estimator.cpp fftwrapper.cpp low_level_fft.cpp)

# specify the executable target to be built
# add_executable(${PROJECT_NAME}-app main.cpp sharpness.cpp s3tools.cpp s3estimator.cpp fftwrapper.cpp low_level_fft.cpp)

# tell it to link the executable target against OpenCV
target_link_libraries(${PROJECT_NAME}  ${OpenCV_LIBS} )

# install(TARGETS ${PROJECT_NAME} ${PROJECT_NAME}-app)


