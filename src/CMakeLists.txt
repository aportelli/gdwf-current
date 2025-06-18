find_package(Grid REQUIRED)

add_executable(${APP_NAME} main.cpp current.hpp)
target_link_libraries(${APP_NAME} PRIVATE Grid::Grid)

install(TARGETS ${APP_NAME} RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
