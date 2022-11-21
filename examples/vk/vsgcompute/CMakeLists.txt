set(SOURCES vsgcompute.cpp)

add_executable(vsgcompute ${SOURCES})

target_link_libraries(vsgcompute vsg::vsg)

install(TARGETS vsgcompute RUNTIME DESTINATION ${VSGEXAMPLES_INSTALL_BINDIR})
