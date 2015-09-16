# This function adds a (serial) unit test executable to be built using cmockery.
function(add_polywog_executable exe)
  add_executable(${exe} ${ARGN})
  target_link_libraries(${exe} ${POLYWOG_LIBRARIES})
endfunction(add_polywog_executable)

