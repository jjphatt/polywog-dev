# Library.
add_library(polywog polywog.c 
            partition_point_cloud_with_neighbors.c
            shape_function.c shepard_shape_function.c mls_shape_function.c
            gmls_functional.c gmls_matrix.c mlpg_quadrature.c fvpm_quadrature.c
            fvpm_interparticle_area.c sph_dynamics.c
            interpreter_register_meshless_functions.c)

set(POLYWOG_LIBRARIES polywog;${POLYMEC_LIBRARIES} PARENT_SCOPE)
set(POLYWOG_LIBRARIES polywog;${POLYMEC_LIBRARIES})

add_subdirectory(tests)

# Installation targets.
install(TARGETS polywog DESTINATION lib)

# Install headers, including generated files.
file(GLOB HEADERS *.h)
install(FILES ${HEADERS} DESTINATION include/polymec/polywog)

