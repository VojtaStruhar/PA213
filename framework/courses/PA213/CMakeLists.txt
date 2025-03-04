################################################################################
# Common Framework for Computer Graphics Courses at FI MUNI.
#
# Copyright (c) Visitlab (https://visitlab.fi.muni.cz)
# All rights reserved.
#
# Course: PA213
################################################################################

# Adds the available project skeletons.
visitlab_add_directory(fluid)
visitlab_add_directory(lsystem)
visitlab_add_directory(path_tracer)
visitlab_add_directory(nurbs_surfaces)

# Adds the available final solutions.
visitlab_add_directory(fluid_final)
visitlab_add_directory(lsystem_final)
visitlab_add_directory(path_tracer_final)
visitlab_add_directory(nurbs_surfaces_final)

# Adds the experimental projects.
visitlab_add_directory(sampling_2D)

# A helper function for generating projects when correcting homework assignments.
#get_filename_component(DIR_NAME "${CMAKE_CURRENT_SOURCE_DIR}" NAME)
#visitlab_generate_lecture(PA213 ${DIR_NAME}) 