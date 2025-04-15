// ################################################################################
// Common Framework for Computer Graphics Courses at FI MUNI.
//
// Copyright (c) Visitlab (https://visitlab.fi.muni.cz)
// All rights reserved.
// ################################################################################

#pragma once

#include "pa213_application.h"
#include "surface.hpp"
#include "surface_config.hpp"
#include "camera.h"
#include <memory>

using GeometryPtr = std::shared_ptr<Geometry>;

class Application : public PA213Application {

    struct DrawFlags
    {
        bool control_points{ true };
        int surface{ 0 };
        bool derivatives_u{ false };
        bool derivatives_v{ false };
        bool normals{ false };
        bool axes{ true };
    };

    SurfaceGenerationParameters surface_gen_params;
    SurfaceResolution surface_res;
    DrawFlags draw_flags;

    nurbs::SurfacePtr surface;

    GeometryPtr geometry_control_points;
    GeometryPtr geometry_surface;
    GeometryPtr geometry_surface_pretty;
    GeometryPtr geometry_derivatives_u;
    GeometryPtr geometry_derivatives_v;
    GeometryPtr geometry_normals;

    Geometry geometry_axes;

    Camera camera;
    ShaderProgram shader_simple;
    ShaderProgram shader_pretty;
    GLuint texture;

    void build_control_grid_geometry();
    void build_surface_geometries();

public:

    Application(int initial_width, int initial_height, std::vector<std::string> arguments = {});
    ~Application();

    // ----------------------------------------------------------------------------
    // Shaders
    // ----------------------------------------------------------------------------
    /**
     * {@copydoc PA213Application::compile_shaders}
     */
    void compile_shaders() override;

    // ----------------------------------------------------------------------------
    // Update
    // ----------------------------------------------------------------------------
    /**
     * {@copydoc PA213Application::update}
     */
    void update(float delta) override;

    // ----------------------------------------------------------------------------
    // Render
    // ----------------------------------------------------------------------------
    /** @copydoc PA213Application::render */
    void render() override;

    // ----------------------------------------------------------------------------
    // Input Events
    // ----------------------------------------------------------------------------

    void on_resize(int width, int height) override;
    void on_mouse_move(double x, double y) override;
    void on_mouse_button(int button, int action, int mods) override;

    // ----------------------------------------------------------------------------
    // GUI
    // ----------------------------------------------------------------------------
    /** @copydoc PA213Application::render_ui */
    void render_ui() override;
};
