// ################################################################################
// Common Framework for Computer Graphics Courses at FI MUNI.
//
// Copyright (c) Visitlab (https://visitlab.fi.muni.cz)
// All rights reserved.
// ################################################################################

#include "application.hpp"
#include "utils.hpp"
#include "nurbs.hpp"
#include "tests.hpp"

#pragma warning(disable : 4996)

static void geometry_push_back_point(std::vector<float>& output, glm::vec3 const& data)
{
    output.push_back(data.x);
    output.push_back(data.y);
    output.push_back(data.z);
}

static void geometry_push_back_line(std::vector<float>& output, glm::vec3 const& a, glm::vec3 const& b)
{
    geometry_push_back_point(output, a);
    geometry_push_back_point(output, b);
}

static void geometry_push_back_triangle(std::vector<float>& output, glm::vec3 const& a, glm::vec3 const& b, glm::vec3 const& c)
{
    geometry_push_back_point(output, a);
    geometry_push_back_point(output, b);
    geometry_push_back_point(output, c);
}

static void geometry_push_back_uvs(std::vector<float>& output, glm::vec2 const& data)
{
    output.push_back(data.x);
    output.push_back(data.y);
}

static void geometry_push_back_uv_triangle(std::vector<float>& output, glm::vec2 const& a, glm::vec2 const& b, glm::vec2 const& c)
{
    geometry_push_back_uvs(output, a);
    geometry_push_back_uvs(output, b);
    geometry_push_back_uvs(output, c);
}

static GeometryPtr make_pretty_geometry_from_grid(
    std::vector<std::vector<glm::vec4> > const& P,
    std::vector<std::vector<glm::vec4> > const& N,
    std::vector<std::vector<glm::vec2> > const& T
    )
{
    return GeometryPtr{ new Geometry{
        GL_TRIANGLES,
        [&P]() -> std::vector<float> // vertices
        {
            std::vector<float> vertices;
            for (std::size_t i = 1ULL; i < P.size(); ++i)
                for (std::size_t j = 1ULL; j < P[i].size(); ++j)
                {
                    geometry_push_back_triangle(vertices, P[i-1][j-1], P[i-1][j], P[i][j-1]);
                    geometry_push_back_triangle(vertices, P[i-1][j], P[i][j-1], P[i][j]);
                }
            return vertices;
        }(),
        {}, // indices
        [&N]() -> std::vector<float> // normals
        {
            std::vector<float> normals;
            for (std::size_t i = 1ULL; i < N.size(); ++i)
                for (std::size_t j = 1ULL; j < N[i].size(); ++j)
                {
                    geometry_push_back_triangle(normals, N[i-1][j-1], N[i-1][j], N[i][j-1]);
                    geometry_push_back_triangle(normals, N[i-1][j], N[i][j-1], N[i][j]);
                }
            return normals;
        }(),
        {}, // colors
        [&T]() -> std::vector<float> // tex_coords
        {
            std::vector<float> uvs;
            for (std::size_t i = 1ULL; i < T.size(); ++i)
                for (std::size_t j = 1ULL; j < T[i].size(); ++j)
                {
                    geometry_push_back_uv_triangle(uvs, T[i-1][j-1], T[i-1][j], T[i][j-1]);
                    geometry_push_back_uv_triangle(uvs, T[i-1][j], T[i][j-1], T[i][j]);
                }
            return uvs;
        }()
    }};
}

static GeometryPtr make_geometry_from_grid(std::vector<std::vector<glm::vec4> > const& G, glm::vec3 const& color)
{
    return GeometryPtr{ new Geometry{
        GL_LINES,
        [&G]() -> std::vector<float> // vertices
        {
            std::vector<float> vertices;
            for (std::size_t i = 0ULL; i < G.size(); ++i)
                for (std::size_t j = 1ULL; j < G[i].size(); ++j)
                    geometry_push_back_line(vertices, G[i][j-1], G[i][j]);
            for (std::size_t i = 1ULL; i < G.size(); ++i)
                for (std::size_t j = 0ULL; j < G[i].size(); ++j)
                    geometry_push_back_line(vertices, G[i-1][j], G[i][j]);
            return vertices;
        }(),
        {}, // indices
        {}, // normals
        [&G, &color]() -> std::vector<float> // colors
        {
            std::vector<float> colors;
            for (std::size_t i = 0ULL; i < G.size(); ++i)
                for (std::size_t j = 1ULL; j < G[i].size(); ++j)
                    geometry_push_back_line(colors, color, color);
            for (std::size_t i = 1ULL; i < G.size(); ++i)
                for (std::size_t j = 0ULL; j < G[i].size(); ++j)
                    geometry_push_back_line(colors, color, color);
            return colors;
        }()
    }};
}

static GeometryPtr make_geometry_from_vector_field(
    std::vector<std::vector<glm::vec4> > const& S,
    std::vector<std::vector<glm::vec4> > const& dS,
    glm::vec3 const& color)
{
    return GeometryPtr{ new Geometry{
        GL_LINES,
        [&S, &dS]() -> std::vector<float> // vertices
        {
            std::vector<float> vertices;
            for (std::size_t i = 0ULL; i < S.size(); ++i)
                for (std::size_t j = 0ULL; j < S[i].size(); ++j)
                    geometry_push_back_line(vertices, S[i][j], S[i][j] + dS[i][j]);
            return vertices;
        }(),
        {}, // indices
        {}, // normals
        [&S, &color]() -> std::vector<float> // colors
        {
            std::vector<float> colors;
            for (std::size_t i = 0ULL; i < S.size(); ++i)
                for (std::size_t j = 0ULL; j < S[i].size(); ++j)
                    geometry_push_back_line(colors, color, color);
            return colors;
        }()
    }};
}

void Application::build_control_grid_geometry()
{
    std::vector<std::vector<glm::vec4> > G{ surface->control_points() };
    for (auto& row : G)
        for (auto& Q : row)
            Q = nurbs::transform_point_from_homogeneous_space(Q);
    geometry_control_points = make_geometry_from_grid(G, glm::vec3{ 1.0f, 1.0f, 0.5f });
}

void Application::build_surface_geometries()
{
    std::vector<std::vector<glm::vec4> > P;
    nurbs::grid_of_surface_points(surface_res.points_u, surface_res.points_v, P,
        surface->control_points(),
        surface->know_vector_u(), surface->degree_u(),
        surface->know_vector_v(), surface->degree_v()
        );

    std::vector<std::vector<glm::vec4> > dSu, dSv;
    nurbs::grid_of_surface_derivatives(
        surface_res.points_u, surface_res.points_v,
        dSu, dSv,
        surface->control_points(),
        surface->know_vector_u(), surface->degree_u(),
        surface->know_vector_v(), surface->degree_v()
        );

    std::vector<std::vector<glm::vec4> > N;
    N.resize(surface_res.points_u, {});
    for (std::uint32_t i = 0ULL; i != surface_res.points_u; ++i)
        for (std::uint32_t j = 0ULL; j != surface_res.points_v; ++j)
            N[i].push_back(-glm::vec4(glm::normalize(glm::cross(glm::vec3(dSu[i][j]), glm::vec3(dSv[i][j]))), 0.0f));

    std::vector<std::vector<glm::vec2> > T;
    T.resize(surface_res.points_u, {});
    for (std::uint32_t i = 0ULL; i != surface_res.points_u; ++i)
        for (std::uint32_t j = 0ULL; j != surface_res.points_v; ++j)
            T[i].push_back(glm::vec2((float)i / (surface_res.points_u - 1U),(float)j / (surface_res.points_v - 1U)));

    geometry_surface = make_geometry_from_grid(P, glm::vec3{ 0.75f, 0.75f, 0.75f });
    geometry_surface_pretty = make_pretty_geometry_from_grid(P,N,T);
    geometry_derivatives_u = make_geometry_from_vector_field(P, dSu, glm::vec3{ 0.75f, 0.25f, 0.25f });
    geometry_derivatives_v = make_geometry_from_vector_field(P, dSv, glm::vec3{ 0.25f, 0.75f, 0.25f });
    geometry_normals = make_geometry_from_vector_field(P, N, glm::vec3{ 0.25f, 0.25f, 0.75f });
}

Application::Application(int initial_width, int initial_height, std::vector<std::string> arguments)
    : PA213Application(initial_width, initial_height, arguments)
    , surface_gen_params{}
    , surface_res{}
    , draw_flags{}
    , surface{ nullptr }
    , geometry_control_points{ nullptr }
    , geometry_surface{ nullptr }
    , geometry_surface_pretty{ nullptr }
    , geometry_derivatives_u{ nullptr }
    , geometry_derivatives_v{ nullptr }
    , geometry_normals{ nullptr }
    , geometry_axes{
        GL_LINES,
        []() -> std::vector<float> // vertices
        {
            float const half_size{ 5.0f };
            return std::vector<float>{
                // x-axis
                -half_size, 0.0f, 0.0f,
                half_size, 0.0f, 0.0f,
                // y-axis
                0.0f, -half_size, 0.0f,
                0.0f,  half_size, 0.0f,
                // z-axis
                0.0f, 0.0f, -half_size,
                0.0f, 0.0f,  half_size,
                // diagonal xy
                1.0f, 0.0f, 0.0f,
                0.0f,  1.0f, 0.0f,
                // diagonal xz
                1.0f, 0.0f, 0.0f,
                0.0f, 0.0f,  1.0f,
                // diagonal yz
                0.0f,  1.0f, 0.0f,
                0.0f, 0.0f,  1.0f,
            };
        }(),
        {}, // indices
        {}, // normals
        []() -> std::vector<float> // colors
        {
            return std::vector<float>{
                // x-axis (red)
                1.0f, 0.0f, 0.0f,
                1.0f, 0.0f, 0.0f,
                // y-axis (green)
                0.0f,  1.0f, 0.0f,
                0.0f,  1.0f, 0.0f,
                // z-axis (blue)
                0.0f, 0.0f,  1.0f,
                0.0f, 0.0f,  1.0f,
                // diagonal xy
                1.0f, 0.0f, 0.0f,
                0.0f,  1.0f, 0.0f,
                // diagonal xz
                1.0f, 0.0f, 0.0f,
                0.0f, 0.0f,  1.0f,
                // diagonal yz
                0.0f,  1.0f, 0.0f,
                0.0f, 0.0f,  1.0f,
            };
        }()
    }
    , camera{}
    , shader_simple{ lecture_shaders_path / "simple.vert", lecture_shaders_path / "simple.frag" }
    , shader_pretty{  lecture_shaders_path / "pretty.vert", lecture_shaders_path / "pretty.frag"  }
    , texture{ TextureUtils::load_texture_2d(lecture_textures_path / "uv_grid.png") }
{
    camera.set_eye_position(-3.14f / 4.0f, 3.14f / 4.0f, 15);
    TextureUtils::set_texture_2d_parameters(texture, GL_REPEAT, GL_REPEAT, GL_LINEAR, GL_LINEAR);
}

Application::~Application()
{
    glDeleteTextures(1, &texture);
}

// ----------------------------------------------------------------------------
// Shaderes
// ----------------------------------------------------------------------------
void Application::compile_shaders() {
    shader_simple = ShaderProgram(lecture_shaders_path / "simple.vert", lecture_shaders_path / "simple.frag");
    shader_pretty = ShaderProgram(lecture_shaders_path / "pretty.vert", lecture_shaders_path / "pretty.frag");
    std::cout << "Shaders are reloaded." << std::endl;
}

// ----------------------------------------------------------------------------
// Update
// ----------------------------------------------------------------------------
void Application::update(float delta) {
    PA213Application::update(delta);

    if (surface_gen_params.changed) {
        surface = nurbs::Surface::make_generic(
            (std::uint32_t)surface_gen_params.points_u, (std::uint32_t)surface_gen_params.degree_u,
            (std::uint32_t)surface_gen_params.points_v, (std::uint32_t)surface_gen_params.degree_v,
            surface_gen_params.bbox_size,
            surface_gen_params.frequency_u, surface_gen_params.phase_u,
            surface_gen_params.frequency_v, surface_gen_params.phase_v,
            (std::uint32_t)(surface_gen_params.weight_u * (surface_gen_params.points_u - 1)),
            (std::uint32_t)(surface_gen_params.weight_v * (surface_gen_params.points_v - 1)),
            surface_gen_params.weight_invert ? 1.0f / surface_gen_params.weight : surface_gen_params.weight
            );
        if (surface != nullptr) {
            build_control_grid_geometry();
            build_surface_geometries();
        } else {
            geometry_control_points = nullptr;
            geometry_surface = nullptr;
            geometry_surface_pretty = nullptr;
            geometry_derivatives_u = nullptr;
            geometry_derivatives_v = nullptr;
            geometry_normals = nullptr;
        }
        surface_gen_params.changed = false;
        surface_res.changed = false;
    }
    if (surface_res.changed) {
        if (surface != nullptr)
            build_surface_geometries();
        surface_res.changed = false;
    }
}

// ----------------------------------------------------------------------------
// Render
// ----------------------------------------------------------------------------
void Application::render() {
    glm::mat4 const projection_matrix = glm::perspective(glm::radians(45.0f), static_cast<float>(width) / static_cast<float>(height), 0.1f, 100.0f);
    glm::mat4 const view_matrix = camera.get_view_matrix();
    glm::mat4 const model_matrix = glm::identity<glm::mat4>();

    glClearColor(0.25f, 0.25f, 0.25f, 1.0f);
    glClearDepth(1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);

    shader_simple.use();
    glUniformMatrix4fv(shader_simple.get_uniform_location("projection"), 1, GL_FALSE, value_ptr(projection_matrix));
    glUniformMatrix4fv(shader_simple.get_uniform_location("view"), 1, GL_FALSE, value_ptr(view_matrix));
    glUniformMatrix4fv(shader_simple.get_uniform_location("model"), 1, GL_FALSE, value_ptr(model_matrix));
    if (draw_flags.control_points && geometry_control_points != nullptr) geometry_control_points->draw();
    if (draw_flags.surface == 0 && geometry_surface != nullptr) geometry_surface->draw();
    if (draw_flags.derivatives_u && geometry_derivatives_u != nullptr) geometry_derivatives_u->draw();
    if (draw_flags.derivatives_v && geometry_derivatives_v != nullptr) geometry_derivatives_v->draw();
    if (draw_flags.normals && geometry_normals != nullptr) geometry_normals->draw();
    if (draw_flags.axes) geometry_axes.draw();

    if (draw_flags.surface > 0 && geometry_surface != nullptr)
    {
        shader_pretty.use();
        glUniformMatrix4fv(shader_pretty.get_uniform_location("projection"), 1, GL_FALSE, value_ptr(projection_matrix));
        glUniformMatrix4fv(shader_pretty.get_uniform_location("view"), 1, GL_FALSE, value_ptr(view_matrix));
        glUniformMatrix4fv(shader_pretty.get_uniform_location("model"), 1, GL_FALSE, value_ptr(model_matrix));
        glUniform3fv(shader_pretty.get_uniform_location("eye_position"), 1, value_ptr(camera.get_eye_position()));
        glUniform1i(shader_pretty.get_uniform_location("use_texture"), draw_flags.surface == 2 ? 1 : 0);
        glBindTexture(GL_TEXTURE_2D, texture);
        geometry_surface_pretty->draw();
    }
}

void Application::on_resize(int width, int height) { IApplication::on_resize(width, height); glViewport(0, 0, width, height); }
void Application::on_mouse_move(double x, double y) { if (!ImGui::GetIO().WantCaptureMouse) camera.on_mouse_move(x, y); }
void Application::on_mouse_button(int button, int action, int mods) { if (!ImGui::GetIO().WantCaptureMouse) camera.on_mouse_button(button, action, mods); }

// ----------------------------------------------------------------------------
// GUI
// ----------------------------------------------------------------------------
void Application::render_ui() {
    const float unit = ImGui::GetFontSize();
    ImGui::Begin("Settings", nullptr, ImGuiWindowFlags_NoDecoration);
    ImGui::SetWindowSize(ImVec2(15 * unit, 36.25f * unit));
    ImGui::SetWindowPos(ImVec2(0 * unit, 0 * unit));
    if (ImGui::SliderInt("Grid-u", &surface_gen_params.points_u, 2, 20)) surface_gen_params.changed = true;
    if (ImGui::SliderInt("Degree-u", &surface_gen_params.degree_u, 1, 10)) surface_gen_params.changed = true;
    if (ImGui::SliderInt("Grid-v", &surface_gen_params.points_v, 2, 20)) surface_gen_params.changed = true;
    if (ImGui::SliderInt("Degree-v", &surface_gen_params.degree_v, 1, 10)) surface_gen_params.changed = true;
    if (ImGui::SliderFloat3("Size", glm::value_ptr(surface_gen_params.bbox_size), 0.1f, 10.0f)) surface_gen_params.changed = true;
    if (ImGui::SliderFloat("Freq-u", &surface_gen_params.frequency_u, 0.1f, 10.0f)) surface_gen_params.changed = true;
    if (ImGui::SliderFloat("Phase-u", &surface_gen_params.phase_u, 0.0f, 2.0f * glm::pi<float>())) surface_gen_params.changed = true;
    if (ImGui::SliderFloat("Freq-v", &surface_gen_params.frequency_v, 0.1f, 10.0f)) surface_gen_params.changed = true;
    if (ImGui::SliderFloat("Phase-v", &surface_gen_params.phase_v, 0.0f, 2.0f * glm::pi<float>())) surface_gen_params.changed = true;
    if (ImGui::SliderFloat("Weight-u", &surface_gen_params.weight_u, 0.0f, 1.0f)) surface_gen_params.changed = true;
    if (ImGui::SliderFloat("Weight-v", &surface_gen_params.weight_v, 0.0f, 1.0f)) surface_gen_params.changed = true;
    if (ImGui::SliderFloat("Weight", &surface_gen_params.weight, 1.0f, 100.0f)) surface_gen_params.changed = true;
    if (ImGui::Checkbox("Weight^-1", &surface_gen_params.weight_invert)) surface_gen_params.changed = true;
    if (ImGui::SliderInt("Res-u", &surface_res.points_u, 2, 100)) surface_res.changed = true;
    if (ImGui::SliderInt("Res-v", &surface_res.points_v, 2, 100)) surface_res.changed = true;
    ImGui::Checkbox("Control points", &draw_flags.control_points);
    ImGui::SliderInt("Surface", &draw_flags.surface, 0, 2);
    ImGui::Checkbox("Derivatives-u", &draw_flags.derivatives_u);
    ImGui::Checkbox("Derivatives-v", &draw_flags.derivatives_v);
    ImGui::Checkbox("Normals", &draw_flags.normals);
    ImGui::Checkbox("Axes", &draw_flags.axes);
    if (ImGui::Button("Run tests")) Tests::instance().run();
    ImGui::End();
}
