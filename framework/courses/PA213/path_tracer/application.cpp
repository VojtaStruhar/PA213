#include "application.hpp"
#include "utils.hpp"
#include <string> 

Application::Application(int initial_width, int initial_height, std::vector<std::string> arguments)
    : OpenGL45Application(initial_width, initial_height, arguments) {
    Application::compile_shaders();
    prepare_cameras();
    prepare_materials();
    prepare_scene();
    prepare_framebuffers();
}

Application::~Application() {}

// ----------------------------------------------------------------------------
// Shaderes
// ----------------------------------------------------------------------------
void Application::compile_shaders() {
    path_tracing_program = ShaderProgram(lecture_shaders_path / "full_screen_quad.vert", lecture_shaders_path / "path_tracing.frag");
    display_result_program = ShaderProgram(lecture_shaders_path / "full_screen_quad.vert", lecture_shaders_path / "display_result.frag");
    std::cout << "Shaders are reloaded." << std::endl;

    distribution = std::uniform_real_distribution<float>(0.0, 1.0);
}

// ----------------------------------------------------------------------------
// Initialize Scene
// ----------------------------------------------------------------------------
void Application::prepare_cameras() {
    // Sets the default camera position.
    camera.set_eye_position(glm::radians(0.f), glm::radians(20.f), 20.f);
    // Computes the projection matrix.
    camera_ubo.set_projection(
        glm::perspective(glm::radians(45.f), static_cast<float>(this->width) / static_cast<float>(this->height), 0.1f, 1000.0f));
    camera_ubo.update_opengl_data();
}

void Application::prepare_materials() {
    // Initializes the materials - the actual values are uploaded in the update function.
    std::vector<PathTracerMaterialData> materials(5);
    materials_ubo = PathTracerMaterialUBO(materials);
}

void Application::prepare_scene() {
    // Initializes the vector with spheres and sets the first static sphere - the rest is set in update function.
    spheres = std::vector<glm::vec4>(4);
    spheres[0] = glm::vec4(0, 2, 0, 4);
    spheres[1] = glm::vec4(-10, 2, 2, 2);
    spheres[2] = glm::vec4(5, 0, 10, 2);

    // Creates a UBO for these random positions.
    // The spheres positions are dynamic so we update them later in update scene method.
    glCreateBuffers(1, &spheres_ubo);
    glNamedBufferStorage(spheres_ubo, sizeof(glm::vec4) * 4 * spheres.size(), nullptr, GL_DYNAMIC_STORAGE_BIT);
}

void Application::prepare_framebuffers() {
    // Creates the framebuffers for rendering.
    glCreateFramebuffers(2, accumulation_fbos);
    glNamedFramebufferDrawBuffers(accumulation_fbos[0], 1, FBOUtils::draw_buffers_constants);
    glNamedFramebufferDrawBuffers(accumulation_fbos[1], 1, FBOUtils::draw_buffers_constants);

    // Creates and binds the required textures.
    resize_fullscreen_textures();
}

void Application::resize_fullscreen_textures() {
    // Removes the previously allocated textures (if any).
    glDeleteTextures(2, accumulation_textures);
    // Creates new textures for G-buffer and set their basic parameters.
    glCreateTextures(GL_TEXTURE_2D, 2, accumulation_textures);
    for (int i = 0; i < 2; ++i) {
        // Initializes the immutable storage.
        glTextureStorage2D(accumulation_textures[i], 1, GL_RGBA32F, width, height);
        // Sets the texture parameters.
        TextureUtils::set_texture_2d_parameters(accumulation_textures[i], GL_CLAMP_TO_EDGE, GL_CLAMP_TO_EDGE, GL_NEAREST, GL_NEAREST);
        // Binds the texture to the framebuffer.
        glNamedFramebufferTexture(accumulation_fbos[i], GL_COLOR_ATTACHMENT0, accumulation_textures[i], 0);
        FBOUtils::check_framebuffer_status(accumulation_fbos[i], ("accumulation buffer [" + std::to_string(i) + "]"));
    }
}

// ----------------------------------------------------------------------------
// Update
// ----------------------------------------------------------------------------
void Application::update(float delta) {
    OpenGL45Application::update(delta);

    if (!accumulate) {
        clear();
    }

    // Updates the main camera.
    const glm::vec3 eye_position = camera.get_eye_position();
    camera_ubo.set_view(lookAt(eye_position, glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(0.0f, 1.0f, 0.0f)));
    camera_ubo.update_opengl_data();

    // Updates the movable sphere(s).
    const glm::vec3 last_sphere_position =
        glm::vec3(15, 15, -12) * glm::vec3(cosf(gui_last_sphere_position / 6.0f) * sinf(gui_last_sphere_position),
                                           sinf(gui_last_sphere_position / 6.0f),
                                           cosf(gui_last_sphere_position / 6.0f) * cosf(gui_last_sphere_position));

    spheres[3] = glm::vec4(last_sphere_position, 1.0);
    glNamedBufferSubData(spheres_ubo, 0, sizeof(glm::vec4) * spheres.size(), spheres.data());

    // Updates the materials.
    materials_ubo.set_material(0, materials[0]);
    materials_ubo.set_material(1, materials[1]);
    materials_ubo.set_material(2, materials[2]);
    materials_ubo.set_material(3, materials[3]);
    materials_ubo.set_material(4, materials[4]);
    materials_ubo.update_opengl_data();
}

// ----------------------------------------------------------------------------
// Render
// ----------------------------------------------------------------------------
void Application::render() {
    // Starts measuring the elapsed time.
    glBeginQuery(GL_TIME_ELAPSED, render_time_query);

    // Increases the number of iterations.
    if(!limit_max_iterations || iterations < max_iterations){
        iterations++;
        path_tracing();
        std::swap(current_read, current_write);
    }
    display_result();

    // Resets the VAO and the program.
    glBindVertexArray(0);
    glUseProgram(0);

    // Stops measuring the elapsed time.
    glEndQuery(GL_TIME_ELAPSED);

    // Waits for OpenGL - don't forget OpenGL is asynchronous.
    glFinish();

    // Evaluates the query.
    GLuint64 render_time;
    glGetQueryObjectui64v(render_time_query, GL_QUERY_RESULT, &render_time);
    fps_gpu = 1000.f / (static_cast<float>(render_time) * 1e-6f);
}

void Application::path_tracing() {
    // Binds one of the accumulation buffers.
    glBindFramebuffer(GL_FRAMEBUFFER, accumulation_fbos[current_write]);
    glViewport(0, 0, width, height);

    // Clears the framebuffer color.
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);
    // We do not need depth and depth test.
    glDisable(GL_DEPTH_TEST);

    // Uses the proper program and sets uniform values.
    path_tracing_program.use();
    path_tracing_program.uniform("resolution", glm::vec2(width, height));
    path_tracing_program.uniform("sampling_strategy", current_sampling);
    path_tracing_program.uniform("current_brdf", current_brdf);
    path_tracing_program.uniform("spheres_count", static_cast<int>(spheres.size()));
    path_tracing_program.uniform("bounces", bounces);
    path_tracing_program.uniform("elapsed_time", static_cast<float>(elapsed_time / 1000.f));
    path_tracing_program.uniform("canonical_random_value", distribution(generator));

    // Binds the texture with accumulated color.
    glBindTextureUnit(0, accumulation_textures[current_read]);

    // Binds the data with the camera and the lights.
    camera_ubo.bind_buffer_base(CameraUBO::DEFAULT_CAMERA_BINDING);

    // Binds the buffers containing the information about the spheres (positions + radii and materials).
    glBindBufferBase(GL_UNIFORM_BUFFER, 4, spheres_ubo);
    materials_ubo.bind_buffer_base(PhongMaterialUBO::DEFAULT_MATERIAL_BINDING);

    // Renders the full screen quad to evaluate every pixel.
    // Binds an empty VAO as we do not need any state.
    glBindVertexArray(empty_vao);
    // Calls a draw command with 3 vertices that are generated in vertex shader.
    glDrawArrays(GL_TRIANGLES, 0, 3);
}

void Application::display_result() {
    // Binds the main framebuffer and clears it.
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glClear(GL_COLOR_BUFFER_BIT);

    // Disables depth test.
    glDisable(GL_DEPTH_TEST);

    // Uses the proper program and sets uniform values.
    display_result_program.use();
    display_result_program.uniform("iterations", iterations);
    display_result_program.uniform("gamma_correction", gamma_correction);
    

    // Binds the texture.
    glBindTextureUnit(0, accumulation_textures[current_read]);

    // Renders the full screen quad to evaluate every pixel.
    // Binds an empty VAO as we do not need any state.
    glBindVertexArray(empty_vao);
    // Calls a draw command with 3 vertices that are generated in vertex shader.
    glDrawArrays(GL_TRIANGLES, 0, 3);
}

void Application::clear() {
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glBindFramebuffer(GL_FRAMEBUFFER, accumulation_fbos[0]);
    glClear(GL_COLOR_BUFFER_BIT);
    glBindFramebuffer(GL_FRAMEBUFFER, accumulation_fbos[1]);
    glClear(GL_COLOR_BUFFER_BIT);

    iterations = 0;
}

// ----------------------------------------------------------------------------
// GUI
// ----------------------------------------------------------------------------
void Application::render_ui() {
    const float unit = ImGui::GetFontSize();

    ImGui::Begin("Settings", nullptr);
    ImGui::PushItemWidth(150);
    std::string fps_cpu_string = "FPS (CPU): " + std::to_string(fps_cpu);
    ImGui::Text(fps_cpu_string.c_str());
    std::string fps_string = "FPS (GPU): " + std::to_string(fps_cpu);
    ImGui::Text(fps_cpu_string.c_str());
    std::string iterations_string = "Iterations: " + std::to_string(iterations);
    ImGui::Text(iterations_string.c_str());
    ImGui::Checkbox("Show Material Box", &show_material_window[0]);
    ImGui::Checkbox("Show Material Sph#1", &show_material_window[1]);
    ImGui::Checkbox("Show Material Sph#2", &show_material_window[2]);
    ImGui::Checkbox("Show Material Sph#3", &show_material_window[3]);
    ImGui::Checkbox("Show Material Sph#4", &show_material_window[4]);
    
    ImGui::Checkbox("Accumulate", &accumulate);
    ImGui::Checkbox("Gama Correction", &gamma_correction);
    if (ImGui::Checkbox("Limit Iterations", &limit_max_iterations)){
        clear();
    }
    if (ImGui::SliderInt("Max Iterations", &max_iterations, 1, 1000)) {
        clear();
    }

    if (ImGui::Combo("Sampling", &current_sampling, SAMPLING_LABELS, IM_ARRAYSIZE(SAMPLING_LABELS))) {
        clear();
    }
    if (ImGui::Combo("BRDF", &current_brdf, BRDF_LABELS, IM_ARRAYSIZE(BRDF_LABELS))) {
        clear();
    }
    if (ImGui::SliderInt("Bounces", &bounces, 1, 15)) {
        clear();
    }
    ImGui::End();

    for (int m = 0; m < spheres.size() + 1; ++m) {
        if (show_material_window[m]) {
            render_material_window(m, unit);
        }
    }
}

void Application::render_material_window(int index, float unit) {
    std::string name;
    if (index == 0) {
        name = "Material Box";
    } else {
        std::ostringstream ss;
        ss << "Material Sph#" << index;
        name = ss.str();
    }

    ImGui::Begin(name.c_str(), nullptr);
    ImGui::PushItemWidth(150);

    if (ImGui::ColorEdit3("Albedo", reinterpret_cast<float*>(&materials[index].albedo))) {
        clear();
    }
    if (index != 0 && ImGui::ColorEdit3("Emission", reinterpret_cast<float*>(&materials[index].emission))) {
        clear();
    }
    if (ImGui::SliderFloat("Roughness", reinterpret_cast<float*>(&materials[index].roughness), 0.001f, 1.0f)) {
        clear();
    }
    if (index == 4 && ImGui::SliderAngle("Position", &gui_last_sphere_position, 0)) {
        clear();
    }
    ImGui::End();
}

// ----------------------------------------------------------------------------
// Input Events
// ----------------------------------------------------------------------------
void Application::on_mouse_button(int button, int action, int mods) {
    OpenGL45Application::on_mouse_button(button, action, mods);
    if (button == GLFW_MOUSE_BUTTON_LEFT || button == GLFW_MOUSE_BUTTON_RIGHT) {
        if (action == GLFW_PRESS) {
            mouse_pressed = true;
        } else {
            mouse_pressed = false;
        }
    }
}

void Application::on_mouse_move(double x, double y) {
    OpenGL45Application::on_mouse_move(x, y);
    if (mouse_pressed && !ImGui::GetIO().WantCaptureMouse) {
        clear();
    }
}

void Application::on_resize(int width, int height) {
    OpenGL45Application::on_resize(width, height);
    /** We make sure we recreate the textures. */
    resize_fullscreen_textures();
    clear();
}
