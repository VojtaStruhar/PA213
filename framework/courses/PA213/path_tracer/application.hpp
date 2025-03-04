#pragma once

#include "camera_ubo.hpp"
#include "light_ubo.hpp"
#include "opengl_4_5_application.hpp"
#include "path_tracer_material.hpp"
#include <random>

class Application : public OpenGL45Application {
    // ----------------------------------------------------------------------------
    // Variables (Geometry)
    // ----------------------------------------------------------------------------
  protected:
    /** The spheres described by its position (xyz) and radius (w). */
    std::vector<glm::vec4> spheres;
    // The UBO with information about individual spheres.
    GLuint spheres_ubo;
    // ----------------------------------------------------------------------------
    // Variables (Materials)
    // ----------------------------------------------------------------------------
  protected:
    /** The UBO storing the data about the materials for each sphere. */
    PathTracerMaterialUBO materials_ubo;
    // ----------------------------------------------------------------------------
    // Variables (Camera)
    // ----------------------------------------------------------------------------
  protected:
    /** The UBO storing the information about camera. */
    CameraUBO camera_ubo;

    // ----------------------------------------------------------------------------
    // Variables (Shaders)
    // ----------------------------------------------------------------------------
  protected:
    /** The path tracing program used for accumulating the radiance in the scene. */
    ShaderProgram path_tracing_program;
    /** The program for rendering the final result. */
    ShaderProgram display_result_program;

    /** The random value generator. */
    std::default_random_engine generator;
    /** Random number distribution that produces floating-point values according to a uniform distribution. */
    std::uniform_real_distribution<float> distribution;

    // ----------------------------------------------------------------------------
    // Variables (Frame Buffers)
    // ----------------------------------------------------------------------------
  protected:
    /** The framebuffer objects that is used for accumulation. */
    GLuint accumulation_fbos[2];
    /** The texture used in @link accumulation_fbos to store the results of the path tracing. */
    GLuint accumulation_textures[2];
    /* The index of the current accumulation buffer to read from into. */
    int current_read = 0;
    /* The index of the current accumulation buffer to write into. */
    int current_write = 1;
    /** The number of iterations */
    int iterations = 0;

    // ----------------------------------------------------------------------------
    // Variables (GUI)
    // ----------------------------------------------------------------------------
  protected:
    /** The light position set in the GUI. */
    float gui_last_sphere_position = glm::radians(360.f);
    /** The flag determining weather a mouse button was pressed and thus we may need to clear the accumulation buffers. */
    bool mouse_pressed = false;
    /** The number of bounces. */
    int bounces = 2;
    /** The flag determining wether to apply gamma correction. */
    bool gamma_correction = true;
    /** The flag determining wether we accumulate the values. */
    bool accumulate = true;
    /** The flag determining wether the GUI for chaning the material should be displayed. */
    bool show_material_window[5] = {false,false,false,false,false};
    /** The maximum number of iterations. */
    bool limit_max_iterations = false;
    /** The maximum number of iterations. */
    int max_iterations = 100;

    /** The current sampling strategy. */
    int current_sampling = 0;
    /** The labels for the sampling strategies. */
    const char* SAMPLING_LABELS[4] = {"Uniform", "Cosine", "GGX", "Mix"};

    /** The current BRDF function. */
    int current_brdf = 0;
    /** The labels for the BRDFs. */
    const char* BRDF_LABELS[4] = {"Lambert", "LambertWithAlbedo", "GGX", "Mix"};

    /** The materials for each sphere */
    PathTracerMaterialData materials[5] = {
        PathTracerMaterialData(glm::vec3(1.0f), glm::vec3(0.0f), 1.0f, false),
        PathTracerMaterialData(glm::vec3(1.0f), glm::vec3(1.0f), 1.0f, false),
        PathTracerMaterialData(glm::vec3(1.0f, 0.0f, 0.0f), glm::vec3(0.0f), 1.0f, false),
        PathTracerMaterialData(glm::vec3(0.0f, 1.0f, 0.0f), glm::vec3(0.0f), 1.0f, false),
        PathTracerMaterialData(glm::vec3(0.0f, 0.0f, 1.0f), glm::vec3(0.0f), 1.0f, false)
    };

    // ----------------------------------------------------------------------------
    // Constructors
    // ----------------------------------------------------------------------------
  public:
    Application(int initial_width, int initial_height, std::vector<std::string> arguments = {});

    /** Destroys the {@link Application} and releases the allocated resources. */
    virtual ~Application();

    // ----------------------------------------------------------------------------
    // Shaders
    // ----------------------------------------------------------------------------
    /**
     * {@copydoc OpenGL45Application::compile_shaders}
     */
    void compile_shaders() override;

    // ----------------------------------------------------------------------------
    // Initialize Scene
    // ----------------------------------------------------------------------------
  public:
    /** Prepares the required cameras. */
    void prepare_cameras();
    /** Prepares the required materials. */
    void prepare_materials();
    /** Prepares the scene objects. */
    void prepare_scene();
    /** Prepares the frame buffer objects. */
    void prepare_framebuffers();
    /** Resizes the full screen textures match the window. */
    void resize_fullscreen_textures();

    // ----------------------------------------------------------------------------
    // Update
    // ----------------------------------------------------------------------------
    /**
     * {@copydoc OpenGL45Application::update}
     */
    void update(float delta) override;

    // ----------------------------------------------------------------------------
    // Render
    // ----------------------------------------------------------------------------
  public:
    /** @copydoc OpenGL45Application::render */
    void render() override;

    /** Uses path tracing on GPU to render the scene. */
    void path_tracing();

    /** Renders the averaged result onto screen. */
    void display_result();

    /** Clears the accumulated buffers. */
    void clear();

    // ----------------------------------------------------------------------------
    // GUI
    // ----------------------------------------------------------------------------
  public:
    /** @copydoc OpenGL45Application::render_ui */
    void render_ui() override;

    /** Renders a window for setting material. */
    void render_material_window(int index, float unit);

    // ----------------------------------------------------------------------------
    // Input Events
    // ----------------------------------------------------------------------------
  public:
    /** @copydoc OpenGL45Application::on_mouse_button */
    void on_mouse_button(int button, int action, int mods) override;

    /** @copydoc OpenGL45Application::on_mouse_move */
    void on_mouse_move(double x, double y) override;

    /** @copydoc OpenGL45Application::on_resize */
    void on_resize(int width, int height) override;
};
