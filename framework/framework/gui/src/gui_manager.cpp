// ################################################################################
// Common Framework for Computer Graphics Courses at FI MUNI.
//
// Copyright (c) Visitlab (https://visitlab.fi.muni.cz)
// All rights reserved.
// ################################################################################

#include "gui_manager.h"

#include "gui_application.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

// ----------------------------------------------------------------------------
// Variables
// ----------------------------------------------------------------------------
float ImGuiManager::xscale = 1.0;
float ImGuiManager::yscale = 1.0;

// ----------------------------------------------------------------------------
// Methods
// ----------------------------------------------------------------------------
void ImGuiManager::init(int width = 1280, int height = 720, std::string title = "Application", int major = 4, int minor = 5) {
    ApplicationManager::init(width, height, title, major, minor);
}

void ImGuiManager::pre_render_loop(IApplication& application) {
    // Setup ImGui context.
    ImGui::CreateContext();
    ImGui::StyleColorsLight();

    // Setup Platform/Renderer backends.
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init();

    // Retrieve DPI scale
    // More at: https://github.com/ocornut/imgui/blob/master/docs/FAQ.md#q-how-should-i-handle-dpi-in-my-application

    glfwGetWindowContentScale(window, &xscale, &yscale);

    ImGuiIO& io = ImGui::GetIO();
    std::filesystem::path font_path = application.get_framework_folder_path() / "gui" / "fonts" / "SourceCodePro-Regular.ttf";
    font_path.make_preferred();
    ImGui::GetStyle().ScaleAllSizes(xscale);
    ImFont* font = io.Fonts->AddFontFromFileTTF(font_path.generic_string().c_str(), xscale * 16);

    if (GUIApplication* gui_application = dynamic_cast<GUIApplication*>(&application)) {
        gui_application->set_scale(1.f / xscale);
    }
}

void ImGuiManager::pre_frame_render() {
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();
}

void ImGuiManager::post_frame_render() {
    ImGui::Render();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}