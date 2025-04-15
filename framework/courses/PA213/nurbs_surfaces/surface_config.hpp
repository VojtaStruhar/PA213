// ################################################################################
// Common Framework for Computer Graphics Courses at FI MUNI.
//
// Copyright (c) Visitlab (https://visitlab.fi.muni.cz)
// All rights reserved.
// ################################################################################

#pragma once

#include <glm/glm.hpp>
#include <cstdint>

struct SurfaceGenerationParameters
{
    int points_u{ 5 };
    int degree_u{ 3 };
    int points_v{ 4 };
    int degree_v{ 2 };
    glm::vec3 bbox_size{ 7.5f, 5.0f, 7.5f };
    float frequency_u{ 2.0f };
    float phase_u{ 0.0f };
    float frequency_v{ 2.0f };
    float phase_v{ 0.0f };
    float weight_u{ 0.5f };
    float weight_v{ 0.5f };
    float weight{ 10.0f };
    bool weight_invert{ false };

    bool changed{ true };
};

struct SurfaceResolution
{
    int points_u{ 50 };
    int points_v{ 50 };

    bool changed{ true };
};
