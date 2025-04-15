// ################################################################################
// Common Framework for Computer Graphics Courses at FI MUNI.
//
// Copyright (c) Visitlab (https://visitlab.fi.muni.cz)
// All rights reserved.
// ################################################################################

#pragma once

#include <glm/glm.hpp>
#include <cstdint>
#include <vector>
#include <memory>

namespace nurbs {

class Surface;
using SurfacePtr = std::shared_ptr<Surface>;

class Surface
{
    std::vector<std::vector<glm::vec4> > P;
    std::vector<float> U;
    std::vector<float> V;
    std::uint32_t p;
    std::uint32_t q;

public:

    Surface(
        std::vector<std::vector<glm::vec4> > const& control_points,
        std::vector<float> const& know_vector_u,
        std::vector<float> const& know_vector_v,
        std::uint32_t degree_u,
        std::uint32_t degree_v
        );

    std::vector<std::vector<glm::vec4> > const& control_points() const { return P; }
    std::vector<float> const& know_vector_u() const { return U; }
    std::vector<float> const& know_vector_v() const { return V; }
    std::uint32_t degree_u() const { return p; }
    std::uint32_t degree_v() const { return q; }

    static SurfacePtr make_generic(
        std::uint32_t points_u, std::uint32_t degree_u,
        std::uint32_t points_v, std::uint32_t degree_v,
        glm::vec3 const& bbox_size,
        float frequency_u, float phase_u,
        float frequency_v, float phase_v,
        std::uint32_t weight_u, std::uint32_t weight_v, float weight
        );
};

}
