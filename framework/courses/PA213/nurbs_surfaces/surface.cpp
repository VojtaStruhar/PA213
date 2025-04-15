// ################################################################################
// Common Framework for Computer Graphics Courses at FI MUNI.
//
// Copyright (c) Visitlab (https://visitlab.fi.muni.cz)
// All rights reserved.
// ################################################################################

#include "surface.hpp"
#include <glm/ext.hpp>

namespace nurbs {

Surface::Surface(
    std::vector<std::vector<glm::vec4> > const& control_points,
    std::vector<float> const& know_vector_u,
    std::vector<float> const& know_vector_v,
    std::uint32_t degree_u,
    std::uint32_t degree_v
    )
    : P{ control_points }
    , U{ know_vector_u }
    , V{ know_vector_v }
    , p{ degree_u }
    , q{ degree_v }
{}

SurfacePtr Surface::make_generic(
    std::uint32_t const points_u, std::uint32_t const degree_u,
    std::uint32_t const points_v, std::uint32_t const degree_v,
    glm::vec3 const& bbox_size,
    float const frequency_u, float const phase_u,
    float const frequency_v, float const phase_v,
    std::uint32_t const weight_u, std::uint32_t const weight_v, float const weight
    )
{
    if (!(
        degree_u > 0U && points_u > degree_u &&
        degree_v > 0U && points_v > degree_v &&
        bbox_size.x > 0.01f && bbox_size.y >= 0.0f && bbox_size.z > 0.01f &&
        frequency_u >= 0.0f && frequency_v >= 0.0f &&
        weight_u < points_u && weight_v < points_v &&
        weight > 0.001f
        ))
        return nullptr;

    std::vector<std::vector<glm::vec4> > P;
    {
        P.reserve(points_u);
        glm::vec3 Qu{ -0.5f * bbox_size.x, 0.0f, -0.5f * bbox_size.z };
        glm::vec3 const D{ (1.0f / (points_u - 1U)) * bbox_size.x, 0.5f * bbox_size.y, (1.0f / (points_v - 1U)) * bbox_size.z };
        for (std::uint32_t i = 0U; i < points_u; ++i)
        {
            P.push_back({});
            P.back().reserve(points_v);
            glm::vec3 Qv{ Qu };
            for (std::uint32_t j = 0U; j < points_v; ++j)
            {
                glm::vec4 Q{ Qv, 1.0f };
                Q.y = (0.5f * (glm::sin(glm::pi<float>() * (phase_u + frequency_u * i / (points_u - 1U))) + 
                               glm::sin(glm::pi<float>() * (phase_v + frequency_v * j / (points_v - 1U)))))
                      * D.y;
                P.back().push_back(Q);
                Qv.z += D.z;
            }
            Qu.x += D.x;
        }
        P[weight_u][weight_v] *= weight;
    }

    std::vector<float> U;
    std::vector<float> V;
    {
        auto const& init_knot_vector = [](std::vector<float>& U, std::uint32_t const p, std::uint32_t const n) {
            std::uint32_t const m{ n + p + 1U };
            U.resize(m + 1U);
            for (std::uint32_t i = 0U; i <= p; ++i)
                U[i] = 0.0f;
            float value{ 1.0f };
            for (std::uint32_t i = p + 1U; i <= m - p - 1U; ++i, value += 1.0f)
                U[i] = value;
            for (std::uint32_t i = 0U; i <= p; ++i)
                U[m - p + i] = value;
        };
        init_knot_vector(U, degree_u, points_u - 1U);
        init_knot_vector(V, degree_v, points_v - 1U);
    }
    return std::make_shared<Surface>(P, U, V, degree_u, degree_v);
}

}
