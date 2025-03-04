// ###############################################################################
// Common Framework for Computer Graphics Courses at FI MUNI.
//
// Copyright (c) Visitlab (https://visitlab.fi.muni.cz)
// All rights reserved.
// ################################################################################

#pragma once

#include "glm/glm.hpp"
#include "material_ubo.hpp"
#include "ubo_impl.hpp" // required for the UBO implementation

/**
 * The structure holding the information about a material data used in a PA213 PathTracer.
 * TODO: Move into framerwok after the course is finished.
 *
 * @author	<a href="mailto:jan.byska@gmail.com">Jan Byška</a>
 */
struct PathTracerMaterialData {
    /** The albedo color of the material (w<1 indicates transparent/glass material). */
    glm::vec4 albedo;
    /** The emission part of the material. */
    glm::vec3 emission;
    /** The roughness of the material. */
    float roughness;

    // ----------------------------------------------------------------------------
    // Constructors
    // ----------------------------------------------------------------------------
    /**
     * Constructs a @link PathTracerMaterialData with default values representing a white plastic material.
     */
    PathTracerMaterialData() : PathTracerMaterialData(glm::vec3(1.0), glm::vec3(0.0f), 1.0f, false) {}

    /**
     * Constructs a @link PathTracerMaterialData with specified values.
     * 
     * @param 	albedo	 	The material albedo color.
     * @param 	emission 	The material emission color.
     * @param 	roughness	The roughness of the material.
     */
    PathTracerMaterialData(const glm::vec3& albedo, const glm::vec3& emission, float roughness, bool glass)
        : albedo(glm::vec4(albedo.x, albedo.y, albedo.z, glass ? 0.0 : 1.0)), roughness(roughness), emission(emission) {}
};

/**
 * Contains the information about a material used in the PBR illumination models.
 *
 * Use this code in shaders for a single material:
 * <code>
 * layout(std140, binding = 3) uniform PathTracerMaterialData
 * {
 *    vec3 albedo;      // The albedo color of the material (w<1 indicates transparent/glass material).
 *    vec3 emission;    // The material emission color.
 *    float roughness;  // The roughness of the material.
 * } material;
 *
 * or a set of multiple materials
 *
 * struct PathTracerMaterialData
 * {
 *    vec3 albedo;      // The albedo color of the material (w<1 indicates transparent/glass material).
 *    vec3 emission;    // The material emission color.
 *    float roughness;  // The roughness of the material.
 * };
 * layout(std140, binding = 3) uniform PathTracerMaterialData
 * {
 *    PBRMaterialData materials[#count];
 * };
 * </code>
 */
class PathTracerMaterialUBO : public MaterialUBO<PathTracerMaterialData> {
    // ----------------------------------------------------------------------------
    // Layout Asserts
    // ----------------------------------------------------------------------------
    static_assert(offsetof(PathTracerMaterialData, albedo) == 0, "Incorrect PathTracerMaterialData layout.");
    static_assert(offsetof(PathTracerMaterialData, emission) == 16, "Incorrect PathTracerMaterialData layout.");
    static_assert(offsetof(PathTracerMaterialData, roughness) == 28, "Incorrect PathTracerMaterialData layout.");
    static_assert(sizeof(PathTracerMaterialData) == 32, "Incorrect PathTracerMaterialData layout.");

    // ----------------------------------------------------------------------------
    // Constructors
    // ----------------------------------------------------------------------------
  public:
    /** Inherits the constructors. */
    using MaterialUBO<PathTracerMaterialData>::MaterialUBO;
};
