// ################################################################################
// Common Framework for Computer Graphics Courses at FI MUNI.
//
// Copyright (c) Visitlab (https://visitlab.fi.muni.cz)
// All rights reserved.
// ################################################################################

#pragma once

#include <cstdint>
#include <glm/glm.hpp>
#include <vector>

namespace nurbs
{
	//////////////////////////////////////////////////////////////////////////////////////
	/// Core functions
	//////////////////////////////////////////////////////////////////////////////////////

	/** See the implementation in the nurbs.cpp file for details. */
	std::uint32_t find_span(float t, std::vector<float> const& U, std::uint32_t p);

	/** See the implementation in the nurbs.cpp file for details. */
	void evaluate_basis_functions(std::vector<float>& N, float t, std::uint32_t i, std::vector<float> const& U,
	                              std::uint32_t p);

	/** See the implementation in the nurbs.cpp file for details. */
	glm::vec4 point_on_curve_in_homogeneous_space(std::uint32_t const k, std::vector<float> const& N,
	                                              std::vector<glm::vec4> const& P, std::uint32_t const p);

	/** See the implementation in the nurbs.cpp file for details. */
	glm::vec4 point_on_surface_in_homogeneous_space(std::uint32_t const k_u, std::vector<float> const& N_u,
	                                                std::uint32_t const k_v, std::vector<float> const& N_v,
	                                                std::vector<std::vector<glm::vec4>> const& P,
	                                                std::uint32_t const p_u, std::uint32_t const p_v);

	/** See the implementation in the nurbs.cpp file for details. */
	glm::vec4 point_on_surface_in_homogeneous_space(float u, float v, std::vector<std::vector<glm::vec4>> const& P,
	                                                std::vector<float> const& U, std::uint32_t p_u,
	                                                std::vector<float> const& V, std::uint32_t p_v);

	/** See the implementation in the nurbs.cpp file for details. */
	void derivative_knot_vector(std::vector<float>& dU, std::vector<float> const& U);

	/** See the implementation in the nurbs.cpp file for details. */
	void derivative_control_grid_u(std::vector<std::vector<glm::vec4>>& dPu,
	                               std::vector<std::vector<glm::vec4>> const& P, std::vector<float> const& U,
	                               std::uint32_t const p);

	/** See the implementation in the nurbs.cpp file for details. */
	void derivative_control_grid_v(std::vector<std::vector<glm::vec4>>& dPv,
	                               std::vector<std::vector<glm::vec4>> const& P, std::vector<float> const& V,
	                               std::uint32_t const q);

	/** See the implementation in the nurbs.cpp file for details. */
	glm::vec4 derivative_using_A(glm::vec4 const& A, glm::vec4 const& dA);

	/** See the implementation in the nurbs.cpp file for details. */
	glm::vec4 derivative_of_surface_u(float const u, float const v, std::vector<std::vector<glm::vec4>> const& P,
	                                  std::vector<float> const& U, std::uint32_t const p, std::vector<float> const& V,
	                                  std::uint32_t const q,
	                                  std::vector<std::vector<glm::vec4>> const& dPu, std::vector<float> const& dU);

	/** See the implementation in the nurbs.cpp file for details. */
	glm::vec4 derivative_of_surface_v(float const u, float const v, std::vector<std::vector<glm::vec4>> const& P,
	                                  std::vector<float> const& U, std::uint32_t const p, std::vector<float> const& V,
	                                  std::uint32_t const q,
	                                  std::vector<std::vector<glm::vec4>> const& dPv, std::vector<float> const& dV);

	//////////////////////////////////////////////////////////////////////////////////////
	/// Utility functions (DO NOT CHANGE THEM)
	//////////////////////////////////////////////////////////////////////////////////////

	/** See the implementation in the nurbs.cpp file for details. */
	glm::vec4 transform_point_from_homogeneous_space(glm::vec4 const& P);

	/** See the implementation in the nurbs.cpp file for details. */
	void grid_of_surface_points(
		std::uint32_t const points_u, std::uint32_t const points_v, std::vector<std::vector<glm::vec4>>& G,
		std::vector<std::vector<glm::vec4>> const& P, std::vector<float> const& U, std::uint32_t p,
		std::vector<float> const& V, std::uint32_t q);

	/** See the implementation in the nurbs.cpp file for details. */
	void grid_of_surface_derivatives(std::uint32_t const vectors_u, std::uint32_t const vectors_v,
	                                 std::vector<std::vector<glm::vec4>>& dSu, std::vector<std::vector<glm::vec4>>& dSv,
	                                 std::vector<std::vector<glm::vec4>> const& P,
	                                 std::vector<float> const& U, std::uint32_t p, std::vector<float> const& V,
	                                 std::uint32_t q);

	float divide_safe(float a, float b);
}
