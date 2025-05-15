// ################################################################################
// Common Framework for Computer Graphics Courses at FI MUNI.
//
// Copyright (c) Visitlab (https://visitlab.fi.muni.cz)
// All rights reserved.
// ################################################################################

#include <iostream>
#include "nurbs.hpp"

namespace nurbs
{
    //////////////////////////////////////////////////////////////////////////////////////
    /// Core functions
    //////////////////////////////////////////////////////////////////////////////////////

    /// @brief Finds the knot span for the parameter t in the knot vector U.
    /// @param t A value for the parameter of the basis functions.
    /// @param U The knot vector.
    /// @param p The degree of the basis functions.
    /// @return An index i to U s.t. [U[i],U[i+1]) is the span for t.
    /// IMPLEMENTATION:
    ///     - The search is done using the binary search in U.
    ///     - Note that the knot vector is non-periodic.
    ///     - See lecture slide 20 (part 1).
    std::uint32_t find_span(float const t, std::vector<float> const& U, std::uint32_t const p)
    {
        uint32_t low_i = p + 1;
        uint32_t high_i = U.size() - (p + 1) - 1;
        uint32_t mid_i = high_i + low_i / 2;
        float mid = U[mid_i];

        while (high_i - low_i > 1) {
            if (mid <= t && U[mid_i + 1] > t) break; // found it

            if (t < mid) {
                high_i = mid_i;
                mid_i = (low_i + high_i) / 2; // or >> 1 ?;
            } else {
                low_i = mid_i;
                mid_i = (low_i + high_i) / 2;
            }
            mid = U[mid_i];
        }

        return mid_i;
    }

    /// @brief Evaluates all basis functions which may be nonzero for the parameter
    ///        t and stores them to the output vector N.
    /// @param N The output vector. Its elements N[0],...,N[p] will hold values of the basis functions.
    /// @param t A value for the parameter of the basis functions.
    /// @param i The index of the know span in U for the parameter t (it is the result from @ref find_span function).
    /// @param U The knot vector.
    /// @param p The degree of the basis functions.
    /// IMPLEMENTATION:
    ///     - Provide an efficient implementation of the algorithm presented in the lecture (part 1).
    ///     - Lecture slides 23 and 31 (part 1) are the most important.
    ///     - Use N.reserve() to reserve (preallocate) required memory in the vector.
    ///     - Use N.push_back() to append an number to the end of the vector.
    void evaluate_basis_functions(std::vector<float>& N,
                                  float const t,
                                  std::uint32_t const i,
                                  std::vector<float> const& U,
                                  std::uint32_t const p)
    {
        std::uint32_t result_size = p + 1U; // The length of the output array (green areas in slides)
        N.reserve(result_size);
    }

    /// @brief Computes the point on the polynomial curve for given values of basis functions.
    ///        The resulting point will stay in the homogeneous space.
    /// @param k Index of a knot span; computed by the @brief find_span.
    /// @param N A vector of evaluated basis functions; computed by the @brief evaluate_basis_functions.
    /// @param P The control points of the curve.
    /// @param p The degree of the curve.
    /// @return The point on the curve in the homogeneous space.
    /// IMPLEMENTATION:
    ///     - See assignment slide 6.
    ///     - Avoid iterating over all values of 'i', i.e. avoid N_i,p that are zero.
    ///         - Hint: The relevant values of 'i' are defined by the knot span and degree of the curve.
    glm::vec4 point_on_curve_in_homogeneous_space(std::uint32_t const k, std::vector<float> const& N,
                                                  std::vector<glm::vec4> const& P, std::uint32_t const p)
    {
        return glm::vec4{0.0};
    }

    /// @brief Computes the point on the surface for given values of basis functions.
    ///        The resulting point will stay in the homogeneous space.
    /// @param k Index of a knot span in the u-direction; computed by the @brief find_span.
    /// @param N A vector of evaluated basis functions in the u-direction; computed by the @brief evaluate_basis_functions.
    /// @param l Index of a knot span in the v-direction; computed by the @brief find_span.
    /// @param M A vector of evaluated basis functions in the v-direction; computed by the @brief evaluate_basis_functions.
    /// @param P The grid of control points of the surface.
    ///          The i-the row (v-direction) consists of control points P[i][0],...,P[i][P[i].size()-1].
    ///          The j-the column (u-direction) consists of control points P[0][j],...,P[P.size()-1][j].
    /// @param p The degree of the surface in the u-direction.
    /// @param q The degree of the surface in the v-direction.
    /// @return The point on the surface in the homogeneous space.
    /// IMPLEMENTATION:
    ///     - See lecture slide 20 (part 2) and also assignment slide 6.
    ///     - There should be a single loop in the function's body.
    ///     - Inside the loop you can use the @brief point_on_curve_in_homogeneous_space function.
    ///     - Avoid iterating over all values of 'i', i.e. avoid N_{i,p} that are zero.
    ///         - Hint: The relevant values of 'i' are defined by the knot span and degree of the curve.
    glm::vec4 point_on_surface_in_homogeneous_space(std::uint32_t const k, std::vector<float> const& N,
                                                    std::uint32_t const l, std::vector<float> const& M,
                                                    std::vector<std::vector<glm::vec4>> const& P, std::uint32_t const p,
                                                    std::uint32_t const q)
    {
        return glm::vec4{0.0};
    }

    /// @brief Computes the point on the surface for given values of parameters u,v.
    ///        The resulting point will stay in the homogeneous space.
    /// @param u The parameter of the surface.
    /// @param v The parameter of the surface.
    /// @param P The grid of control points of the surface.
    ///          The i-the row (v-direction) consists of control points P[i][0],...,P[i][P[i].size()-1].
    ///          The j-the column (u-direction) consists of control points P[0][j],...,P[P.size()-1][j].
    /// @param U The knot vector in the u-direction.
    /// @param p The degree of the surface in the u-direction.
    /// @param V The knot vector in the u-direction.
    /// @param q The degree of the surface in the v-direction.
    /// @return The point on the surface in the homogeneous space.
    /// IMPLEMENTATION:
    ///     - Compute both knot spans and all possibly non-zero basis functions.
    ///          - Use @brief find_span and @brief evaluate_basis_functions.
    ///     - Then use the function @brief point_on_surface_in_homogeneous_space above.
    ///     - See lecture slide 20 (part 2).
    glm::vec4 point_on_surface_in_homogeneous_space(float u, float v, std::vector<std::vector<glm::vec4>> const& P,
                                                    std::vector<float> const& U, std::uint32_t p,
                                                    std::vector<float> const& V, std::uint32_t q)
    {
        return glm::vec4{0.0};
    }

    /// @brief Given the knot vector U it computes the know vector U' used for
    ///        for representation of derivative curve/surface.
    /// @param dU The resulting knot vector (for the derivative).
    /// @param U The knot vector of the curve/surface.
    /// IMPLEMENTATION:
    ///     - See lecture slide 31 (part 2).
    ///     - Use function dU.assign(...) to set the result.
    void derivative_knot_vector(std::vector<float>& dU, std::vector<float> const& U) {}

    /// @brief Given the control grid points P_{i,j}^w of a surface it computes the control grid
    ///        points Q_{i,j}^w of its partial u-derivative.
    /// @param dPu The resulting control grid.
    /// @param P The grid of control points of the surface.
    ///          The i-the row (v-direction) consists of control points P[i][0],...,P[i][P[i].size()-1].
    ///          The j-the column (u-direction) consists of control points P[0][j],...,P[P.size()-1][j].
    /// @param U The knot vector in the u-direction.
    /// @param p The degree of the surface in the u-direction.
    /// IMPLEMENTATION:
    ///     - See lecture slide 34 (part 2).
    ///     - Use dPu.resize(..., {}) to initialize the resulting grid with right number of empty rows.
    ///     - Use dPu.push_back() to append Q_{i,j}^w to the end of the vector.
    ///     - The implementation should likely contain a nested loop (inner and outer 'for' loop) for rows and columns.
    void derivative_control_grid_u(std::vector<std::vector<glm::vec4>>& dPu,
                                   std::vector<std::vector<glm::vec4>> const& P, std::vector<float> const& U,
                                   std::uint32_t const p) {}

    /// @brief Given the control grid points P_{i,j}^w of a surface it computes the control grid
    ///        points Q_{i,j}^w of its partial v-derivative.
    /// @param dPv The resulting control grid.
    /// @param P The grid of control points of the surface.
    ///          The i-the row (v-direction) consists of control points P[i][0],...,P[i][P[i].size()-1].
    ///          The j-the column (u-direction) consists of control points P[0][j],...,P[P.size()-1][j].
    /// @param V The knot vector in the v-direction.
    /// @param q The degree of the surface in the v-direction.
    /// IMPLEMENTATION:
    ///     - See lecture slide 33 (part 2).
    ///     - The function is quite similar to the @ref derivative_control_grid_u function, however you must compute the derivatives along the columns.
    ///     - Use dPv.resize(..., {}) to initialize the resulting grid with right number of empty rows.
    ///     - Fore each row use dPv[i].reserve() to reserve (preallocate) required memory in the column vector.
    ///     - Use dPu.push_back() to append Q_{i,j}^w to the end of the vector.
    ///     - The implementation should likely contain a nested loop (inner and outer 'for' loop) for rows and columns
    void derivative_control_grid_v(std::vector<std::vector<glm::vec4>>& dPv,
                                   std::vector<std::vector<glm::vec4>> const& P, std::vector<float> const& V,
                                   std::uint32_t const q) {}

    /// @brief Computes surface derivative S_\alpha(u,v) from the point A(u,v) and derivative
    ///        A_\alpha(u,v) in homogeneous coordinates.
    /// @param A Point of the surface in homogeneous space, at which we compute the derivative.
    /// @param dA The derivative in homogeneous space of the surface at the point A.
    /// @return The derivative of the surface at point A projected back to the affine space.
    /// IMPLEMENTATION:
    ///     - See lecture slide 32 and 34 (part 2).
    ///     - Tip: To get w(u,v) think about the relationship between rational curves represented in affine and homogeneous spaces (slide 14, part 2).
    glm::vec4 derivative_using_A(glm::vec4 const& A, glm::vec4 const& dA)
    {
        return glm::vec4{0.0};
    }

    /// @brief Computes partial u-derivative of the surface at the point corresponding
    ///        to parameters u,v. The resulting vector is already in the affine space.
    /// @param u The parameter of the surface.
    /// @param v The parameter of the surface.
    /// @param P The grid of control points of the surface.
    ///          The i-the row (v-direction) consists of control points P[i][0],...,P[i][P[i].size()-1].
    ///          The j-the column (u-direction) consists of control points P[0][j],...,P[P.size()-1][j].
    /// @param U The knot vector in the u-direction.
    /// @param p The degree of the surface in the u-direction.
    /// @param V The knot vector in the v-direction.
    /// @param q The degree of the surface in the v-direction.
    /// @param dPu The control grid of the partial u-derivative; computed by @brief derivative_control_grid_u.
    /// @param dU The knot vector of the derivative; computed by @brief derivative_knot_vector.
    /// @return The partial u-derivative of the surface in the affine space.
    /// IMPLEMENTATION:
    ///     - See lecture slide 34 (part 2) and assignment slide 6.
    ///     - There is no loop in the implementation.
    ///     - Use function @brief derivative_using_A and also @brief point_on_surface_in_homogeneous_space (twice).
    glm::vec4 derivative_of_surface_u(float const u, float const v, std::vector<std::vector<glm::vec4>> const& P,
                                      std::vector<float> const& U, std::uint32_t const p, std::vector<float> const& V,
                                      std::uint32_t const q,
                                      std::vector<std::vector<glm::vec4>> const& dPu, std::vector<float> const& dU)
    {
        return glm::vec4{0.0};
    }

    /// @brief Computes partial u-derivative of the surface at the point corresponding
    ///        to parameters u,v. The resulting vector is already in the affine space.
    /// @param u The parameter of the surface.
    /// @param v The parameter of the surface.
    /// @param P The grid of control points of the surface.
    ///          The i-the row (v-direction) consists of control points P[i][0],...,P[i][P[i].size()-1].
    ///          The j-the column (u-direction) consists of control points P[0][j],...,P[P.size()-1][j].
    /// @param U The knot vector in the u-direction.
    /// @param p The degree of the surface in the u-direction.
    /// @param V The knot vector in the v-direction.
    /// @param q The degree of the surface in the v-direction.
    /// @param dPv The control grid of the partial v-derivative; computed by @brief derivative_control_grid_v.
    /// @param dV The knot vector of the derivative; computed by @brief derivative_knot_vector.
    /// @return The partial v-derivative of the surface in the affine space.
    /// IMPLEMENTATION:
    ///     - See lecture slide 33 (part 2) and assignment slide 6.
    ///     - There is no loop in the implementation.
    ///     - Use function @brief derivative_using_A and also @brief point_on_surface_in_homogeneous_space (twice).
    glm::vec4 derivative_of_surface_v(float const u, float const v, std::vector<std::vector<glm::vec4>> const& P,
                                      std::vector<float> const& U, std::uint32_t const p, std::vector<float> const& V,
                                      std::uint32_t const q,
                                      std::vector<std::vector<glm::vec4>> const& dPv, std::vector<float> const& dV)
    {
        return glm::vec4{0.0};
    }

    //////////////////////////////////////////////////////////////////////////////////////
    /// Utility functions (aleady implemented)
    //////////////////////////////////////////////////////////////////////////////////////

    /// @brief Transforms the passed point P from the homogeneous space to the affine space.
    /// @param P A point in the homogeneous space
    /// @return The point P transformed to the affine space.
    glm::vec4 transform_point_from_homogeneous_space(glm::vec4 const& P)
    {
        return (1.0f / P.w) * P;
    }

    /// @brief Computes a grid of "points_u x points_v" points on the surface.
    ///        The computed points will be stored to the output grid G.
    /// @param points_u The number of grid points to be generated in the u direction.
    /// @param points_v The number of grid points to be generated in the v direction.
    /// @param G The output grid to be filled by computed points on the surface.
    ///          The grid is organized in the same way as the grid P (see below).
    /// @param P The grid of control points of the surface.
    ///          The i-the row (v-direction) consists of control points P[i][0],...,P[i][P[i].size()-1].
    ///          The j-the column (u-direction) consists of control points P[0][j],...,P[P.size()-1][j].
    /// @param U The knot vector.
    /// @param p The degree of the basis functions used on the know vector U.
    /// @param V The knot vector.
    /// @param q The degree of the basis functions used on the know vector V.
    /// NOTES:
    ///     - Surface parameters u,v are regularly distributed over the entire domain of the surface.
    void grid_of_surface_points(
        std::uint32_t const points_u, std::uint32_t const points_v, std::vector<std::vector<glm::vec4>>& G,
        std::vector<std::vector<glm::vec4>> const& P, std::vector<float> const& U, std::uint32_t p,
        std::vector<float> const& V, std::uint32_t q)
    {
        float const u_min{U[p]};
        float const u_max{U[U.size() - p - 1U]};
        float const v_min{V[q]};
        float const v_max{V[V.size() - q - 1U]};
        G.reserve(points_u);
        for (std::uint32_t i = 0U; i < points_u; ++i) {
            float const u{u_min + (i / (float)(points_u - 1U)) * (u_max - u_min)};
            G.push_back({});
            G.back().reserve(points_v);
            for (std::uint32_t j = 0U; j < points_v; ++j) {
                float const v{v_min + (j / (float)(points_v - 1U)) * (v_max - v_min)};
                G.back().push_back(
                    transform_point_from_homogeneous_space(point_on_surface_in_homogeneous_space(u, v, P, U, p, V, q)));
            }
        }
    }

    /// @brief Computes two grids of "vectors_u x vectors_v" partial derivatives on the surface.
    ///        The computed vectors will be stored to the output vector dSu,dSv.
    /// @param vectors_u The number of grid vectors to be generated in the u direction.
    /// @param vectors_v The number of grid vectors to be generated in the v direction.
    /// @param dSu The output grid filled in by partial u-derivatives of the surface.
    ///            The grid is organized in the same way as the grid P (see below).
    /// @param dSv The output grid filled in by partial v-derivatives of the surface.
    ///            The grid is organized in the same way as the grid P (see below).
    /// @param P The grid of control points of the surface.
    ///          The i-the row (v-direction) consists of control points P[i][0],...,P[i][P[i].size()-1].
    ///          The j-the column (u-direction) consists of control points P[0][j],...,P[P.size()-1][j].
    /// @param U The knot vector.
    /// @param p The degree of the basis functions used on the know vector U.
    /// @param V The knot vector.
    /// @param q The degree of the basis functions used on the know vector V.
    /// NOTES:
    ///     - Surface parameters u,v are regularly distributed over the entire domain of the surface.
    void grid_of_surface_derivatives(std::uint32_t const vectors_u, std::uint32_t const vectors_v,
                                     std::vector<std::vector<glm::vec4>>& dSu, std::vector<std::vector<glm::vec4>>& dSv,
                                     std::vector<std::vector<glm::vec4>> const& P,
                                     std::vector<float> const& U, std::uint32_t p, std::vector<float> const& V,
                                     std::uint32_t q)
    {
        std::vector<std::vector<glm::vec4>> dPu, dPv;
        derivative_control_grid_u(dPu, P, U, p);
        derivative_control_grid_v(dPv, P, V, q);
        std::vector<float> dU, dV;
        derivative_knot_vector(dU, U);
        derivative_knot_vector(dV, V);

        float const u_min{U[p]};
        float const u_max{U[U.size() - p - 1U]};
        float const v_min{V[q]};
        float const v_max{V[V.size() - q - 1U]};
        dSu.reserve(vectors_u);
        dSv.reserve(vectors_u);
        for (std::uint32_t i = 0U; i < vectors_u; ++i) {
            float const u{u_min + (i / (float)(vectors_u - 1U)) * (u_max - u_min)};
            dSu.push_back({});
            dSu.back().reserve(vectors_v);
            dSv.push_back({});
            dSv.back().reserve(vectors_v);
            for (std::uint32_t j = 0U; j < vectors_v; ++j) {
                float const v{v_min + (j / (float)(vectors_v - 1U)) * (v_max - v_min)};
                dSu.back().push_back(derivative_of_surface_u(u, v, P, U, p, V, q, dPu, dU));
                dSv.back().push_back(derivative_of_surface_v(u, v, P, U, p, V, q, dPv, dV));
            }
        }
    }
}
