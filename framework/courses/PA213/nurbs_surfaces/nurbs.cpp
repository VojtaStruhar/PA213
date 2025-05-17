// ################################################################################
// Common Framework for Computer Graphics Courses at FI MUNI.
//
// Copyright (c) Visitlab (https://visitlab.fi.muni.cz)
// All rights reserved.
// ################################################################################

#include <iostream>
#include "nurbs.hpp"

#define PRINT_VECTOR(vec) \
    std::cout << "vec{"; \
    for (int j = 0; j < vec.size(); j++) { \
        std::cout << N_u[j] << (j < vec.size() - 1 ? ", " : "");\
    }\
    std::cout << "}" std::endl;


#define PRINT(stuff) std::cout << "[DEBUG] " << stuff << std::endl;

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
        uint32_t low_i = p;
        uint32_t high_i = U.size() - p - 1;

        if (t == U[high_i]) {
            PRINT(high_i - 1)
            return high_i - 1; // return the lower index of the two
        }

        while (high_i - low_i > 1) {
            uint32_t mid_i = low_i + (high_i - low_i) / 2;

            if (t < U[mid_i]) {
                high_i = mid_i;
            } else {
                low_i = mid_i;
            }
        }
        return low_i;
    }

    /// @brief Evaluates all basis functions which may be nonzero for the parameter
    ///        t and stores them to the output vector N.
    /// @param N The output vector. Its elements N[0],...,N[p] will hold values of the basis functions.
    /// @param t A value for the parameter of the basis functions.
    /// @param i The index of the knot span in U for the parameter t (it is the result from @related find_span function).
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
        N.clear();
        N.reserve(p + 1);
        N.push_back(1.0f);

        std::vector<float> left(p + 1);
        std::vector<float> right(p + 1);

        for (std::uint32_t j = 1; j <= p; ++j) {
            left[j] = t - U[i + 1 - j];
            right[j] = U[i + j] - t;
            float saved = 0.0f;
            for (std::uint32_t r = 0; r < j; ++r) {
                float temp = N[r] / (right[r + 1] + left[j - r]);
                N[r] = saved + right[r + 1] * temp;
                saved = left[j - r] * temp;
            }
            N.push_back(saved);
        }
    }

    /// @brief Computes the point on the polynomial curve for given values of basis functions.
    ///        The resulting point will stay in the homogeneous space.
    /// @param k Index of a knot span; computed by the @related find_span.
    /// @param N A vector of evaluated basis functions; computed by the @related evaluate_basis_functions.
    /// @param P The control points of the curve.
    /// @param p The degree of the curve.
    /// @return The point on the curve in the homogeneous space.
    /// IMPLEMENTATION:
    ///     - See assignment slide 6.
    ///     - Avoid iterating over all values of 'i', i.e. avoid N_i,p that are zero.
    ///         - Hint: The relevant values of 'i' are defined by the knot span and degree of the curve.
    glm::vec4 point_on_curve_in_homogeneous_space(std::uint32_t const k,
                                                  std::vector<float> const& N,
                                                  std::vector<glm::vec4> const& P,
                                                  std::uint32_t const p)
    {
        glm::vec4 point(0.0f, 0.0f, 0.0f, 0.0f);

        for (std::uint32_t i = 0; i <= p; ++i) {
            std::uint32_t point_index = k - p + i;

            point += N[i] * P[point_index];
        }

        return point;
    }

    /// @brief Computes the point on the surface for given values of basis functions.
    ///        The resulting point will stay in the homogeneous space.
    /// @param k_u Index of a knot span in the u-direction; computed by the @ref find_span .
    /// @param N_u A vector of evaluated basis functions in the u-direction; computed by the @ref evaluate_basis_functions.
    /// @param k_v Index of a knot span in the v-direction; computed by the @ref find_span.
    /// @param N_v A vector of evaluated basis functions in the v-direction; computed by the @ref evaluate_basis_functions.
    /// @param P The grid of control points of the surface.
    ///          The i-the row (v-direction) consists of control points P[i][0],...,P[i][P[i].size()-1].
    ///          The j-the column (u-direction) consists of control points P[0][j],...,P[P.size()-1][j].
    /// @param p_u The degree of the surface in the u-direction.
    /// @param p_v The degree of the surface in the v-direction.
    /// @return The point on the surface in the homogeneous space.
    /// IMPLEMENTATION:
    ///     - See lecture slide 20 (part 2) and also assignment slide 6.
    ///     - There should be a single loop in the function's body.
    ///     - Inside the loop you can use the @related point_on_curve_in_homogeneous_space function.
    ///     - Avoid iterating over all values of 'i', i.e. avoid N_{i,p} that are zero.
    ///         - Hint: The relevant values of 'i' are defined by the knot span and degree of the curve.
    glm::vec4 point_on_surface_in_homogeneous_space(std::uint32_t const k_u,
                                                    std::vector<float> const& N_u,
                                                    std::uint32_t const k_v,
                                                    std::vector<float> const& N_v,
                                                    std::vector<std::vector<glm::vec4>> const& P,
                                                    std::uint32_t const p_u,
                                                    std::uint32_t const p_v)
    {
        glm::vec4 result(0.0f, 0.0f, 0.0f, 0.0f);

        for (std::uint32_t i = 0; i <= p_u; ++i) {
            // same as previous
            std::uint32_t point_index = k_u - p_u + i;
            glm::vec4 curve_point = point_on_curve_in_homogeneous_space(k_v, N_v, P[point_index], p_v);

            result += N_u[i] * curve_point;
        }
        return result;
    }

    /// Computes the point on the surface for given values of parameters u,v.
    /// The resulting point will stay in the homogeneous space.
    /// /// IMPLEMENTATION:
    ///     - Compute both knot spans and all possibly non-zero basis functions.
    ///          - Use @related find_span and @related evaluate_basis_functions.
    ///     - Then use the function @related point_on_surface_in_homogeneous_space above.
    ///     - See lecture slide 20 (part 2).
    /// @param u The parameter of the surface.
    /// @param v The parameter of the surface.
    /// @param P The grid of control points of the surface.
    ///          The i-th row (v-direction) consists of control points P[i][0],...,P[i][P[i].size()-1].
    ///          The j-th column (u-direction) consists of control points P[0][j],...,P[P.size()-1][j].
    /// @param U The knot vector in the u-direction.
    /// @param p_u The degree of the surface in the u-direction.
    /// @param V The knot vector in the u-direction.
    /// @param p_v The degree of the surface in the v-direction.
    /// @return The point on the surface in the homogeneous space.
    glm::vec4 point_on_surface_in_homogeneous_space(float u,
                                                    float v,
                                                    std::vector<std::vector<glm::vec4>> const& P,
                                                    std::vector<float> const& U,
                                                    std::uint32_t p_u,
                                                    std::vector<float> const& V,
                                                    std::uint32_t p_v)
    {
        auto k_u = find_span(u, U, p_u);
        auto k_v = find_span(v, V, p_v);

        std::vector<float> N_u;
        evaluate_basis_functions(N_u, u, k_u, U, p_u);

        std::vector<float> N_v;
        evaluate_basis_functions(N_v, v, k_v, V, p_v);

        return point_on_surface_in_homogeneous_space(k_u,
                                                     N_u,
                                                     k_v,
                                                     N_v,
                                                     P,
                                                     p_u,
                                                     p_v);
    }

    /// @brief Given the knot vector U it computes the know vector U' used for
    ///        for representation of derivative curve/surface.
    /// @param dU The resulting knot vector (for the derivative).
    /// @param U The knot vector of the curve/surface.
    /// IMPLEMENTATION:
    ///     - See lecture slide 31 (part 2).
    ///     - Use function dU.assign(...) to set the result.
    void derivative_knot_vector(std::vector<float>& dU, std::vector<float> const& U)
    {
        dU.assign(U.begin() + 1, U.end() - 1);
    }

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
                                   std::uint32_t const p)
    {
        std::size_t n_rows = P.size();
        std::size_t n_cols = P[0].size();

        dPu.resize(n_rows);

        // TODO: <= ???, also rename to Y
        for (std::size_t i = 0; i < n_rows - 1; ++i) {
            // TODO: <= ???, also rename to X
            for (std::size_t j = 0; j < n_cols; ++j) {
                std::size_t u_idx1 = j + 1;
                std::size_t u_idx2 = j + p + 1;

                glm::vec4 d_point = static_cast<float>(p) * (P[i + 1][j] - P[i][j]);
                float denom = U[u_idx2] - U[u_idx1];

                if (std::abs(denom) < 1e-6f) denom = 1e-6f;

                d_point /= denom;

                dPu[i].push_back(d_point);
            }
        }
    }

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
    ///     - The function is quite similar to the @related derivative_control_grid_u function, however you must compute the derivatives along the columns.
    ///     - Use dPv.resize(..., {}) to initialize the resulting grid with right number of empty rows.
    ///     - Fore each row use dPv[i].reserve() to reserve (preallocate) required memory in the column vector.
    ///     - Use dPu.push_back() to append Q_{i,j}^w to the end of the vector.
    ///     - The implementation should likely contain a nested loop (inner and outer 'for' loop) for rows and columns
    void derivative_control_grid_v(std::vector<std::vector<glm::vec4>>& dPv,
                                   std::vector<std::vector<glm::vec4>> const& P,
                                   std::vector<float> const& V,
                                   std::uint32_t const q)
    {
        std::size_t n_rows = P.size();
        std::size_t n_cols = P[0].size();

        dPv.resize(n_rows);

        // TODO: <= ???, also rename to Y
        for (std::size_t i = 0; i < n_rows; ++i) {
            // TODO: <= ???, also rename to X
            for (std::size_t j = 0; j < n_cols - 1; ++j) {
                std::size_t u_idx1 = j + 1;
                std::size_t u_idx2 = j + q + 1;

                glm::vec4 d_point = static_cast<float>(q) * (P[i][j + 1] - P[i][j]);
                float denom = V[u_idx2] - V[u_idx1];

                if (std::abs(denom) < 1e-6f) denom = 1e-6f;

                d_point /= denom;

                dPv[i].push_back(d_point);
            }
        }
    }

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
        float w_uv = A.w;
        float dw_uv = dA.w;


        return (1.0f / w_uv) * dA - (dw_uv / (w_uv * w_uv)) * A;
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
    /// @param dPu The control grid of the partial u-derivative; computed by @related derivative_control_grid_u.
    /// @param dU The knot vector of the derivative; computed by @related derivative_knot_vector.
    /// @return The partial u-derivative of the surface in the affine space.
    /// IMPLEMENTATION:
    ///     - See lecture slide 34 (part 2) and assignment slide 6.
    ///     - There is no loop in the implementation.
    ///     - Use function @related derivative_using_A and also @related point_on_surface_in_homogeneous_space (twice).
    glm::vec4 derivative_of_surface_u(float const u,
                                      float const v,
                                      std::vector<std::vector<glm::vec4>> const& P,
                                      std::vector<float> const& U,
                                      std::uint32_t const p,
                                      std::vector<float> const& V,
                                      std::uint32_t const q,
                                      std::vector<std::vector<glm::vec4>> const& dPu,
                                      std::vector<float> const& dU)
    {
        const glm::vec4 A = point_on_surface_in_homogeneous_space(u, v, P, U, p, V, q);
        const glm::vec4 dA = point_on_surface_in_homogeneous_space(u, v, dPu, dU, p - 1, V, q);
        return derivative_using_A(A, dA);
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
    /// @param dPv The control grid of the partial v-derivative; computed by @related derivative_control_grid_v.
    /// @param dV The knot vector of the derivative; computed by @related derivative_knot_vector.
    /// @return The partial v-derivative of the surface in the affine space.
    /// IMPLEMENTATION:
    ///     - See lecture slide 33 (part 2) and assignment slide 6.
    ///     - There is no loop in the implementation.
    ///     - Use function @related derivative_using_A and also @related point_on_surface_in_homogeneous_space (twice).
    glm::vec4 derivative_of_surface_v(float const u,
                                      float const v,
                                      std::vector<std::vector<glm::vec4>> const& P,
                                      std::vector<float> const& U,
                                      std::uint32_t const p,
                                      std::vector<float> const& V,
                                      std::uint32_t const q,
                                      std::vector<std::vector<glm::vec4>> const& dPv,
                                      std::vector<float> const& dV)
    {
        // same as above
        const glm::vec4 A = point_on_surface_in_homogeneous_space(u, v, P, U, p, V, q);
        const glm::vec4 dA = point_on_surface_in_homogeneous_space(u, v, dPv, U, p, dV, q - 1);
        return derivative_using_A(A, dA);
    }

    //////////////////////////////////////////////////////////////////////////////////////
    /// TODO: Utility functions (aleady implemented)
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

    inline float divide_safe(float a, float b)
    {
        if (b == 0.0f) return 0.f;
        return a / b;
    }
}
