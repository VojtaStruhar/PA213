// ################################################################################
// Common Framework for Computer Graphics Courses at FI MUNI.
//
// Copyright (c) Visitlab (https://visitlab.fi.muni.cz)
// All rights reserved.
// ################################################################################

#include "tests.hpp"
#include "surface.hpp"
#include "nurbs.hpp"
#include <iostream>
#include <fstream>

static void compute_current_results(Tests::CurrentResults& current, Tests::TestConfig const& config)
{
    auto const surface = nurbs::Surface::make_generic(
        (std::uint32_t)config.params.points_u, (std::uint32_t)config.params.degree_u,
        (std::uint32_t)config.params.points_v, (std::uint32_t)config.params.degree_v,
        config.params.bbox_size,
        config.params.frequency_u, config.params.phase_u,
        config.params.frequency_v, config.params.phase_v,
        (std::uint32_t)(config.params.weight_u * (config.params.points_u - 1)),
        (std::uint32_t)(config.params.weight_v * (config.params.points_v - 1)),
        config.params.weight_invert ? 1.0f / config.params.weight : config.params.weight
        );

    nurbs::grid_of_surface_points(config.resolution.points_u, config.resolution.points_v, current.points,
        surface->control_points(),
        surface->know_vector_u(), surface->degree_u(),
        surface->know_vector_v(), surface->degree_v()
        );

    nurbs::grid_of_surface_derivatives(
        config.resolution.points_u, config.resolution.points_v,
        current.derivatives_u, current.derivatives_v,
        surface->control_points(),
        surface->know_vector_u(), surface->degree_u(),
        surface->know_vector_v(), surface->degree_v()
        );
}

static bool compare_float(float const& current, float const& expected)
{
    return glm::abs(expected - current) <= 0.001f;
}

static bool compare_vec(glm::vec3 const& current, glm::vec3 const& expected)
{
    if (!compare_float(current.x, expected.x))
        return false;
    if (!compare_float(current.y, expected.y))
        return false;
    if (!compare_float(current.z, expected.z))
        return false;
    return true;
}

static bool compare_results(std::vector<std::vector<glm::vec4> > const& current, std::vector<std::vector<glm::vec4> > const& expected)
{
    if (current.size() != expected.size())
        return false;
    for (std::size_t i = 0ULL; i < current.size(); ++i)
    {
        auto const& crow{ current.at(i) };
        auto const& erow{ expected.at(i) };
        if (crow.size() != erow.size())
            return false;
        for (std::size_t j = 0ULL; j < crow.size(); ++j)
            if (!compare_vec(crow.at(j), erow.at(j)))
                return false;
    }
    return true;
}

static bool compare_results(Tests::CurrentResults const& current, Tests::CurrentResults const& expected)
{
    if (!compare_results(current.points, expected.points))
        return false;
    if (!compare_results(current.derivatives_u, expected.derivatives_u))
        return false;
    if (!compare_results(current.derivatives_v, expected.derivatives_v))
        return false;
    return true;
}

Tests& Tests::instance()
{
    static Tests tests;
    return tests;
}

Tests::Tests()
    : m_tests{
#include "./test_cases/test_case_0.inc"
#include "./test_cases/test_case_1.inc"
#include "./test_cases/test_case_2.inc"
#include "./test_cases/test_case_3.inc"
#include "./test_cases/test_case_4.inc"
#include "./test_cases/test_case_5.inc"
#include "./test_cases/test_case_6.inc"
#include "./test_cases/test_case_7.inc"
    }
{}

bool Tests::run() const
{
    std::cout << "*** Starting tests ********************************\n";
    std::cout.flush();

    bool result{ true };
    for (std::size_t i = 0ULL; i < tests().size(); ++i)
    {
        std::cout << "Test " << i << ": ";
        std::cout.flush();

        bool const res{ run_one(tests().at(i)) };

        std::cout << (res ? "OK" : "FAILED") << "\n";
        std::cout.flush();

        result = result && res;
    }
    std::cout << "SUMMARY: " << (result ? "All tests PASSED." : "Some test(s) FAILED!") << " \n";
    std::cout.flush();

    return result;
}

bool Tests::run_one(Test const& test) const
{
    CurrentResults current{};
    compute_current_results(current, test.config);
    return compare_results(current, test.expected);
}


void __generate_test_case_includes(std::filesystem::path const& output_dir, Tests::TestConfig const& config)
{
    std::ofstream ostr{ (output_dir / "test_case.inc").string().c_str() };
    
    Tests::CurrentResults current;
    compute_current_results(current, config);

    auto const& dump_grid = [&ostr](std::vector<std::vector<glm::vec4> > const& G) {
        auto const& dump_point = [&ostr](glm::vec4 const& v) {
            ostr << "glm::vec4{ " << v.x << ", " << v.y << ", " << v.z << ", " << v.w << " }";
        };
        for (auto const& R : G)
        {
            ostr << "            std::vector<glm::vec4>{\n                ";
            for (auto const& v : R)
            {
                dump_point(v);
                ostr << ", ";
            }
            ostr << "},\n";
        }
    };

    ostr << "Test{\n";
    ostr << "    TestConfig{\n";
    ostr << "        SurfaceGenerationParameters{ "
         << config.params.points_u << ", "
         << config.params.degree_u << ", "
         << config.params.points_v << ", "
         << config.params.degree_v << ", "
         << "glm::vec3{ " << config.params.bbox_size.x << ", "
                 << config.params.bbox_size.y << ", "
                 << config.params.bbox_size.z << "}, "
         << "(float)" << config.params.frequency_u << ", "
         << "(float)" << config.params.phase_u << ", "
         << "(float)" << config.params.frequency_v << ", "
         << "(float)" << config.params.phase_v << ", "
         << "(float)" << config.params.weight_u << ", "
         << "(float)" << config.params.weight_v << ", "
         << "(float)" << config.params.weight << ", "
         << "(float)" << config.params.weight_invert << ", "
         << "false },\n";
    ostr << "        SurfaceResolution{ "
         << config.resolution.points_u << ", "
         << config.resolution.points_v << ", "
         << "false }\n";
    ostr << "    },\n";
    ostr << "    ExpectedResults{\n";
    ostr << "        std::vector<std::vector<glm::vec4> >{ // points\n";
    dump_grid(current.points);
    ostr << "        },\n";
    ostr << "        std::vector<std::vector<glm::vec4> >{ // derivatives_u\n";
    dump_grid(current.derivatives_u);
    ostr << "        },\n";
    ostr << "        std::vector<std::vector<glm::vec4> >{ // derivatives_v\n";
    dump_grid(current.derivatives_v);
    ostr << "        }\n";
    ostr << "    }\n";
    ostr << "},\n";
}
