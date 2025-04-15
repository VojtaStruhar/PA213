// ################################################################################
// Common Framework for Computer Graphics Courses at FI MUNI.
//
// Copyright (c) Visitlab (https://visitlab.fi.muni.cz)
// All rights reserved.
// ################################################################################

#pragma once

#include "surface_config.hpp"
#include <glm/glm.hpp>
#include <filesystem>
#include <cstdint>
#include <vector>

struct Tests
{
    struct TestConfig
    {
        SurfaceGenerationParameters params{};
        SurfaceResolution resolution{};
    };

    struct CurrentResults
    {
        std::vector<std::vector<glm::vec4> > points{};
        std::vector<std::vector<glm::vec4> > derivatives_u{};
        std::vector<std::vector<glm::vec4> > derivatives_v{};
    };

    using ExpectedResults = CurrentResults;

    struct Test
    {
        TestConfig config{};
        ExpectedResults expected{};
    };

    static Tests& instance();

    bool run() const;
    bool run_one(Test const& test) const;

    std::vector<Test> const& tests() const { return m_tests; }

private:

    Tests();

    Tests(Tests const&);
    Tests(Tests&&);
    Tests& operator=(Tests const&);
    Tests& operator=(Tests&&);

    std::vector<Test> m_tests;
};

void __generate_test_case_includes(std::filesystem::path const& output_dir, Tests::TestConfig const& config);
