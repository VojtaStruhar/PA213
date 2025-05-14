Test{
    TestConfig{
        SurfaceGenerationParameters{ 4, 2, 7, 3, glm::vec3{ 7.5, 5, 7.5}, (float)2, (float)0, (float)2, (float)0, (float)0.5, (float)0.5, (float)1, false, false },
        SurfaceResolution{ 5, 5, false }
    },
    ExpectedResults{
        std::vector<std::vector<glm::vec4> >{ // points
            std::vector<glm::vec4>{
                glm::vec4{ -3.75, 0, -3.75, 1 }, glm::vec4{ -3.75, 0.90211, -1.35417, 1 }, glm::vec4{ -3.75, -1.04308e-07, 0, 1 }, glm::vec4{ -3.75, -0.90211, 1.35417, 1 }, glm::vec4{ -3.75, 2.18557e-07, 3.75, 1 }, },
            std::vector<glm::vec4>{
                glm::vec4{ -1.5625, 0.541266, -3.75, 1 }, glm::vec4{ -1.5625, 1.44338, -1.35417, 1 }, glm::vec4{ -1.5625, 0.541266, 0, 1 }, glm::vec4{ -1.5625, -0.360844, 1.35417, 1 }, glm::vec4{ -1.5625, 0.541266, 3.75, 1 }, },
            std::vector<glm::vec4>{
                glm::vec4{ 0, -5.96046e-08, -3.75, 1 }, glm::vec4{ 0, 0.90211, -1.35417, 1 }, glm::vec4{ 0, -1.19209e-07, 0, 1 }, glm::vec4{ 0, -0.90211, 1.35417, 1 }, glm::vec4{ 0, 1.78814e-07, 3.75, 1 }, },
            std::vector<glm::vec4>{
                glm::vec4{ 1.5625, -0.541266, -3.75, 1 }, glm::vec4{ 1.5625, 0.360844, -1.35417, 1 }, glm::vec4{ 1.5625, -0.541266, 0, 1 }, glm::vec4{ 1.5625, -1.44338, 1.35417, 1 }, glm::vec4{ 1.5625, -0.541266, 3.75, 1 }, },
            std::vector<glm::vec4>{
                glm::vec4{ 3.75, 2.18557e-07, -3.75, 1 }, glm::vec4{ 3.75, 0.90211, -1.35417, 1 }, glm::vec4{ 3.75, 1.3411e-07, 0, 1 }, glm::vec4{ 3.75, -0.90211, 1.35417, 1 }, glm::vec4{ 3.75, 4.37114e-07, 3.75, 1 }, },
        },
        std::vector<std::vector<glm::vec4> >{ // derivatives_u
            std::vector<glm::vec4>{
                glm::vec4{ 5, 2.16506, 0, 0 }, glm::vec4{ 5, 2.16506, 0, 0 }, glm::vec4{ 5, 2.16506, 0, 0 }, glm::vec4{ 5, 2.16506, 0, 0 }, glm::vec4{ 5, 2.16506, 0, 0 }, },
            std::vector<glm::vec4>{
                glm::vec4{ 3.75, 0, 0, 0 }, glm::vec4{ 3.75, 0, 0, 0 }, glm::vec4{ 3.75, 0, 0, 0 }, glm::vec4{ 3.75, -2.38419e-07, 0, 0 }, glm::vec4{ 3.75, 0, 0, 0 }, },
            std::vector<glm::vec4>{
                glm::vec4{ 2.5, -2.16506, 0, 0 }, glm::vec4{ 2.5, -2.16506, 0, 0 }, glm::vec4{ 2.5, -2.16506, 0, 0 }, glm::vec4{ 2.5, -2.16506, 0, 0 }, glm::vec4{ 2.5, -2.16506, 0, 0 }, },
            std::vector<glm::vec4>{
                glm::vec4{ 3.75, 3.57628e-07, 0, 0 }, glm::vec4{ 3.75, 3.57628e-07, 0, 0 }, glm::vec4{ 3.75, 2.38419e-07, 0, 0 }, glm::vec4{ 3.75, 1.19209e-07, 0, 0 }, glm::vec4{ 3.75, 3.57628e-07, 0, 0 }, },
            std::vector<glm::vec4>{
                glm::vec4{ 5, 2.16506, 0, 0 }, glm::vec4{ 5, 2.16506, 0, 0 }, glm::vec4{ 5, 2.16506, 0, 0 }, glm::vec4{ 5, 2.16506, 0, 0 }, glm::vec4{ 5, 2.16506, 0, 0 }, },
        },
        std::vector<std::vector<glm::vec4> >{ // derivatives_v
            std::vector<glm::vec4>{
                glm::vec4{ 0, 3.2476, 3.75, 0 }, glm::vec4{ 0, -0.541266, 1.5625, 0 }, glm::vec4{ 0, -1.08253, 1.25, 0 }, glm::vec4{ 0, -0.541266, 1.5625, 0 }, glm::vec4{ 0, 3.2476, 3.75, 0 }, },
            std::vector<glm::vec4>{
                glm::vec4{ 0, 3.2476, 3.75, 0 }, glm::vec4{ 0, -0.541266, 1.5625, 0 }, glm::vec4{ 0, -1.08253, 1.25, 0 }, glm::vec4{ 0, -0.541266, 1.5625, 0 }, glm::vec4{ 0, 3.2476, 3.75, 0 }, },
            std::vector<glm::vec4>{
                glm::vec4{ 0, 3.2476, 3.75, 0 }, glm::vec4{ 0, -0.541266, 1.5625, 0 }, glm::vec4{ 0, -1.08253, 1.25, 0 }, glm::vec4{ 0, -0.541266, 1.5625, 0 }, glm::vec4{ 0, 3.2476, 3.75, 0 }, },
            std::vector<glm::vec4>{
                glm::vec4{ 0, 3.2476, 3.75, 0 }, glm::vec4{ 0, -0.541266, 1.5625, 0 }, glm::vec4{ 0, -1.08253, 1.25, 0 }, glm::vec4{ 0, -0.541266, 1.5625, 0 }, glm::vec4{ 0, 3.2476, 3.75, 0 }, },
            std::vector<glm::vec4>{
                glm::vec4{ 0, 3.2476, 3.75, 0 }, glm::vec4{ 0, -0.541266, 1.5625, 0 }, glm::vec4{ 0, -1.08253, 1.25, 0 }, glm::vec4{ 0, -0.541266, 1.5625, 0 }, glm::vec4{ 0, 3.2476, 3.75, 0 }, },
        }
    }
},
