Test{
    TestConfig{
        SurfaceGenerationParameters{ 5, 1, 4, 1, glm::vec3{ 7.5, 5, 7.5}, (float)2, (float)0, (float)2, (float)0, (float)0.5, (float)0.5, (float)10, false, false },
        SurfaceResolution{ 5, 5, false }
    },
    ExpectedResults{
        std::vector<std::vector<glm::vec4> >{ // points
            std::vector<glm::vec4>{
                glm::vec4{ -3.75, 0, -3.75, 1 }, glm::vec4{ -3.75, 0.811899, -1.875, 1 }, glm::vec4{ -3.75, -5.96046e-08, 0, 1 }, glm::vec4{ -3.75, -0.811899, 1.875, 1 }, glm::vec4{ -3.75, 2.18557e-07, 3.75, 1 }, },
            std::vector<glm::vec4>{
                glm::vec4{ -1.875, 1.25, -3.75, 1 }, glm::vec4{ -1.875, 2.0619, -1.875, 1 }, glm::vec4{ -1.875, 1.25, 0, 1 }, glm::vec4{ -1.875, 0.438101, 1.875, 1 }, glm::vec4{ -1.875, 1.25, 3.75, 1 }, },
            std::vector<glm::vec4>{
                glm::vec4{ 0, -1.09278e-07, -3.75, 1 }, glm::vec4{ 0, 1.04761, -1.33065, 1 }, glm::vec4{ 0, 0.885708, -1.02273, 1 }, glm::vec4{ 0, -0.811899, 1.875, 1 }, glm::vec4{ 0, 1.09278e-07, 3.75, 1 }, },
            std::vector<glm::vec4>{
                glm::vec4{ 1.875, -1.25, -3.75, 1 }, glm::vec4{ 1.875, -0.438101, -1.875, 1 }, glm::vec4{ 1.875, -1.25, 0, 1 }, glm::vec4{ 1.875, -2.0619, 1.875, 1 }, glm::vec4{ 1.875, -1.25, 3.75, 1 }, },
            std::vector<glm::vec4>{
                glm::vec4{ 3.75, 2.18557e-07, -3.75, 1 }, glm::vec4{ 3.75, 0.811899, -1.875, 1 }, glm::vec4{ 3.75, 1.78814e-07, 0, 1 }, glm::vec4{ 3.75, -0.811899, 1.875, 1 }, glm::vec4{ 3.75, 4.37114e-07, 3.75, 1 }, },
        },
        std::vector<std::vector<glm::vec4> >{ // derivatives_u
            std::vector<glm::vec4>{
                glm::vec4{ 1.875, 1.25, 0, 0 }, glm::vec4{ 1.875, 1.25, 0, 0 }, glm::vec4{ 1.875, 1.25, 0, 0 }, glm::vec4{ 1.875, 1.25, 0, 0 }, glm::vec4{ 1.875, 1.25, 0, 0 }, },
            std::vector<glm::vec4>{
                glm::vec4{ 1.875, -1.25, 0, 0 }, glm::vec4{ 14.5312, -7.86073, 4.21875, 0 }, glm::vec4{ 10.3125, -2.00361, -5.625, 0 }, glm::vec4{ 1.875, -1.25, 0, 0 }, glm::vec4{ 1.875, -1.25, 0, 0 }, },
            std::vector<glm::vec4>{
                glm::vec4{ 1.875, -1.25, 0, 0 }, glm::vec4{ 0.241935, -0.191705, -0.0702394, 0 }, glm::vec4{ 0.340909, -0.38831, 0.18595, 0 }, glm::vec4{ 1.875, -1.25, 0, 0 }, glm::vec4{ 1.875, -1.25, 0, 0 }, },
            std::vector<glm::vec4>{
                glm::vec4{ 1.875, 1.25, 0, 0 }, glm::vec4{ 1.875, 1.25, 0, 0 }, glm::vec4{ 1.875, 1.25, 0, 0 }, glm::vec4{ 1.875, 1.25, 0, 0 }, glm::vec4{ 1.875, 1.25, 0, 0 }, },
            std::vector<glm::vec4>{
                glm::vec4{ 1.875, 1.25, 0, 0 }, glm::vec4{ 1.875, 1.25, 0, 0 }, glm::vec4{ 1.875, 1.25, 0, 0 }, glm::vec4{ 1.875, 1.25, 0, 0 }, glm::vec4{ 1.875, 1.25, 0, 0 }, },
        },
        std::vector<std::vector<glm::vec4> >{ // derivatives_v
            std::vector<glm::vec4>{
                glm::vec4{ 0, 1.08253, 2.5, 0 }, glm::vec4{ 0, 1.08253, 2.5, 0 }, glm::vec4{ 0, -2.16506, 2.5, 0 }, glm::vec4{ 0, 1.08253, 2.5, 0 }, glm::vec4{ 0, 1.08253, 2.5, 0 }, },
            std::vector<glm::vec4>{
                glm::vec4{ 0, 1.08253, 2.5, 0 }, glm::vec4{ 0, 1.08253, 2.5, 0 }, glm::vec4{ 0, -2.16506, 2.5, 0 }, glm::vec4{ 0, 1.08253, 2.5, 0 }, glm::vec4{ 0, 1.08253, 2.5, 0 }, },
            std::vector<glm::vec4>{
                glm::vec4{ 0, 10.8253, 25, 0 }, glm::vec4{ 0, 0.180234, 0.416233, 0 }, glm::vec4{ 0, -0.715724, 0.826446, 0 }, glm::vec4{ 0, 1.08253, 2.5, 0 }, glm::vec4{ 0, 1.08253, 2.5, 0 }, },
            std::vector<glm::vec4>{
                glm::vec4{ 0, 1.08253, 2.5, 0 }, glm::vec4{ 0, 1.08253, 2.5, 0 }, glm::vec4{ 0, -2.16506, 2.5, 0 }, glm::vec4{ 0, 1.08253, 2.5, 0 }, glm::vec4{ 0, 1.08253, 2.5, 0 }, },
            std::vector<glm::vec4>{
                glm::vec4{ 0, 1.08253, 2.5, 0 }, glm::vec4{ 0, 1.08253, 2.5, 0 }, glm::vec4{ 0, -2.16506, 2.5, 0 }, glm::vec4{ 0, 1.08253, 2.5, 0 }, glm::vec4{ 0, 1.08253, 2.5, 0 }, },
        }
    }
},
