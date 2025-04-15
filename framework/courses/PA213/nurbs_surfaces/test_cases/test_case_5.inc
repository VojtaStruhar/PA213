Test{
    TestConfig{
        SurfaceGenerationParameters{ 4, 2, 7, 3, glm::vec3{ 7.5, 5, 7.5}, (float)0.1, (float)0, (float)0.1, (float)0, (float)0.5, (float)0.5, (float)1, (float)0, false },
        SurfaceResolution{ 5, 5, false }
    },
    ExpectedResults{
        std::vector<std::vector<glm::vec4> >{ // points
            std::vector<glm::vec4>{
                glm::vec4{ -3.75, 0, -3.75, 1 }, glm::vec4{ -3.75, 0.125164, -1.35417, 1 }, glm::vec4{ -3.75, 0.195454, 0, 1 }, glm::vec4{ -3.75, 0.265074, 1.35417, 1 }, glm::vec4{ -3.75, 0.386271, 3.75, 1 }, },
            std::vector<glm::vec4>{
                glm::vec4{ -1.5625, 0.114149, -3.75, 1 }, glm::vec4{ -1.5625, 0.239313, -1.35417, 1 }, glm::vec4{ -1.5625, 0.309603, 0, 1 }, glm::vec4{ -1.5625, 0.379223, 1.35417, 1 }, glm::vec4{ -1.5625, 0.50042, 3.75, 1 }, },
            std::vector<glm::vec4>{
                glm::vec4{ 0, 0.195275, -3.75, 1 }, glm::vec4{ 0, 0.320439, -1.35417, 1 }, glm::vec4{ 0, 0.390729, 0, 1 }, glm::vec4{ 0, 0.460349, 1.35417, 1 }, glm::vec4{ 0, 0.581546, 3.75, 1 }, },
            std::vector<glm::vec4>{
                glm::vec4{ 1.5625, 0.275331, -3.75, 1 }, glm::vec4{ 1.5625, 0.400496, -1.35417, 1 }, glm::vec4{ 1.5625, 0.470785, 0, 1 }, glm::vec4{ 1.5625, 0.540405, 1.35417, 1 }, glm::vec4{ 1.5625, 0.661603, 3.75, 1 }, },
            std::vector<glm::vec4>{
                glm::vec4{ 3.75, 0.386271, -3.75, 1 }, glm::vec4{ 3.75, 0.511436, -1.35417, 1 }, glm::vec4{ 3.75, 0.581725, 0, 1 }, glm::vec4{ 3.75, 0.651345, 1.35417, 1 }, glm::vec4{ 3.75, 0.772542, 3.75, 1 }, },
        },
        std::vector<std::vector<glm::vec4> >{ // derivatives_u
            std::vector<glm::vec4>{
                glm::vec4{ 5, 0.261321, 0, 0 }, glm::vec4{ 5, 0.261321, 0, 0 }, glm::vec4{ 5, 0.261321, 0, 0 }, glm::vec4{ 5, 0.261321, 0, 0 }, glm::vec4{ 5, 0.261321, 0, 0 }, },
            std::vector<glm::vec4>{
                glm::vec4{ 3.75, 0.195275, 0, 0 }, glm::vec4{ 3.75, 0.195275, 0, 0 }, glm::vec4{ 3.75, 0.195275, 0, 0 }, glm::vec4{ 3.75, 0.195275, 0, 0 }, glm::vec4{ 3.75, 0.195275, 0, 0 }, },
            std::vector<glm::vec4>{
                glm::vec4{ 2.5, 0.129229, 0, 0 }, glm::vec4{ 2.5, 0.129229, 0, 0 }, glm::vec4{ 2.5, 0.129229, 0, 0 }, glm::vec4{ 2.5, 0.129229, 0, 0 }, glm::vec4{ 2.5, 0.129229, 0, 0 }, },
            std::vector<glm::vec4>{
                glm::vec4{ 3.75, 0.190996, 0, 0 }, glm::vec4{ 3.75, 0.190996, 0, 0 }, glm::vec4{ 3.75, 0.190996, 0, 0 }, glm::vec4{ 3.75, 0.190996, 0, 0 }, glm::vec4{ 3.75, 0.190996, 0, 0 }, },
            std::vector<glm::vec4>{
                glm::vec4{ 5, 0.252763, 0, 0 }, glm::vec4{ 5, 0.252763, 0, 0 }, glm::vec4{ 5, 0.252763, 0, 0 }, glm::vec4{ 5, 0.252763, 0, 0 }, glm::vec4{ 5, 0.252763, 0, 0 }, },
        },
        std::vector<std::vector<glm::vec4> >{ // derivatives_v
            std::vector<glm::vec4>{
                glm::vec4{ 0, 0.19626, 3.75, 0 }, glm::vec4{ 0, 0.0813717, 1.5625, 0 }, glm::vec4{ 0, 0.0646145, 1.25, 0 }, glm::vec4{ 0, 0.0798989, 1.5625, 0 }, glm::vec4{ 0, 0.188242, 3.75, 0 }, },
            std::vector<glm::vec4>{
                glm::vec4{ 0, 0.19626, 3.75, 0 }, glm::vec4{ 0, 0.0813717, 1.5625, 0 }, glm::vec4{ 0, 0.0646145, 1.25, 0 }, glm::vec4{ 0, 0.0798989, 1.5625, 0 }, glm::vec4{ 0, 0.188242, 3.75, 0 }, },
            std::vector<glm::vec4>{
                glm::vec4{ 0, 0.19626, 3.75, 0 }, glm::vec4{ 0, 0.0813717, 1.5625, 0 }, glm::vec4{ 0, 0.0646145, 1.25, 0 }, glm::vec4{ 0, 0.0798989, 1.5625, 0 }, glm::vec4{ 0, 0.188242, 3.75, 0 }, },
            std::vector<glm::vec4>{
                glm::vec4{ 0, 0.19626, 3.75, 0 }, glm::vec4{ 0, 0.0813717, 1.5625, 0 }, glm::vec4{ 0, 0.0646145, 1.25, 0 }, glm::vec4{ 0, 0.0798989, 1.5625, 0 }, glm::vec4{ 0, 0.188242, 3.75, 0 }, },
            std::vector<glm::vec4>{
                glm::vec4{ 0, 0.19626, 3.75, 0 }, glm::vec4{ 0, 0.0813717, 1.5625, 0 }, glm::vec4{ 0, 0.0646145, 1.25, 0 }, glm::vec4{ 0, 0.0798989, 1.5625, 0 }, glm::vec4{ 0, 0.188242, 3.75, 0 }, },
        }
    }
},
