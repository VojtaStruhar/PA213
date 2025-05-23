Test{
    TestConfig{
        SurfaceGenerationParameters{ 5, 4, 4, 3, glm::vec3{ 7.5, 5, 7.5}, (float)2, (float)0, (float)2, (float)0, (float)0.5, (float)0.5, (float)10, false, false },
        SurfaceResolution{ 5, 5, false }
    },
    ExpectedResults{
        std::vector<std::vector<glm::vec4> >{ // points
            std::vector<glm::vec4>{
                glm::vec4{ -3.75, 0, -3.75, 1 }, glm::vec4{ -3.75, 0.304462, -1.875, 1 }, glm::vec4{ -3.75, -2.4827e-09, 0, 1 }, glm::vec4{ -3.75, -0.304462, 1.875, 1 }, glm::vec4{ -3.75, 2.18557e-07, 3.75, 1 }, },
            std::vector<glm::vec4>{
                glm::vec4{ -1.875, 0.46875, -3.75, 1 }, glm::vec4{ -1.04114, 0.910774, -1.59705, 1 }, glm::vec4{ -1.09527, 0.723996, -0.519823, 1 }, glm::vec4{ -1.47991, 0.357775, 1.21652, 1 }, glm::vec4{ -1.875, 0.46875, 3.75, 1 }, },
            std::vector<glm::vec4>{
                glm::vec4{ 0, -1.61425e-08, -3.75, 1 }, glm::vec4{ 0, 0.761523, -1.50786, 1 }, glm::vec4{ 0, 0.604725, -0.698276, 1 }, glm::vec4{ 0, 0.141948, 0.869205, 1 }, glm::vec4{ 0, 1.46529e-07, 3.75, 1 }, },
            std::vector<glm::vec4>{
                glm::vec4{ 1.875, -0.46875, -3.75, 1 }, glm::vec4{ 1.04114, 0.390202, -1.59705, 1 }, glm::vec4{ 1.09527, 0.176364, -0.519823, 1 }, glm::vec4{ 1.47991, -0.382181, 1.21652, 1 }, glm::vec4{ 1.875, -0.46875, 3.75, 1 }, },
            std::vector<glm::vec4>{
                glm::vec4{ 3.75, 2.18557e-07, -3.75, 1 }, glm::vec4{ 3.75, 0.304462, -1.875, 1 }, glm::vec4{ 3.75, 2.03651e-07, 0, 1 }, glm::vec4{ 3.75, -0.304462, 1.875, 1 }, glm::vec4{ 3.75, 4.37114e-07, 3.75, 1 }, },
        },
        std::vector<std::vector<glm::vec4> >{ // derivatives_u
            std::vector<glm::vec4>{
                glm::vec4{ 7.5, 5, 0, 0 }, glm::vec4{ 7.5, 5, 0, 0 }, glm::vec4{ 7.5, 5, 0, 0 }, glm::vec4{ 7.5, 5, 0, 0 }, glm::vec4{ 7.5, 5, 0, 0 }, },
            std::vector<glm::vec4>{
                glm::vec4{ 7.5, -0.625, 0, 0 }, glm::vec4{ 6.63402, 0.0603378, 0.823148, 0 }, glm::vec4{ 6.81026, 0.430112, -1.61947, 0 }, glm::vec4{ 7.58278, 0.321184, -2.77189, 0 }, glm::vec4{ 7.5, -0.625, 0, 0 }, },
            std::vector<glm::vec4>{
                glm::vec4{ 7.5, -2.5, 0, 0 }, glm::vec4{ 3.09428, -1.03143, 0, 0 }, glm::vec4{ 3.31034, -1.10345, 0, 0 }, glm::vec4{ 5.08609, -1.69536, 0, 0 }, glm::vec4{ 7.5, -2.5, 0, 0 }, },
            std::vector<glm::vec4>{
                glm::vec4{ 7.5, -0.624999, 0, 0 }, glm::vec4{ 6.63402, -1.98916, -0.823148, 0 }, glm::vec4{ 6.81026, -2.37489, 1.61947, 0 }, glm::vec4{ 7.58278, -2.13936, 2.77189, 0 }, glm::vec4{ 7.5, -0.624999, 0, 0 }, },
            std::vector<glm::vec4>{
                glm::vec4{ 7.5, 5, 0, 0 }, glm::vec4{ 7.5, 5, 0, 0 }, glm::vec4{ 7.5, 5, 0, 0 }, glm::vec4{ 7.5, 5, 0, 0 }, glm::vec4{ 7.5, 5, 0, 0 }, },
        },
        std::vector<std::vector<glm::vec4> >{ // derivatives_v
            std::vector<glm::vec4>{
                glm::vec4{ 0, 3.2476, 7.5, 0 }, glm::vec4{ 0, -0.405949, 7.5, 0 }, glm::vec4{ 0, -1.6238, 7.5, 0 }, glm::vec4{ 0, -0.405949, 7.5, 0 }, glm::vec4{ 0, 3.2476, 7.5, 0 }, },
            std::vector<glm::vec4>{
                glm::vec4{ 10.6787, 6.74327, 21.7383, 0 }, glm::vec4{ 0.617361, -0.123568, 4.37036, 0 }, glm::vec4{ -0.910951, -1.24673, 4.98836, 0 }, glm::vec4{ -2.07892, -1.33852, 9.38451, 0 }, glm::vec4{ 0, 3.2476, 7.5, 0 }, },
            std::vector<glm::vec4>{
                glm::vec4{ 0, 14.2082, 32.8125, 0 }, glm::vec4{ 0, 0.0839437, 3.29624, 0 }, glm::vec4{ 0, -1.25054, 3.92675, 0 }, glm::vec4{ 0, -2.2935, 9.63326, 0 }, glm::vec4{ 0, 3.2476, 7.5, 0 }, },
            std::vector<glm::vec4>{
                glm::vec4{ -10.6787, 12.0826, 21.7383, 0 }, glm::vec4{ -0.617361, 0.185112, 4.37036, 0 }, glm::vec4{ 0.910951, -1.7022, 4.98836, 0 }, glm::vec4{ 2.07892, -2.37798, 9.38451, 0 }, glm::vec4{ 0, 3.2476, 7.5, 0 }, },
            std::vector<glm::vec4>{
                glm::vec4{ 0, 3.2476, 7.5, 0 }, glm::vec4{ 0, -0.405949, 7.5, 0 }, glm::vec4{ 0, -1.6238, 7.5, 0 }, glm::vec4{ 0, -0.405949, 7.5, 0 }, glm::vec4{ 0, 3.2476, 7.5, 0 }, },
        }
    }
},
