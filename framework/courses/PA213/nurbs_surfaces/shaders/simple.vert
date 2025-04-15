#version 450

// ----------------------------------------------------------------------------
// Input Variables
// ----------------------------------------------------------------------------

layout(location = 0) uniform mat4 model;
layout(location = 1) uniform mat4 view;
layout(location = 2) uniform mat4 projection;

layout(location = 0) in vec3 position;
layout(location = 5) in vec3 color;

// ----------------------------------------------------------------------------
// Output Variables
// ----------------------------------------------------------------------------

layout(location = 0) out vec3 fs_color;

// ----------------------------------------------------------------------------
// Main Method
// ----------------------------------------------------------------------------

void main()
{
	fs_color = color;

    gl_Position = projection * view * model * vec4(position, 1.0);
}
