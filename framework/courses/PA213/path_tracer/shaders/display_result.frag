#version 450 core

// ----------------------------------------------------------------------------
// Input Variables
// ----------------------------------------------------------------------------
in VertexData
{
	vec2 tex_coord;  // The vertex texture coordinates.
} in_data;

// The texture with accumulated color.
layout (binding = 0) uniform sampler2D accumulated_color;
// The number of iterations accumulated in the texture.
uniform int iterations;
// The flag determinig wether we apply gamma correction to the color.
uniform bool gamma_correction;

// ----------------------------------------------------------------------------
// Output Variables
// ----------------------------------------------------------------------------
// The final output color.
layout (location = 0) out vec4 final_color;

// ----------------------------------------------------------------------------
// Main Method
// ----------------------------------------------------------------------------
void main()
{
	// We applie the transformation to the color from the texture using the uniform matrix.
	vec3 accumulated = texture(accumulated_color, in_data.tex_coord).rgb / iterations;
	// We apply gamma correction to the color.
	if(gamma_correction){
		accumulated = pow(accumulated, vec3(1.0/5.2));
	}
	final_color = vec4(accumulated,1);
}