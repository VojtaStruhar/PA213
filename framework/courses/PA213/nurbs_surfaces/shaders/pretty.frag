#version 450

// ----------------------------------------------------------------------------
// Input Variables
// ----------------------------------------------------------------------------

layout(location = 0) in vec3 fs_position;
layout(location = 1) in vec3 fs_normal;
layout(location = 2) in vec2 fs_texture_coordinate;

layout(location = 3) uniform vec3 eye_position;
layout(location = 4) uniform bool use_texture;

layout(binding = 0) uniform sampler2D diffuse_texture;

// ----------------------------------------------------------------------------
// Output Variables
// ----------------------------------------------------------------------------

layout(location = 0) out vec4 final_color;

// ----------------------------------------------------------------------------
// Main Method
// ----------------------------------------------------------------------------

void main()
{
	vec3 light_vector = vec3(1,1,1);
	vec3 ambient = use_texture ? vec3(0) : vec3(0.1);
	vec3 diffuse = use_texture ? texture(diffuse_texture, fs_texture_coordinate).rgb : vec3(0.75,0.75,0.75);
	vec3 specular = vec3(0.5);
	float shininess = 100.0;

	vec3 L = normalize(light_vector);
	vec3 N = normalize(fs_normal); 
	vec3 E = normalize(eye_position - fs_position); 	
	vec3 H = normalize(L + E); 

	float NdotL = max(dot(N, L), use_texture ? 0.25 : 0.0);
	float NdotH = max(dot(N, H), 0.0);
	
	final_color = vec4(ambient + NdotL * diffuse + pow(NdotH, shininess) * specular, 1.0);
}
