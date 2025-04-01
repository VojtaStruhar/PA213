#version 450 core

// ----------------------------------------------------------------------------
// Global Variables
// ----------------------------------------------------------------------------

// Due to floating-point precision errors, when a ray intersects geometry at a surface, the point of intersection could possibly be just below the surface.
// The subsequent reflection and shadow rays would then bounce off the *inside* wall of the surface. This is known as self-intersection.
// We, therefore, use a small EPSILON to offset the subsequent rays.
const float EPSILON = 1e-3;
// The size of the box in the scene.
float BOX_SIZE=15;
// The PI value.
float PI = 3.14159265359;
// The number of lights in the scene.
int LIGHT_COUNT = 0;

// ----------------------------------------------------------------------------
// Input Variables
// ----------------------------------------------------------------------------
// The interpolated data from the vertex shader.
in VertexData
{
	vec2 tex_coord;  // The texture coordinates on the screen (we are using full-screen quad).
} in_data;

// The UBO with camera data.	
layout (std140, binding = 0) uniform CameraBuffer
{
	mat4 projection;	  // The projection matrix.
	mat4 projection_inv;  // The inverse of the projection matrix.
	mat4 view;			  // The view matrix
	mat4 view_inv;		  // The inverse of the view matrix.
	mat3 view_it;		  // The inverse of the transpose of the top-left part 3x3 of the view matrix
	vec3 eye_position;	  // The position of the eye in world space.
};

// The structure defining a material.
struct PathTracerMaterialData
{
    vec4 albedo;      // The albedo color of the material (w<1 indicates transparent/glass material).
    vec3 emission;    // The material emission color.
    float roughness;  // The roughness of the material.
};

 // The buffer with individual materials (material[0] contains material for the box).
layout (std140, binding = 3) uniform PBRMaterialBuffer
{
   PathTracerMaterialData materials[5];
};

// The buffer with spheres.
layout (std140, binding = 4) uniform SpheresBuffer
{
	vec4 spheres[4];
};

// The windows size.
uniform vec2 resolution;
// The current sampling strategy.
uniform int sampling_strategy;
// The current brdf.
uniform int current_brdf;
// The number of spheres to render.
uniform int spheres_count;
// The number of bounces.
uniform int bounces;
// The elapsed time from the begining.
uniform float elapsed_time;
// The random value from the canonical uniform distribution.
uniform float canonical_random_value;
// The texture with already accumulated radiance.
layout (binding = 0) uniform sampler2D accumulated_radiance; 

// ----------------------------------------------------------------------------
// Output Variables
// ----------------------------------------------------------------------------
// The final output color.
layout (location = 0) out vec4 final_color;

// ----------------------------------------------------------------------------
// Path Tracing Structures
// ----------------------------------------------------------------------------
// The definition of a ray.
struct Ray {
    vec3 origin;     // The ray origin.
    vec3 direction;  // The ray direction.
};
// The definition of a ray.
struct BSDFSample {
    vec3 direction;  // The sample direction.
	float pdf;       // The probability density function associated with this sample.
	bool diffuse;	 // The flag indicating whether the sample is diffuse (or specular).
	bool transmission;  // The flag indicating whether the sample is transmission (or reflection).
};
// The definition of an intersection.
struct Hit {
    float t;				  // The distance between the ray origin and the intersection points along the ray. 
	vec3 intersection;        // The intersection point.
    vec3 normal;              // The surface normal at the interesection point.
	vec3 ffnormal;			  // The front face normal (useful for transmission).
	bool emissive;			  // The flag indicating whether the object is emissive.
	bool glass;				  // The flag indicating whether the object is glass.
	PathTracerMaterialData material; // The material of the object at the intersection point.
};
const Hit miss = Hit(1e20, vec3(0.0), vec3(0.0), vec3(0.0), false, false, PathTracerMaterialData(vec4(0),vec3(0),0));

// Checks wether the material is emissive.
bool IsEmissive(PathTracerMaterialData material){
	return material.emission != vec3(0);
}

// Checks wether the material is transparent.
bool IsGlass(PathTracerMaterialData material){
	return material.albedo.w < 1.0;
}

// ----------------------------------------------------------------------------
// Intersection Functions
// ----------------------------------------------------------------------------

// Computes an intersection between a ray and a sphere defined by its center and radius.
Hit RaySphereIntersection(Ray ray, vec3 center, float radius, int i) {
	vec3 oc = ray.origin - center;
	float b = dot(ray.direction, oc);
	float c = dot(oc, oc) - (radius*radius);

	float det = b*b - c;
	if (det < 0.0) return miss;

	float t = -b - sqrt(det);
	if (t < 0.0) t = -b + sqrt(det);
	if (t < 0.0) return miss;

	vec3 intersection = ray.origin + t * ray.direction;
	vec3 normal = normalize(intersection - center);

	vec3 ffnormal = dot(normal, ray.direction) <= 0.0 ? normal : normal * -1.0;

    return Hit(t, intersection, normal, ffnormal, IsEmissive(materials[i]), IsGlass(materials[i]), materials[i]);
}

// Computes an intersection between a ray and a plane defined by its normal and one point inside the plane.
Hit RayPlaneIntersection(Ray ray, vec3 normal, vec3 point) {
	float nd = dot(normal, ray.direction);
	if(nd > 0) return miss;
 	vec3 sp = point - ray.origin;
	float t = dot(sp, normal) / nd;
    if (t < 0.0) return miss;

	vec3 intersection = ray.origin + t * ray.direction;
	
	
	float size_with_epsilon = BOX_SIZE+EPSILON;
	float yoffset = 2.0 + EPSILON;
	if(intersection.x > size_with_epsilon || intersection.x < -size_with_epsilon || intersection.y > size_with_epsilon || intersection.y < -yoffset || intersection.z > size_with_epsilon || intersection.z < -size_with_epsilon) return miss;

	vec3 ffnormal = dot(normal, ray.direction) <= 0.0 ? normal : normal * -1.0;
    return Hit(t, intersection, normal, ffnormal, IsEmissive(materials[0]), IsGlass(materials[0]), materials[0]);
}

// Evaluates the intersections of the ray with the scene objects and returns the closes hit.
Hit ClosestHit(Ray ray){
	// Sets the closes hit either to miss or to an intersection with one of the plane representing the box.
	Hit planes[6] = {miss, miss, miss, miss, miss, miss};
	planes[0] = RayPlaneIntersection(ray, vec3(0, 1, 0), vec3(0,/*-BOX_SIZE*/-2,0)); // bottom
	planes[1] = RayPlaneIntersection(ray, vec3(-1, 0, 0), vec3(BOX_SIZE,0,0)); // right
	planes[2] = RayPlaneIntersection(ray, vec3(1, 0, 0), vec3(-BOX_SIZE,0,0)); // left
	planes[3] = RayPlaneIntersection(ray, vec3(0, -1, 0), vec3(0,BOX_SIZE,0)); // top
	planes[4] = RayPlaneIntersection(ray, vec3(0, 0, 1), vec3(0,0,-BOX_SIZE)); // back
	planes[5] = RayPlaneIntersection(ray, vec3(0, 0, -1), vec3(0,0,BOX_SIZE)); // front
	Hit closest_hit = miss;
	for(int i = 0; i < 6; i++){
		if(planes[i].t < closest_hit.t){
			closest_hit = planes[i];
		}
	}
	
	// Iterates over all spheres.  
	for(int i = 0; i < spheres_count; i++){
		vec3 center = spheres[i].xyz;
		Hit intersection = RaySphereIntersection(ray, center, spheres[i].w, i+1);
		if(intersection.t < closest_hit.t){
			closest_hit = intersection;
		}
	}

	// Returns the closest hit.
    return closest_hit;
}

// ----------------------------------------------------------------------------
// Random
// ----------------------------------------------------------------------------

// The seed we change every time we call Rand method.
float random_seed=0;
// The random number generator based on the vector from a uniform distribution.
float Rand()
{
	random_seed += canonical_random_value;
	return fract(sin(dot(vec2(random_seed*0.1) * gl_FragCoord.xy, vec2(12.9898, 78.233))) * 43758.5453);
}

vec2 Rand2() {
	return vec2(Rand(), Rand());
}

vec3 Rand3() {
	return vec3(Rand(), Rand(), Rand());
}

/// Slide #56
vec3 LocalToWorldCoords(vec3 local, vec3 n) {
	local = normalize(local);
	vec3 U = abs(n.z) < 0.999 ? vec3(0, 0, 1) : vec3(0, 1, 0);
	vec3 T = cross(U, n);
	vec3 B = cross(n, T);

	vec3 result = T * local.x + B * local.y + n * local.z;
	return normalize(result);
}

// ----------------------------------------------------------------------------
// Assignment: Uniform Lambertian BRDF
// ----------------------------------------------------------------------------

BSDFSample SampleUniformLambert(Hit hit, Ray ray) {
	// TIP: Implement a function for the world-to-local transition. It will be used several times.

	// TODO: 1a: Implement an algorithm to sample a unit hemisphere using uniform distribution (see slide #54).
	//          You will also need to tranform the sample from the local to the world coordinates (see slide #56).
	//   Hints: The surface normal is stored in 'hit.normal'. 
	//         	Use Rand() to generate uniformly distributed values between [0,1].
	//          Do not forget to normalize all the vectors.
	vec3 dir = vec3(1.0);

	vec2 xi = Rand2();
	float theta = 2 * PI * xi.x;
	float r = sqrt(xi.y);

	// Slide #54
	float x = sqrt(1 - xi.x * xi.x) * cos(2 * PI * xi.y);
	float y = sqrt(1 - xi.x * xi.x) * sin(2 * PI * xi.y);
	float z = xi.x;

	dir = vec3(x, y, z);
	dir = LocalToWorldCoords(dir, hit.normal);


	// TODO: 1b: Set a correct PDF value for uniform samples (see slide #54).
	float pdf = 1.0 / (2.0 * PI);

	// DO NOT MODIFY
    return BSDFSample(dir, pdf, true, false);
}

vec3 EvaluateLambert(Hit hit, BSDFSample Wi, vec3 Wr) {

	// TODO: 2: Return a correct value for Lambertian BRDF function (for now use 1.0 for albedo).
	//         See slide #67.
	return vec3(1.0 / PI);
}

// ----------------------------------------------------------------------------
// Assignment: Cosine-Weigted Lambertian BRDF
// ----------------------------------------------------------------------------

BSDFSample SampleCosineWeightedLambert(Hit hit, Ray ray){
	// TODO: 3a: Implement an algorithm to sample a unit hemisphere using cosine-weighted distribution (see slide #63).
	//          You will also need to tranform the sample from the local to the world coordinates (see slide #56).
	//   Hints: The surface normal is stored in 'hit.normal'. 
	//         	Use Rand() to generate uniformly distributed values between [0,1]. 
	vec3 dir = vec3(1.0);

	vec2 xi = Rand2();
	float theta = 2 * PI * xi.x;
	float r = sqrt(xi.y);

	// Slide #63
	float x = sqrt(xi.x) * cos(2 * PI * xi.y);
	float y = sqrt(xi.x) * sin(2 * PI * xi.y);
	float z = sqrt(1 - xi.x);

	dir = vec3(x, y, z);
	dir = LocalToWorldCoords(dir, hit.normal);

	// TODO: 3b: Set a correct PDF value for cosine-weighted samples (see slide #63).
	float pdf = cos(theta) / PI;

	// DO NOT MODIFY
    return BSDFSample(dir, pdf, true, false);
}

vec3 EvaluateLambertWithAlbedo(Hit hit, BSDFSample Wi, vec3 Wr) {
    // TODO: 4: Return a correct value for Lambertian BRDF function (using actuall albedo).
	//  Hints: The surface albedo is stored in 'hit.material.albedo'. 
	//         See slide #67.
	return hit.material.albedo.rgb / PI;
}

// ----------------------------------------------------------------------------
// Assignment: Sampling Smith GGX
// ----------------------------------------------------------------------------

BSDFSample SampleSmithGGX(Hit hit, Ray ray) {
	// WARNING: Probably the most complex task (?)

	// TODO: 5a: Implement an algorithm for sampling a unit hemisphere according to SmithGGX
	//          distribution (see slides #89 and #90).
	//          You will also need to tranform the sample from the local to the world coordinates (see slide #56).
	//   Hints: The 'alpha' term is stored in 'hit.material.roughness'.
	//          The surface normal is stored in 'hit.normal'.
	//         	Use Rand() to generate uniformly distributed values between [0,1].
	//          You can reflect the vector Wh along the surface normal using reflect() function.
	//          You can compute arccos using acos() function.
	//          The ray.direction stores the '-omega_r' vector;
	//          Avoid dividing by 0.
	vec3 dir = vec3(1.0);

	vec2 xi = Rand2();
	float a = hit.material.roughness;
	float theta = acos(sqrt((1 - xi.x) / (xi.x * (a * a - 1) + 1)));
	float r = sqrt(xi.y);

	// Slide #63 - same as Cosine-Weighted Sampling
	float x = sqrt(xi.x) * cos(2 * PI * xi.y);
	float y = sqrt(xi.x) * sin(2 * PI * xi.y);
	float z = sqrt(1 - xi.x);

	dir = vec3(x, y, z);
	dir = LocalToWorldCoords(dir, hit.normal);

	// TODO: 5b: Set a correct PDF value for samples based on SmithGGX distribution (see slides #89 and #90).
	float a2 = (a * a);
	float term1 = (a2 - 1) * cos(theta) * cos(theta) + 1;
	float pdf = (a2 * cos(theta)) / (PI * term1 * term1);

	// DO NOT MODIFY
    return BSDFSample(dir, pdf, true, false);
}

float FresnelSchlick(float R0, vec3 Wi, vec3 N){
	// TODO: 6: Implement Schlick approximation of the Fresnel equations (see slide #79).
	//  Hints: Do not forget to clamp the dot product value to range [0,1]; 
	//         useful functions: pow, clamp
	float result = R0 + (1 - R0) * pow(1 - clamp(dot(Wi, N), 0, 1), 5.0);
	return result;
}

/// AKA Geometric Attenuation for SmithGGX
float SmithGGXMaskingShadowing(Hit hit, vec3 Wi, vec3 Wr){
	// TODO: 7: Implement SmithGGX masking and shadowing term (see slide #83)
	//  Hints: The 'alpha' term is stored in 'hit.material.roughness'
	//         The surface normal is stored in 'hit.normal'
	//         The angles theta_i and theta_r are angles between the vectors Wi and Wr and normal respectively.
	//         Avoid dividing by 0.
	float theta_i = dot(Wi, hit.normal);
	float theta_r = dot(Wr, hit.normal);
	float alpha = hit.material.roughness;

	float top = 2 * cos(theta_i) * cos(theta_r);
	float bottom_1 = cos(theta_r) * sqrt(alpha * alpha  + (1 - alpha * alpha) * cos(theta_i) * cos(theta_i));
	float bottom_2 = cos(theta_i) * sqrt(alpha * alpha  + (1 - alpha * alpha) * cos(theta_r) * cos(theta_r));


	return top / (bottom_1 + bottom_2);
}

float nonzero(float value) {
//	if (abs(value) < 0.001) {
//		return value < 0 ? -0.001 : 0.001;
//	}
	return value;
}

vec3 EvaluateSmithGGX_BRDF(Hit hit, BSDFSample Wi, vec3 Wr){
	// TODO: 8: Implement Microfacet BRDF (see slide #81)
	//  Hints: The 'alpha' term is stored in 'hit.material.roughness'.
	//         The surface normal is stored in 'hit.normal'.
	//         To get the Wi vector use 'Wi.direction'.
	//         When computing FresnelSchlick use R0 = 1.0;
	//         You want to multiply the final value by the albedo which is store in 'hit.material.albedo'.
	//         Avoid dividing by 0.

	float alpha = hit.material.roughness;
	vec3 Wh = (Wr + Wi.direction) / abs(Wr + Wi.direction);

	vec3 N = hit.normal;
	float theta_i = dot(Wi.direction, N);
	float theta_r = dot(Wr, N);
	float theta_h = dot(Wh, N);

	float F = FresnelSchlick(1.0, Wi.direction, Wh);
	float G = SmithGGXMaskingShadowing(hit, Wi.direction, Wr);

	float NDF_top = alpha * alpha;
	float NDF_bottom = ((NDF_top - 1) * cos(theta_h) * cos(theta_h) + 1);
	float D = NDF_top / nonzero(PI * NDF_bottom * NDF_bottom);

	float result_bottom = 4 * cos(theta_i) * cos(theta_r);

	return hit.material.albedo.rgb * (F * G * D / nonzero(result_bottom));
}

// ----------------------------------------------------------------------------
// Assignment: Combined Diffuse-Specular Sampling
// ----------------------------------------------------------------------------

BSDFSample SampleCombinedDiffuseSpecular(Hit hit, Ray ray){
	// TODO 9a: Use roughness stored in 'hit.material.roughness' and Rand() function to
	//          randomly perform either sampling based on the specular (SampleSmithGGX) 
	//          or diffuse (SampleCosineWeightedLambert) properties (see slide #95). 
	//   Hints: You can call the methods you have implemented earlier directly.
	//        
	vec3 dir = vec3(1.0);

	// TODO: 9b: The variable diffuse_sample should be set to true for diffuse samples,
	//          and to false for specular sample.
	bool diffuse_sample = true;

	// TODO: 9c: Set a correct PDF value weigted based on the surface roughness (see slide #95).
	float pdf = 1.0;
	
	// DO NOT MODIFY
	return BSDFSample(dir, pdf, diffuse_sample, false);
}

vec3 EvaluateCombinedDiffuseSpecular(Hit hit, BSDFSample Wi, vec3 Wr){
	// TODO: 10: Return Lambert BRDF (with albedo) for diffuse samples and Microfacet BRDF for specular samples.
	//   Hints: You can use the Wi.diffuse variable (you have set in TASK 9b, by setting diffuse_sample).
	return vec3(1.0);
}

// ----------------------------------------------------------------------------
// Path Tracing Evaluation
// ----------------------------------------------------------------------------

BSDFSample GetNextSample(Hit hit, Ray ray){
	switch(sampling_strategy){
		case 0:
			return SampleUniformLambert(hit, ray);
		case 1:	
			return SampleCosineWeightedLambert(hit, ray);
		case 2:
			return SampleSmithGGX(hit, ray);
		case 3:
			return SampleCombinedDiffuseSpecular(hit, ray);
	}
	return BSDFSample(vec3(0), 0, false, false);
}

vec3 ComputeBRDF (Hit hit, BSDFSample Wi, vec3 Wr){
	switch(current_brdf){
		case 0:
			return EvaluateLambert(hit, Wi, Wr);
		case 1:	
			return EvaluateLambertWithAlbedo(hit, Wi, Wr);
		case 2: 
			return EvaluateSmithGGX_BRDF(hit, Wi, Wr);
		case 3:
			return EvaluateCombinedDiffuseSpecular(hit, Wi, Wr);
	}
	return vec3(0);
}

// Traces the ray trough the scene and accumulates the radiance.
vec3 Trace(Ray ray) {
    // The accumulated radiance and attenuation.
	vec3 radiance = vec3(0.0);
    vec3 throughput = vec3(1.0);

    for (int i = 0; i < bounces; ++i) {
        Hit hit = ClosestHit(ray);
        if(hit == miss) return radiance;

		if(hit.emissive){
			radiance += hit.material.emission * throughput;
			return radiance;
		}

		BSDFSample bsdf_sample = GetNextSample(hit, ray);
		vec3 bsdf = ComputeBRDF(hit, bsdf_sample, -ray.direction);

		float cosTheta = max(dot(hit.normal, bsdf_sample.direction),0.0);
		if(bsdf_sample.pdf > EPSILON) {
			throughput *= bsdf * cosTheta / bsdf_sample.pdf;
			ray = Ray(hit.intersection +  (hit.glass ? -EPSILON * hit.ffnormal : EPSILON * hit.ffnormal), bsdf_sample.direction);
		} else {
			return radiance;
		}
	}

    return radiance;
}

// ----------------------------------------------------------------------------
// Main Method
// ----------------------------------------------------------------------------
void main()
{
	// The aspect ratio tells us how many times the window is wider (or narrower) w.r.t. its height.
	float aspect_ratio = resolution.x/resolution.y;

	// Computes the nubmer of lights in the scene
	for(int i = 0; i < spheres_count; i++){
		if(IsEmissive(materials[i+1])){
			LIGHT_COUNT++;
		}
	}

	// We use the texture coordinates and the aspect ratio to map the coordinates to the screen space.
	vec2 uv = (2.0*in_data.tex_coord - 1.0) * vec2(aspect_ratio, 1.0);

	// Computes the ray origin and ray direction using the view matrix.
	vec3 P = vec3(view_inv * vec4(uv, -1.0, 1.0));
	vec3 direction = normalize(P - eye_position);
	Ray ray = Ray(eye_position, direction);

	// We pass the ray to the trace function.
	// TODO: Clamp the returned radiance to avoid any numberical instabilities we may have introduced?
	vec3 radiance = Trace(ray); 

	// Adds the new radiance to the already accumulated radiance.
	final_color = vec4(radiance, 1.0) + texture(accumulated_radiance, in_data.tex_coord);
}