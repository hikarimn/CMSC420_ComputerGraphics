#version 130

uniform vec3 uLight, uLight2, uColor;
uniform int uUseTexture;
uniform sampler2D uTexUnit;
uniform sampler2D uTexUnit1;
uniform int uSphere;

in vec3 vNormal;
in vec3 vPosition;
in vec2 vTexCoord;
in vec3 vTangent; //new!
in vec3 vBinormal; //new!
in mat4 vNTMat;

out vec4 fragColor;

void main() {
	vec3 tolight = normalize(uLight - vPosition);
	vec3 tolight2 = normalize(uLight2 - vPosition);
	vec3 normal = normalize(vNormal);

	if (uUseTexture == 2){
		vec3 nvec = texture2D(uTexUnit1,vTexCoord).xyz;
		nvec = nvec * 2 - vec3(1,1,1);
		normal = normalize(vNTMat * mat4(mat3(vTangent, vBinormal, vNormal)) * vec4(nvec,0)).xyz;
	}
	

	float diffuse = max(0.0, dot(normal, tolight));
	diffuse += max(0.0, dot(normal, tolight2));
	  
	if(uUseTexture == 0) {
	  vec3 intensity = uColor * diffuse;
	  if (uSphere == 1) intensity = uColor;

	  fragColor = vec4(intensity, 1.0);
  } else {
	fragColor = texture(uTexUnit, vTexCoord)*diffuse;
  }
}
