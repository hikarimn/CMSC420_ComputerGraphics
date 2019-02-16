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
in mat4 vNTMat;  //new!

out vec4 fragColor;

void main() {
	vec3 tolight = normalize(uLight - vPosition);
	vec3 tolight2 = normalize(uLight2 - vPosition);
	vec3 normal = normalize(vNormal);
	vec3 halfVector = normalize(tolight - vPosition);
	vec3 halfVector2 = normalize(tolight2 - vPosition);
	
	//you will multiply , the normal data fetched from the texture
	//vec3 nvec = texture2D(uTexUnit,vTexCoord);  //new!
	//vec4 n = vec4(nvec, 0.0));  //new!
	//vec3 color = texture2D(uTexUnit,vTexCoord);  //new!  how do i use this?

	//vec4 newvec = n * 

	if (uUseTexture == 2){
		vec3 nvec = texture2D(uTexUnit1,vTexCoord).xyz;
		nvec = nvec * 2 - vec3(1,1,1);
		normal = normalize(vNTMat * mat4(mat3(vTangent, vBinormal, vNormal)) * vec4(nvec,0)).xyz;
	}

	float diffuse = max(0.0, dot(normal, tolight));
	//diffuse += max(0.0, dot(normal, tolight2));
	  

	if(uUseTexture == 0) {
	  vec3 intensity = uColor * diffuse;
	  if (uSphere == 1) intensity = uColor;
	  fragColor = vec4(intensity, 1.0);
  } else {
	fragColor = texture(uTexUnit, vTexCoord)*diffuse;
  }

  float specular = pow(max(0.0,dot(normal,halfVector)),12);
  //specular += pow(max(0.0,dot(normal,halfVector2)),12);
  fragColor += 0.3*vec4(1.0)*specular;
}