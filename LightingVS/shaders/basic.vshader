#version 130

uniform mat4 uProjMatrix;
uniform mat4 uModelViewMatrix;
uniform mat4 uNormalMatrix;

/*
  Cvec3f p, n;  //point coordinates, normal coordinates
  Cvec2f t; //texture coordinates
  Cvec3f tan, bin; //tangent vector, binormal vector

  VertexPNT(float x, float y, float z,
		   float nx, float ny, float nz,
		   float tu, float tv,
		   float tx, float ty, float tz,
		   float bx, float by, float bz)
	: p(x,y,z), n(nx, ny, nz), t(tu,tv), tan(tx,ty,tz),bin(bx,by,bz)
*/


in vec3 aPosition;
in vec3 aNormal;
in vec2 aTexCoord;
in vec3 aTangent;
in vec3 aBinormal;


out vec3 vNormal;
out vec3 vPosition;
out vec2 vTexCoord;
out vec3 vTangent;
out vec3 vBinormal;
out mat4 vNTMat;   

void main() {
  vNormal = vec3(uNormalMatrix * vec4(aNormal, 0.0));
  vTangent = vec3(uNormalMatrix * vec4(aTangent, 0.0));
  vBinormal = vec3(uNormalMatrix * vec4(aBinormal, 0.0));

  mat4 T = mat4(mat3(aTangent, aBinormal, aNormal));

  // send position (eye coordinates) to fragment shader
  vec4 tPosition = uModelViewMatrix * vec4(aPosition, 1.0);

  vPosition = vec3(tPosition);
  vTexCoord = aTexCoord;
  gl_Position = uProjMatrix * tPosition;
  vNTMat = (uNormalMatrix)*T;
}