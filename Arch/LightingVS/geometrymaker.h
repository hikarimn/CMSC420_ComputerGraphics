#ifndef GEOMETRYMAKER_H
#define GEOMETRYMAKER_H

#include <cmath>
#include <iostream>
#include <stdio.h>      /* printf */
#include "cvec.h"

//--------------------------------------------------------------------------------
// Helpers for creating some special geometries such as plane, cubes, and spheres
//--------------------------------------------------------------------------------


// A generic vertex structure containing position, normal, and texture information
// Used by make* functions to pass vertex information to the caller
struct GenericVertex {
  Cvec3f pos;
  Cvec3f normal; //normal vector
  Cvec2f tex;
  Cvec3f tangent, binormal;  // tangent vector, binormal vector

  GenericVertex(
    float x, float y, float z,
    float nx, float ny, float nz,
    float tu, float tv,
    float tx, float ty, float tz,
    float bx, float by, float bz)
    : pos(x, y ,z), normal(nx, ny, nz), tex(tu, tv), tangent(tx, ty, tz), binormal(bx, by, bz)
  {}
};


struct SmallVertex {
	Cvec3f pos;
	Cvec3f normal;

	SmallVertex(
		const Cvec3f& p, const Cvec3f& n)
		: pos(p), normal(n)
	{}
};

inline void getPlaneVbIbLen(int& vbLen, int& ibLen) {
  vbLen = 4;
  ibLen = 6;
}

template<typename VtxOutIter, typename IdxOutIter>
void makePlane(float size, VtxOutIter vtxIter, IdxOutIter idxIter) {
  float h = size / 2.0;
  *vtxIter = GenericVertex(    -h, 0, -h, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, -1);
  *(++vtxIter) = GenericVertex(-h, 0,  h, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, -1);
  *(++vtxIter) = GenericVertex( h, 0,  h, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, -1);
  *(++vtxIter) = GenericVertex( h, 0, -h, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, -1);
  *idxIter = 0;
  *(++idxIter) = 1;
  *(++idxIter) = 2;
  *(++idxIter) = 0;
  *(++idxIter) = 2;
  *(++idxIter) = 3;
}
inline void getTriangleVbIbLen(int& vbLen, int& ibLen) {
	vbLen = 3;
	ibLen = 3;
}

/*
GenericVertex(
    float x, float y, float z,
    float nx, float ny, float nz,
    float tu, float tv,
    float tx, float ty, float tz,
    float bx, float by, float bz)
    : pos(x,y,z), normal(nx,ny,nz), tex(tu, tv), tangent(tx, ty, tz), binormal(bx, by, bz)
	*/
template<typename VtxOutIter, typename IdxOutIter>
void makeTriangle(float size, VtxOutIter vtxIter, IdxOutIter idxIter) {
	float h = size / 2.0;
	*vtxIter = GenericVertex(0, 0, h-(2*h*0.86602), 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, -1);
	*(++vtxIter) = GenericVertex(-h, 0, h, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, -1);
	*(++vtxIter) = GenericVertex(h, 0, h, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, -1);
	
	*idxIter = 0;
	*(++idxIter) = 1;
	*(++idxIter) = 2;
}

inline void getCubeVbIbLen(int& vbLen, int& ibLen) {
  vbLen = 24;
  ibLen = 36;
}

template<typename VtxOutIter, typename IdxOutIter>
void makeCube(float size, VtxOutIter vtxIter, IdxOutIter idxIter) {
  float h = size / 2.0;
#define DEFV(x, y, z, nx, ny, nz, tu, tv) { \
    *vtxIter = GenericVertex(x h, y h, z h, \
                             nx, ny, nz, tu, tv, \
                             tan[0], tan[1], tan[2], \
                             bin[0], bin[1], bin[2]); \
    ++vtxIter; \
}
  Cvec3f tan(0, 1, 0), bin(0, 0, 1);
  DEFV(+, -, -, 1, 0, 0, 0, 0); // facing +X
  DEFV(+, +, -, 1, 0, 0, 1, 0);
  DEFV(+, +, +, 1, 0, 0, 1, 1);
  DEFV(+, -, +, 1, 0, 0, 0, 1);

  tan = Cvec3f(0, 0, 1);
  bin = Cvec3f(0, 1, 0);
  DEFV(-, -, -, -1, 0, 0, 0, 0); // facing -X
  DEFV(-, -, +, -1, 0, 0, 1, 0);
  DEFV(-, +, +, -1, 0, 0, 1, 1);
  DEFV(-, +, -, -1, 0, 0, 0, 1);

  tan = Cvec3f(0, 0, 1);
  bin = Cvec3f(1, 0, 0);
  DEFV(-, +, -, 0, 1, 0, 0, 0); // facing +Y
  DEFV(-, +, +, 0, 1, 0, 1, 0);
  DEFV(+, +, +, 0, 1, 0, 1, 1);
  DEFV(+, +, -, 0, 1, 0, 0, 1);

  tan = Cvec3f(1, 0, 0);
  bin = Cvec3f(0, 0, 1);
  DEFV(-, -, -, 0, -1, 0, 0, 0); // facing -Y
  DEFV(+, -, -, 0, -1, 0, 1, 0);
  DEFV(+, -, +, 0, -1, 0, 1, 1);
  DEFV(-, -, +, 0, -1, 0, 0, 1);

  tan = Cvec3f(1, 0, 0);
  bin = Cvec3f(0, 1, 0);
  DEFV(-, -, +, 0, 0, 1, 0, 0); // facing +Z
  DEFV(+, -, +, 0, 0, 1, 1, 0);
  DEFV(+, +, +, 0, 0, 1, 1, 1);
  DEFV(-, +, +, 0, 0, 1, 0, 1);

  tan = Cvec3f(0, 1, 0);
  bin = Cvec3f(1, 0, 0);
  DEFV(-, -, -, 0, 0, -1, 0, 0); // facing -Z
  DEFV(-, +, -, 0, 0, -1, 1, 0);
  DEFV(+, +, -, 0, 0, -1, 1, 1);
  DEFV(+, -, -, 0, 0, -1, 0, 1);
#undef DEFV

  for (int v = 0; v < 24; v +=4) {
    *idxIter = v;
    *++idxIter = v + 1;
    *++idxIter = v + 2;
    *++idxIter = v;
    *++idxIter = v + 2;
    *++idxIter = v + 3;
    ++idxIter;
  }
}


inline void getSurfaceVbIbLen(int s_steps, bool wrap_s, int t_steps, bool wrap_t, int& vbLen, int& ibLen) {
	assert(s_steps > 1);
	vbLen = (s_steps + (wrap_s ? 0 : 1)) * (t_steps + (wrap_t ? 0 : 1));
	ibLen = (s_steps * 2 + 2)*t_steps;
}

template<typename vMaker, typename nMaker, typename VtxOutIter, typename IdxOutIter>
void makeSurface(float start_s, float step_s, unsigned short s_steps, bool wrap_s,
	float start_t, float step_t, unsigned short t_steps, bool wrap_t, vMaker makeV, nMaker makeN,
	VtxOutIter vtxIter, IdxOutIter idxIter) {
	/* start_s = 0, step_s = 9, s_steps = 20, start_t = 0, step_t = 0.5, t_steps = 2 */

	unsigned short t, s, next_t;

	unsigned short sCount = s_steps + (wrap_s ? 0 : 1); // 20 + 0 = 20
	unsigned short tCount = t_steps + (wrap_t ? 0 : 1); // 2 + 1 = 3

	// Make the vertices
	for (t = 0; t < tCount; t++) { //t: 0 - 3
		for (s = 0; s < sCount; s++) { //s: 0 - 20
			float paramS = start_s + s * step_s; // 0 + s * 9 -> paramS: 0 - 180
			float paramT = start_t + t * step_t; // 0 + t * 2 -> paramT: 0 -  6
			Cvec3f vertex = makeV(paramS, paramT);
			Cvec3f normal = makeN(paramS, paramT);
			*vtxIter = SmallVertex(vertex, normal);
			++vtxIter;
		}
	}

	// Construct a series of triangle strips, one per t step
	// We use an indexed representation for the triangle strips
	for (t = 0; t < t_steps - 1; t++) {
		next_t = t + 1;
		for (s = 0; s < s_steps; s++) {
			*idxIter = t * sCount + s;
			++idxIter;
			*idxIter = next_t * sCount + s;
			++idxIter;
		}
		if (wrap_s)
			s = 0;
		else
			s++;
		*idxIter = t * sCount + s;
		++idxIter;
		*idxIter = next_t * sCount + s;
		++idxIter;
	}
	if (wrap_t)
		next_t = 0;
	else
		next_t = t + 1;
	for (s = 0; s < s_steps; s++) {
		*idxIter = t * sCount + s;
		++idxIter;
		*idxIter = next_t * sCount + s;
		++idxIter;
	}
	if (wrap_s)
		s = 0;
	else
		s++;
	*idxIter = t * sCount + s;
	++idxIter;
	*idxIter = next_t * sCount + s;
	++idxIter;
}

/*tring to make an arch************************************/
inline void getArchVbIbLen(int steps, int& vbLen, int& ibLen) {
	vbLen = (steps + 1) * 2;
	ibLen = (steps + 1) * 2;
}

template<typename vMaker,  typename VtxOutIter, typename IdxOutIter>
void makeArch(float a1, float a2,  unsigned short steps, Cvec3f point1, Cvec3f point2, vMaker makeV,
	VtxOutIter vtxIter, IdxOutIter idxIter) {

	VtxOutIter a, b;

	unsigned short t, next;
	float step1 = ( -point1[0])*2/steps;
	float step2 = ( -point2[0]) * 2 /steps;
	// Make the vertices
	for (t = 0; t <= steps; t++) {
		float paramT1 = point1[0] + t * step1;
		float paramT2 = point2[0] + t * step2;
		Cvec3f vertex1 = makeV(paramT1, a1, point1);
		Cvec3f vertex2 = makeV(paramT2, a2, point2);
		Cvec3f normal;
		if (t != 0) {
			Cvec3f edge1 = b->p - a->p;
			Cvec3f edge2 = a->p - vertex1;
			normal = cross(edge1, edge2).normalize();
			a->n = normal;
			b->n = normal;
		}
		//Cvec3f normal = cross(b-a, vertex1 - a).normalize();
		*vtxIter = SmallVertex(vertex1, normal);
		a = vtxIter++;
		*vtxIter = SmallVertex(vertex2, normal);
		b = vtxIter++;
	}

	// Construct a series of triangle strips, one per t step
	// We use an indexed representation for the triangle strips
	for (t = 0; t < steps*2+2; t++) {
		*idxIter = t;
		++idxIter;
	}

}
/*************************************************************/

//http://www.songho.ca/opengl/gl_sphere.html
inline void getSphereVbIbLen(int slices, int stacks, int& vbLen, int& ibLen) {
  assert(slices > 1);
  assert(stacks >= 2);
  vbLen = (slices + 1) * (stacks + 1);
  ibLen = slices * stacks * 6;
}

template<typename VtxOutIter, typename IdxOutIter>
void makeSphere(float radius, int slices, int stacks, VtxOutIter vtxIter, IdxOutIter idxIter) {
  using namespace std;
  assert(slices > 1);
  assert(stacks >= 2);

  const double radPerSlice = 2 * CS175_PI / slices;
  const double radPerStack = CS175_PI / stacks;

  vector<double> longSin(slices+1), longCos(slices+1);
  vector<double> latSin(stacks+1), latCos(stacks+1);
  for (int i = 0; i < slices + 1; ++i) {
    longSin[i] = sin(radPerSlice * i);
    longCos[i] = cos(radPerSlice * i);
  }
  for (int i = 0; i < stacks + 1; ++i) {
    latSin[i] = sin(radPerStack * i);
    latCos[i] = cos(radPerStack * i);
  }

  for (int i = 0; i < slices + 1; ++i) {
    for (int j = 0; j < stacks + 1; ++j) {
      float x = longCos[i] * latSin[j];
      float y = longSin[i] * latSin[j];
      float z = latCos[j];

      Cvec3f n(x, y, z);
      Cvec3f t(-longSin[i], longCos[i], 0);
      Cvec3f b = cross(n, t);

      *vtxIter = GenericVertex(
        x * radius, y * radius, z * radius,
        x, y, z,
        1.0/slices*i, 1.0/stacks*j,
        t[0], t[1], t[2],
        b[0], b[1], b[2]);
      ++vtxIter;

      if (i < slices && j < stacks ) {
        *idxIter = (stacks+1) * i + j;
        *++idxIter = (stacks+1) * i + j + 1;
        *++idxIter = (stacks+1) * (i + 1) + j + 1;

        *++idxIter = (stacks+1) * i + j;
        *++idxIter = (stacks+1) * (i + 1) + j + 1;
        *++idxIter = (stacks+1) * (i + 1) + j;
        ++idxIter;
      }
    }
  }
}


#endif
