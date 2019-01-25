#ifndef RIGTFORM_H
#define RIGTFORM_H

#include <iostream>
#include <cassert>

#include "matrix4.h"
#include "quat.h"

class RigTForm {
  Cvec3 t_; // translation component
  Quat r_;  // rotation component represented as a quaternion

public:
  RigTForm() : t_(0) {
    assert(norm2(Quat(1,0,0,0) - r_) < CS175_EPS2);
  }

  RigTForm(const Cvec3& t, const Quat& r) {
	  t_ = t;
	  r_ = r;
  }

  explicit RigTForm(const Cvec3& t) {
	  t_ = t;
	  r_ = Quat(1, 0, 0, 0);
  }

  explicit RigTForm(const Quat& r) {
	  t_ = Cvec3(0,0,0);
	  r_ = r;
  }

  Cvec3 getTranslation() const {
    return t_;
  }

  Quat getRotation() const {
    return r_;
  }

  RigTForm& setTranslation(const Cvec3& t) {
    t_ = t;
    return *this;
  }

  RigTForm& setRotation(const Quat& r) {
    r_ = r;
    return *this;
  }

  Cvec4 operator * (const Cvec4& a) const {
	  Matrix4 t;
	  t(0, 3) = this->getTranslation()[0];
	  t(1, 3) = this->getTranslation()[1];
	  t(2, 3) = this->getTranslation()[2];

	  return (t*quatToMatrix(this->getRotation()))*a;
  }

  RigTForm operator * (const RigTForm& a) const {
	  Quat rot1 = this->getRotation();
	  Quat rot2 = a.getRotation();
	  Cvec3 tra1 = this->getTranslation();
	  Cvec3 tra2 = a.getTranslation();
	  Cvec4 vec4 = Cvec4(tra1, 1) + (rot1 * Cvec4(tra2, 1));
	  Cvec3 vec3 = Cvec3(vec4[0], vec4[1], vec4[2]);
	  RigTForm result;
	  result.setTranslation(vec3);
	  result.setRotation(rot1 * rot2);
	  return result;
  }
};

inline RigTForm inv(const RigTForm& tform) {
	Cvec4 vec4 = inv(tform.getRotation()) * Cvec4(tform.getTranslation(),1);
	Cvec3 vec3 = Cvec3(-vec4[0], -vec4[1], -vec4[2]);
	RigTForm result;
	result.setTranslation(vec3);
	result.setRotation(inv(tform.getRotation()));
	return result;
}

inline RigTForm transFact(const RigTForm& tform) {
  return RigTForm(tform.getTranslation());
}

inline RigTForm linFact(const RigTForm& tform) {
  return RigTForm(tform.getRotation());
}

inline Matrix4 rigTFormToMatrix(const RigTForm& tform) {
	Matrix4 mat4;
	Cvec3 vec3 = tform.getTranslation();
	mat4(0, 3) = vec3[0];
	mat4(1, 3) = vec3[1];
	mat4(2, 3) = vec3[2];
	mat4 = mat4 * quatToMatrix(tform.getRotation());
	return mat4;
}

#endif
