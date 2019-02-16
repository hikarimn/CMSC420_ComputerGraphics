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

  RigTForm(const Cvec3& t, const Quat& r): t_(t),r_(r) {
  }

  explicit RigTForm(const Cvec3& t) : t_(t) {
  }

  explicit RigTForm(const Quat& r) : r_(r) {
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
	  return Matrix4::makeTranslation(t_)*quatToMatrix(r_)*a;
  }

  RigTForm operator * (const RigTForm& a) const {
	  return RigTForm(t_ + Cvec3(r_*Cvec4(a.getTranslation())),r_*a.getRotation());
  }
};

inline RigTForm inv(const RigTForm& tform) {
	Quat i = inv(tform.getRotation());
	return RigTForm(-Cvec3(i*Cvec4(tform.getTranslation())), i);
}

inline RigTForm transFact(const RigTForm& tform) {
  return RigTForm(tform.getTranslation());
}

inline RigTForm linFact(const RigTForm& tform) {
  return RigTForm(tform.getRotation());
}

inline Matrix4 rigTFormToMatrix(const RigTForm& tform) {
	return Matrix4::makeTranslation(tform.getTranslation())*quatToMatrix(tform.getRotation());
}

#endif
