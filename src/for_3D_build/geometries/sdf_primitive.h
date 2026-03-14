/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for    *
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH   *
 * (smoothed particle hydrodynamics), a meshless computational method using  *
 * particle discretization.                                                  *
 *                                                                           *
 * SPHinXsys is partially funded by German Research Foundation               *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,            *
 *  HU1527/12-1 and HU1527/12-4.                                             *
 *                                                                           *
 * Portions copyright (c) 2017-2025 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file sdf_primitive.h
 * @brief All signed distance function (sdf) primitives use local coordinates.
 * All rotation is around the x-axis.
 * @details Here, we only give popular primitives,
 * For more signed distance function, please check the website:
 * https://iquilezles.org/articles/.
 * @author	Xiangyu Hu
 */

#ifndef SDF_PRIMITIVE_H
#define SDF_PRIMITIVE_H

#include "base_data_type.h"
#include "geometric_primitive.h"
#include "scalar_functions.h"

namespace SPH
{
//----------------------------------------------------------------------
// 3D geometric primitives for signed distance function definition
//----------------------------------------------------------------------
class SDFBall
{
    Real radius_;

  public:
    explicit SDFBall(Real radius) : radius_(radius) {}
    virtual ~SDFBall() {}
    void setParameters(Real radius) { radius_ = radius; }
    Real operator()(const Vec3d &point) const { return point.norm() - radius_; }
    BoundingBox3d findBounds() const { return BoundingBox3d(Vec3d::Constant(radius_)); }
};

class SDFBox
{
    Vec3d halfsize_;

  public:
    explicit SDFBox(const Vec3d &halfsize) : halfsize_(halfsize) {}
    void setParameters(const Vec3d &halfsize) { halfsize_ = halfsize; }
    Real operator()(const Vec3d &point) const;
    BoundingBox3d findBounds() const { return BoundingBox3d(halfsize_); }
};

class SDFCylinder
{
    Real halflength_, radius_;

  public:
    explicit SDFCylinder(Real halflength, Real radius) : halflength_(halflength), radius_(radius) {}
    void setParameters(Real halflength, Real radius);
    Real operator()(const Vec3d &point) const;
    BoundingBox3d findBounds() const;
};

class SDFCapsule
{
    Real halflength_, radius_;

  public:
    explicit SDFCapsule(Real halflength, Real radius) : halflength_(halflength), radius_(radius) {}
    void setParameters(Real halflength, Real radius);
    Real operator()(const Vec3d &point) const;
    BoundingBox3d findBounds() const;
};

class SDFCone
{
    Real halfheight_, radius_;

  public:
    explicit SDFCone(Real halfheight, Real radius) : halfheight_(halfheight), radius_(radius) {}
    void setParameters(Real halfheight, Real radius);
    Real operator()(const Vec3d &point) const;
    BoundingBox3d findBounds() const;
};

class SDFCappedCone
{
    Real halfheight_, radius1_, radius2_;

  public:
    explicit SDFCappedCone(Real halfheight, Real radius1, Real radius2)
        : halfheight_(halfheight), radius1_(radius1), radius2_(radius2) {}
    void setParameters(Real halfheight, Real radius1, Real radius2);
    Real operator()(const Vec3d &point) const;
    BoundingBox3d findBounds() const;
};
//----------------------------------------------------------------------
// Extended geometric primitives
//----------------------------------------------------------------------
template <typename InputType, typename ExtensionType>
class SDFExtension
{
    InputType input_;
    ExtensionType extend_;

  public:
    explicit SDFExtension(const InputType &input, const ExtensionType &extension);
    template <typename... InputArgs, typename... ExtendArgs>
    void setParameters(const InputArgs &...inputArgs, const ExtendArgs &...extensionArgs);
    Real operator()(const Vec3d &point) const { return extend_(input_, point); }
    BoundingBox3d findBounds() const { return extend_.findBounds(input_); }
};

class SDFRound
{
    Real radius_;

  public:
    explicit SDFRound(Real radius) : radius_(radius) {}
    void setParameters(Real radius) { radius_ = radius; }
    template <typename InputType>
    Real operator()(const InputType &input, const Vec3d &point) const { return input(point) - radius_; }
    template <typename InputType>
    auto findBounds(const InputType &input) const { return input.findBounds().expand(radius_); }
};

class SDFOnion
{
    Real radius_;

  public:
    explicit SDFOnion(Real radius) : radius_(radius) {}
    void setParameters(Real radius) { radius_ = radius; }
    template <typename InputType>
    Real operator()(const InputType &input, const Vec3d &point) const { return ABS(input(point)) - radius_; }
    template <typename InputType>
    auto findBounds(const InputType &input) const { return input.findBounds().expand(SMAX(radius_, 0.0)); }
};

class SDFChamfer
{
    Real chamfer_size_;

  public:
    explicit SDFChamfer(Real chamfer_size) : chamfer_size_(chamfer_size) {}
    void setParameters(Real chamfer_size) { chamfer_size_ = chamfer_size; }
    template <typename InputType>
    Real operator()(const InputType &input, const Vec3d &point) const;
    template <typename InputType>
    auto findBounds(const InputType &input) const;
};

class SDFScale
{
    Real scale_factor_;

  public:
    explicit SDFScale(Real scale_factor) : scale_factor_(scale_factor) {}
    Real getParameters() const { return scale_factor_; }
    void setParameters(Real scale_factor) { scale_factor_ = scale_factor; }
    template <typename InputType>
    Real operator()(const InputType &input, const Vec3d &point) const;
    template <typename InputType>
    auto findBounds(const InputType &input) const;
};

class SDFTransform
{
    Transform3d transform_;

  public:
    template <typename... Args>
    explicit SDFTransform(Args &&...args) : transform_(std::forward<Args>(args)...) {}
    template <typename... Args>
    void setParameters(Args &&...args);
    template <typename InputType>
    Real operator()(const InputType &input, const Vec3d &point) const;
    template <typename InputType>
    auto findBounds(const InputType &input) const;
};

template <typename OperationType, typename InputType1, typename InputType2>
class SDFOperation
{
    OperationType operation_;
    InputType1 input1_;
    InputType2 input2_;

  public:
    explicit SDFOperation(const OperationType &operation, const InputType1 &input1, const InputType2 &input2)
        : operation_(operation), input1_(input1), input2_(input2) {}
    template <typename... OperationTypeArgs, typename... InputType1Args, typename... InputType2Args>
    void setParameters(OperationTypeArgs &&...operationArgs, InputType1Args &&...input1Args, InputType2Args &&...input2Args);
    Real operator()(const Vec3d &point) const { return operation_(point, input1_, input2_); }
    auto findBounds() const { return operation_.findBounds(input1_, input2_); }
};

struct SDFAddition
{
    template <typename Input1, typename Input2>
    Real operator()(const Vec3d &point, const Input1 &input1, const Input2 &input2) const;
    template <typename Input1, typename Input2>
    auto findBounds(const Input1 &input1, const Input2 &input2) const;
};

struct SDFSubtraction
{
    template <typename Input1, typename Input2>
    Real operator()(const Vec3d &point, const Input1 &input1, const Input2 &input2) const;
    template <typename Input1, typename Input2>
    auto findBounds(const Input1 &input1, const Input2 &input2) const { return input1.findBounds(); }
};

struct SDFIntersection
{
    template <typename Input1, typename Input2>
    Real operator()(const Vec3d &point, const Input1 &input1, const Input2 &input2) const;
    template <typename Input1, typename Input2>
    auto findBounds(const Input1 &input1, const Input2 &input2) const;
};

class SDFSmoothAddition
{
    Real smoothing_length_; // scaled with mesh size

  public:
    explicit SDFSmoothAddition(Real smoothing_length) : smoothing_length_(smoothing_length) {}
    void setParameters(Real smoothing_length) { smoothing_length_ = smoothing_length; }
    template <typename Input1, typename Input2>
    Real operator()(const Vec3d &point, const Input1 &input1, const Input2 &input2) const;
    template <typename Input1, typename Input2>
    auto findBounds(const Input1 &input1, const Input2 &input2) const;
};

class SDFSmoothSubtraction
{
    Real smoothing_length_; // scaled with mesh size

  public:
    explicit SDFSmoothSubtraction(Real smoothing_length) : smoothing_length_(smoothing_length) {}
    void setParameters(Real smoothing_length) { smoothing_length_ = smoothing_length; }
    template <typename Input1, typename Input2>
    Real operator()(const Vec3d &point, const Input1 &input1, const Input2 &input2) const;
    template <typename Input1, typename Input2>
    auto findBounds(const Input1 &input1, const Input2 &input2) const;
};

class SDFSmoothIntersection
{
    Real smoothing_length_; // scaled with mesh size

  public:
    explicit SDFSmoothIntersection(Real smoothing_length) : smoothing_length_(smoothing_length) {}
    void setParameters(Real smoothing_length) { smoothing_length_ = smoothing_length; }
    template <typename Input1, typename Input2>
    Real operator()(const Vec3d &point, const Input1 &input1, const Input2 &input2) const;
    template <typename Input1, typename Input2>
    auto findBounds(const Input1 &input1, const Input2 &input2) const;
};
//----------------------------------------------------------------------
// 3D geometric primitives derived from 2D primitives
//----------------------------------------------------------------------
class SDFExtrusion
{
    Real height_;

  public:
    explicit SDFExtrusion(Real height) : height_(height) {}
    void setParameters(Real height) { height_ = height; }
    template <typename Input2D>
    Real operator()(const Input2D &input, const Vec3d &point) const;
};

class SDFRotation
{
    Real angle_;

  public:
    explicit SDFRotation(Real angle) : angle_(angle) {}
    void setParameters(Real angle) { angle_ = angle; }
    template <typename Input2D>
    Real operator()(const Input2D &input, const Vec3d &point) const;
};
//----------------------------------------------------------------------
// 3D geometric primitives derived from 3D primitives
//----------------------------------------------------------------------
class SDFElongation
{
    Real elongation_factor_;

  public:
    explicit SDFElongation(Real elongation_factor) : elongation_factor_(elongation_factor) {}
    Real getParameters() const { return elongation_factor_; }
    void setParameters(Real elongation_factor) { elongation_factor_ = elongation_factor; }
    template <typename Input3D>
    Real operator()(const Input3D &input, const Vec3d &point) const;
};
} // namespace SPH

#endif // SDF_PRIMITIVE_H
