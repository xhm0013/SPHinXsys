#ifndef SDF_PRIMITIVE_HPP
#define SDF_PRIMITIVE_HPP

#include "sdf_primitive.h"

namespace SPH
{
//=================================================================================================//
template <typename InputType, typename ExtensionType>
SDFExtension<InputType, ExtensionType>::SDFExtension(const InputType &input, const ExtensionType &extension)
    : input_(input), extend_(extension) {}
//=================================================================================================//
template <typename InputType, typename ExtensionType>
template <typename... InputArgs, typename... ExtendArgs>
void SDFExtension<InputType, ExtensionType>::setParameters(
    const InputArgs &...inputArgs, const ExtendArgs &...extensionArgs)
{
    input_.setParameters(inputArgs...);
    extend_.setParameters(extensionArgs...);
}
//=================================================================================================//
template <typename InputType>
Real SDFChamfer::operator()(const InputType &input, const Vec3d &point) const
{
    Real sd = input(point);
    return sd + chamfer_size_ - sqrt(chamfer_size_ * chamfer_size_ + sd * sd);
}
//=================================================================================================//
template <typename InputType>
auto SDFChamfer::findBounds(const InputType &input) const
{
    return input.findBounds().expand(chamfer_size_);
}
//=================================================================================================//
template <typename InputType>
Real SDFScale::operator()(const InputType &input, const Vec3d &point) const
{
    return input(point / scale_factor_) * scale_factor_;
}
//=================================================================================================//
template <typename InputType>
auto SDFScale::findBounds(const InputType &input) const
{
    return input.findBounds().scale(scale_factor_);
}
//=================================================================================================//
template <typename... Args>
void SDFTransform::setParameters(Args &&...args)
{
    transform_ = Transform3d(std::forward<Args>(args)...);
}
//=================================================================================================//
template <typename InputType>
Real SDFTransform::operator()(const InputType &input, const Vec3d &point) const
{
    Vec3d transformed_point = transform_.shiftBaseStationToFrame(point);
    return input(transformed_point);
}
//=================================================================================================//
template <typename InputType>
auto SDFTransform::findBounds(const InputType &input) const
{
    BoundingBox3d original_bound = input.findBounds();
    Vec3d bb_min = Vec3d::Constant(MaxReal);
    Vec3d bb_max = Vec3d::Constant(-MaxReal);
    for (auto x : {original_bound.lower_.x(), original_bound.upper_.x()})
    {
        for (auto y : {original_bound.lower_.y(), original_bound.upper_.y()})
        {
            for (auto z : {original_bound.lower_.z(), original_bound.upper_.z()})
            {
                bb_min = bb_min.cwiseMin(transform_.shiftFrameStationToBase(Vec3d(x, y, z)));
                bb_max = bb_max.cwiseMax(transform_.shiftFrameStationToBase(Vec3d(x, y, z)));
            }
        }
    }
    return BoundingBox3d(bb_min, bb_max);
}
//=================================================================================================//
template <typename OperationType, typename InputType1, typename InputType2>
template <typename... OperationTypeArgs, typename... InputType1Args, typename... InputType2Args>
void SDFOperation<OperationType, InputType1, InputType2>::setParameters(
    OperationTypeArgs &&...operationArgs, InputType1Args &&...input1Args, InputType2Args &&...input2Args)
{
    operation_.setParameters(std::forward<OperationTypeArgs>(operationArgs)...);
    input1_.setParameters(std::forward<InputType1Args>(input1Args)...);
    input2_.setParameters(std::forward<InputType2Args>(input2Args)...);
}
//=================================================================================================//
template <typename Input1, typename Input2>
Real SDFAddition::operator()(const Vec3d &point, const Input1 &input1, const Input2 &input2) const
{
    return SMIN(input1(point), input2(point));
}
//=================================================================================================//
template <typename Input1, typename Input2>
auto SDFAddition::findBounds(const Input1 &input1, const Input2 &input2) const
{
    return input1.findBounds().add(input2.findBounds());
}
//=================================================================================================//
template <typename Input1, typename Input2>
Real SDFSubtraction::operator()(const Vec3d &point, const Input1 &input1, const Input2 &input2) const
{
    return SMAX(input1(point), -input2(point));
}
//=================================================================================================//
template <typename Input1, typename Input2>
Real SDFIntersection::operator()(const Vec3d &point, const Input1 &input1, const Input2 &input2) const
{
    return SMAX(input1(point), input2(point));
}
//=================================================================================================//
template <typename Input1, typename Input2>
auto SDFIntersection::findBounds(const Input1 &input1, const Input2 &input2) const
{
    return input1.findBounds().intersect(input2.findBounds());
}
//=================================================================================================//
template <typename Input1, typename Input2>
Real SDFSmoothAddition::operator()(
    const Vec3d &point, const Input1 &input1, const Input2 &input2) const
{
    Real d1 = input1(point);
    Real d2 = input2(point);
    Real k = 4.0 * smoothing_length_;
    Real h = SMAX(k - ABS(d1 - d2), 0.0);
    return SMIN(d1, d2) - h * h * 0.25 / k;
}
//=================================================================================================//
template <typename Input1, typename Input2>
auto SDFSmoothAddition::findBounds(const Input1 &input1, const Input2 &input2) const
{
    return input1.findBounds().add(input2.findBounds());
}
//=================================================================================================//
template <typename Input1, typename Input2>
Real SDFSmoothSubtraction::operator()(
    const Vec3d &point, const Input1 &input1, const Input2 &input2) const
{
    Real d1 = -input1(point);
    Real d2 = input2(point);
    Real k = 4.0 * smoothing_length_;
    Real h = SMAX(k - ABS(d1 - d2), 0.0);
    return -(SMIN(d1, d2) - h * h * 0.25 / k);
}
//=================================================================================================//
template <typename Input1, typename Input2>
auto SDFSmoothSubtraction::findBounds(const Input1 &input1, const Input2 &input2) const
{
    return input1.findBounds();
}
//=================================================================================================//
template <typename Input1, typename Input2>
Real SDFSmoothIntersection::operator()(
    const Vec3d &point, const Input1 &input1, const Input2 &input2) const
{
    Real d1 = -input1(point);
    Real d2 = -input2(point);
    Real k = 4.0 * smoothing_length_;
    Real h = SMAX(k - ABS(d1 - d2), 0.0);
    return -(SMIN(d1, d2) - h * h * 0.25 / k);
}
//=================================================================================================//
template <typename Input1, typename Input2>
auto SDFSmoothIntersection::findBounds(const Input1 &input1, const Input2 &input2) const
{
    return input1.findBounds().intersect(input2.findBounds());
}
//=================================================================================================//
template <typename Input2D>
Real SDFExtrusion::operator()(const Input2D &input, const Vec3d &point) const
{
    Vec2d point_2d = point.tail(2); // Project to 2D plane
    Real axial = point[0];
    if (axial < 0.0 || axial > height_)
        return MaxReal; // Outside the extrusion's height
    return input(point_2d);
}
//=================================================================================================//
template <typename Input2D>
Real SDFRotation::operator()(const Input2D &input, const Vec3d &point) const
{
    Real cos_angle = std::cos(angle_);
    Real sin_angle = std::sin(angle_);
    Vec2d rotated_point;
    rotated_point[0] = cos_angle * point[0] - sin_angle * point[1];
    rotated_point[1] = sin_angle * point[0] + cos_angle * point[1];
    return input(rotated_point);
}
//=================================================================================================//
template <typename Input3D>
Real SDFElongation::operator()(const Input3D &input, const Vec3d &point) const
{
    Vec3d elongated_point = point;
    elongated_point[0] *= elongation_factor_; // Elongate along x-axis
    return input(elongated_point);
}
//=================================================================================================//
} // namespace SPH
#endif // SDF_PRIMITIVE_HPP