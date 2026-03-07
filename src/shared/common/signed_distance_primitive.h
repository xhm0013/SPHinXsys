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
 * @file 	signed_distance_primitive.h
 * @brief 	All primitives use local coordinates. All rotation is around the x-axis.
 * @details Here, we only give popular primitives,
 * For more signed distance function, please check the website: https://iquilezles.org/articles/.
 * @author	Xiangyu Hu
 */

#ifndef SIGNED_DISTANCE_PRIMITIVE_H
#define SIGNED_DISTANCE_PRIMITIVE_H

#include <base_data_type.h>

namespace SPH
{
//----------------------------------------------------------------------
// Shared geometric primitives for signed distance function definition
//----------------------------------------------------------------------
class SDBall
{
    Real radius_;

  public:
    explicit SDBall(Real radius) : radius_(radius) {}
    Real getParameters() const { return radius_; }
    void setParameters(Real radius) { radius_ = radius; }
    template <typename VecType>
    Real operator()(const VecType &point) const
    {
        return point.norm() - radius_;
    }
};

template <int N>
class SDBox
{
    Eigen::Matrix<Real, N, 1> halfsize_;

  public:
    explicit SDBox(const Eigen::Matrix<Real, N, 1> &halfsize) : halfsize_(halfsize) {}
    Eigen::Matrix<Real, N, 1> getParameters() const { return halfsize_; }
    void setParameters(const Eigen::Matrix<Real, N, 1> &halfsize) { halfsize_ = halfsize; }
    Real operator()(const Eigen::Matrix<Real, N, 1> &point) const
    {
        Eigen::Matrix<Real, N, 1> d = point.cwiseAbs() - halfsize_;
        return d.cwiseMax(Eigen::Matrix<Real, N, 1>::Zero()).norm() + SMIN(d.maxCoeff(), 0.0);
    }
};

class SDRound
{
    Real radius_;

  public:
    explicit SDRound(Real radius) : radius_(radius) {}
    Real getParameters() const { return radius_; }
    void setParameters(Real radius) { radius_ = radius; }
    template <typename VecType, typename Input>
    Real operator()(const VecType &point, const Input &input) const
    {
        return input(point) - radius_;
    }
};

class SDOnion
{
    Real radius_;

  public:
    explicit SDOnion(Real radius) : radius_(radius) {}
    Real getParameters() const { return radius_; }
    void setParameters(Real radius) { radius_ = radius; }
    template <typename VecType, typename Input>
    Real operator()(const VecType &point, const Input &input) const
    {
        return ABS(input(point)) - radius_;
    }
};

class SDChamfer
{
    Real chamfer_size_;

  public:
    explicit SDChamfer(Real chamfer_size) : chamfer_size_(chamfer_size) {}
    Real getParameters() const { return chamfer_size_; }
    void setParameters(Real chamfer_size) { chamfer_size_ = chamfer_size; }
    template <typename VecType, typename Input>
    Real operator()(const VecType &point, const Input &input) const
    {
        return input(point) + chamfer_size_ - sqrt(chamfer_size_ * chamfer_size_ + input(point) * input(point));
    }
};

class SDScale
{
    Real scale_factor_;

  public:
    explicit SDScale(Real scale_factor) : scale_factor_(scale_factor) {}
    Real getParameters() const { return scale_factor_; }
    void setParameters(Real scale_factor) { scale_factor_ = scale_factor; }
    template <typename VecType, typename Input>
    Real operator()(const VecType &point, const Input &input) const
    {
        return input(point / scale_factor_) * scale_factor_;
    }
};

struct SDUnion
{
    template <typename VecType, typename Input1, typename Input2>
    Real operator()(const VecType &point, const Input1 &input1, const Input2 &input2) const
    {
        return SMIN(input1(point), input2(point));
    }
};

struct SDsubtraction
{
    template <typename VecType, typename Input1, typename Input2>
    Real operator()(const VecType &point, const Input1 &input1, const Input2 &input2) const
    {
        return SMAX(input1(point), -input2(point));
    }
};

struct SDIntersection
{
    template <typename Input1, typename Input2>
    Real operator()(const Vec2d &point, const Input1 &input1, const Input2 &input2) const
    {
        return SMAX(input1(point), input2(point));
    }
};

class SDSmoothUnion
{
    Real blend_factor_;

  public:
    explicit SDSmoothUnion(Real blend_factor) : blend_factor_(blend_factor) {}
    Real getParameters() const { return blend_factor_; }
    void setParameters(Real blend_factor) { blend_factor_ = blend_factor; }
    template <typename VecType, typename Input1, typename Input2>
    Real operator()(const VecType &point, const Input1 &input1, const Input2 &input2) const
    {
        Real d1 = input1(point);
        Real d2 = input2(point);
        Real h = std::clamp(0.5 + 0.5 * (d2 - d1) / blend_factor_, 0.0, 1.0);
        return d1 * (1.0 - h) + d2 * h - blend_factor_ * h * (1.0 - h);
    }
};

class SDSmoothSubtraction
{
    Real blend_factor_;

  public:
    explicit SDSmoothSubtraction(Real blend_factor) : blend_factor_(blend_factor) {}
    Real getParameters() const { return blend_factor_; }
    void setParameters(Real blend_factor) { blend_factor_ = blend_factor; }
    template <typename VecType, typename Input1, typename Input2>
    Real operator()(const VecType &point, const Input1 &input1, const Input2 &input2) const
    {
        Real d1 = input1(point);
        Real d2 = -input2(point);
        Real h = std::clamp(0.5 + 0.5 * (d2 - d1) / blend_factor_, 0.0, 1.0);
        return d1 * (1.0 - h) + d2 * h - blend_factor_ * h * (1.0 - h);
    }
};

class SDSmoothIntersection
{
    Real blend_factor_;

  public:
    explicit SDSmoothIntersection(Real blend_factor) : blend_factor_(blend_factor) {}
    Real getParameters() const { return blend_factor_; }
    void setParameters(Real blend_factor) { blend_factor_ = blend_factor; }
    template <typename VecType, typename Input1, typename Input2>
    Real operator()(const VecType &point, const Input1 &input1, const Input2 &input2) const
    {
        Real d1 = input1(point);
        Real d2 = input2(point);
        Real h = std::clamp(0.5 + 0.5 * (d2 - d1) / blend_factor_, 0.0, 1.0);
        return d1 * (1.0 - h) + d2 * h + blend_factor_ * h * (1.0 - h);
    }
};

class SDSymmetry
{
    int axis_; // 0 for x-axis, 1 for y-axis, 2 for z-axis

  public:
    explicit SDSymmetry(int axis) : axis_(axis) {}
    int getParameters() const { return axis_; }
    void setParameters(int axis) { axis_ = axis; }
    template <typename VecType, typename Input>
    Real operator()(const VecType &point, const Input &input) const
    {
        VecType symmetric_point = point;
        symmetric_point[axis_] = ABS(point[axis_]);
        return input(symmetric_point);
    }
};

class SDRepeat
{
    Vec3d period_;

  public:
    explicit SDRepeat(const Vec3d &period) : period_(period) {}
    Vec3d getParameters() const { return period_; }
    void setParameters(const Vec3d &period) { period_ = period; }
    template <typename VecType, typename Input>
    Real operator()(const VecType &point, const Input &input) const
    {
        VecType repeated_point = point;
        for (int i = 0; i < 3; ++i)
        {
            if (period_[i] > 0.0)
                repeated_point[i] = std::fmod(point[i] + 0.5 * period_[i], period_[i]) - 0.5 * period_[i];
        }
        return input(repeated_point);
    }
};
//----------------------------------------------------------------------
// 3D geometric primitives for signed distance function definition
//----------------------------------------------------------------------
class SDRoundedBox
{
    template <typename VecType>
    Real operator()(const VecType &point, const VecType &halfsize, Real radius) const
    {
        VecType d = point.cwiseAbs() - halfsize + VecType::Constant(radius);
        return d.cwiseMax(VecType::Zero()).norm() - radius + SMIN(d.maxCoeff(), 0.0) - radius;
    }
};

class SDCylinder
{
    Real operator()(const Vec3d &point, Real halflength, Real radius) const
    {
        Real axial = point[0];
        Real radial_distance = point.tail(2).norm();
        Real dh = ABS(axial) - halflength;
        Real dr = radial_distance - radius;
        return SMAX(dh, dr);
    }
};

class SDRoundedCylinder
{
    Real operator()(const Vec3d &point, Real halflength, Real radius) const
    {
        Real axial = point[0];
        Real radial_distance = point.tail(2).norm();
        Real dh = ABS(axial) - halflength + radius;
        Real dr = radial_distance - radius;
        return SMAX(dh, dr) - radius;
    }
};

class SDCapsule
{
    Real operator()(const Vec3d &point, Real halflength, Real radius) const
    {
        Real axial = point[0];
        Real radial_distance = point.tail(2).norm();
        if (axial < 0.0)
            return (point - Vec3d(0.0, 0.0, 0.0)).norm() - radius; // bottom hemisphere
        else if (axial > halflength)
            return (point - Vec3d(halflength, 0.0, 0.0)).norm() - radius; // top hemisphere
        else
            return radial_distance - radius; // cylindrical part
    }
};

class SDCone
{
    Real operator()(const Vec3d &point, Real halflength, Real radius) const
    {
        Real axial = point[0];
        Real radial_distance = point.tail(2).norm();
        if (axial < 0.0 || axial > halflength)
            return std::numeric_limits<Real>::max(); // outside the cone's height
        Real local_radius = (halflength - axial) / halflength * radius;
        return radial_distance - local_radius;
    }
};

class SDRoundedCone
{
    Real operator()(const Vec3d &point, Real halflength, Real radius) const
    {
        Real axial = point[0];
        Real radial_distance = point.tail(2).norm();
        if (axial < 0.0 || axial > halflength)
            return std::numeric_limits<Real>::max(); // outside the cone's height
        Real local_radius = (halflength - axial) / halflength * radius;
        return radial_distance - local_radius + radius;
    }
};

class SDCappedCone
{
    Real operator()(const Vec3d &point, Real halflength, Real radius) const
    {
        Real axial = point[0];
        Real radial_distance = point.tail(2).norm();
        if (axial < 0.0)
            return (point - Vec3d(0.0, 0.0, 0.0)).norm() - radius; // bottom hemisphere
        else if (axial > halflength)
            return (point - Vec3d(halflength, 0.0, 0.0)).norm() - radius; // top hemisphere
        else
        {
            Real local_radius = (halflength - axial) / halflength * radius;
            return radial_distance - local_radius; // conical part
        }
    }
};
//----------------------------------------------------------------------
// 2D geometric primitives for signed distance function definition
//----------------------------------------------------------------------
class SDTrapezoid
{
    Real operator()(const Vec2d &point, Real halflength, Real top_halfwidth, Real bottom_halfwidth) const
    {
        Real axial = point[0];
        Real lateral = point[1];
        if (axial < 0.0 || axial > halflength)
            return std::numeric_limits<Real>::max(); // outside the trapezoid's height
        Real local_halfwidth = (halflength - axial) / halflength * top_halfwidth + axial / halflength * bottom_halfwidth;
        return ABS(lateral) - local_halfwidth;
    }
};

class SDParallelogram
{
    Real operator()(const Vec2d &point, Real halflength, Real halfwidth, Real skew) const
    {
        Real axial = point[0];
        Real lateral = point[1];
        if (axial < 0.0 || axial > halflength)
            return std::numeric_limits<Real>::max(); // outside the parallelogram's height
        Real local_halfwidth = halfwidth + skew * (axial / halflength - 0.5);
        return ABS(lateral) - local_halfwidth;
    }
};

class SDEquilateralTriangle
{
    Real operator()(const Vec2d &point, Real halflength) const
    {
        Real axial = point[0];
        Real lateral = point[1];
        if (axial < 0.0 || axial > halflength)
            return std::numeric_limits<Real>::max(); // outside the triangle's height
        Real local_halfwidth = (halflength - axial) / halflength * 0.5 * halflength;
        return ABS(lateral) - local_halfwidth;
    }
};

class SDIsoscelesTriangle
{
    Real operator()(const Vec2d &point, Real halflength, Real top_halfwidth, Real bottom_halfwidth) const
    {
        Real axial = point[0];
        Real lateral = point[1];
        if (axial < 0.0 || axial > halflength)
            return std::numeric_limits<Real>::max(); // outside the triangle's height
        Real local_halfwidth = (halflength - axial) / halflength * top_halfwidth + axial / halflength * bottom_halfwidth;
        return ABS(lateral) - local_halfwidth;
    }
};

class SDTriangle
{
    Real operator()(const Vec2d &point, const Vec2d &vertex1, const Vec2d &vertex2, const Vec2d &vertex3) const
    {
        // Barycentric technique for signed distance to triangle
        Vec2d v0 = vertex2 - vertex1;
        Vec2d v1 = vertex3 - vertex1;
        Vec2d v2 = point - vertex1;

        Real d00 = v0.dot(v0);
        Real d01 = v0.dot(v1);
        Real d11 = v1.dot(v1);
        Real d20 = v2.dot(v0);
        Real d21 = v2.dot(v1);

        Real denom = d00 * d11 - d01 * d01;
        if (denom == 0.0)
            return std::numeric_limits<Real>::max(); // Degenerate triangle

        Real inv_denom = 1.0 / denom;
        Real u = (d11 * d20 - d01 * d21) * inv_denom;
        Real v = (d00 * d21 - d01 * d20) * inv_denom;

        if (u >= 0.0 && v >= 0.0 && u + v <= 1.0)
            return -std::sqrt((v2 - u * v0 - v * v1).squaredNorm()); // Inside triangle
        else
            return std::sqrt((v2 - u * v0 - v * v1).squaredNorm()); // Outside triangle
    }
};

class SDQuadrilateral
{
    Real operator()(const Vec2d &point, const Vec2d &vertex1, const Vec2d &vertex2, const Vec2d &vertex3, const Vec2d &vertex4) const
    {
        // Approximate signed distance to quadrilateral by taking the minimum distance to its edges
        SDEquilateralTriangle tri1;
        SDEquilateralTriangle tri2;
        Real d1 = tri1(point, vertex1, vertex2, vertex3);
        Real d2 = tri2(point, vertex1, vertex3, vertex4);
        return SMIN(d1, d2);
    }
};

class SDPolygon
{
    Real operator()(const Vec2d &point, const std::vector<Vec2d> &vertices) const
    {
        // Approximate signed distance to polygon by taking the minimum distance to its edges
        Real min_distance = std::numeric_limits<Real>::max();
        size_t n = vertices.size();
        for (size_t i = 0; i < n; ++i)
        {
            Vec2d v1 = vertices[i];
            Vec2d v2 = vertices[(i + 1) % n];
            SDEquilateralTriangle tri;
            Real d = tri(point, v1, v2, point); // Treat edge as a degenerate triangle
            min_distance = SMIN(min_distance, d);
        }
        return min_distance;
    }
};

class SDPie
{
    Real operator()(const Vec2d &point, Real radius, Real start_angle, Real end_angle) const
    {
        Real angle = std::atan2(point[1], point[0]);
        if (angle < 0.0)
            angle += 2.0 * M_PI; // Normalize angle to [0, 2π]
        if (angle >= start_angle && angle <= end_angle)
            return point.norm() - radius; // Inside the pie slice
        else
            return std::numeric_limits<Real>::max(); // Outside the pie slice
    }
};

class SDAnnulus
{
    Real operator()(const Vec2d &point, Real inner_radius, Real outer_radius) const
    {
        Real r = point.norm();
        if (r < inner_radius)
            return inner_radius - r; // Inside the inner circle
        else if (r > outer_radius)
            return r - outer_radius; // Outside the outer circle
        else
            return -SMIN(r - inner_radius, outer_radius - r); // Inside the annulus
    }
};

class SDCutDisk
{
    Real operator()(const Vec2d &point, Real radius, Real cut_angle) const
    {
        Real angle = std::atan2(point[1], point[0]);
        if (angle < 0.0)
            angle += 2.0 * M_PI; // Normalize angle to [0, 2π]
        if (angle <= cut_angle)
            return point.norm() - radius; // Inside the cut disk
        else
            return std::numeric_limits<Real>::max(); // Outside the cut disk
    }
};

class SDWedge
{
    Real operator()(const Vec2d &point, Real radius, Real start_angle, Real end_angle) const
    {
        Real angle = std::atan2(point[1], point[0]);
        if (angle < 0.0)
            angle += 2.0 * M_PI; // Normalize angle to [0, 2π]
        if (angle >= start_angle && angle <= end_angle)
            return point.norm() - radius; // Inside the wedge
        else
            return std::numeric_limits<Real>::max(); // Outside the wedge
    }
};

class SDArc
{
    Real operator()(const Vec2d &point, Real radius, Real start_angle, Real end_angle) const
    {
        Real angle = std::atan2(point[1], point[0]);
        if (angle < 0.0)
            angle += 2.0 * M_PI; // Normalize angle to [0, 2π]
        if (angle >= start_angle && angle <= end_angle)
            return point.norm() - radius; // Inside the arc
        else
            return std::numeric_limits<Real>::max(); // Outside the arc
    }
};

class SDEllipse
{
    Real operator()(const Vec2d &point, Real a, Real b) const
    {
        return std::sqrt((point[0] * point[0]) / (a * a) + (point[1] * point[1]) / (b * b)) - 1.0;
    }
};
//----------------------------------------------------------------------
// 3D geometric primitives derived from 2D primitives
//----------------------------------------------------------------------
class SDExtrusion
{
    Real height_;

  public:
    explicit SDExtrusion(Real height) : height_(height) {}
    Real getParameters() const { return height_; }
    void setParameters(Real height) { height_ = height; }
    template <typename Input2D>
    Real operator()(const Vec3d &point, const Input2D &input) const
    {
        Vec2d point_2d = point.tail(2); // Project to 2D plane
        Real axial = point[0];
        if (axial < 0.0 || axial > height_)
            return MaxReal; // Outside the extrusion's height
        return input(point_2d);
    }
};

class SDRotation
{
    Real angle_;

  public:
    explicit SDRotation(Real angle) : angle_(angle) {}
    Real getParameters() const { return angle_; }
    void setParameters(Real angle) { angle_ = angle; }
    template <typename Input2D>
    Real operator()(const Vec3d &point, const Input2D &input) const
    {
        Real cos_angle = std::cos(angle_);
        Real sin_angle = std::sin(angle_);
        Vec2d rotated_point;
        rotated_point[0] = cos_angle * point[0] - sin_angle * point[1];
        rotated_point[1] = sin_angle * point[0] + cos_angle * point[1];
        return input(rotated_point);
    }
};
//----------------------------------------------------------------------
// 3D geometric primitives derived from 3D primitives
//----------------------------------------------------------------------
class SDElongation
{
    Real elongation_factor_;

  public:
    explicit SDElongation(Real elongation_factor) : elongation_factor_(elongation_factor) {}
    Real getParameters() const { return elongation_factor_; }
    void setParameters(Real elongation_factor) { elongation_factor_ = elongation_factor; }
    template <typename Input3D>
    Real operator()(const Vec3d &point, const Input3D &input) const
    {
        Vec3d elongated_point = point;
        elongated_point[0] *= elongation_factor_; // Elongate along x-axis
        return input(elongated_point);
    }
};
} // namespace SPH

#endif // SIGNED_DISTANCE_PRIMITIVE_H
