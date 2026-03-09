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
// Shared geometric primitives for signed distance function definition
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
    Real operator()(const Vec3d &point) const
    {
        Vec3d d = point.cwiseAbs() - halfsize_;
        return d.cwiseMax(Vec3d::Zero()).norm() + SMIN(d.maxCoeff(), 0.0);
    }
    BoundingBox3d findBounds() const { return BoundingBox3d(halfsize_); }
};

class SDFCylinder
{
    Real halflength_, radius_;

  public:
    explicit SDFCylinder(Real halflength, Real radius)
        : halflength_(halflength), radius_(radius) {}

    void setParameters(Real halflength, Real radius)
    {
        halflength_ = halflength;
        radius_ = radius;
    }

    Real operator()(const Vec3d &point) const
    {
        Real axial = point[0];
        Real radial_distance = point.tail(2).norm();
        Real dh = ABS(axial) - halflength_;
        Real dr = radial_distance - radius_;
        return SMAX(dh, dr);
    }

    BoundingBox3d findBounds() const
    {
        return BoundingBox3d(Vec3d(halflength_, radius_, radius_));
    }
};

class SDFCapsule
{
    Real halflength_, radius_;

  public:
    explicit SDFCapsule(Real halflength, Real radius)
        : halflength_(halflength), radius_(radius) {}
    void setParameters(Real halflength, Real radius)
    {
        halflength_ = halflength;
        radius_ = radius;
    }
    Real operator()(const Vec3d &point) const
    {
        Real axial = point[0];
        Real radial_distance = point.tail(2).norm();
        if (axial < 0.0)
            return (point - Vec3d(0.0, 0.0, 0.0)).norm() - radius_; // bottom hemisphere
        else if (axial > halflength_)
            return (point - Vec3d(halflength_, 0.0, 0.0)).norm() - radius_; // top hemisphere
        else
            return radial_distance - radius_; // cylindrical part
    }

    BoundingBox3d findBounds() const
    {
        return BoundingBox3d(Vec3d(halflength_ + radius_, radius_, radius_));
    }
};

class SDFCone
{
    Real height_, theta_;
    Real sin_t_, cos_t_, tan_t_;

  public:
    explicit SDFCone(Real height, Real theta) : height_(height), theta_(theta)
    {
        sin_t_ = std::sin(theta_);
        cos_t_ = std::cos(theta_);
        tan_t_ = std::tan(theta_);
    }
    void setParameters(Real height, Real theta)
    {
        height_ = height;
        theta_ = theta;
        sin_t_ = std::sin(theta_);
        cos_t_ = std::cos(theta_);
        tan_t_ = std::tan(theta_);
    }
    Real operator()(const Vec3d &point) const
    {
        // c is the unit vector along the cone slope in the (x, radial) plane
        Vec2d c(cos_t_, sin_t_);

        // q is the point projected onto the (x, radial) plane
        // x is the axial distance, q_rad is the distance from the x-axis
        Real q_rad = point.segment<2>(1).norm(); // norm of (y, z)
        Vec2d q(point.x(), q_rad);
        // Project q onto the line defined by the cone slope
        // dot product gives the projection length
        Real dot_qc = q.dot(c);
        Vec2d w = q - c * std::clamp(dot_qc, 0.0, height_ / cos_t_);
        // Calculate distance to the base (cap)
        // The cap is at x = h, with radial distance <= h * tan(theta)
        Vec2d base_proj(q.x() - height_, q_rad - std::clamp(q_rad, 0.0, height_ * tan_t_));
        // Determine if the point is inside or outside for the sign
        // s1: height check, s2: slope check
        Real s = SMAX(q.dot(Vec2d(c.y(), -c.x())), q.x() - height_);
        return std::sqrt(SMIN(w.squaredNorm(), base_proj.squaredNorm())) * SGN(s);
    }

    BoundingBox3d findBounds() const
    {
        Real max_radius = height_ * tan_t_;
        return BoundingBox3d(Vec3d(height_ + max_radius, max_radius, max_radius));
    }
};

class SDFRoundedCone
{
    Real halflength_, radius_;

  public:
    explicit SDFRoundedCone(Real halflength, Real radius)
        : halflength_(halflength), radius_(radius) {}
    void setParameters(Real halflength, Real radius)
    {
        halflength_ = halflength;
        radius_ = radius;
    }
    Real operator()(const Vec3d &point) const
    {
        Real axial = point[0];
        Real radial_distance = point.tail(2).norm();
        if (axial < 0.0 || axial > halflength_)
            return MaxReal; // outside the cone's height
        Real local_radius = (halflength_ - axial) / halflength_ * radius_;
        return radial_distance - local_radius + radius_;
    }
    BoundingBox3d findBounds() const
    {
        return BoundingBox3d(Vec3d(halflength_ + radius_, radius_, radius_));
    }
};

class SDFCappedCone
{
    Real halflength_, radius_;

  public:
    explicit SDFCappedCone(Real halflength, Real radius)
        : halflength_(halflength), radius_(radius) {}
    void setParameters(Real halflength, Real radius)
    {
        halflength_ = halflength;
        radius_ = radius;
    }
    Real operator()(const Vec3d &point) const
    {
        Real axial = point[0];
        Real radial_distance = point.tail(2).norm();
        if (axial < 0.0)
            return (point - Vec3d(0.0, 0.0, 0.0)).norm() - radius_; // bottom hemisphere
        else if (axial > halflength_)
            return (point - Vec3d(halflength_, 0.0, 0.0)).norm() - radius_; // top hemisphere
        else
        {
            Real local_radius = (halflength_ - axial) / halflength_ * radius_;
            return radial_distance - local_radius; // conical part
        }
    }
    BoundingBox3d findBounds() const
    {
        return BoundingBox3d(Vec3d(halflength_ + radius_, radius_, radius_));
    }
};
//----------------------------------------------------------------------
// Extended geometric primitives
//----------------------------------------------------------------------
template <typename InputType, typename ExtendType>
class SDFExtend
{
    InputType input_;
    ExtendType extend_;

  public:
    explicit SDFExtend(const InputType &input, const ExtendType &extended) : input_(input), extend_(extended) {}
    template <typename... InputArgs, typename... ExtendArgs>
    void setParameters(const InputArgs &...inputArgs, const ExtendArgs &...extendedArgs)
    {
        input_.setParameters(inputArgs...);
        extend_.setParameters(extendedArgs...);
    }
    template <typename VecType>
    Real operator()(const VecType &point) const { return extend_(input_, point); }
    BoundingBox3d findBounds() const { return extend_.findBounds(input_); }
};

class SDFRound
{
    Real radius_;

  public:
    explicit SDFRound(Real radius) : radius_(radius) {}
    void setParameters(Real radius) { radius_ = radius; }
    template <typename InputType, typename VecType>
    Real operator()(const InputType &input, const VecType &point) const { return input(point) - radius_; }

    template <typename InputType>
    auto findBounds(const InputType &input) const
    {
        return input.findBounds().expand(radius_);
    }
};

class SDFOnion
{
    Real radius_;

  public:
    explicit SDFOnion(Real radius) : radius_(radius) {}
    void setParameters(Real radius) { radius_ = radius; }
    template <typename InputType, typename VecType>
    Real operator()(const InputType &input, const VecType &point) const { return ABS(input(point)) - radius_; }

    template <typename InputType>
    auto findBounds(const InputType &input) const
    {
        return input.findBounds().expand(SMAX(radius_, 0.0));
    }
};

class SDFChamfer
{
    Real chamfer_size_;

  public:
    explicit SDFChamfer(Real chamfer_size) : chamfer_size_(chamfer_size) {}
    void setParameters(Real chamfer_size) { chamfer_size_ = chamfer_size; }

    template <typename InputType, typename VecType>
    Real operator()(const InputType &input, const VecType &point) const
    {
        Real sd = input(point);
        return sd + chamfer_size_ - sqrt(chamfer_size_ * chamfer_size_ + sd * sd);
    }

    template <typename InputType>
    auto findBounds(const InputType &input) const
    {
        return input.findBounds().expand(chamfer_size_);
    }
};

class SDFScale
{
    Real scale_factor_;

  public:
    explicit SDFScale(Real scale_factor) : scale_factor_(scale_factor) {}
    Real getParameters() const { return scale_factor_; }
    void setParameters(Real scale_factor) { scale_factor_ = scale_factor; }
    template <typename InputType, typename VecType>
    Real operator()(const InputType &input, const VecType &point) const
    {
        return input(point / scale_factor_) * scale_factor_;
    }

    template <typename InputType>
    auto findBounds(const InputType &input) const
    {
        return input.findBounds().scale(scale_factor_);
    }
};

class SDFTransform
{
    Transform3d transform_;

  public:
    template <typename... Args>
    explicit SDFTransform(Args &&...args) : transform_(std::forward<Args>(args)...) {}

    template <typename... Args>
    void setParameters(Args &&...args)
    {
        transform_ = Transform3d(std::forward<Args>(args)...);
    }

    template <typename InputType, typename VecType>
    Real operator()(const InputType &input, const VecType &point) const
    {
        Vecd transformed_point = transform_.shiftBaseStationToFrame(point);
        return input(transformed_point);
    }

    template <typename InputType>
    auto findBounds(const InputType &input) const
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
};

struct SDFAddition
{
    template <typename VecType, typename Input1, typename Input2>
    Real operator()(const VecType &point, const Input1 &input1, const Input2 &input2) const
    {
        return SMIN(input1(point), input2(point));
    }
};

struct SDFSubtraction
{
    template <typename VecType, typename Input1, typename Input2>
    Real operator()(const VecType &point, const Input1 &input1, const Input2 &input2) const
    {
        return SMAX(input1(point), -input2(point));
    }
};

struct SDFIntersection
{
    template <typename Input1, typename Input2>
    Real operator()(const Vec2d &point, const Input1 &input1, const Input2 &input2) const
    {
        return SMAX(input1(point), input2(point));
    }
};

class SDFSmoothAddition
{
    Real blend_factor_;

  public:
    explicit SDFSmoothAddition(Real blend_factor) : blend_factor_(blend_factor) {}
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

class SDFSmoothSubtraction
{
    Real blend_factor_;

  public:
    explicit SDFSmoothSubtraction(Real blend_factor) : blend_factor_(blend_factor) {}
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

class SDFSmoothIntersection
{
    Real blend_factor_;

  public:
    explicit SDFSmoothIntersection(Real blend_factor) : blend_factor_(blend_factor) {}
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

class SDFSymmetry
{
    int axis_; // 0 for x-axis, 1 for y-axis, 2 for z-axis

  public:
    explicit SDFSymmetry(int axis) : axis_(axis) {}
    void setParameters(int axis) { axis_ = axis; }
    template <typename VecType, typename Input>
    Real operator()(const Input &input, const VecType &point) const
    {
        VecType symmetric_point = point;
        symmetric_point[axis_] = ABS(point[axis_]);
        return input(symmetric_point);
    }
};

class SDFRepeat
{
    Vec3d period_;

  public:
    explicit SDFRepeat(const Vec3d &period) : period_(period) {}
    void setParameters(const Vec3d &period) { period_ = period; }
    template <typename VecType, typename Input>
    Real operator()(const Input &input, const VecType &point) const
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
// 2D geometric primitives for signed distance function definition
//----------------------------------------------------------------------
class SDFTrapezoid
{
    Real halflength_, top_halfwidth_, bottom_halfwidth_;

  public:
    explicit SDFTrapezoid(Real halflength, Real top_halfwidth, Real bottom_halfwidth)
        : halflength_(halflength), top_halfwidth_(top_halfwidth), bottom_halfwidth_(bottom_halfwidth) {}
    void setParameters(Real halflength, Real top_halfwidth, Real bottom_halfwidth)
    {
        halflength_ = halflength;
        top_halfwidth_ = top_halfwidth;
        bottom_halfwidth_ = bottom_halfwidth;
    }
    Real operator()(const Vec2d &point) const
    {
        Real axial = point[0];
        Real lateral = point[1];
        if (axial < 0.0 || axial > halflength_)
            return MaxReal; // outside the trapezoid's height
        Real local_halfwidth = (halflength_ - axial) / halflength_ * top_halfwidth_ + axial / halflength_ * bottom_halfwidth_;
        return ABS(lateral) - local_halfwidth;
    }
};

class SDFParallelogram
{
    Real halflength_, halfwidth_, skew_;

  public:
    explicit SDFParallelogram(Real halflength, Real halfwidth, Real skew)
        : halflength_(halflength), halfwidth_(halfwidth), skew_(skew) {}
    Real operator()(const Vec2d &point) const
    {
        Real axial = point[0];
        Real lateral = point[1];
        if (axial < 0.0 || axial > halflength_)
            return MaxReal; // outside the parallelogram's height
        Real local_halfwidth = halfwidth_ + skew_ * (axial / halflength_ - 0.5);
        return ABS(lateral) - local_halfwidth;
    }
};

class SDFEquilateralTriangle
{
    Real halflength_;

  public:
    explicit SDFEquilateralTriangle(Real halflength) : halflength_(halflength) {}
    void setParameters(Real halflength) { halflength_ = halflength; }
    Real operator()(const Vec2d &point) const
    {
        Real axial = point[0];
        Real lateral = point[1];
        if (axial < 0.0 || axial > halflength_)
            return MaxReal; // outside the triangle's height
        Real local_halfwidth = (halflength_ - axial) / halflength_ * 0.5 * halflength_ * std::sqrt(3.0);
        return ABS(lateral) - local_halfwidth;
    }
};

class SDFIsoscelesTriangle
{
    Real halflength_, top_halfwidth_, bottom_halfwidth_;

  public:
    explicit SDFIsoscelesTriangle(Real halflength, Real top_halfwidth, Real bottom_halfwidth)
        : halflength_(halflength), top_halfwidth_(top_halfwidth), bottom_halfwidth_(bottom_halfwidth) {}
    void setParameters(Real halflength, Real top_halfwidth, Real bottom_halfwidth)
    {
        halflength_ = halflength;
        top_halfwidth_ = top_halfwidth;
        bottom_halfwidth_ = bottom_halfwidth;
    }

    Real operator()(const Vec2d &point) const
    {
        Real axial = point[0];
        Real lateral = point[1];
        if (axial < 0.0 || axial > halflength_)
            return MaxReal; // outside the triangle's height
        Real local_halfwidth = (halflength_ - axial) / halflength_ * top_halfwidth_ + axial / halflength_ * bottom_halfwidth_;
        return ABS(lateral) - local_halfwidth;
    }
};

class SDFTriangle
{
    Vec2d vertex1_, vertex2_, vertex3_;

  public:
    explicit SDFTriangle(const Vec2d &v1, const Vec2d &v2, const Vec2d &v3)
        : vertex1_(v1), vertex2_(v2), vertex3_(v3) {}
    void setParameters(const Vec2d &v1, const Vec2d &v2, const Vec2d &v3)
    {
        vertex1_ = v1;
        vertex2_ = v2;
        vertex3_ = v3;
    }

    Real operator()(const Vec2d &point) const
    {
        // Barycentric technique for signed distance to triangle
        Vec2d v0 = vertex2_ - vertex1_;
        Vec2d v1 = vertex3_ - vertex1_;
        Vec2d v2 = point - vertex1_;

        Real d00 = v0.dot(v0);
        Real d01 = v0.dot(v1);
        Real d11 = v1.dot(v1);
        Real d20 = v2.dot(v0);
        Real d21 = v2.dot(v1);

        Real denom = d00 * d11 - d01 * d01;
        if (denom == 0.0)
            return MaxReal; // Degenerate triangle

        Real inv_denom = 1.0 / denom;
        Real u = (d11 * d20 - d01 * d21) * inv_denom;
        Real v = (d00 * d21 - d01 * d20) * inv_denom;

        if (u >= 0.0 && v >= 0.0 && u + v <= 1.0)
            return -std::sqrt((v2 - u * v0 - v * v1).squaredNorm()); // Inside triangle
        else
            return std::sqrt((v2 - u * v0 - v * v1).squaredNorm()); // Outside triangle
    }
};

class SDFQuadrilateral
{
    Vec2d vertex1_, vertex2_, vertex3_, vertex4_;

  public:
    SDFQuadrilateral(const Vec2d &v1, const Vec2d &v2, const Vec2d &v3, const Vec2d &v4)
        : vertex1_(v1), vertex2_(v2), vertex3_(v3), vertex4_(v4) {}
    void setParameters(const Vec2d &v1, const Vec2d &v2, const Vec2d &v3, const Vec2d &v4)
    {
        vertex1_ = v1;
        vertex2_ = v2;
        vertex3_ = v3;
        vertex4_ = v4;
    }
    Real operator()(const Vec2d &point) const
    {
        SDFTriangle tri1(vertex1_, vertex2_, vertex3_);
        SDFTriangle tri2(vertex1_, vertex3_, vertex4_);
        Real d1 = tri1(point);
        Real d2 = tri2(point);
        return SMIN(d1, d2);
    }
};

class SDFPolygon
{
    std::vector<Vec2d> vertices_;

  public:
    explicit SDFPolygon(const std::vector<Vec2d> &vertices) : vertices_(vertices) {}
    void setParameters(const std::vector<Vec2d> &vertices)
    {
        vertices_ = vertices;
    }

    Real operator()(const Vec2d &point) const
    {
        // Approximate signed distance to polygon by taking the minimum distance to its edges
        Real min_distance = MaxReal;
        size_t n = vertices_.size();
        for (size_t i = 0; i < n; ++i)
        {
            Vec2d v1 = vertices_[i];
            Vec2d v2 = vertices_[(i + 1) % n];
            SDFTriangle tri(v1, v2, point);
            Real d = tri(point);
            min_distance = SMIN(min_distance, d);
        }
        return min_distance;
    }
};

class SDFPie
{
    Real radius_, start_angle_, end_angle_;

  public:
    explicit SDFPie(Real radius, Real start_angle, Real end_angle)
        : radius_(radius), start_angle_(start_angle), end_angle_(end_angle) {}
    void setParameters(Real radius, Real start_angle, Real end_angle)
    {
        radius_ = radius;
        start_angle_ = start_angle;
        end_angle_ = end_angle;
    }
    Real operator()(const Vec2d &point) const
    {
        Real angle = std::atan2(point[1], point[0]);
        if (angle < 0.0)
            angle += 2.0 * M_PI; // Normalize angle to [0, 2π]
        if (angle >= start_angle_ && angle <= end_angle_)
            return point.norm() - radius_; // Inside the pie slice
        else
            return MaxReal; // Outside the pie slice
    }
};

class SDFAnnulus
{
    Real inner_radius_, outer_radius_;

  public:
    explicit SDFAnnulus(Real inner_radius, Real outer_radius)
        : inner_radius_(inner_radius), outer_radius_(outer_radius) {}
    void setParameters(Real inner_radius, Real outer_radius)
    {
        inner_radius_ = inner_radius;
        outer_radius_ = outer_radius;
    }
    Real operator()(const Vec2d &point) const
    {
        Real r = point.norm();
        if (r < inner_radius_)
            return inner_radius_ - r; // Inside the inner circle
        else if (r > outer_radius_)
            return r - outer_radius_; // Outside the outer circle
        else
            return -SMIN(r - inner_radius_, outer_radius_ - r); // Inside the annulus
    }
};

class SDFCutDisk
{
    Real radius_, cut_angle_;

  public:
    explicit SDFCutDisk(Real radius, Real cut_angle)
        : radius_(radius), cut_angle_(cut_angle) {}
    void setParameters(Real radius, Real cut_angle)
    {
        radius_ = radius;
        cut_angle_ = cut_angle;
    }
    Real operator()(const Vec2d &point) const
    {
        Real angle = std::atan2(point[1], point[0]);
        if (angle < 0.0)
            angle += 2.0 * M_PI; // Normalize angle to [0, 2π]
        if (angle <= cut_angle_)
            return point.norm() - radius_; // Inside the cut disk
        else
            return MaxReal; // Outside the cut disk
    }
};

class SDFWedge
{
    Real radius_, start_angle_, end_angle_;

  public:
    explicit SDFWedge(Real radius, Real start_angle, Real end_angle)
        : radius_(radius), start_angle_(start_angle), end_angle_(end_angle) {}
    void setParameters(Real radius, Real start_angle, Real end_angle)
    {
        radius_ = radius;
        start_angle_ = start_angle;
        end_angle_ = end_angle;
    }
    Real operator()(const Vec2d &point) const
    {
        Real angle = std::atan2(point[1], point[0]);
        if (angle < 0.0)
            angle += 2.0 * M_PI; // Normalize angle to [0, 2π]
        if (angle >= start_angle_ && angle <= end_angle_)
            return point.norm() - radius_; // Inside the wedge
        else
            return MaxReal; // Outside the wedge
    }
};

class SDFArc
{
    Real radius_, start_angle_, end_angle_;

  public:
    explicit SDFArc(Real radius, Real start_angle, Real end_angle)
        : radius_(radius), start_angle_(start_angle), end_angle_(end_angle) {}
    void setParameters(Real radius, Real start_angle, Real end_angle)
    {
        radius_ = radius;
        start_angle_ = start_angle;
        end_angle_ = end_angle;
    }
    Real operator()(const Vec2d &point) const
    {
        Real angle = std::atan2(point[1], point[0]);
        if (angle < 0.0)
            angle += 2.0 * M_PI; // Normalize angle to [0, 2π]
        if (angle >= start_angle_ && angle <= end_angle_)
            return point.norm() - radius_; // Inside the arc
        else
            return MaxReal; // Outside the arc
    }
};

class SDFEllipse
{
    Real a_, b_;

  public:
    explicit SDFEllipse(Real a, Real b) : a_(a), b_(b) {}
    void setParameters(Real a, Real b)
    {
        a_ = a;
        b_ = b;
    }
    Real operator()(const Vec2d &point) const
    {
        return std::sqrt((point[0] * point[0]) / (a_ * a_) + (point[1] * point[1]) / (b_ * b_)) - 1.0;
    }
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
    Real operator()(const Input2D &input, const Vec3d &point) const
    {
        Vec2d point_2d = point.tail(2); // Project to 2D plane
        Real axial = point[0];
        if (axial < 0.0 || axial > height_)
            return MaxReal; // Outside the extrusion's height
        return input(point_2d);
    }
};

class SDFRotation
{
    Real angle_;

  public:
    explicit SDFRotation(Real angle) : angle_(angle) {}
    void setParameters(Real angle) { angle_ = angle; }
    template <typename Input2D>
    Real operator()(const Input2D &input, const Vec3d &point) const
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
class SDFElongation
{
    Real elongation_factor_;

  public:
    explicit SDFElongation(Real elongation_factor) : elongation_factor_(elongation_factor) {}
    Real getParameters() const { return elongation_factor_; }
    void setParameters(Real elongation_factor) { elongation_factor_ = elongation_factor; }
    template <typename Input3D>
    Real operator()(const Input3D &input, const Vec3d &point) const
    {
        Vec3d elongated_point = point;
        elongated_point[0] *= elongation_factor_; // Elongate along x-axis
        return input(elongated_point);
    }
};
} // namespace SPH

#endif // SDF_PRIMITIVE_H
