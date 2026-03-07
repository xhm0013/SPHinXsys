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
 * @file 	geometric_element.h
 * @brief tbd.
 * @author Xiangyu Hu
 */

#ifndef GEOMETRIC_ELEMENT_H
#define GEOMETRIC_ELEMENT_H

#include "base_data_type_package.h"
#include "signed_distance_primitive.h"

namespace SPH
{

class GeometricBox
{
  public:
    explicit GeometricBox(const Vecd &halfsize);
    ~GeometricBox() {};

    bool checkContain(const Vecd &probe_point)
    {
        bool is_contained = true;
        for (int i = 0; i != Dimensions; ++i)
        {
            if (ABS(probe_point[i]) > halfsize_[i]) // outside the box
            {
                is_contained = false;
                break;
            }
        }
        return is_contained;
    };

    Vecd findClosestPoint(const Vecd &probe_point);
    BoundingBoxd findBounds();

  protected:
    Vecd halfsize_;
};

class GeometricBall
{
  public:
    explicit GeometricBall(Real radius);
    ~GeometricBall() {};

    bool checkContain(const Vecd &probe_point);
    Vecd findClosestPoint(const Vecd &probe_point);
    BoundingBoxd findBounds();

  protected:
    Real radius_;
};

class GeometricCylinder
{
  public:
    explicit GeometricCylinder(Real radius, Real halflength);
    ~GeometricCylinder() {};

    bool checkContain(const Vecd &probe_point)
    {
        if (ABS(probe_point[0]) > halflength_)
            return false;
        return probe_point.tail(Dimensions - 1).norm() <= radius_;
    };

    Vecd findClosestPoint(const Vecd &probe_point);
    BoundingBoxd findBounds();

  protected:
    Real radius_;
    Real halflength_;
};

template <typename SDPrimitiveType>
class GeometricElementSD
{
    SDPrimitiveType sd_primitive_;

    Vecd probeNormalDirection(const Vecd &probe_point)
    {
        Real eps = 1e-4;
        Vecd normal_direction;
        for (int i = 0; i != Dimensions; ++i)
        {
            Vecd probe_point_plus = probe_point;
            probe_point_plus[i] += eps;
            Vecd probe_point_minus = probe_point;
            probe_point_minus[i] -= eps;
            normal_direction[i] = (sd_primitive_(probe_point_plus) - sd_primitive_(probe_point_minus)) / (2.0 * eps);
        }
        return normal_direction.normalized();
    }

  public:
    template <typename... Args>
    explicit GeometricElementSD(Args &&...args) : sd_primitive_(std::forward<Args>(args)...) {}
    ~GeometricElementSD() {};

    Real checkContain(const Vecd &probe_point)
    {
        return sd_primitive_(probe_point) > 0.0 ? sd_primitive_(probe_point) : -sd_primitive_(probe_point);
    }

    Vecd findClosestPoint(const Vecd &probe_point)
    {
        Real signed_distance = sd_primitive_(probe_point);
        Vecd normal_direction = probeNormalDirection(probe_point);
        return probe_point - signed_distance * normal_direction;
    }
};
} // namespace SPH
#endif // GEOMETRIC_ELEMENT_H
