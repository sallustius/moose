//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

class FVFaceInterface
{
public:
  /// This codifies a set of available ways to interpolate with elem+neighbor
  /// solution information to calculate values (e.g. solution, material
  /// properties, etc.) at the face (centroid).  These methods are used in the
  /// class's interpolate functions.  Some interpolation methods are only meant
  /// to be used with advective terms (e.g. upwind), others are more
  /// generically applicable.
  enum class InterpMethod
  {
    /// (elem+neighbor)/2
    Average,
    /// weighted
    Upwind,
    // Rhie-Chow
    RhieChow
  };

  FVFaceInterface() = default;

protected:
  /// Provides interpolation of face values for non-advection-specific purposes
  /// (although it can/will still be used by advective kernels sometimes).  The
  /// interpolated value is stored in result.  This should be called when a
  /// face value needs to be computed using elem and neighbor information (e.g. a
  /// material property, solution value, etc.).  elem and neighbor represent the
  /// property/value to compute the face value for.
  template <typename T, typename T2, typename T3>
  void interpolate(InterpMethod m, T & result, const T2 & elem, const T3 & neighbor)
  {
    switch (m)
    {
      case InterpMethod::Average:
        result = (elem + neighbor) * 0.5;
        break;
      default:
        mooseError("unsupported interpolation method for FVFaceInterface::interpolate");
    }
  }

  /// Provides interpolation of face values for advective flux kernels.  This
  /// should be called by advective kernels when a u_face value is needed from
  /// u_elem and u_neighbor.  The interpolated value is stored in result.  elem
  /// and neighbor represent the property/value being advected in the elem and
  /// neighbor elements respectively.  advector represents the vector quantity at
  /// the face that is doing the advecting (e.g. the flow velocity at the
  /// face); this value often will have been computed using a call to the
  /// non-advective interpolate function.
  template <typename T, typename T2, typename T3, typename Vector>
  void interpolate(
      InterpMethod m, T & result, const T2 & elem, const T3 & neighbor, const Vector & advector)
  {
    switch (m)
    {
      case InterpMethod::Average:
        result = (elem + neighbor) * 0.5;
        break;
      case InterpMethod::Upwind:
        if (advector * normal() > 0)
          result = elem;
        else
          result = neighbor;
        break;
      default:
        mooseError("unsupported interpolation method for FVFaceInterface::interpolate");
    }
  }

  /**
   * Return the normal at the face
   */
  virtual const ADRealVectorValue & normal() const = 0;
};
