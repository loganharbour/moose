#pragma once

#include "RayTracingStudy.h"

class BasicPICRayStudy : public RayTracingStudy
{
public:
  BasicPICRayStudy(const InputParameters & parameters);

  static InputParameters validParams();

  virtual void generateRays() override;

  virtual void postOnSegment(const THREAD_ID tid, const std::shared_ptr<Ray> & ray) override;

private:
  /**
   * @return The max distance the given ray should travel at the current time
   *
   * Uses the current ray position to sample the velocity field in \p velocity_function
   */
  Real maxDistance(const Ray & ray) const;

  /// The starting points
  const std::vector<Point> & _start_points;

  /// The starting directions
  const std::vector<Point> & _start_directions;

  /// The function that represents the velocity field
  const Function & _velocity_function;

  /// Whether or not we've generated rays yet
  bool _has_generated;

  /// The banked rays to be used on the next timestep
  std::vector<std::shared_ptr<Ray>> _banked_rays;
};
