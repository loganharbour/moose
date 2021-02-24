//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "RayTracingStudy.h"
#include "CoupleableMooseVariableDependencyIntermediateInterface.h"

// Local includes
#include "ClaimRays.h"

class MPMRayStudyBase : public RayTracingStudy,
                        public CoupleableMooseVariableDependencyIntermediateInterface
{
public:
  MPMRayStudyBase(const InputParameters & parameters);

  static InputParameters validParams();

  void onCompleteRay(const std::shared_ptr<Ray> & ray) override;

protected:
  /**
   * Pure virtual to be overridden by derived classes to fill the initial particles.
   *
   * Fill the initial particles (rays) into _initial_particles.
   */
  virtual void generateInitialParticles() = 0;

  /**
   * Generates an initial particle at the given point. To be called within derived
   * classes during generateInitialParticles().
   *
   * Also returns a reference to the generated Ray, for the purpose of setting data.
   */
  const std::shared_ptr<Ray> & generateInitialParticle(const Point & point,
                                                       const Elem * elem = nullptr);

  /// Index into the Ray aux data for the particle mass
  const RayDataIndex _particle_mass_index;

  /// Whether or not the initial particles generated are replicated
  const bool _initial_particles_replicated;

private:
  /**
   * Generates and claims the initial particles (if needed), and moves the particles (rays)
   * from the bank to be executed with the correct displacement
   */
  void generateRays() override;

  /**
   * Swaps the particle banks and transfers the particle masses to the mass field
   */
  void postExecuteStudy() override;

  /// The displacement variables - 1 for each dimension
  const std::vector<const VariableValue *> _displacements;
  /// The mass variable
  const VariableValue & _mass;

  /// The object used to claim Rays
  ClaimRays _claim_rays;

  /// The mapping from element -> rays (particles)
  std::unique_ptr<std::unordered_map<const Elem *, std::vector<std::shared_ptr<Ray>>>> _particles;
  /// The old mapping from element -> rays (particles)
  std::unique_ptr<std::unordered_map<const Elem *, std::vector<std::shared_ptr<Ray>>>>
      _old_particles;

  /// Filler for the initial particles; derived classes will fill this within generateInitialParticles()
  std::vector<std::shared_ptr<Ray>> _initial_particles;
  /// Filler for the initial rays that are claimed
  std::vector<std::shared_ptr<Ray>> _initial_local_particles;

  /// The MOOSE variables that contain the displacements: used for setting active vars in reinit
  std::set<MooseVariableFEBase *> _displacement_vars;
  /// The MOOSE variable that contains the mass: used for setting active vars in reinit
  MooseVariableFEBase * _mass_var;

  /// Whether or not the initial particles have been generated
  bool _generated_initial_particles;
  /// Whether or not we are currently generating initial particles
  bool _generating_initial_particles;
};
