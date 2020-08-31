//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "AQViewFactorRayStudy.h"

// MOOSE includes
#include "Conversion.h"
#include "TimedPrint.h"
#include "MooseUtils.h"

// Local includes
#include "AQViewFactorRayBC.h"

// libMesh includes
#include "libmesh/parallel_algebra.h"
#include "libmesh/parallel_sync.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/fe_base.h"
#include "libmesh/quadrature.h"

registerMooseObject("HeatConductionApp", AQViewFactorRayStudy);

void
AQViewFactorRayStudy::getLegHRQ(unsigned int order, std::vector<Real> & x, std::vector<Real> & w)
{
  switch (order)
  {
    case 2:
      x = {0.2113248654051872, 0.7886751345948129};
      w = {0.5, 0.5};
      break;
    case 4:
      x = {0.06943184420297366, 0.3300094782075719, 0.6699905217924281, 0.9305681557970262};
      w = {0.1739274225687273, 0.3260725774312733, 0.3260725774312729, 0.1739274225687271};
      break;
    case 6:
      x = {0.03376524289842392,
           0.1693953067668676,
           0.3806904069584015,
           0.6193095930415985,
           0.8306046932331324,
           0.9662347571015761};
      w = {0.08566224618958523,
           0.1803807865240692,
           0.2339569672863452,
           0.2339569672863455,
           0.1803807865240691,
           0.0856622461895854};
      break;
    case 8:
      x = {0.01985507175123175,
           0.1016667612931865,
           0.2372337950418355,
           0.4082826787521751,
           0.5917173212478248,
           0.7627662049581646,
           0.8983332387068135,
           0.9801449282487682};
      w = {0.050614268145188,
           0.1111905172266872,
           0.1568533229389437,
           0.1813418916891811,
           0.181341891689181,
           0.1568533229389435,
           0.1111905172266871,
           0.05061426814518794};
      break;
    case 10:
      x = {0.01304673574141424,
           0.06746831665550779,
           0.1602952158504878,
           0.2833023029353764,
           0.4255628305091844,
           0.5744371694908157,
           0.7166976970646236,
           0.8397047841495122,
           0.9325316833444921,
           0.9869532642585856};
      w = {0.03333567215434423,
           0.07472567457529,
           0.1095431812579912,
           0.134633359654998,
           0.1477621123573761,
           0.1477621123573765,
           0.1346333596549973,
           0.109543181257991,
           0.07472567457529056,
           0.03333567215434426};
      break;
    case 12:
      x = {0.009219682876640378,
           0.04794137181476255,
           0.1150486629028475,
           0.2063410228566914,
           0.3160842505009099,
           0.4373832957442654,
           0.5626167042557346,
           0.68391574949909,
           0.7936589771433087,
           0.8849513370971525,
           0.9520586281852376,
           0.9907803171233598};
      w = {0.02358766819325601,
           0.05346966299765912,
           0.08003916427167318,
           0.101583713361533,
           0.1167462682691775,
           0.1245735229067015,
           0.1245735229067015,
           0.1167462682691773,
           0.101583713361533,
           0.08003916427167339,
           0.05346966299765927,
           0.02358766819325585};
      break;
    case 16:
      x = {0.005299532504174975,
           0.02771248846338353,
           0.06718439880608396,
           0.1222977958224984,
           0.1910618777986781,
           0.2709916111713863,
           0.3591982246103707,
           0.4524937450811815,
           0.5475062549188188,
           0.6408017753896296,
           0.7290083888286136,
           0.808938122201322,
           0.8777022041775016,
           0.9328156011939159,
           0.9722875115366161,
           0.994700467495825};
      w = {0.01357622970587707,
           0.03112676196932375,
           0.04757925584124661,
           0.06231448562776682,
           0.07479799440828847,
           0.08457825969750156,
           0.09130170752246193,
           0.09472530522753399,
           0.09472530522753478,
           0.09130170752246126,
           0.08457825969750159,
           0.0747979944082882,
           0.06231448562776774,
           0.04757925584124636,
           0.03112676196932393,
           0.01357622970587707};
      break;
    case 20:
      x = {0.003435700407452502, 0.01801403636104321, 0.04388278587433714, 0.08044151408889055,
           0.1268340467699246,   0.1819731596367425,  0.2445664990245865,  0.3131469556422902,
           0.3861070744291775,   0.4617367394332513,  0.5382632605667487,  0.6138929255708228,
           0.6868530443577098,   0.7554335009754136,  0.8180268403632576,  0.8731659532300755,
           0.9195584859111094,   0.9561172141256626,  0.9819859636389568,  0.9965642995925473};
      w = {0.008807003569576043, 0.0203007149001937,  0.03133602416705427, 0.04163837078835239,
           0.0509650599086201,   0.05909726598075941, 0.06584431922458812, 0.0710480546591913,
           0.07458649323630168,  0.0763766935653629,  0.07637669356536299, 0.07458649323630215,
           0.07104805465919138,  0.06584431922458857, 0.05909726598075964, 0.05096505990861998,
           0.04163837078835197,  0.03133602416705467, 0.02030071490019304, 0.008807003569576187};
      break;
    case 30:
      x = {0.001553257962675192, 0.008165938360126357, 0.01998906751584634, 0.03689997628536273,
           0.05871973210397358,  0.08521711880861577,  0.1161112839475868,  0.1510747526033421,
           0.1897369085053784,   0.23168792592899,     0.2764831152309553,  0.3236476372345609,
           0.372681536916055,    0.4230650431957081,   0.474264078722341,   0.5257359212776591,
           0.5769349568042921,   0.627318463083945,    0.6763523627654391,  0.7235168847690446,
           0.7683120740710101,   0.8102630914946214,   0.8489252473966579,  0.8838887160524131,
           0.9147828811913843,   0.9412802678960264,   0.963100023714637,   0.9800109324841538,
           0.9918340616398736,   0.9984467420373246};
      w = {0.003984096248083338, 0.009233234155545484, 0.01439235394166194, 0.01939959628481315,
           0.02420133641529709,  0.0287465781088094,   0.03298711494109029, 0.03687798736885261,
           0.04037794761471006,  0.04344989360054158,  0.04606126111889319, 0.04818436858732199,
           0.04979671029339737,  0.05088119487420244,  0.05142632644677967, 0.05142632644677915,
           0.05088119487420278,  0.04979671029339732,  0.0481843685873224,  0.04606126111889252,
           0.04344989360054204,  0.0403779476147104,   0.03687798736885246, 0.0329871149410903,
           0.02874657810880945,  0.02420133641529774,  0.01939959628481316, 0.01439235394166161,
           0.009233234155545658, 0.003984096248083273};
      break;
    case 40:
      x = {0.0008811451447204299, 0.004636880650271624, 0.01137002500811307, 0.02104159039310433,
           0.03359359586066191,   0.04895059651556288,  0.06702024839387016, 0.08769388458334409,
           0.1108471742867403,    0.1363408724050365,   0.1640216576929103,  0.1937230551660097,
           0.2252664374524357,    0.2584620991569107,   0.2931103978141975,  0.3290029545871207,
           0.3659239074963729,    0.4036512096493143,   0.4419579646623724,  0.4806137912469748,
           0.5193862087530253,    0.5580420353376276,   0.5963487903506857,  0.6340760925036268,
           0.6709970454128792,    0.7068896021858027,   0.7415379008430895,  0.774733562547564,
           0.80627694483399,      0.8359783423070898,   0.8636591275949634,  0.8891528257132597,
           0.9123061154166558,    0.9329797516061296,   0.9510494034844372,  0.9664064041393384,
           0.978958409606896,     0.9886299749918872,   0.9953631193497287,  0.9991188548552798};
      w = {0.002260638549266528, 0.005249142265576765, 0.008210529190954059, 0.01112292459708328,
           0.01396850349001133,  0.01673009764127386,  0.01939108398723628,  0.02193545409283664,
           0.02434790381753616,  0.02661392349196839,  0.02871988454969542,  0.03065312124646457,
           0.03240200672830051,  0.03395602290761735,  0.03530582369564322,  0.03644329119790184,
           0.03736158452898407,  0.03805518095031343,  0.03851990908212326,  0.03875297398921271,
           0.03875297398921225,  0.03851990908212417,  0.0380551809503133,   0.03736158452898412,
           0.03644329119790179,  0.03530582369564304,  0.03395602290761734,  0.0324020067283003,
           0.0306531212464645,   0.02871988454969549,  0.02661392349196847,  0.02434790381753585,
           0.02193545409283717,  0.01939108398723586,  0.01673009764127338,  0.01396850349001175,
           0.01112292459708374,  0.008210529190954169, 0.005249142265576059, 0.002260638549266752};
      break;
    case 80:
      x = {0.000223088674184635, 0.001175067800881335, 0.002886229517155892, 0.005354348750122362,
           0.008575713630685655, 0.01254542970713612,  0.01725745547810054,  0.0227046168281827,
           0.02887861934506397,  0.03577006141377725,  0.04336844871412099,  0.05166221028061474,
           0.06063871616089328,  0.07028429666844438,  0.08058426320987233,  0.0915229306592682,
           0.1030836412476974,   0.115248789932479,    0.1279998512082015,   0.1413174073189502,
           0.1551811778289863,   0.16957005050694,     0.1844621134765639,   0.199834688585124,
           0.2156643659386451,   0.2319270395514338,   0.2485979440556073,   0.2656516924147276,
           0.2830623145841219,   0.3008032970590154,   0.3188476232502564,   0.3371678146261489,
           0.3557359725577439,   0.374523820803864,    0.3935027485711671,   0.412643854083677,
           0.431917988595428,    0.4512958007792078,   0.4707477814237899,   0.4902443083716032,
           0.5097556916283972,   0.5292522185762104,   0.5487041992207923,   0.568082011404572,
           0.5873561459163236,   0.6064972514288333,   0.6254761791961362,   0.644264027442256,
           0.6628321853738509,   0.6811523767497436,   0.6991967029409848,   0.716937685415878,
           0.7343483075852724,   0.7514020559443926,   0.7680729604485659,   0.7843356340613548,
           0.8001653114148759,   0.815537886523436,    0.8304299494930596,   0.8448188221710137,
           0.85868259268105,     0.8720001487917985,   0.8847512100675206,   0.8969163587523026,
           0.9084770693407318,   0.9194157367901274,   0.9297157033315555,   0.9393612838391068,
           0.9483377897193852,   0.9566315512858787,   0.9642299385862227,   0.9711213806549364,
           0.9772953831718172,   0.9827425445218994,   0.9874545702928637,   0.9914242863693143,
           0.9946456512498776,   0.9971137704828441,   0.9988249321991187,   0.9997769113258153};
      w = {
          0.0005724750015935829, 0.001331766794756083, 0.00209015656234747,  0.002845461225701455,
          0.003596452384058864,  0.004341972634630063, 0.005080883020551938, 0.005812057060399283,
          0.006534380796200345,  0.007246754020254332, 0.00794809179186279,  0.008637326028134728,
          0.009313407104149343,  0.009975305439070716, 0.01062201305789127,  0.01125254512316596,
          0.01186594143296515,   0.01246126788205799,  0.01303761788378269,  0.01359411375024307,
          0.01412990802863836,   0.01464418479163394,  0.01513616087977912,  0.01560508709405745,
          0.01605024933674352,   0.01647096969882277,  0.01686660749230598,  0.01723656022587697,
          0.01758026452237389,   0.0178971969767078,   0.01818687495291785,  0.01844885731913827,
          0.01868274511936515,   0.018888182181001,    0.01906485565723848,  0.01921249650347986,
          0.01933087988703829,   0.01941982552952643,  0.01947919798138537,  0.01950890682815304,
          0.01950890682815289,   0.01947919798138433,  0.01941982552952594,  0.01933087988703873,
          0.01921249650347963,   0.01906485565723888,  0.01888818218100091,  0.01868274511936495,
          0.0184488573191383,    0.01818687495291781,  0.01789719697670791,  0.01758026452237373,
          0.01723656022587688,   0.01686660749230561,  0.01647096969882261,  0.01605024933674383,
          0.01560508709405768,   0.01513616087977928,  0.01464418479163459,  0.01412990802863829,
          0.0135941137502433,    0.01303761788378273,  0.01246126788205774,  0.0118659414329652,
          0.01125254512316672,   0.01062201305789101,  0.009975305439071061, 0.009313407104150202,
          0.008637326028134436,  0.007948091791862205, 0.007246754020254351, 0.006534380796200811,
          0.005812057060399114,  0.0050808830205515,   0.004341972634630382, 0.003596452384058713,
          0.00284546122570141,   0.002090156562346913, 0.001331766794756246, 0.0005724750015935458};
      break;
    default:
      ::mooseError("Gauss-Legendre order ", order, " is not supported");
  }
}

DenseMatrix<Real>
AQViewFactorRayStudy::aqRoationMatrix(const Point & normal) const
{
  DenseMatrix<Real> rot(3, 3);
  if (_mesh.dimension() == 2)
  {
    mooseError("2D damnit");
  }
  else if (_mesh.dimension() == 3)
  {
    // Create a local coordinate system around normal
    Point tx(getOrthonormalVector(normal, _mesh.dimension()));
    Point ty = normal.cross(tx);

    // Create rotation matrix and rotate vector omega
    for (unsigned int j = 0; j < 3; ++j)
    {
      rot(j, 0) = tx(j);
      rot(j, 1) = ty(j);
      rot(j, 2) = normal(j);
    }
  }
  return rot;
}

void
AQViewFactorRayStudy::getLegChebHRQ(unsigned int cheb_order,
                                    unsigned int leg_order,
                                    std::vector<std::pair<Real, Real>> & x,
                                    std::vector<Real> & w)
{
  std::vector<Real> cheb_x(cheb_order);
  std::vector<Real> cheb_w(cheb_order);

  for (unsigned int j = 0; j < cheb_order; ++j)
  {
    unsigned int k = j + 1;
    cheb_x[j] = (2 * k - 1) * M_PI / cheb_order;
    cheb_w[j] = 2 * M_PI / cheb_order;
  }

  // Legendre quadrature
  std::vector<Real> gleg_hr_x;
  std::vector<Real> gleg_hr_w;

  getLegHRQ(leg_order, gleg_hr_x, gleg_hr_w);

  // product quadrature
  unsigned int prod_order = cheb_order * leg_order;
  x.resize(prod_order);
  w.resize(prod_order);

  unsigned int l = 0;
  for (unsigned int i = 0; i < cheb_order; ++i)
    for (unsigned int j = 0; j < leg_order; ++j)
    {
      x[l].first = cheb_x[i];
      x[l].second = gleg_hr_x[j];
      w[l] = gleg_hr_w[j] * cheb_w[i];
      ++l;
    }
}

Point
AQViewFactorRayStudy::getOrthonormalVector(const Point & v, unsigned int dim)
{
  if (MooseUtils::absoluteFuzzyLessEqual(v.norm(), 0))
    ::mooseError("Vector v has norm close to 0");
  if (dim != 2 && dim != 3)
    ::mooseError("Dimension must be 2 or 3 but is provided as ", dim);

  if (v(0) == 0)
    return Point(1, 0, 0);
  else if (v(1) == 0)
    return Point(0, 1, 0);
  else if (v(2) == 0 && dim == 3)
    return Point(0, 0, 1);

  Point t(-v(1), v(0), 0);
  return t / t.norm();
}

Point
AQViewFactorRayStudy::getAngularDirection(unsigned int l) const
{
  Point direction;
  if (_mesh.dimension() == 2)
    mooseError("2D damnit");
  else if (_mesh.dimension() == 3)
  {
    DenseVector<Real> omega(3);
    Real phi = _aq_angles[l].first;
    Real mu = _aq_angles[l].second;
    omega(0) = sqrt(1 - mu * mu) * cos(phi);
    omega(1) = sqrt(1 - mu * mu) * sin(phi);
    omega(2) = mu;
    DenseVector<Real> omega_p(3);
    _rotation_matrix.vector_mult(omega_p, omega);
    direction = Point(omega_p(0), omega_p(1), omega_p(2));
  }
  return direction;
}

InputParameters
AQViewFactorRayStudy::validParams()
{
  auto params = RayTracingStudy::validParams();

  params.addRequiredParam<std::vector<BoundaryName>>(
      "boundary", "The list of boundary IDs from the mesh where this boundary condition applies");

  MooseEnum qorders("CONSTANT FIRST SECOND THIRD FOURTH FIFTH SIXTH SEVENTH EIGHTH NINTH TENTH "
                    "ELEVENTH TWELFTH THIRTEENTH FOURTEENTH FIFTEENTH SIXTEENTH SEVENTEENTH "
                    "EIGHTTEENTH NINTEENTH TWENTIETH",
                    "CONSTANT");
  params.addParam<MooseEnum>("face_order", qorders, "The face quadrature rule order");

  MooseEnum qtypes("GAUSS GRID", "GRID");
  params.addParam<MooseEnum>("face_type", qtypes, "The face quadrature type");

  params.addParam<unsigned int>(
      "polar_quad_order",
      16,
      "Order of the polar quadrature [polar angle is between ray and normal]. Must be even.");
  params.addParam<unsigned int>(
      "azimuthal_quad_order",
      8,
      "Order of the azimuthal quadrature per quadrant [azimuthal angle is measured in "
      "a plan perpendicular to the normal].");

  MooseEnum convention("positive=0 negative=1", "positive");
  params.addParam<MooseEnum>(
      "internal_convention",
      convention,
      "The convention for spawning rays from internal sidesets; denotes the sign of the dot "
      "product between a ray and the internal sideset side normal");

  // Shouldn't ever need RayKernels for view factors
  params.set<bool>("ray_kernel_coverage_check") = false;
  params.suppressParameter<bool>("ray_kernel_coverage_check");

  // So that the study executes before the RayTracingViewFactor
  params.set<bool>("force_preaux") = true;
  params.suppressParameter<bool>("force_preaux");

  // Need to use internal sidesets
  params.set<bool>("use_internal_sidesets") = true;
  params.suppressParameter<bool>("use_internal_sidesets");

  // May have duplicate external boundaries on one boundary so ignore
  params.set<bool>("external_ray_bc_coverage_check") = false;
  params.suppressParameter<bool>("external_ray_bc_coverage_check");

  // No need to use Ray registration
  params.set<bool>("_use_ray_registration") = false;
  // Do not need to bank Rays on completion
  params.set<bool>("_bank_rays_on_completion") = false;

  return params;
}

AQViewFactorRayStudy::AQViewFactorRayStudy(const InputParameters & parameters)
  : RayTracingStudy(parameters),
    _bnd_ids(_mesh.getBoundaryIDs(getParam<std::vector<BoundaryName>>("boundary"))),
    _internal_convention(getParam<MooseEnum>("internal_convention")),
    _ray_index_start_bnd_id(registerRayAuxData("start_bnd_id")),
    _ray_index_start_total_weight(registerRayAuxData("start_total_weight")),
    _fe_face(FEBase::build(_mesh.dimension(), FEType(CONSTANT, MONOMIAL))),
    _q_face(QBase::build(Moose::stringToEnum<QuadratureType>(getParam<MooseEnum>("face_type")),
                         _mesh.dimension() - 1,
                         Moose::stringToEnum<Order>(getParam<MooseEnum>("face_order")))),
    _threaded_vf_info(libMesh::n_threads()),
    _threaded_cached_ray_bcs(libMesh::n_threads())
{
  _fe_face->attach_quadrature_rule(_q_face.get());
  _fe_face->get_normals();
  _fe_face->get_xyz();

  // create angular quadrature
  if (_mesh.dimension() == 2)
    mooseError("2D damnit");
  else if (_mesh.dimension() == 3)
    getLegChebHRQ(4 * getParam<unsigned int>("azimuthal_quad_order"),
                  getParam<unsigned int>("polar_quad_order"),
                  _aq_angles,
                  _aq_weights);
  _naq = _aq_weights.size();
}

void
AQViewFactorRayStudy::initialSetup()
{
  RayTracingStudy::initialSetup();

  // We optimized away RayKernels, so don't allow them
  std::vector<RayKernelBase *> ray_kernels;
  RayTracingStudy::getRayKernels(ray_kernels, 0);
  if (!ray_kernels.empty())
    mooseError("The AQViewFactorRayStudy '", name(), "' is not compatible with RayKernels");

  // Make sure we have only one RayBC: a ViewFactorRayBC with the same boundaries
  std::vector<RayBC *> ray_bcs;
  RayTracingStudy::getRayBCs(ray_bcs, 0);
  if (ray_bcs.size() != 1)
    mooseError("The AQViewFactorRayStudy '",
               name(),
               "' requires one and only one RayBC, a ViewFactorRayBC");
  AQViewFactorRayBC * vf_ray_bc = dynamic_cast<AQViewFactorRayBC *>(ray_bcs[0]);
  if (!vf_ray_bc)
    mooseError("The AQViewFactorRayStudy '",
               name(),
               "' requires one and only one RayBC, a ViewFactorRayBC");
  if (!vf_ray_bc->hasBoundary(_bnd_ids))
    mooseError("The AQViewFactorRayBC '",
               vf_ray_bc->name(),
               "' must be applied to the same boundaries as the AQViewFactorRayStudy '",
               name(),
               "'");

  // Cache our one RayBC per thread so that we don't spend unnecessary time querying for RayBCs on
  // every boundary
  for (THREAD_ID tid = 0; tid < libMesh::n_threads(); ++tid)
  {
    RayTracingStudy::getRayBCs(ray_bcs, tid);
    _threaded_cached_ray_bcs[tid] = ray_bcs;
  }

  // find domain length scale == length of diagonal across domain
  _domain_length_scale = 0.0;
  for (unsigned int j = 0; j < _mesh.dimension(); ++j)
  {
    Real d = _mesh.getMaxInDimension(j) - _mesh.getMinInDimension(j);
    _domain_length_scale += d * d;
  }
  _domain_length_scale = std::sqrt(_domain_length_scale);
}

void
AQViewFactorRayStudy::preExecuteStudy()
{
  RayTracingStudy::preExecuteStudy();

  // Clear and zero the view factor maps we're about to accumulate into for each thread
  for (THREAD_ID tid = 0; tid < libMesh::n_threads(); ++tid)
  {
    _threaded_vf_info[tid].clear();
    for (const BoundaryID from_id : _bnd_ids)
      for (const BoundaryID to_id : _bnd_ids)
        _threaded_vf_info[tid][from_id][to_id] = 0;
  }

  generatePoints();

  generateRayIDs();
}

void
AQViewFactorRayStudy::postExecuteStudy()
{
  RayTracingStudy::postExecuteStudy();

  // Finalize the cumulative _vf_info];
  _vf_info.clear();
  for (const BoundaryID from_id : _bnd_ids)
    for (const BoundaryID to_id : _bnd_ids)
    {
      Real & entry = _vf_info[from_id][to_id];

      // Zero before summing
      entry = 0;

      // Sum over threads
      for (THREAD_ID tid = 0; tid < libMesh::n_threads(); ++tid)
        entry += _threaded_vf_info[tid][from_id][to_id];

      // Sum over processors
      _communicator.sum(entry);

      _vf_info[from_id][to_id] = entry;
    }
}

void
AQViewFactorRayStudy::addToViewFactorInfo(Real value,
                                          const BoundaryID from_id,
                                          const BoundaryID to_id,
                                          const THREAD_ID tid)
{
  if (!currentlyPropagating() && !currentlyGenerating())
    mooseError("addToViewFactorInfo() can only be called during Ray generation and propagation");
  mooseAssert(_threaded_vf_info[tid].count(from_id),
              "Threaded view factor info does not have from boundary");
  mooseAssert(_threaded_vf_info[tid][from_id].count(to_id),
              "Threaded view factor info does not have from -> to boundary");

  _threaded_vf_info[tid][from_id][to_id] += value;
}

Real
AQViewFactorRayStudy::viewFactorInfo(const BoundaryID from_id, const BoundaryID to_id) const
{
  auto it = _vf_info.find(from_id);
  if (it == _vf_info.end())
    mooseError("From boundary id ", from_id, " not in view factor map.");

  auto itt = it->second.find(to_id);
  if (itt == it->second.end())
    mooseError("From boundary id ", from_id, " to boundary_id ", to_id, " not in view factor map.");
  return itt->second;
}

AQViewFactorRayStudy &
AQViewFactorRayStudy::castFromStudy(const InputParameters & params)
{
  RayTracingStudy * study = params.getCheckedPointerParam<RayTracingStudy *>("_ray_tracing_study");
  if (!study)
    ::mooseError("Fails earlier");
  AQViewFactorRayStudy * vf_study = dynamic_cast<AQViewFactorRayStudy *>(study);
  if (!vf_study)
    ::mooseError(params.get<std::string>("_type"),
                 " '",
                 params.get<std::string>("_object_name"),
                 "' must be paired with a AQViewFactorRayStudy");
  return *vf_study;
}

void
AQViewFactorRayStudy::generatePoints()
{
  _internal_bnd_ids.clear();

  const auto & normals = _fe_face->get_normals();
  const auto & points = _fe_face->get_xyz();
  const auto & weights = _fe_face->get_JxW();

  _start_info.clear();

  // Starting elements we have that are on the wrong side of an internal boundary
  std::vector<StartElem> send_start_info;

  // Get all possible points on the user defined boundaries on this proc
  for (const BndElement * belem : *_mesh.getBoundaryElementRange())
  {
    const Elem * elem = belem->_elem;
    const auto side = belem->_side;
    const auto bnd_id = belem->_bnd_id;

    // Skip if we don't own you
    if (elem->processor_id() != processor_id())
      continue;

    // Skip if the boundary id isn't one we're looking for
    if (std::find(_bnd_ids.begin(), _bnd_ids.end(), bnd_id) == _bnd_ids.end())
      continue;

    // Sanity check on QGRID not working on some types
    if (_q_face->type() == QGRID && elem->type() == TET4)
      mooseError(
          "Cannot use GRID quadrature type with tetrahedral elements in AQViewFactorRayStudy '",
          _name,
          "'");

    // The elem/side that we will actually start the trace from
    // (this may change on internal sidesets)
    const Elem * start_elem = elem;
    auto start_side = side;

    // Reinit this face
    _fe_face->reinit(elem, side);
    const auto & normal = normals[0];

    // See if this boundary is internal
    const Elem * neighbor = elem->neighbor_ptr(side);
    if (neighbor)
    {
      if (!neighbor->active())
        mooseError("AQViewFactorRayStudy does not work with adaptivity");

      // Mark if this is an internal sideset for future use
      _internal_bnd_ids.insert(bnd_id);

      // With the positive convention, the Rays that we want to spawn from internal boundaries
      // have positive dot products with the outward normal on the side. This is not efficient
      // for the ray tracer, so set them up to be spawned from the other element on this
      // side instead
      if (_internal_convention == 0)
      {
        start_elem = neighbor;
        start_side = neighbor->which_neighbor_am_i(elem);
      }
    }

    // If the start elem is local, add to _start_info, otherwise add to send_start_info
    // so that we can ship it to who will actually start it
    auto & add_start_info =
        start_elem->processor_id() == processor_id() ? _start_info : send_start_info;
    add_start_info.emplace_back(
        elem, side, start_elem, start_side, bnd_id, normal, points, weights);
  }

  // Decide on the internal boundaries
  _communicator.set_union(_internal_bnd_ids);

  // Send the additional start info if any
  if (_internal_convention == 0)
  {
    // Fill the elements to send to each proc based on who owns the elem
    std::unordered_map<processor_id_type, std::vector<StartElem *>> send_start_map;
    for (StartElem & start_elem : send_start_info)
      send_start_map[start_elem._elem->processor_id()].push_back(&start_elem);

    // Communicate and fill into our _start_info
    auto start_functor = [this](processor_id_type, const std::vector<StartElem *> & start_elems) {
      _start_info.reserve(_start_info.size() + start_elems.size());
      for (StartElem * start_elem : start_elems)
        _start_info.emplace_back(std::move(*start_elem));
    };
    Parallel::push_parallel_packed_range(_communicator, send_start_map, this, start_functor);
  }

  // Decide on the processors that have boundary elements so that we know who to send to
  std::vector<dof_id_type> have_boundary(_communicator.size());
  _communicator.allgather((dof_id_type)_start_info.size(), have_boundary);
}

void
AQViewFactorRayStudy::generateRayIDs()
{
  std::size_t num_local_rays = 0;
  std::size_t num_local_start_points = 0;
  for (const auto & start_elem : _start_info)
    for (unsigned int j = 0; j < start_elem._points.size(); ++j)
    {
      num_local_start_points += 1;
      num_local_rays += _naq;
    }

  // Decide on the starting Ray IDs for each rank. Each rank will have a contiguous set of IDs for
  // its rays and there will be no holes in the IDs between ranks Number of rays this proc has
  // Vector of sizes across all processors (for rank 0 only)
  std::vector<std::size_t> proc_sizes;
  // Send this procesor's local size to rank 0
  _comm.gather(0, num_local_rays, proc_sizes);
  // Rank 0 has the proc sizes and will decide on a starting ID for every rank
  std::vector<RayID> proc_starting_id;
  if (processor_id() == 0)
  {
    proc_starting_id.resize(_comm.size());
    RayID current_starting_id = 0;
    for (processor_id_type pid = 0; pid < _comm.size(); ++pid)
    {
      proc_starting_id[pid] = current_starting_id;
      current_starting_id += proc_sizes[pid];
    }
  }

  // Send the starting ID to each rank
  _comm.scatter(proc_starting_id, _current_starting_id, 0);
  if (!num_local_rays)
    _current_starting_id = Ray::INVALID_RAY_ID;

  // Print out totals while we're here
  std::size_t num_total_points = num_local_start_points;
  std::size_t num_total_rays = num_local_rays;
  _communicator.sum(num_total_points);
  _communicator.sum(num_total_rays);
  _console << "View factor study generated " << num_total_points
           << " points with an angular quadrature of " << _naq << " directions per point requiring "
           << num_total_rays << " rays" << std::endl;
}

void
AQViewFactorRayStudy::generateRays()
{
  CONSOLE_TIMED_PRINT("View factor study generating rays");

  // loop through all starting points and spawn rays from each
  for (const auto & start_elem : _start_info)
  {
    const bool start_is_external = start_elem._elem->neighbor_ptr(start_elem._side) == nullptr;
    auto normal = start_elem._normal;

    // if the boundary is external, flip the normal;
    if (start_is_external && _internal_convention == 0)
      normal *= -1;

    // Create rotation matrix and rotate vector omega
    _rotation_matrix = aqRoationMatrix(normal);

    for (std::size_t start_i = 0; start_i < start_elem._points.size(); ++start_i)
    {
      const auto & start_point = start_elem._points[start_i];

      // loop through all directions in quadrature set
      for (unsigned int l = 0; l < _naq; ++l)
      {
        const auto direction = getAngularDirection(l);
        const auto start_weight =
            start_elem._weights[start_i] * _aq_weights[l] * _aq_angles[l].second;
        const Point mock_end_point = start_point + _domain_length_scale * direction;

        // Create a Ray and add it to the buffer for future tracing
        std::shared_ptr<Ray> ray = _ray_pool.acquire(start_point,
                                                     mock_end_point,
                                                     rayDataSize(),
                                                     rayAuxDataSize(),
                                                     start_elem._start_elem,
                                                     start_elem._start_side);
        addToWorkingBuffer(ray);
        ray->setID(_current_starting_id++);
        ray->setAuxData(_ray_index_start_bnd_id, start_elem._bnd_id);
        ray->setAuxData(_ray_index_start_total_weight, start_weight);

        ray.reset();

        chunkyTraceAndBuffer();
      }
    }
  }
}

namespace libMesh
{
namespace Parallel
{

unsigned int
Packing<AQViewFactorRayStudy::StartElem *>::packed_size(
    typename std::vector<Real>::const_iterator in)
{
  unsigned int total_size = 0;

  // Number of points
  unsigned int num_points = *in++;

  // Number of points, elem, side, start_elem, start_side, bnd_id, normal
  total_size += 9;
  // Points and weights
  total_size += num_points * 4;

  return total_size;
}

unsigned int
Packing<AQViewFactorRayStudy::StartElem *>::packable_size(
    const AQViewFactorRayStudy::StartElem * const start_elem, const void *)
{
  unsigned int total_size = 0;

  // Number of points, elem, side, start_elem, start_side, bnd_id, normal
  total_size += 9;
  // Points and weights
  total_size += start_elem->_points.size() * 4;

  return total_size;
}

template <>
AQViewFactorRayStudy::StartElem *
Packing<AQViewFactorRayStudy::StartElem *>::unpack(std::vector<Real>::const_iterator in,
                                                   AQViewFactorRayStudy * study)
{
  // Number of points
  const std::size_t num_points = *in++;

  // Elem id
  const dof_id_type elem_id = *in++;
  const Elem * elem =
      elem_id != DofObject::invalid_id ? study->meshBase().query_elem_ptr(elem_id) : nullptr;

  // Side
  const unsigned short side = *in++;

  // Start elem id
  const dof_id_type start_elem_id = *in++;
  const Elem * start_elem = start_elem_id != DofObject::invalid_id
                                ? study->meshBase().query_elem_ptr(start_elem_id)
                                : nullptr;

  // Start side
  const unsigned short start_side = *in++;

  // Boundary ID
  const BoundaryID bnd_id = *in++;

  // Normal
  Point normal;
  normal(0) = *in++;
  normal(1) = *in++;
  normal(2) = *in++;

  // Points
  std::vector<Point> points(num_points);
  for (std::size_t i = 0; i < num_points; ++i)
  {
    points[i](0) = *in++;
    points[i](1) = *in++;
    points[i](2) = *in++;
  }

  // Weights
  std::vector<Real> weights(num_points);
  for (std::size_t i = 0; i < num_points; ++i)
    weights[i] = *in++;

  // We store this with std::move, so this is safe if we do it right!
  return new AQViewFactorRayStudy::StartElem(elem,
                                             side,
                                             start_elem,
                                             start_side,
                                             bnd_id,
                                             std::move(normal),
                                             std::move(points),
                                             std::move(weights));
}

template <>
void
Packing<AQViewFactorRayStudy::StartElem *>::pack(
    const AQViewFactorRayStudy::StartElem * const start_elem,
    std::back_insert_iterator<std::vector<Real>> data_out,
    const AQViewFactorRayStudy *)
{
  // Number of points
  data_out = start_elem->_points.size();

  // Elem id
  data_out = start_elem->_elem->id();

  // Side
  data_out = start_elem->_side;

  // Start elem id
  data_out = start_elem->_start_elem->id();

  // Start side
  data_out = start_elem->_start_side;

  // Boundary id
  data_out = start_elem->_bnd_id;

  // Normal
  data_out = start_elem->_normal(0);
  data_out = start_elem->_normal(1);
  data_out = start_elem->_normal(2);

  // Points
  for (const auto & point : start_elem->_points)
  {
    data_out = point(0);
    data_out = point(1);
    data_out = point(2);
  }

  // Weights
  std::copy(start_elem->_weights.begin(), start_elem->_weights.end(), data_out);
}

} // namespace Parallel

} // namespace libMesh
