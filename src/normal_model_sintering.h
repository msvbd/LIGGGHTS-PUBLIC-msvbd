/* ----------------------------------------------------------------------
    This is the

    ██╗     ██╗ ██████╗  ██████╗  ██████╗ ██╗  ██╗████████╗███████╗
    ██║     ██║██╔════╝ ██╔════╝ ██╔════╝ ██║  ██║╚══██╔══╝██╔════╝
    ██║     ██║██║  ███╗██║  ███╗██║  ███╗███████║   ██║   ███████╗
    ██║     ██║██║   ██║██║   ██║██║   ██║██╔══██║   ██║   ╚════██║
    ███████╗██║╚██████╔╝╚██████╔╝╚██████╔╝██║  ██║   ██║   ███████║
    ╚══════╝╚═╝ ╚═════╝  ╚═════╝  ╚═════╝ ╚═╝  ╚═╝   ╚═╝   ╚══════╝®

    DEM simulation engine, released by
    DCS Computing Gmbh, Linz, Austria
    http://www.dcs-computing.com, office@dcs-computing.com

    LIGGGHTS® is part of CFDEM®project:
    http://www.liggghts.com | http://www.cfdem.com

    Core developer and main author:
    Christoph Kloss, christoph.kloss@dcs-computing.com

    LIGGGHTS® is open-source, distributed under the terms of the GNU Public
    License, version 2 or later. It is distributed in the hope that it will
    be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. You should have
    received a copy of the GNU General Public License along with LIGGGHTS®.
    If not, see http://www.gnu.org/licenses . See also top-level README
    and LICENSE files.

    LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
    the producer of the LIGGGHTS® software and the CFDEM®coupling software
    See http://www.cfdem.com/terms-trademark-policy for details.

-------------------------------------------------------------------------
    Contributing author and copyright for this file:

    Christoph Kloss (DCS Computing GmbH, Linz)
    Christoph Kloss (JKU Linz)
    Richard Berger (JKU Linz)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#ifdef NORMAL_MODEL
NORMAL_MODEL(SINTERING,sintering,13)
#else
#ifndef NORMAL_MODEL_SINTERING_H_
#define NORMAL_MODEL_SINTERING_H_
#include "global_properties.h"
#include "fix_property_atom.h"
#include <cmath>
#include "normal_model_base.h"
#include "fix_mesh_surface.h"

namespace LIGGGHTS {

namespace ContactModels
{
  class ContactModelBase;

  template<>
  class NormalModel<SINTERING> : public NormalModelBase
  {
  public:
    NormalModel(LAMMPS * lmp, IContactHistorySetup* hsetup, class ContactModelBase *c) :
      NormalModelBase(lmp, hsetup, c),
      Yeff(NULL),
      Geff(NULL),
      betaeff(NULL),
      limitForce(false),
      displayedSettings(false),
      heating(false),
      heating_track(false),
      elastic_potential_offset_(-1),
      elasticpotflag_(false),
      fix_dissipated_(NULL),
      dissipatedflag_(false),
      overlap_offset_(0.0),
      disable_when_bonded_(false),
      bond_history_offset_(-1),
      dissipation_history_offset_(-1),
      cmb(c)
    {
      
    }

    void registerSettings(Settings & settings)
    {
      

      //TODO error->one(FLERR,"TODO here also check if right surface model used");
    }

    inline void postSettings(IContactHistorySetup * hsetup, ContactModelBase *cmb)
    {
        
    }

    void connectToProperties(PropertyRegistry & registry)
    {
      
      
    }

    // effective exponent for stress-strain relationship
    
    inline double stressStrainExponent()
    {
      return 1.5;
    }

    void dissipateElasticPotential(SurfacesCloseData &scdata)
    {
        
    }

    inline void surfacesIntersect(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces)
    {
      
    }

    void surfacesClose(SurfacesCloseData &scdata, ForceData&, ForceData&)
    {
        
    }

    void beginPass(SurfacesIntersectData&, ForceData&, ForceData&){}
    void endPass(SurfacesIntersectData&, ForceData&, ForceData&){}

  protected:
    double ** Yeff;
    double ** Geff;
    double ** betaeff;

    bool tangential_damping;
    bool limitForce;
    bool displayedSettings;
    bool heating;
    bool heating_track;
    int elastic_potential_offset_;
    bool elasticpotflag_;
    FixPropertyAtom *fix_dissipated_;
    bool dissipatedflag_;
    int overlap_offset_;
    bool disable_when_bonded_;
    int bond_history_offset_;
    int dissipation_history_offset_;
    class ContactModelBase *cmb;

  };

}

}
#endif
#endif
