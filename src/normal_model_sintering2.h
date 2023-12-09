// /* ----------------------------------------------------------------------
//     This is the

//     ██╗     ██╗ ██████╗  ██████╗  ██████╗ ██╗  ██╗████████╗███████╗
//     ██║     ██║██╔════╝ ██╔════╝ ██╔════╝ ██║  ██║╚══██╔══╝██╔════╝
//     ██║     ██║██║  ███╗██║  ███╗██║  ███╗███████║   ██║   ███████╗
//     ██║     ██║██║   ██║██║   ██║██║   ██║██╔══██║   ██║   ╚════██║
//     ███████╗██║╚██████╔╝╚██████╔╝╚██████╔╝██║  ██║   ██║   ███████║
//     ╚══════╝╚═╝ ╚═════╝  ╚═════╝  ╚═════╝ ╚═╝  ╚═╝   ╚═╝   ╚══════╝®

//     DEM simulation engine, released by
//     DCS Computing Gmbh, Linz, Austria
//     http://www.dcs-computing.com, office@dcs-computing.com

//     LIGGGHTS® is part of CFDEM®project:
//     http://www.liggghts.com | http://www.cfdem.com

//     Core developer and main author:
//     Christoph Kloss, christoph.kloss@dcs-computing.com

//     LIGGGHTS® is open-source, distributed under the terms of the GNU Public
//     License, version 2 or later. It is distributed in the hope that it will
//     be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
//     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. You should have
//     received a copy of the GNU General Public License along with LIGGGHTS®.
//     If not, see http://www.gnu.org/licenses . See also top-level README
//     and LICENSE files.

//     LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
//     the producer of the LIGGGHTS® software and the CFDEM®coupling software
//     See http://www.cfdem.com/terms-trademark-policy for details.

// -------------------------------------------------------------------------
//     Contributing author and copyright for this file:

//     Christoph Kloss (DCS Computing GmbH, Linz)
//     Christoph Kloss (JKU Linz)
//     Richard Berger (JKU Linz)

//     Copyright 2012-     DCS Computing GmbH, Linz
//     Copyright 2009-2012 JKU Linz
// ------------------------------------------------------------------------- */

// #ifdef NORMAL_MODEL
// NORMAL_MODEL(SINTERING,sintering,13)
// #else
// #ifndef NORMAL_MODEL_SINTERING_H_
// #define NORMAL_MODEL_SINTERING_H_
// #include "global_properties.h"
// #include "fix_property_atom.h"
// #include <cmath>
// #include "normal_model_base.h"
// #include "fix_mesh_surface.h"

// namespace LIGGGHTS {

// namespace ContactModels
// {
//   class ContactModelBase;

//   template<>
//   class NormalModel<SINTERING> : public NormalModelBase
//   {
//   public:
//     NormalModel(LAMMPS * lmp, IContactHistorySetup* hsetup, class ContactModelBase *c) :
//       NormalModelBase(lmp, hsetup, c),
//       displayedSettings(false),
//       cmb(c)
//     {

//     }

//     void registerSettings(Settings & settings)
//     {
//       //TODO error->one(FLERR,"TODO here also check if right surface model used");
//     }

//     inline void postSettings(IContactHistorySetup * hsetup, ContactModelBase *cmb)
//     {
        
  
//     }

//     void connectToProperties(PropertyRegistry & registry)
//     {

//     }

//     // effective exponent for stress-strain relationship
    
//     inline double stressStrainExponent()
//     {
//       return 1.5;
//     }

//     inline void surfacesIntersect(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces)
//     {
//       if (sidata.contact_flags) *sidata.contact_flags |= CONTACT_NORMAL_MODEL;

//       const double radi = sidata.radi;
//       const double radj = sidata.radj;
//       double reff = sidata.is_wall ? radi : (radi*radj/(radi+radj)); 
//       const double meff=sidata.meff;

//       if(sidata.deltan < 0)
//         error->one(FLERR, "sidata.deltan < 0!");
//       const double sqrtval = sqrt(reff*sidata.deltan);
//       #ifdef LIGGGHTS_DEBUG
//         if(std::isnan(sqrtval))
//           error->one(FLERR, "sqrtval is NaN!");
//       #endif

//       const double deltan = sidata.deltan;
//       double a_s = 4.0*reff*deltan;   // a_s dopočítáno
//       const double sqrtFiveOverSix = 0.91287092917527685576161630466800355658790782499663875;
//       const double surfaceEnergy = 1.0;
//       const double difusionParam = 0.000001;

//       // relative translational velocity
//       const double vr1 = sidata.v_i[0] - sidata.v_j[0];
//       const double vr2 = sidata.v_i[1] - sidata.v_j[1];
//       const double vr3 = sidata.v_i[2] - sidata.v_j[2];

//       // normal component
//       const double vn = vr1 * sidata.en[0] + vr2 * sidata.en[1] + vr3 * sidata.en[2];
//       const double vn1 = vn * sidata.en[0];
//       const double vn2 = vn * sidata.en[1];
//       const double vn3 = vn * sidata.en[2];
      
//       // if vr is not equal sidata.vn
//       // if (vn != sidata.vn) {
//       //   printf("vr = %e, sidata.vn = %e\n", vn, sidata.vn);
//       //   error->one(FLERR, "vr != sidata.vn!");
//       // }
      
//       if(!displayedSettings)
//       {
//         displayedSettings = true;
//       }

//       const double Fn_sintering = 1.125*2.0*M_PI*reff*surfaceEnergy;
//       const double Fn_viscous = vn * M_PI* pow(a_s,2) / 8.0 / difusionParam;
//       double Fn = Fn_sintering - Fn_viscous;

//       sidata.Fn = Fn;

//       // apply normal force
//       if(sidata.is_wall) {
//         const double Fn_sintering = Fn_sintering * sidata.area_ratio;
//         const double Fn_viscous = Fn_viscous * sidata.area_ratio;
//         i_forces.delta_F[0] += Fn_sintering * sidata.en[0] - Fn_viscous * vn1;
//         i_forces.delta_F[1] += Fn_sintering * sidata.en[1] - Fn_viscous * vn2;
//         i_forces.delta_F[2] += Fn_sintering * sidata.en[2] - Fn_viscous * vn3;
//       } else {
//         i_forces.delta_F[0] += Fn_sintering * sidata.en[0] - Fn_viscous * vn1;
//         i_forces.delta_F[1] += Fn_sintering * sidata.en[1] - Fn_viscous * vn2;
//         i_forces.delta_F[2] += Fn_sintering * sidata.en[2] - Fn_viscous * vn3;

//         j_forces.delta_F[0] += -i_forces.delta_F[0];
//         j_forces.delta_F[1] += -i_forces.delta_F[1];
//         j_forces.delta_F[2] += -i_forces.delta_F[2];
//       }
      
//     }

//     void surfacesClose(SurfacesCloseData &scdata, ForceData&, ForceData&)
//     {
//         if (scdata.contact_flags)
//             *scdata.contact_flags |= CONTACT_NORMAL_MODEL;
//     }

//     void beginPass(SurfacesIntersectData&, ForceData&, ForceData&){}
//     void endPass(SurfacesIntersectData&, ForceData&, ForceData&){}

//   protected:
//     bool displayedSettings;
//     FixPropertyAtom *fix_dissipated_;
//     class ContactModelBase *cmb;
//   };

// }

// }
// #endif
// #endif
