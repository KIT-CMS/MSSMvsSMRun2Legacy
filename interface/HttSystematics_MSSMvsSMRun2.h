#ifndef MSSMvsSMRun2Legacy_HttSystematics_MSSMvsSMRun2_h
#define MSSMvsSMRun2Legacy_HttSystematics_MSSMvsSMRun2_h
#include "CombineHarvester/CombineTools/interface/CombineHarvester.h"

namespace ch {
// Run2 MSSM (with SM categories) analysis systematics
// Implemented in src/HttSystematics_MSSMvsSMRun2.cc
void AddMSSMvsSMRun2Systematics(CombineHarvester& cb, bool jetfakes, bool embedding, bool regional_jec, bool ggh_wg1, bool qqh_wg1, int era, bool mva=false, bool sm=false);
}

#endif
