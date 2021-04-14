#include "CombineHarvester/MSSMvsSMRun2Legacy/interface/HttSystematics_MSSMvsSMRun2.h"
#include "CombineHarvester/CombineTools/interface/Process.h"
#include "CombineHarvester/CombineTools/interface/Systematics.h"
#include "CombineHarvester/CombineTools/interface/Utilities.h"
#include <string>
#include <vector>

using namespace std;

namespace ch {

using ch::syst::SystMap;
using ch::syst::SystMapAsymm;
using ch::syst::era;
using ch::syst::channel;
using ch::syst::bin_id;
using ch::syst::mass;
using ch::syst::process;
using ch::syst::bin;
using ch::JoinStr;

void AddMSSMvsSMRun2Systematics(CombineHarvester &cb, bool jetfakes, bool embedding, bool regional_jec, bool ggh_wg1, bool qqh_wg1, int era, bool mva, bool sm) {

  // ##########################################################################
  // Define groups of signal processes
  // ##########################################################################

  std::vector<std::string> signals_ggH = {
      // STXS stage 0
      "ggH_htt125",
      // STXS stage 1.1
      "ggH_FWDH_htt",
      "ggH_PTH_200_300_htt",
      "ggH_PTH_300_450_htt",
      "ggH_PTH_450_650_htt",
      "ggH_PTH_GT650_htt",
      "ggH_0J_PTH_0_10_htt",
      "ggH_0J_PTH_GT10_htt",
      "ggH_1J_PTH_0_60_htt",
      "ggH_1J_PTH_60_120_htt",
      "ggH_1J_PTH_120_200_htt",
      "ggH_GE2J_MJJ_0_350_PTH_0_60_htt",
      "ggH_GE2J_MJJ_0_350_PTH_60_120_htt",
      "ggH_GE2J_MJJ_0_350_PTH_120_200_htt",
      "ggH_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_0_25_htt",
      "ggH_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_GT25_htt",
      "ggH_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_0_25_htt",
      "ggH_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_GT25_htt",
      };
  std::vector<std::string> signals_qqH = {
      // STXS stage 0
      "qqH_htt125",
      // STXS stage 1
      "qqH_FWDH_htt",
      "qqH_0J_htt",
      "qqH_1J_htt",
      "qqH_GE2J_MJJ_0_60_htt",
      "qqH_GE2J_MJJ_60_120_htt",
      "qqH_GE2J_MJJ_120_350_htt",
      "qqH_GE2J_MJJ_GT350_PTH_GT200_htt",
      "qqH_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_0_25_htt",
      "qqH_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_GT25_htt",
      "qqH_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_0_25_htt",
      "qqH_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_GT25_htt"
      };
  std::vector<std::string> signals_VH = {
      // STXS stage 0
      "WH125", "ZH125", "ttH125"};
  std::vector<std::string> signals_ggHToWW = {
     // STXS stage 0
     "ggHWW125"};
  std::vector<std::string> signals_qqHToWW = {
     // STXS stage 0
     "qqHWW125"};
  std::vector<std::string> signals_VHToWW = {
      // STXS stage 0
      "WHWW125", "ZHWW125"};
  std::vector<std::string> signals = JoinStr({signals_ggH, signals_qqH, signals_VH, {"qqh", "ggh"}});
  std::vector<std::string> signals_HWW = JoinStr({signals_ggHToWW, signals_qqHToWW, signals_VHToWW});

  std::vector<std::string> mssm_ggH_signals = {"ggH_t", "ggH_b", "ggH_i", "ggh_t", "ggh_b", "ggh_i", "ggA_t", "ggA_b", "ggA_i"};
  std::vector<std::string> mssm_bbH_signals = {"bbA", "bbH", "bbh", "bbH_500", "bbH_1400"};
  std::vector<std::string> mssm_signals = JoinStr({mssm_ggH_signals, mssm_bbH_signals});
  std::vector<std::string> jetFakes = {"jetFakes"};
  if(sm == true){
    jetFakes = {"jetFakesSM"};
    std::cout << "[INFO] Using jetFakesSM instead of jetFakes" << std::endl;
  }
  // All processes being taken from simulation
  std::vector<std::string> mc_processes =
      JoinStr({
              signals,
              signals_HWW,
              mssm_signals,
              {"ZTT", "TT", "TTT", "TTL", "TTJ", "W", "ZJ", "ZL", "VV", "VVT", "VVL", "VVJ", "ST"}
              });

  std::vector<int> nobtag_catagories = {2, 3, 4, 5, 6, 7, 8, 9,10,
                                    11,12,13,14,15,16,17,18,19,20,
                                    21,22,23,24,25,26,27,28,29,30,
                                    31,32,33,34}; // SM and MSSM no-btag categories
  std::vector<int> btag_catagories = {35,36,37};

  std::vector<int> mssm_categories = {300,2,32,33,34,35,36,37}; // Useful in we need to use different treatment of some uncertainties for 

  std::vector<int> mssm_nobtag_catagories = {32,33,34};

   // ##########################################################################
   // Uncertainty: b tagging acceptance uncertainties for pdf and scale and hdamp variations.
   // References:
   // - "Talk in MSSM HTT Meeting by Danny"
   //   (https://indico.cern.ch/event/990440/contributions/4167707/attachments/2167087/3659678/MSSM_signal.pdf)
   //   (file source: https://cernbox.cern.ch/index.php/s/9cgjdpaTeqYEaFZ
   // Notes:
   // ##########################################################################

   cb.cp().process(mssm_bbH_signals).AddSyst(cb, "pdf_bbH_ACCEPT", "lnN", SystMap<channel,ch::syst::era,bin_id,mass>::init
     ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"60"}, 0.996)
     ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"80"}, 0.996)
     ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"100"}, 0.996)
     ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"120"}, 0.996)
     ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"125"}, 0.996)
     ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"130"}, 0.996)
     ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"140"}, 0.995)
     ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"160"}, 0.995)
     ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"180"}, 0.995)
     ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"200"}, 0.995)
     ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"250"}, 0.994)
     ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"300"}, 0.994)
     ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"350"}, 0.993)
     ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"400"}, 0.993)
     ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"450"}, 0.993)
     ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"500"}, 0.993)
     ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"600"}, 0.993)
     ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"700"}, 0.992)
     ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"800"}, 0.990)
     ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"900"}, 0.988)
     ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"1000"}, 0.987)
     ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"1200"}, 0.986)
     ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"1400"}, 0.985)
     ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"1600"}, 0.984)
     ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"1800"}, 0.984)
     ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"2000"}, 0.983)
     ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"2300"}, 0.982)
     ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"2600"}, 0.980)
     ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"2900"}, 0.979)
     ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"3200"}, 0.976)
     ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"3500"}, 0.974)
     ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {"60"}, 0.996)
     ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {"80"}, 0.996)
     ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {"100"}, 0.996)
     ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {"120"}, 0.995)
     ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {"125"}, 0.995)
     ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {"130"}, 0.995)
     ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {"140"}, 0.995)
     ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {"160"}, 0.994)
     ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {"180"}, 0.994)
     ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {"200"}, 0.994)
     ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {"250"}, 0.994)
     ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {"300"}, 0.993)
     ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {"350"}, 0.993)
     ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {"400"}, 0.993)
     ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {"450"}, 0.993)
     ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {"500"}, 0.993)
     ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {"600"}, 0.993)
     ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {"700"}, 0.992)
     ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {"800"}, 0.990)
     ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {"900"}, 0.988)
     ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {"1000"}, 0.987)
     ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {"1200"}, 0.986)
     ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {"1400"}, 0.985)
     ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {"1600"}, 0.984)
     ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {"1800"}, 0.984)
     ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {"2000"}, 0.983)
     ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {"2300"}, 0.981)
     ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {"2600"}, 0.980)
     ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {"2900"}, 0.979)
     ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {"3200"}, 0.976)
     ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {"3500"}, 0.974)
     ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {"60"}, 0.996)
     ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {"80"}, 0.996)
     ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {"100"}, 0.996)
     ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {"120"}, 0.995)
     ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {"125"}, 0.995)
     ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {"130"}, 0.995)
     ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {"140"}, 0.995)
     ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {"160"}, 0.994)
     ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {"180"}, 0.994)
     ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {"200"}, 0.994)
     ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {"250"}, 0.994)
     ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {"300"}, 0.993)
     ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {"350"}, 0.993)
     ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {"400"}, 0.993)
     ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {"450"}, 0.993)
     ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {"500"}, 0.993)
     ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {"600"}, 0.993)
     ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {"700"}, 0.991)
     ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {"800"}, 0.990)
     ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {"900"}, 0.988)
     ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {"1000"}, 0.986)
     ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {"1200"}, 0.986)
     ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {"1400"}, 0.985)
     ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {"1600"}, 0.984)
     ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {"1800"}, 0.984)
     ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {"2000"}, 0.983)
     ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {"2300"}, 0.981)
     ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {"2600"}, 0.980)
     ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {"2900"}, 0.978)
     ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {"3200"}, 0.976)
     ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {"3500"}, 0.973)
     ({"et","mt","tt"}, {"2016"}, btag_catagories, {"60"}, 1.023)
     ({"et","mt","tt"}, {"2016"}, btag_catagories, {"80"}, 1.020)
     ({"et","mt","tt"}, {"2016"}, btag_catagories, {"100"}, 1.017)
     ({"et","mt","tt"}, {"2016"}, btag_catagories, {"120"}, 1.015)
     ({"et","mt","tt"}, {"2016"}, btag_catagories, {"125"}, 1.015)
     ({"et","mt","tt"}, {"2016"}, btag_catagories, {"130"}, 1.015)
     ({"et","mt","tt"}, {"2016"}, btag_catagories, {"140"}, 1.016)
     ({"et","mt","tt"}, {"2016"}, btag_catagories, {"160"}, 1.017)
     ({"et","mt","tt"}, {"2016"}, btag_catagories, {"180"}, 1.017)
     ({"et","mt","tt"}, {"2016"}, btag_catagories, {"200"}, 1.016)
     ({"et","mt","tt"}, {"2016"}, btag_catagories, {"250"}, 1.015)
     ({"et","mt","tt"}, {"2016"}, btag_catagories, {"300"}, 1.014)
     ({"et","mt","tt"}, {"2016"}, btag_catagories, {"350"}, 1.013)
     ({"et","mt","tt"}, {"2016"}, btag_catagories, {"400"}, 1.012)
     ({"et","mt","tt"}, {"2016"}, btag_catagories, {"450"}, 1.012)
     ({"et","mt","tt"}, {"2016"}, btag_catagories, {"500"}, 1.011)
     ({"et","mt","tt"}, {"2016"}, btag_catagories, {"600"}, 1.010)
     ({"et","mt","tt"}, {"2016"}, btag_catagories, {"700"}, 1.012)
     ({"et","mt","tt"}, {"2016"}, btag_catagories, {"800"}, 1.013)
     ({"et","mt","tt"}, {"2016"}, btag_catagories, {"900"}, 1.014)
     ({"et","mt","tt"}, {"2016"}, btag_catagories, {"1000"}, 1.016)
     ({"et","mt","tt"}, {"2016"}, btag_catagories, {"1200"}, 1.017)
     ({"et","mt","tt"}, {"2016"}, btag_catagories, {"1400"}, 1.018)
     ({"et","mt","tt"}, {"2016"}, btag_catagories, {"1600"}, 1.018)
     ({"et","mt","tt"}, {"2016"}, btag_catagories, {"1800"}, 1.019)
     ({"et","mt","tt"}, {"2016"}, btag_catagories, {"2000"}, 1.019)
     ({"et","mt","tt"}, {"2016"}, btag_catagories, {"2300"}, 1.020)
     ({"et","mt","tt"}, {"2016"}, btag_catagories, {"2600"}, 1.022)
     ({"et","mt","tt"}, {"2016"}, btag_catagories, {"2900"}, 1.023)
     ({"et","mt","tt"}, {"2016"}, btag_catagories, {"3200"}, 1.025)
     ({"et","mt","tt"}, {"2016"}, btag_catagories, {"3500"}, 1.027)
     ({"em"}, {"2016"}, btag_catagories, {"60"}, 1.022)
     ({"em"}, {"2016"}, btag_catagories, {"80"}, 1.019)
     ({"em"}, {"2016"}, btag_catagories, {"100"}, 1.016)
     ({"em"}, {"2016"}, btag_catagories, {"120"}, 1.015)
     ({"em"}, {"2016"}, btag_catagories, {"125"}, 1.015)
     ({"em"}, {"2016"}, btag_catagories, {"130"}, 1.015)
     ({"em"}, {"2016"}, btag_catagories, {"140"}, 1.014)
     ({"em"}, {"2016"}, btag_catagories, {"160"}, 1.016)
     ({"em"}, {"2016"}, btag_catagories, {"180"}, 1.015)
     ({"em"}, {"2016"}, btag_catagories, {"200"}, 1.015)
     ({"em"}, {"2016"}, btag_catagories, {"250"}, 1.013)
     ({"em"}, {"2016"}, btag_catagories, {"300"}, 1.012)
     ({"em"}, {"2016"}, btag_catagories, {"350"}, 1.010)
     ({"em"}, {"2016"}, btag_catagories, {"400"}, 1.009)
     ({"em"}, {"2016"}, btag_catagories, {"450"}, 1.008)
     ({"em"}, {"2016"}, btag_catagories, {"500"}, 1.007)
     ({"em"}, {"2016"}, btag_catagories, {"600"}, 1.006)
     ({"em"}, {"2016"}, btag_catagories, {"700"}, 1.007)
     ({"em"}, {"2016"}, btag_catagories, {"800"}, 1.008)
     ({"em"}, {"2016"}, btag_catagories, {"900"}, 1.010)
     ({"em"}, {"2016"}, btag_catagories, {"1000"}, 1.011)
     ({"em"}, {"2016"}, btag_catagories, {"1200"}, 1.012)
     ({"em"}, {"2016"}, btag_catagories, {"1400"}, 1.014)
     ({"em"}, {"2016"}, btag_catagories, {"1600"}, 1.014)
     ({"em"}, {"2016"}, btag_catagories, {"1800"}, 1.013)
     ({"em"}, {"2016"}, btag_catagories, {"2000"}, 1.013)
     ({"em"}, {"2016"}, btag_catagories, {"2300"}, 1.013)
     ({"em"}, {"2016"}, btag_catagories, {"2600"}, 1.015)
     ({"em"}, {"2016"}, btag_catagories, {"2900"}, 1.017)
     ({"em"}, {"2016"}, btag_catagories, {"3200"}, 1.019)
     ({"em"}, {"2016"}, btag_catagories, {"3500"}, 1.022)
     ({"et","mt","tt"}, {"2017"}, btag_catagories, {"60"}, 1.023)
     ({"et","mt","tt"}, {"2017"}, btag_catagories, {"80"}, 1.019)
     ({"et","mt","tt"}, {"2017"}, btag_catagories, {"100"}, 1.016)
     ({"et","mt","tt"}, {"2017"}, btag_catagories, {"120"}, 1.016)
     ({"et","mt","tt"}, {"2017"}, btag_catagories, {"125"}, 1.016)
     ({"et","mt","tt"}, {"2017"}, btag_catagories, {"130"}, 1.016)
     ({"et","mt","tt"}, {"2017"}, btag_catagories, {"140"}, 1.015)
     ({"et","mt","tt"}, {"2017"}, btag_catagories, {"160"}, 1.016)
     ({"et","mt","tt"}, {"2017"}, btag_catagories, {"180"}, 1.016)
     ({"et","mt","tt"}, {"2017"}, btag_catagories, {"200"}, 1.015)
     ({"et","mt","tt"}, {"2017"}, btag_catagories, {"250"}, 1.014)
     ({"et","mt","tt"}, {"2017"}, btag_catagories, {"300"}, 1.013)
     ({"et","mt","tt"}, {"2017"}, btag_catagories, {"350"}, 1.012)
     ({"et","mt","tt"}, {"2017"}, btag_catagories, {"400"}, 1.011)
     ({"et","mt","tt"}, {"2017"}, btag_catagories, {"450"}, 1.011)
     ({"et","mt","tt"}, {"2017"}, btag_catagories, {"500"}, 1.010)
     ({"et","mt","tt"}, {"2017"}, btag_catagories, {"600"}, 1.010)
     ({"et","mt","tt"}, {"2017"}, btag_catagories, {"700"}, 1.011)
     ({"et","mt","tt"}, {"2017"}, btag_catagories, {"800"}, 1.012)
     ({"et","mt","tt"}, {"2017"}, btag_catagories, {"900"}, 1.014)
     ({"et","mt","tt"}, {"2017"}, btag_catagories, {"1000"}, 1.015)
     ({"et","mt","tt"}, {"2017"}, btag_catagories, {"1200"}, 1.015)
     ({"et","mt","tt"}, {"2017"}, btag_catagories, {"1400"}, 1.016)
     ({"et","mt","tt"}, {"2017"}, btag_catagories, {"1600"}, 1.016)
     ({"et","mt","tt"}, {"2017"}, btag_catagories, {"1800"}, 1.016)
     ({"et","mt","tt"}, {"2017"}, btag_catagories, {"2000"}, 1.016)
     ({"et","mt","tt"}, {"2017"}, btag_catagories, {"2300"}, 1.018)
     ({"et","mt","tt"}, {"2017"}, btag_catagories, {"2600"}, 1.019)
     ({"et","mt","tt"}, {"2017"}, btag_catagories, {"2900"}, 1.020)
     ({"et","mt","tt"}, {"2017"}, btag_catagories, {"3200"}, 1.022)
     ({"et","mt","tt"}, {"2017"}, btag_catagories, {"3500"}, 1.024)
     ({"em"}, {"2017"}, btag_catagories, {"60"}, 1.022)
     ({"em"}, {"2017"}, btag_catagories, {"80"}, 1.019)
     ({"em"}, {"2017"}, btag_catagories, {"100"}, 1.016)
     ({"em"}, {"2017"}, btag_catagories, {"120"}, 1.015)
     ({"em"}, {"2017"}, btag_catagories, {"125"}, 1.015)
     ({"em"}, {"2017"}, btag_catagories, {"130"}, 1.014)
     ({"em"}, {"2017"}, btag_catagories, {"140"}, 1.013)
     ({"em"}, {"2017"}, btag_catagories, {"160"}, 1.014)
     ({"em"}, {"2017"}, btag_catagories, {"180"}, 1.013)
     ({"em"}, {"2017"}, btag_catagories, {"200"}, 1.013)
     ({"em"}, {"2017"}, btag_catagories, {"250"}, 1.011)
     ({"em"}, {"2017"}, btag_catagories, {"300"}, 1.010)
     ({"em"}, {"2017"}, btag_catagories, {"350"}, 1.008)
     ({"em"}, {"2017"}, btag_catagories, {"400"}, 1.007)
     ({"em"}, {"2017"}, btag_catagories, {"450"}, 1.006)
     ({"em"}, {"2017"}, btag_catagories, {"500"}, 1.006)
     ({"em"}, {"2017"}, btag_catagories, {"600"}, 1.005)
     ({"em"}, {"2017"}, btag_catagories, {"700"}, 1.006)
     ({"em"}, {"2017"}, btag_catagories, {"800"}, 1.007)
     ({"em"}, {"2017"}, btag_catagories, {"900"}, 1.009)
     ({"em"}, {"2017"}, btag_catagories, {"1000"}, 1.010)
     ({"em"}, {"2017"}, btag_catagories, {"1200"}, 1.011)
     ({"em"}, {"2017"}, btag_catagories, {"1400"}, 1.012)
     ({"em"}, {"2017"}, btag_catagories, {"1600"}, 1.012)
     ({"em"}, {"2017"}, btag_catagories, {"1800"}, 1.011)
     ({"em"}, {"2017"}, btag_catagories, {"2000"}, 1.011)
     ({"em"}, {"2017"}, btag_catagories, {"2300"}, 1.011)
     ({"em"}, {"2017"}, btag_catagories, {"2600"}, 1.012)
     ({"em"}, {"2017"}, btag_catagories, {"2900"}, 1.014)
     ({"em"}, {"2017"}, btag_catagories, {"3200"}, 1.017)
     ({"em"}, {"2017"}, btag_catagories, {"3500"}, 1.020)
     ({"et","mt","tt"}, {"2018"}, btag_catagories, {"60"}, 1.022)
     ({"et","mt","tt"}, {"2018"}, btag_catagories, {"80"}, 1.018)
     ({"et","mt","tt"}, {"2018"}, btag_catagories, {"100"}, 1.015)
     ({"et","mt","tt"}, {"2018"}, btag_catagories, {"120"}, 1.016)
     ({"et","mt","tt"}, {"2018"}, btag_catagories, {"125"}, 1.016)
     ({"et","mt","tt"}, {"2018"}, btag_catagories, {"130"}, 1.015)
     ({"et","mt","tt"}, {"2018"}, btag_catagories, {"140"}, 1.014)
     ({"et","mt","tt"}, {"2018"}, btag_catagories, {"160"}, 1.016)
     ({"et","mt","tt"}, {"2018"}, btag_catagories, {"180"}, 1.016)
     ({"et","mt","tt"}, {"2018"}, btag_catagories, {"200"}, 1.015)
     ({"et","mt","tt"}, {"2018"}, btag_catagories, {"250"}, 1.014)
     ({"et","mt","tt"}, {"2018"}, btag_catagories, {"300"}, 1.012)
     ({"et","mt","tt"}, {"2018"}, btag_catagories, {"350"}, 1.011)
     ({"et","mt","tt"}, {"2018"}, btag_catagories, {"400"}, 1.010)
     ({"et","mt","tt"}, {"2018"}, btag_catagories, {"450"}, 1.010)
     ({"et","mt","tt"}, {"2018"}, btag_catagories, {"500"}, 1.010)
     ({"et","mt","tt"}, {"2018"}, btag_catagories, {"600"}, 1.009)
     ({"et","mt","tt"}, {"2018"}, btag_catagories, {"700"}, 1.010)
     ({"et","mt","tt"}, {"2018"}, btag_catagories, {"800"}, 1.012)
     ({"et","mt","tt"}, {"2018"}, btag_catagories, {"900"}, 1.013)
     ({"et","mt","tt"}, {"2018"}, btag_catagories, {"1000"}, 1.015)
     ({"et","mt","tt"}, {"2018"}, btag_catagories, {"1200"}, 1.015)
     ({"et","mt","tt"}, {"2018"}, btag_catagories, {"1400"}, 1.016)
     ({"et","mt","tt"}, {"2018"}, btag_catagories, {"1600"}, 1.016)
     ({"et","mt","tt"}, {"2018"}, btag_catagories, {"1800"}, 1.016)
     ({"et","mt","tt"}, {"2018"}, btag_catagories, {"2000"}, 1.016)
     ({"et","mt","tt"}, {"2018"}, btag_catagories, {"2300"}, 1.018)
     ({"et","mt","tt"}, {"2018"}, btag_catagories, {"2600"}, 1.019)
     ({"et","mt","tt"}, {"2018"}, btag_catagories, {"2900"}, 1.020)
     ({"et","mt","tt"}, {"2018"}, btag_catagories, {"3200"}, 1.022)
     ({"et","mt","tt"}, {"2018"}, btag_catagories, {"3500"}, 1.024)
     ({"em"}, {"2018"}, btag_catagories, {"60"}, 1.022)
     ({"em"}, {"2018"}, btag_catagories, {"80"}, 1.018)
     ({"em"}, {"2018"}, btag_catagories, {"100"}, 1.015)
     ({"em"}, {"2018"}, btag_catagories, {"120"}, 1.015)
     ({"em"}, {"2018"}, btag_catagories, {"125"}, 1.015)
     ({"em"}, {"2018"}, btag_catagories, {"130"}, 1.014)
     ({"em"}, {"2018"}, btag_catagories, {"140"}, 1.013)
     ({"em"}, {"2018"}, btag_catagories, {"160"}, 1.014)
     ({"em"}, {"2018"}, btag_catagories, {"180"}, 1.013)
     ({"em"}, {"2018"}, btag_catagories, {"200"}, 1.013)
     ({"em"}, {"2018"}, btag_catagories, {"250"}, 1.011)
     ({"em"}, {"2018"}, btag_catagories, {"300"}, 1.009)
     ({"em"}, {"2018"}, btag_catagories, {"350"}, 1.008)
     ({"em"}, {"2018"}, btag_catagories, {"400"}, 1.006)
     ({"em"}, {"2018"}, btag_catagories, {"450"}, 1.006)
     ({"em"}, {"2018"}, btag_catagories, {"500"}, 1.006)
     ({"em"}, {"2018"}, btag_catagories, {"600"}, 1.005)
     ({"em"}, {"2018"}, btag_catagories, {"700"}, 1.006)
     ({"em"}, {"2018"}, btag_catagories, {"800"}, 1.007)
     ({"em"}, {"2018"}, btag_catagories, {"900"}, 1.008)
     ({"em"}, {"2018"}, btag_catagories, {"1000"}, 1.009)
     ({"em"}, {"2018"}, btag_catagories, {"1200"}, 1.010)
     ({"em"}, {"2018"}, btag_catagories, {"1400"}, 1.011)
     ({"em"}, {"2018"}, btag_catagories, {"1600"}, 1.011)
     ({"em"}, {"2018"}, btag_catagories, {"1800"}, 1.010)
     ({"em"}, {"2018"}, btag_catagories, {"2000"}, 1.010)
     ({"em"}, {"2018"}, btag_catagories, {"2300"}, 1.010)
     ({"em"}, {"2018"}, btag_catagories, {"2600"}, 1.011)
     ({"em"}, {"2018"}, btag_catagories, {"2900"}, 1.013)
     ({"em"}, {"2018"}, btag_catagories, {"3200"}, 1.016)
     ({"em"}, {"2018"}, btag_catagories, {"3500"}, 1.019));

   cb.cp().process(mssm_bbH_signals).AddSyst(cb, "QCDscaleAndHdamp_bbH_ACCEPT", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
    ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"60"}, 1.010, 0.993)
    ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"80"}, 1.013, 0.992)
    ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"100"}, 1.016, 0.992)
    ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"120"}, 1.009, 0.993)
    ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"125"}, 1.007, 0.993)
    ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"130"}, 1.009, 0.994)
    ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"140"}, 1.012, 0.995)
    ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"160"}, 1.013, 0.994)
    ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"180"}, 1.013, 0.994)
    ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"200"}, 1.013, 0.994)
    ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"250"}, 1.013, 0.994)
    ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"300"}, 1.014, 0.993)
    ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"350"}, 1.014, 0.993)
    ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"400"}, 1.014, 0.993)
    ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"450"}, 1.014, 0.990)
    ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"500"}, 1.013, 0.988)
    ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"600"}, 1.012, 0.983)
    ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"700"}, 1.016, 0.984)
    ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"800"}, 1.020, 0.984)
    ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"900"}, 1.024, 0.984)
    ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"1000"}, 1.028, 0.985)
    ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"1200"}, 1.027, 0.980)
    ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"1400"}, 1.027, 0.975)
    ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"1600"}, 1.031, 0.977)
    ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"1800"}, 1.034, 0.979)
    ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"2000"}, 1.038, 0.981)
    ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"2300"}, 1.040, 0.979)
    ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"2600"}, 1.040, 0.970)
    ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"2900"}, 1.041, 0.961)
    ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"3200"}, 1.046, 0.962)
    ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"3500"}, 1.051, 0.962)
    ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {"60"}, 1.010, 0.992)
    ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {"80"}, 1.015, 0.992)
    ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {"100"}, 1.020, 0.993)
    ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {"120"}, 1.010, 0.992)
    ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {"125"}, 1.008, 0.992)
    ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {"130"}, 1.009, 0.993)
    ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {"140"}, 1.012, 0.994)
    ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {"160"}, 1.014, 0.993)
    ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {"180"}, 1.014, 0.993)
    ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {"200"}, 1.014, 0.993)
    ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {"250"}, 1.014, 0.993)
    ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {"300"}, 1.014, 0.993)
    ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {"350"}, 1.014, 0.993)
    ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {"400"}, 1.014, 0.993)
    ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {"450"}, 1.014, 0.990)
    ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {"500"}, 1.013, 0.988)
    ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {"600"}, 1.012, 0.982)
    ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {"700"}, 1.016, 0.983)
    ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {"800"}, 1.020, 0.984)
    ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {"900"}, 1.024, 0.984)
    ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {"1000"}, 1.028, 0.985)
    ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {"1200"}, 1.028, 0.980)
    ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {"1400"}, 1.028, 0.975)
    ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {"1600"}, 1.032, 0.977)
    ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {"1800"}, 1.036, 0.979)
    ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {"2000"}, 1.040, 0.981)
    ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {"2300"}, 1.040, 0.979)
    ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {"2600"}, 1.040, 0.971)
    ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {"2900"}, 1.041, 0.964)
    ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {"3200"}, 1.046, 0.964)
    ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {"3500"}, 1.050, 0.963)
    ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {"60"}, 1.011, 0.992)
    ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {"80"}, 1.016, 0.992)
    ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {"100"}, 1.021, 0.993)
    ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {"120"}, 1.011, 0.992)
    ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {"125"}, 1.009, 0.992)
    ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {"130"}, 1.010, 0.993)
    ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {"140"}, 1.013, 0.994)
    ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {"160"}, 1.015, 0.993)
    ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {"180"}, 1.015, 0.993)
    ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {"200"}, 1.015, 0.993)
    ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {"250"}, 1.015, 0.993)
    ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {"300"}, 1.014, 0.992)
    ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {"350"}, 1.014, 0.992)
    ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {"400"}, 1.014, 0.992)
    ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {"450"}, 1.014, 0.989)
    ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {"500"}, 1.013, 0.986)
    ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {"600"}, 1.012, 0.981)
    ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {"700"}, 1.016, 0.982)
    ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {"800"}, 1.020, 0.982)
    ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {"900"}, 1.024, 0.983)
    ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {"1000"}, 1.028, 0.984)
    ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {"1200"}, 1.028, 0.978)
    ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {"1400"}, 1.029, 0.973)
    ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {"1600"}, 1.033, 0.975)
    ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {"1800"}, 1.036, 0.978)
    ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {"2000"}, 1.040, 0.980)
    ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {"2300"}, 1.042, 0.978)
    ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {"2600"}, 1.042, 0.970)
    ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {"2900"}, 1.043, 0.962)
    ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {"3200"}, 1.048, 0.963)
    ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {"3500"}, 1.053, 0.964)
    ({"et","mt","tt"}, {"2016"}, btag_catagories, {"60"}, 0.938, 1.042)
    ({"et","mt","tt"}, {"2016"}, btag_catagories, {"80"}, 0.932, 1.038)
    ({"et","mt","tt"}, {"2016"}, btag_catagories, {"100"}, 0.926, 1.035)
    ({"et","mt","tt"}, {"2016"}, btag_catagories, {"120"}, 0.964, 1.028)
    ({"et","mt","tt"}, {"2016"}, btag_catagories, {"125"}, 0.974, 1.026)
    ({"et","mt","tt"}, {"2016"}, btag_catagories, {"130"}, 0.969, 1.024)
    ({"et","mt","tt"}, {"2016"}, btag_catagories, {"140"}, 0.958, 1.019)
    ({"et","mt","tt"}, {"2016"}, btag_catagories, {"160"}, 0.958, 1.021)
    ({"et","mt","tt"}, {"2016"}, btag_catagories, {"180"}, 0.959, 1.020)
    ({"et","mt","tt"}, {"2016"}, btag_catagories, {"200"}, 0.961, 1.020)
    ({"et","mt","tt"}, {"2016"}, btag_catagories, {"250"}, 0.964, 1.018)
    ({"et","mt","tt"}, {"2016"}, btag_catagories, {"300"}, 0.967, 1.017)
    ({"et","mt","tt"}, {"2016"}, btag_catagories, {"350"}, 0.971, 1.015)
    ({"et","mt","tt"}, {"2016"}, btag_catagories, {"400"}, 0.974, 1.014)
    ({"et","mt","tt"}, {"2016"}, btag_catagories, {"450"}, 0.976, 1.017)
    ({"et","mt","tt"}, {"2016"}, btag_catagories, {"500"}, 0.978, 1.019)
    ({"et","mt","tt"}, {"2016"}, btag_catagories, {"600"}, 0.983, 1.025)
    ({"et","mt","tt"}, {"2016"}, btag_catagories, {"700"}, 0.978, 1.023)
    ({"et","mt","tt"}, {"2016"}, btag_catagories, {"800"}, 0.974, 1.022)
    ({"et","mt","tt"}, {"2016"}, btag_catagories, {"900"}, 0.970, 1.020)
    ({"et","mt","tt"}, {"2016"}, btag_catagories, {"1000"}, 0.965, 1.019)
    ({"et","mt","tt"}, {"2016"}, btag_catagories, {"1200"}, 0.966, 1.024)
    ({"et","mt","tt"}, {"2016"}, btag_catagories, {"1400"}, 0.967, 1.030)
    ({"et","mt","tt"}, {"2016"}, btag_catagories, {"1600"}, 0.964, 1.027)
    ({"et","mt","tt"}, {"2016"}, btag_catagories, {"1800"}, 0.960, 1.024)
    ({"et","mt","tt"}, {"2016"}, btag_catagories, {"2000"}, 0.957, 1.021)
    ({"et","mt","tt"}, {"2016"}, btag_catagories, {"2300"}, 0.957, 1.022)
    ({"et","mt","tt"}, {"2016"}, btag_catagories, {"2600"}, 0.956, 1.032)
    ({"et","mt","tt"}, {"2016"}, btag_catagories, {"2900"}, 0.955, 1.042)
    ({"et","mt","tt"}, {"2016"}, btag_catagories, {"3200"}, 0.950, 1.042)
    ({"et","mt","tt"}, {"2016"}, btag_catagories, {"3500"}, 0.946, 1.041)
    ({"em"}, {"2016"}, btag_catagories, {"60"}, 0.950, 1.042)
    ({"em"}, {"2016"}, btag_catagories, {"80"}, 0.937, 1.042)
    ({"em"}, {"2016"}, btag_catagories, {"100"}, 0.924, 1.041)
    ({"em"}, {"2016"}, btag_catagories, {"120"}, 0.965, 1.034)
    ({"em"}, {"2016"}, btag_catagories, {"125"}, 0.975, 1.032)
    ({"em"}, {"2016"}, btag_catagories, {"130"}, 0.973, 1.027)
    ({"em"}, {"2016"}, btag_catagories, {"140"}, 0.968, 1.017)
    ({"em"}, {"2016"}, btag_catagories, {"160"}, 0.966, 1.019)
    ({"em"}, {"2016"}, btag_catagories, {"180"}, 0.967, 1.018)
    ({"em"}, {"2016"}, btag_catagories, {"200"}, 0.969, 1.018)
    ({"em"}, {"2016"}, btag_catagories, {"250"}, 0.972, 1.016)
    ({"em"}, {"2016"}, btag_catagories, {"300"}, 0.976, 1.015)
    ({"em"}, {"2016"}, btag_catagories, {"350"}, 0.979, 1.013)
    ({"em"}, {"2016"}, btag_catagories, {"400"}, 0.983, 1.012)
    ({"em"}, {"2016"}, btag_catagories, {"450"}, 0.984, 1.013)
    ({"em"}, {"2016"}, btag_catagories, {"500"}, 0.986, 1.014)
    ({"em"}, {"2016"}, btag_catagories, {"600"}, 0.988, 1.015)
    ({"em"}, {"2016"}, btag_catagories, {"700"}, 0.985, 1.014)
    ({"em"}, {"2016"}, btag_catagories, {"800"}, 0.982, 1.014)
    ({"em"}, {"2016"}, btag_catagories, {"900"}, 0.980, 1.013)
    ({"em"}, {"2016"}, btag_catagories, {"1000"}, 0.977, 1.013)
    ({"em"}, {"2016"}, btag_catagories, {"1200"}, 0.976, 1.024)
    ({"em"}, {"2016"}, btag_catagories, {"1400"}, 0.974, 1.035)
    ({"em"}, {"2016"}, btag_catagories, {"1600"}, 0.972, 1.028)
    ({"em"}, {"2016"}, btag_catagories, {"1800"}, 0.970, 1.022)
    ({"em"}, {"2016"}, btag_catagories, {"2000"}, 0.968, 1.015)
    ({"em"}, {"2016"}, btag_catagories, {"2300"}, 0.971, 1.014)
    ({"em"}, {"2016"}, btag_catagories, {"2600"}, 0.969, 1.023)
    ({"em"}, {"2016"}, btag_catagories, {"2900"}, 0.967, 1.033)
    ({"em"}, {"2016"}, btag_catagories, {"3200"}, 0.962, 1.034)
    ({"em"}, {"2016"}, btag_catagories, {"3500"}, 0.956, 1.035)
    ({"et","mt","tt"}, {"2017"}, btag_catagories, {"60"}, 0.944, 1.043)
    ({"et","mt","tt"}, {"2017"}, btag_catagories, {"80"}, 0.932, 1.036)
    ({"et","mt","tt"}, {"2017"}, btag_catagories, {"100"}, 0.920, 1.029)
    ({"et","mt","tt"}, {"2017"}, btag_catagories, {"120"}, 0.962, 1.026)
    ({"et","mt","tt"}, {"2017"}, btag_catagories, {"125"}, 0.972, 1.025)
    ({"et","mt","tt"}, {"2017"}, btag_catagories, {"130"}, 0.969, 1.023)
    ({"et","mt","tt"}, {"2017"}, btag_catagories, {"140"}, 0.962, 1.018)
    ({"et","mt","tt"}, {"2017"}, btag_catagories, {"160"}, 0.960, 1.019)
    ({"et","mt","tt"}, {"2017"}, btag_catagories, {"180"}, 0.962, 1.018)
    ({"et","mt","tt"}, {"2017"}, btag_catagories, {"200"}, 0.963, 1.018)
    ({"et","mt","tt"}, {"2017"}, btag_catagories, {"250"}, 0.967, 1.016)
    ({"et","mt","tt"}, {"2017"}, btag_catagories, {"300"}, 0.970, 1.015)
    ({"et","mt","tt"}, {"2017"}, btag_catagories, {"350"}, 0.974, 1.013)
    ({"et","mt","tt"}, {"2017"}, btag_catagories, {"400"}, 0.978, 1.012)
    ({"et","mt","tt"}, {"2017"}, btag_catagories, {"450"}, 0.980, 1.015)
    ({"et","mt","tt"}, {"2017"}, btag_catagories, {"500"}, 0.981, 1.018)
    ({"et","mt","tt"}, {"2017"}, btag_catagories, {"600"}, 0.984, 1.023)
    ({"et","mt","tt"}, {"2017"}, btag_catagories, {"700"}, 0.980, 1.021)
    ({"et","mt","tt"}, {"2017"}, btag_catagories, {"800"}, 0.976, 1.020)
    ({"et","mt","tt"}, {"2017"}, btag_catagories, {"900"}, 0.972, 1.018)
    ({"et","mt","tt"}, {"2017"}, btag_catagories, {"1000"}, 0.968, 1.017)
    ({"et","mt","tt"}, {"2017"}, btag_catagories, {"1200"}, 0.969, 1.022)
    ({"et","mt","tt"}, {"2017"}, btag_catagories, {"1400"}, 0.970, 1.027)
    ({"et","mt","tt"}, {"2017"}, btag_catagories, {"1600"}, 0.967, 1.024)
    ({"et","mt","tt"}, {"2017"}, btag_catagories, {"1800"}, 0.964, 1.022)
    ({"et","mt","tt"}, {"2017"}, btag_catagories, {"2000"}, 0.961, 1.019)
    ({"et","mt","tt"}, {"2017"}, btag_catagories, {"2300"}, 0.962, 1.020)
    ({"et","mt","tt"}, {"2017"}, btag_catagories, {"2600"}, 0.961, 1.027)
    ({"et","mt","tt"}, {"2017"}, btag_catagories, {"2900"}, 0.960, 1.035)
    ({"et","mt","tt"}, {"2017"}, btag_catagories, {"3200"}, 0.956, 1.035)
    ({"et","mt","tt"}, {"2017"}, btag_catagories, {"3500"}, 0.952, 1.035)
    ({"em"}, {"2017"}, btag_catagories, {"60"}, 0.952, 1.042)
    ({"em"}, {"2017"}, btag_catagories, {"80"}, 0.935, 1.038)
    ({"em"}, {"2017"}, btag_catagories, {"100"}, 0.919, 1.035)
    ({"em"}, {"2017"}, btag_catagories, {"120"}, 0.962, 1.030)
    ({"em"}, {"2017"}, btag_catagories, {"125"}, 0.973, 1.029)
    ({"em"}, {"2017"}, btag_catagories, {"130"}, 0.972, 1.025)
    ({"em"}, {"2017"}, btag_catagories, {"140"}, 0.969, 1.018)
    ({"em"}, {"2017"}, btag_catagories, {"160"}, 0.963, 1.017)
    ({"em"}, {"2017"}, btag_catagories, {"180"}, 0.965, 1.016)
    ({"em"}, {"2017"}, btag_catagories, {"200"}, 0.967, 1.016)
    ({"em"}, {"2017"}, btag_catagories, {"250"}, 0.972, 1.014)
    ({"em"}, {"2017"}, btag_catagories, {"300"}, 0.978, 1.013)
    ({"em"}, {"2017"}, btag_catagories, {"350"}, 0.983, 1.011)
    ({"em"}, {"2017"}, btag_catagories, {"400"}, 0.988, 1.010)
    ({"em"}, {"2017"}, btag_catagories, {"450"}, 0.988, 1.012)
    ({"em"}, {"2017"}, btag_catagories, {"500"}, 0.988, 1.013)
    ({"em"}, {"2017"}, btag_catagories, {"600"}, 0.989, 1.016)
    ({"em"}, {"2017"}, btag_catagories, {"700"}, 0.986, 1.015)
    ({"em"}, {"2017"}, btag_catagories, {"800"}, 0.984, 1.014)
    ({"em"}, {"2017"}, btag_catagories, {"900"}, 0.982, 1.012)
    ({"em"}, {"2017"}, btag_catagories, {"1000"}, 0.979, 1.011)
    ({"em"}, {"2017"}, btag_catagories, {"1200"}, 0.978, 1.019)
    ({"em"}, {"2017"}, btag_catagories, {"1400"}, 0.978, 1.027)
    ({"em"}, {"2017"}, btag_catagories, {"1600"}, 0.977, 1.022)
    ({"em"}, {"2017"}, btag_catagories, {"1800"}, 0.976, 1.018)
    ({"em"}, {"2017"}, btag_catagories, {"2000"}, 0.975, 1.013)
    ({"em"}, {"2017"}, btag_catagories, {"2300"}, 0.972, 1.012)
    ({"em"}, {"2017"}, btag_catagories, {"2600"}, 0.968, 1.017)
    ({"em"}, {"2017"}, btag_catagories, {"2900"}, 0.965, 1.022)
    ({"em"}, {"2017"}, btag_catagories, {"3200"}, 0.963, 1.027)
    ({"em"}, {"2017"}, btag_catagories, {"3500"}, 0.961, 1.032)
    ({"et","mt","tt"}, {"2018"}, btag_catagories, {"60"}, 0.941, 1.042)
    ({"et","mt","tt"}, {"2018"}, btag_catagories, {"80"}, 0.930, 1.034)
    ({"et","mt","tt"}, {"2018"}, btag_catagories, {"100"}, 0.920, 1.027)
    ({"et","mt","tt"}, {"2018"}, btag_catagories, {"120"}, 0.962, 1.026)
    ({"et","mt","tt"}, {"2018"}, btag_catagories, {"125"}, 0.972, 1.026)
    ({"et","mt","tt"}, {"2018"}, btag_catagories, {"130"}, 0.969, 1.023)
    ({"et","mt","tt"}, {"2018"}, btag_catagories, {"140"}, 0.962, 1.017)
    ({"et","mt","tt"}, {"2018"}, btag_catagories, {"160"}, 0.961, 1.019)
    ({"et","mt","tt"}, {"2018"}, btag_catagories, {"180"}, 0.962, 1.018)
    ({"et","mt","tt"}, {"2018"}, btag_catagories, {"200"}, 0.964, 1.018)
    ({"et","mt","tt"}, {"2018"}, btag_catagories, {"250"}, 0.968, 1.016)
    ({"et","mt","tt"}, {"2018"}, btag_catagories, {"300"}, 0.972, 1.015)
    ({"et","mt","tt"}, {"2018"}, btag_catagories, {"350"}, 0.975, 1.013)
    ({"et","mt","tt"}, {"2018"}, btag_catagories, {"400"}, 0.979, 1.012)
    ({"et","mt","tt"}, {"2018"}, btag_catagories, {"450"}, 0.980, 1.015)
    ({"et","mt","tt"}, {"2018"}, btag_catagories, {"500"}, 0.982, 1.018)
    ({"et","mt","tt"}, {"2018"}, btag_catagories, {"600"}, 0.985, 1.024)
    ({"et","mt","tt"}, {"2018"}, btag_catagories, {"700"}, 0.981, 1.022)
    ({"et","mt","tt"}, {"2018"}, btag_catagories, {"800"}, 0.978, 1.020)
    ({"et","mt","tt"}, {"2018"}, btag_catagories, {"900"}, 0.974, 1.019)
    ({"et","mt","tt"}, {"2018"}, btag_catagories, {"1000"}, 0.970, 1.017)
    ({"et","mt","tt"}, {"2018"}, btag_catagories, {"1200"}, 0.970, 1.022)
    ({"et","mt","tt"}, {"2018"}, btag_catagories, {"1400"}, 0.971, 1.027)
    ({"et","mt","tt"}, {"2018"}, btag_catagories, {"1600"}, 0.968, 1.024)
    ({"et","mt","tt"}, {"2018"}, btag_catagories, {"1800"}, 0.965, 1.021)
    ({"et","mt","tt"}, {"2018"}, btag_catagories, {"2000"}, 0.962, 1.018)
    ({"et","mt","tt"}, {"2018"}, btag_catagories, {"2300"}, 0.962, 1.020)
    ({"et","mt","tt"}, {"2018"}, btag_catagories, {"2600"}, 0.962, 1.027)
    ({"et","mt","tt"}, {"2018"}, btag_catagories, {"2900"}, 0.961, 1.034)
    ({"et","mt","tt"}, {"2018"}, btag_catagories, {"3200"}, 0.957, 1.033)
    ({"et","mt","tt"}, {"2018"}, btag_catagories, {"3500"}, 0.953, 1.032)
    ({"em"}, {"2018"}, btag_catagories, {"60"}, 0.951, 1.041)
    ({"em"}, {"2018"}, btag_catagories, {"80"}, 0.935, 1.037)
    ({"em"}, {"2018"}, btag_catagories, {"100"}, 0.920, 1.033)
    ({"em"}, {"2018"}, btag_catagories, {"120"}, 0.963, 1.026)
    ({"em"}, {"2018"}, btag_catagories, {"125"}, 0.974, 1.024)
    ({"em"}, {"2018"}, btag_catagories, {"130"}, 0.973, 1.021)
    ({"em"}, {"2018"}, btag_catagories, {"140"}, 0.971, 1.015)
    ({"em"}, {"2018"}, btag_catagories, {"160"}, 0.963, 1.016)
    ({"em"}, {"2018"}, btag_catagories, {"180"}, 0.965, 1.015)
    ({"em"}, {"2018"}, btag_catagories, {"200"}, 0.968, 1.015)
    ({"em"}, {"2018"}, btag_catagories, {"250"}, 0.973, 1.013)
    ({"em"}, {"2018"}, btag_catagories, {"300"}, 0.979, 1.012)
    ({"em"}, {"2018"}, btag_catagories, {"350"}, 0.984, 1.010)
    ({"em"}, {"2018"}, btag_catagories, {"400"}, 0.990, 1.009)
    ({"em"}, {"2018"}, btag_catagories, {"450"}, 0.990, 1.012)
    ({"em"}, {"2018"}, btag_catagories, {"500"}, 0.990, 1.015)
    ({"em"}, {"2018"}, btag_catagories, {"600"}, 0.991, 1.021)
    ({"em"}, {"2018"}, btag_catagories, {"700"}, 0.988, 1.018)
    ({"em"}, {"2018"}, btag_catagories, {"800"}, 0.986, 1.015)
    ({"em"}, {"2018"}, btag_catagories, {"900"}, 0.984, 1.013)
    ({"em"}, {"2018"}, btag_catagories, {"1000"}, 0.981, 1.010)
    ({"em"}, {"2018"}, btag_catagories, {"1200"}, 0.980, 1.018)
    ({"em"}, {"2018"}, btag_catagories, {"1400"}, 0.979, 1.027)
    ({"em"}, {"2018"}, btag_catagories, {"1600"}, 0.979, 1.022)
    ({"em"}, {"2018"}, btag_catagories, {"1800"}, 0.979, 1.017)
    ({"em"}, {"2018"}, btag_catagories, {"2000"}, 0.979, 1.012)
    ({"em"}, {"2018"}, btag_catagories, {"2300"}, 0.971, 1.011)
    ({"em"}, {"2018"}, btag_catagories, {"2600"}, 0.968, 1.014)
    ({"em"}, {"2018"}, btag_catagories, {"2900"}, 0.966, 1.018)
    ({"em"}, {"2018"}, btag_catagories, {"3200"}, 0.964, 1.023)
    ({"em"}, {"2018"}, btag_catagories, {"3500"}, 0.963, 1.029));

  // ##########################################################################
  // Uncertainty: Lumi
  // References:
  // - "CMS Luminosity Measurements for the 2016 Data Taking Period"
  //   (PAS, https://cds.cern.ch/record/2257069)
  // - Recommendation twiki
  //    https://twiki.cern.ch/twiki/bin/view/CMS/TWikiLUM#LumiComb  
  // Notes:
  // - FIXME: Adapt for fake factor and embedding
  // ##########################################################################

  float lumi_unc = 1.0;
  float lumi_xy_fact = 1.0;
  float lumi_len_scale = 1.0;
  float lumi_beam_beam = 1.0;
  float lumi_dyn_beta = 1.0;
  float lumi_beam_curr = 1.0;
  float lumi_ghost = 1.0;
  if (era == 2016) {
      lumi_unc = 1.022;
      lumi_xy_fact = 1.009;
      lumi_len_scale = 1.000;
      lumi_beam_beam = 1.004;
      lumi_dyn_beta = 1.005;
      lumi_beam_curr = 1.000;
      lumi_ghost = 1.004;
  } else if (era == 2017) {
      lumi_unc = 1.020;
      lumi_xy_fact = 1.008;
      lumi_len_scale = 1.003;
      lumi_beam_beam = 1.004;
      lumi_dyn_beta = 1.005;
      lumi_beam_curr = 1.003;
      lumi_ghost = 1.001;
  } else if (era == 2018) {
      lumi_unc = 1.015;
      lumi_xy_fact = 1.020;
      lumi_len_scale = 1.002;
      lumi_beam_beam = 1.000;
      lumi_dyn_beta = 1.000;
      lumi_beam_curr = 1.002;
      lumi_ghost = 1.000;
  }
  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process(mc_processes)
      .AddSyst(cb, "lumi_13TeV_$ERA", "lnN", SystMap<>::init(lumi_unc));
  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process(mc_processes)
      .AddSyst(cb, "lumi_13TeV_XY", "lnN", SystMap<>::init(lumi_xy_fact));
  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process(mc_processes)
      .AddSyst(cb, "lumi_13TeV_LS", "lnN", SystMap<>::init(lumi_len_scale));
  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process(mc_processes)
      .AddSyst(cb, "lumi_13TeV_BBD", "lnN", SystMap<>::init(lumi_beam_beam));
  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process(mc_processes)
      .AddSyst(cb, "lumi_13TeV_DB", "lnN", SystMap<>::init(lumi_dyn_beta));
  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process(mc_processes)
      .AddSyst(cb, "lumi_13TeV_BCC", "lnN", SystMap<>::init(lumi_beam_curr));
  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process(mc_processes)
      .AddSyst(cb, "lumi_13TeV_GS", "lnN", SystMap<>::init(lumi_ghost));

  // ##########################################################################
  // Uncertainty: ggH Reweighting Hdamp uncertainty
  // References:
  // - 
  // Notes: for Hdamp scales t, b, and i components are decorrelated
  // ##########################################################################

   cb.cp()
       .channel({"et", "mt", "tt", "em"})
       .process({"ggh_i","ggH_i","ggA_i"})
       .AddSyst(cb, "Hdamp_ggH_i_REWEIGHT", "shape", SystMap<>::init(1.00));

   cb.cp()
       .channel({"et", "mt", "tt", "em"})
       .process({"ggh_t","ggH_t","ggA_t"})
       .AddSyst(cb, "Hdamp_ggH_t_REWEIGHT", "shape", SystMap<>::init(1.00));

   cb.cp()
       .channel({"et", "mt", "tt", "em"})
       .process({"ggh_b","ggH_b","ggA_b"})
       .AddSyst(cb, "Hdamp_ggH_b_REWEIGHT", "shape", SystMap<>::init(1.00));

  // ##########################################################################
  // Uncertainty: ggH Reweighting QCDscale uncertainty
  // References:
  // - 
  // Notes: t,b, and i are correlated in this case
  // ##########################################################################


  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process(mssm_ggH_signals)
      .AddSyst(cb, "QCDscale_ggH_REWEIGHT", "shape", SystMap<>::init(1.00));


  // ##########################################################################
  // Uncertainty: Prefiring
  // References:
  // - "https://twiki.cern.ch/twiki/bin/viewauth/CMS/L1ECALPrefiringWeightRecipe"
  // Notes:
  // - FIXME: assumed as uncorrelated accross the years for now, what is the recommendation?
  // ##########################################################################
  if (era != 2018) {
  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_prefiring", "shape", SystMap<>::init(1.00));
  }

  // ##########################################################################
  // Uncertainty: Trigger efficiency
  // References:
  // Notes:
  // - FIXME: References?
  // ##########################################################################

  cb.cp()
    .channel({"et"})
    .process(mc_processes)
    .AddSyst(cb, "CMS_eff_trigger_et_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
    .channel({"et"})
    .process(mc_processes)
    .AddSyst(cb, "CMS_eff_xtrigger_l_et_$ERA", "shape", SystMap<>::init(1.00));
  // 100% uncorrelated for embedded
  cb.cp()
    .channel({"et"})
    .process({"EMB"})
    .AddSyst(cb, "CMS_eff_trigger_emb_et_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
    .channel({"et"})
    .process({"EMB"})
    .AddSyst(cb, "CMS_eff_xtrigger_l_emb_et_$ERA", "shape", SystMap<>::init(1.00));

  cb.cp()
      .channel({"mt"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_eff_trigger_mt_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"mt"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_eff_xtrigger_l_mt_$ERA", "shape", SystMap<>::init(1.00));

  cb.cp()
      .channel({"em"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_eff_trigger_em_$ERA", "lnN", SystMap<>::init(1.02));

  // 100% uncorrelated for embedded
  cb.cp()
      .channel({"mt"})
      .process({"EMB"})
      .AddSyst(cb, "CMS_eff_trigger_emb_mt_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"mt"})
      .process({"EMB"})
      .AddSyst(cb, "CMS_eff_xtrigger_l_emb_mt_$ERA", "shape", SystMap<>::init(1.00));

  cb.cp()
      .channel({"em"})
      .process({"EMB"})
      .AddSyst(cb, "CMS_eff_trigger_emb_em_$ERA", "lnN", SystMap<>::init(1.02));

  // Tau trigger efficiencies implemented as shape uncertainties in all channels.
  std::string tauTriggerdmbins[4] = {"0", "1", "10", "11"};
  for (auto tauTriggerbin: tauTriggerdmbins)
  {
      // lt cross trigger
      cb.cp()
          .channel({"mt", "et"})
          .process(mc_processes)
          .AddSyst(cb, "CMS_eff_xtrigger_t_$CHANNEL_dm"+tauTriggerbin+"_$ERA", "shape", SystMap<>::init(1.00));
          
      cb.cp()
          .channel({"mt", "et"})
          .process({"EMB"})
          .AddSyst(cb, "CMS_eff_xtrigger_t_emb_$CHANNEL_dm"+tauTriggerbin+"_$ERA", "shape", SystMap<>::init(0.866));

      // Correlated component acting on Embedded
      cb.cp()
          .channel({"mt", "et"})
          .process({"EMB"})
          .AddSyst(cb, "CMS_eff_xtrigger_t_$CHANNEL_dm"+tauTriggerbin+"_$ERA", "shape", SystMap<>::init(0.5));
      

      // di-tau trigger
      cb.cp()
          .channel({"tt"})
          .process(mc_processes)
          .AddSyst(cb, "CMS_eff_xtrigger_t_tt_dm"+tauTriggerbin+"_$ERA", "shape", SystMap<>::init(1.00));

      cb.cp()
          .channel({"tt"})
          .process({"EMB"})
          .AddSyst(cb, "CMS_eff_xtrigger_t_emb_tt_dm"+tauTriggerbin+"_$ERA", "shape", SystMap<>::init(0.866));

      // Correlated component acting on Embedded
      cb.cp()
          .channel({"tt"})
          .process({"EMB"})
          .AddSyst(cb, "CMS_eff_xtrigger_t_tt_dm"+tauTriggerbin+"_$ERA", "shape", SystMap<>::init(0.5));
      if (!sm) {
          // di-tau trigger
          cb.cp()
              .channel({"tt"})
              .process(mc_processes)
              .AddSyst(cb, "CMS_eff_xtrigger_t_tt_dm"+tauTriggerbin+"_highpT_$ERA", "shape", SystMap<>::init(1.00));

          cb.cp()
              .channel({"tt"})
              .process({"EMB"})
              .AddSyst(cb, "CMS_eff_xtrigger_t_emb_tt_dm"+tauTriggerbin+"_highpT_$ERA", "shape", SystMap<>::init(0.866));

          // Correlated component acting on Embedded
          cb.cp()
              .channel({"tt"})
              .process({"EMB"})
              .AddSyst(cb, "CMS_eff_xtrigger_t_tt_dm"+tauTriggerbin+"_highpT_$ERA", "shape", SystMap<>::init(0.5));
      }
  }

 if (!sm) {
  // Single tau trigger
  cb.cp()
      .channel({"mt", "et", "tt"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_eff_trigger_single_t_$ERA", "shape", SystMap<>::init(1.00));

  cb.cp()
      .channel({"mt", "et", "tt"})
      .process({"EMB"})
      .AddSyst(cb, "CMS_eff_trigger_single_t_emb_$ERA", "shape", SystMap<>::init(0.866));
  
  // Correlated component acting on Embedded
  cb.cp()
      .channel({"mt", "et", "tt"})
      .process({"EMB"})
      .AddSyst(cb, "CMS_eff_trigger_single_t_$ERA", "shape", SystMap<>::init(0.5));
 }
  // ##########################################################################
  // Uncertainty: Electron, muon and tau ID efficiency
  // References:
  // Notes:
  // - FIXME: Adapt for fake factor and embedding
  // - FIXME: Handling of ZL in fully-hadronic channel?
  // - FIXME: References?
  // ##########################################################################

  // 3% in Tau ID SF with different anti-l fake WP
  cb.cp()
      .channel({"mt", "tt"})
      .process(JoinStr({signals, {"EMB", "ZTT", "TTT", "TTL", "VVT", "VVL"}}))
      .AddSyst(cb, "CMS_eff_t_wp_$ERA", "lnN", SystMap<>::init(1.03));

  std::string tauIDptbins[5] = {"30-35", "35-40", "40-500", "highpT_100-500", "highpT_500-inf"};
  if (sm) {
    tauIDptbins[3] = "500-1000";
    tauIDptbins[4] = "1000-inf";
  }
  std::string tauIDdmbins[4] = {"0", "1", "10", "11"};
  // Common component acting on MC
  
  // Electron ID
  cb.cp()
      .channel({"et", "em"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_eff_e", "lnN", SystMap<>::init(1.02));

  // Muon ID
  cb.cp()
      .channel({"mt", "em"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_eff_m", "lnN", SystMap<>::init(1.02));

  // Tau ID: et and mt with 1 real tau
      
  for (auto tauIDbin : tauIDptbins){ //first part correlated between channels for IDvsJets
    cb.cp()
        .channel({"et", "mt"})
        .process(JoinStr({signals, {"ZTT", "TTT", "TTL", "VVT", "VVL"}}))
        .AddSyst(cb, "CMS_eff_t_"+tauIDbin+"_$ERA", "shape", SystMap<>::init(1.0));
  }
  cb.cp() //second part uncorrelated between channels for IDvsLep
      .channel({"et", "mt"})
      .process(JoinStr({signals, {"ZTT", "TTT", "TTL", "VVT", "VVL"}}))
      .AddSyst(cb, "CMS_eff_t_$CHANNEL_$ERA", "lnN", SystMap<>::init(1.01));

  // Tau ID: tt with 2 real taus
  for (auto tauIDbin : tauIDdmbins){
    cb.cp()
        .channel({"tt"})
        .process(JoinStr({signals, {"ZTT", "TTT", "TTL", "VVT", "VVL"}}))
        .AddSyst(cb, "CMS_eff_t_dm"+tauIDbin+"_$ERA", "shape", SystMap<>::init(1.0));
  }
  if (!sm){
      cb.cp()
      .channel({"tt"})
      .process(JoinStr({signals, {"ZTT", "TTT", "TTL", "VVT", "VVL"}}))
      .AddSyst(cb, "CMS_eff_t_highpT_100-500_$ERA", "shape", SystMap<>::init(1.0));
      cb.cp()
        .channel({"tt"})
        .process(JoinStr({signals, {"ZTT", "TTT", "TTL", "VVT", "VVL"}}))
        .AddSyst(cb, "CMS_eff_t_highpT_500-inf_$ERA", "shape", SystMap<>::init(1.0));
  }
  cb.cp()
      .channel({"tt"})
      .process(JoinStr({signals, {"ZTT", "TTT", "TTL", "VVT", "VVL"}}))
      .AddSyst(cb, "CMS_eff_t_$CHANNEL_$ERA", "lnN", SystMap<>::init(1.014));

  // Component for EMB only

  // Electron ID
  cb.cp()
      .channel({"et", "em"})
      .process({"EMB"})
      .AddSyst(cb, "CMS_eff_e_emb", "lnN", SystMap<>::init(1.017));

  // Muon ID
  cb.cp()
      .channel({"mt", "em"})
      .process({"EMB"})
      .AddSyst(cb, "CMS_eff_m_emb", "lnN", SystMap<>::init(1.017));

  // Tau ID: et and mt with 1 real tau
  for (auto tauIDbin : tauIDptbins){
    cb.cp()
        .channel({"et", "mt"})
        .process({"EMB"})
        .AddSyst(cb, "CMS_eff_t_emb_"+tauIDbin+"_$ERA", "shape", SystMap<>::init(0.866));
  }
  cb.cp()
      .channel({"et", "mt"})
      .process({"EMB"})
      .AddSyst(cb, "CMS_eff_t_emb_$CHANNEL_$ERA", "lnN", SystMap<>::init(1.0087));

  // Tau ID: tt with 2 real taus
  for (auto tauIDbin : tauIDdmbins){
    cb.cp()
        .channel({"tt"})
        .process({"EMB"})
        .AddSyst(cb, "CMS_eff_t_emb_dm"+tauIDbin+"_$ERA", "shape", SystMap<>::init(0.866));
  }
  if (!sm){
    cb.cp()
        .channel({"tt"})
        .process({"EMB"})
        .AddSyst(cb, "CMS_eff_t_emb_highpT_100-500_$ERA", "shape", SystMap<>::init(1.0));
    cb.cp()
        .channel({"tt"})
        .process({"EMB"})
        .AddSyst(cb, "CMS_eff_t_emb_highpT_500-inf_$ERA", "shape", SystMap<>::init(1.0));
  }
  cb.cp()
      .channel({"tt"})
      .process({"EMB"})
      .AddSyst(cb, "CMS_eff_t_emb_$CHANNEL_$ERA", "lnN", SystMap<>::init(1.012));


  // Common NP acting on EMB

  // Electron ID
  cb.cp()
      .channel({"et", "em"})
      .process({"EMB"})
      .AddSyst(cb, "CMS_eff_e", "lnN", SystMap<>::init(1.01));

  // Muon ID
  cb.cp()
      .channel({"mt", "em"})
      .process({"EMB"})
      .AddSyst(cb, "CMS_eff_m", "lnN", SystMap<>::init(1.01));

  // Tau ID: et and mt with 1 real tau
  for (auto tauIDbin : tauIDptbins){
    cb.cp()
        .channel({"et", "mt"})
        .process({"EMB"})
        .AddSyst(cb, "CMS_eff_t_"+tauIDbin+"_$ERA", "shape", SystMap<>::init(0.5));
  }
  cb.cp()
      .channel({"et", "mt"})
      .process({"EMB"})
      .AddSyst(cb, "CMS_eff_t_$CHANNEL_$ERA", "lnN", SystMap<>::init(1.005));

  // Tau ID: tt with 2 real taus
  for (auto tauIDbin : tauIDdmbins){
    cb.cp()
        .channel({"tt"})
        .process({"EMB"})
        .AddSyst(cb, "CMS_eff_t_dm"+tauIDbin+"_$ERA", "shape", SystMap<>::init(0.5));
  }
  if (!sm){
    cb.cp()
        .channel({"tt"})
        .process({"EMB"})
        .AddSyst(cb, "CMS_eff_t_highpT_100-500_$ERA", "shape", SystMap<>::init(1.0));
    cb.cp()
        .channel({"tt"})
        .process({"EMB"})
        .AddSyst(cb, "CMS_eff_t_highpT_500-inf_$ERA", "shape", SystMap<>::init(1.0));
  }
  cb.cp()
      .channel({"tt"})
      .process({"EMB"})
      .AddSyst(cb, "CMS_eff_t_$CHANNEL_$ERA", "lnN", SystMap<>::init(1.007));



  // ##########################################################################
  // Uncertainty: b-tag and mistag efficiency
  // References:
  // Notes:
  // - FIXME: References?
  // ##########################################################################

  // SM btag uncertainties uses shape systematic 
  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .bin_id(mssm_categories, false)
      .process(mc_processes)
      .AddSyst(cb, "CMS_htt_eff_b_$ERA", "shape", SystMap<>::init(1.00));

  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .bin_id(mssm_categories, false)
      .process(mc_processes)
      .AddSyst(cb, "CMS_htt_mistag_b_$ERA", "shape", SystMap<>::init(1.00));

  // Classic MSSM categories btag uncertainties will be lnN
  cb.cp().process({"W"}).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id>::init
  ({"em"}, {"2016"}, {35}, 1.0, 1.0)
  ({"em"}, {"2017"}, {35}, 1.0, 1.024)
  ({"em"}, {"2018"}, {35}, 1.0, 1.069)
  ({"em"}, {"2016"}, {36}, 1.0, 1.0)
  ({"em"}, {"2017"}, {36}, 1.0, 1.004)
  ({"em"}, {"2018"}, {36}, 0.962, 1.0)
  ({"em"}, {"2016"}, {37}, 1.0, 1.0)
  ({"em"}, {"2017"}, {37}, 0.987, 1.015)
  ({"em"}, {"2018"}, {37}, 0.981, 1.015)
  ({"em"}, {"2016"}, {32}, 1.0, 1.0)
  ({"em"}, {"2017"}, {32}, 1.0, 0.996)
  ({"em"}, {"2018"}, {32}, 1.0, 0.99)
  ({"em"}, {"2016"}, {33}, 1.0, 1.0)
  ({"em"}, {"2017"}, {33}, 1.0, 1.0)
  ({"em"}, {"2018"}, {33}, 1.003, 1.0)
  ({"em"}, {"2016"}, {34}, 0.998, 0.998)
  ({"em"}, {"2017"}, {34}, 1.001, 0.999)
  ({"em"}, {"2018"}, {34}, 1.001, 0.999)
  ({"em"}, {"2016"}, {2}, 1.001, 1.001)
  ({"em"}, {"2017"}, {2}, 1.0, 1.0)
  ({"em"}, {"2018"}, {2}, 1.0, 1.0)
  );
  
  cb.cp().process({"ZL"}).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id>::init
  ({"em"}, {"2016"}, {35}, 1.0, 1.0)
  ({"em"}, {"2017"}, {35}, 0.993, 1.0)
  ({"em"}, {"2018"}, {35}, 0.952, 1.0)
  ({"em"}, {"2016"}, {36}, 0.99, 1.0)
  ({"em"}, {"2017"}, {36}, 0.993, 1.005)
  ({"em"}, {"2018"}, {36}, 0.985, 1.016)
  ({"em"}, {"2016"}, {37}, 0.987, 0.987)
  ({"em"}, {"2017"}, {37}, 0.983, 1.007)
  ({"em"}, {"2018"}, {37}, 0.988, 1.008)
  ({"em"}, {"2016"}, {32}, 1.0, 1.0)
  ({"em"}, {"2017"}, {32}, 1.001, 1.0)
  ({"em"}, {"2018"}, {32}, 1.006, 1.0)
  ({"em"}, {"2016"}, {33}, 0.999, 0.999)
  ({"em"}, {"2017"}, {33}, 1.0, 1.0)
  ({"em"}, {"2018"}, {33}, 1.001, 0.999)
  ({"em"}, {"2016"}, {34}, 1.001, 1.001)
  ({"em"}, {"2017"}, {34}, 1.002, 1.0)
  ({"em"}, {"2018"}, {34}, 1.0, 0.999)
  ({"em"}, {"2016"}, {2}, 1.002, 1.002)
  ({"em"}, {"2017"}, {2}, 0.999, 0.999)
  ({"em"}, {"2018"}, {2}, 1.001, 1.001)
  ({"et"}, {"2016"}, {35}, 0.992, 1.033)
  ({"et"}, {"2017"}, {35}, 0.966, 1.04)
  ({"et"}, {"2018"}, {35}, 0.968, 1.105)
  ({"et"}, {"2016"}, {36}, 1.0, 1.02)
  ({"et"}, {"2017"}, {36}, 0.96, 1.069)
  ({"et"}, {"2018"}, {36}, 0.983, 1.069)
  ({"et"}, {"2016"}, {32}, 1.0, 0.999)
  ({"et"}, {"2017"}, {32}, 1.003, 1.001)
  ({"et"}, {"2018"}, {32}, 1.001, 0.996)
  ({"et"}, {"2016"}, {33}, 1.0, 0.999)
  ({"et"}, {"2017"}, {33}, 1.006, 1.001)
  ({"et"}, {"2018"}, {33}, 1.001, 0.996)
  ({"mt"}, {"2016"}, {35}, 0.984, 1.008)
  ({"mt"}, {"2017"}, {35}, 0.97, 1.033)
  ({"mt"}, {"2018"}, {35}, 0.978, 1.04)
  ({"mt"}, {"2016"}, {36}, 0.978, 1.037)
  ({"mt"}, {"2017"}, {36}, 0.972, 1.062)
  ({"mt"}, {"2018"}, {36}, 0.99, 1.022)
  ({"mt"}, {"2016"}, {32}, 1.0, 1.0)
  ({"mt"}, {"2017"}, {32}, 1.003, 1.001)
  ({"mt"}, {"2018"}, {32}, 1.001, 0.999)
  ({"mt"}, {"2016"}, {33}, 1.001, 0.999)
  ({"mt"}, {"2017"}, {33}, 1.004, 0.999)
  ({"mt"}, {"2018"}, {33}, 1.001, 0.999)
  ({"tt"}, {"2016"}, {32}, 0.971, 1.017)
  ({"tt"}, {"2017"}, {32}, 1.004, 1.062)
  ({"tt"}, {"2018"}, {32}, 0.948, 1.038)
  ({"tt"}, {"2016"}, {35}, 1.002, 0.999)
  ({"tt"}, {"2017"}, {35}, 1.003, 1.0)
  ({"tt"}, {"2018"}, {35}, 1.004, 0.997)
  );
  
  cb.cp().process({"TTL"}).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id>::init
  ({"em"}, {"2016"}, {35}, 0.99, 1.009)
  ({"em"}, {"2017"}, {35}, 0.982, 1.016)
  ({"em"}, {"2018"}, {35}, 0.99, 1.009)
  ({"em"}, {"2016"}, {36}, 0.99, 1.009)
  ({"em"}, {"2017"}, {36}, 0.982, 1.016)
  ({"em"}, {"2018"}, {36}, 0.99, 1.009)
  ({"em"}, {"2016"}, {37}, 0.99, 1.009)
  ({"em"}, {"2017"}, {37}, 0.981, 1.017)
  ({"em"}, {"2018"}, {37}, 0.99, 1.009)
  ({"em"}, {"2016"}, {32}, 1.057, 0.944)
  ({"em"}, {"2017"}, {32}, 1.11, 0.89)
  ({"em"}, {"2018"}, {32}, 1.073, 0.926)
  ({"em"}, {"2016"}, {33}, 1.053, 0.947)
  ({"em"}, {"2017"}, {33}, 1.11, 0.895)
  ({"em"}, {"2018"}, {33}, 1.072, 0.926)
  ({"em"}, {"2016"}, {34}, 1.052, 0.947)
  ({"em"}, {"2017"}, {34}, 1.108, 0.892)
  ({"em"}, {"2018"}, {34}, 1.068, 0.931)
  ({"em"}, {"2016"}, {2}, 1.0, 1.0)
  ({"em"}, {"2017"}, {2}, 1.0, 1.0)
  ({"em"}, {"2018"}, {2}, 1.0, 1.0)
  ({"et"}, {"2016"}, {35}, 0.964, 0.984)
  ({"et"}, {"2017"}, {35}, 0.983, 1.018)
  ({"et"}, {"2018"}, {35}, 0.991, 1.009)
  ({"et"}, {"2016"}, {36}, 0.961, 0.98)
  ({"et"}, {"2017"}, {36}, 0.98, 1.018)
  ({"et"}, {"2018"}, {36}, 0.991, 1.009)
  ({"et"}, {"2016"}, {32}, 1.027, 0.917)
  ({"et"}, {"2017"}, {32}, 1.116, 0.882)
  ({"et"}, {"2018"}, {32}, 1.073, 0.929)
  ({"et"}, {"2016"}, {33}, 1.024, 0.925)
  ({"et"}, {"2017"}, {33}, 1.13, 0.883)
  ({"et"}, {"2018"}, {33}, 1.075, 0.926)
  ({"mt"}, {"2016"}, {35}, 0.965, 0.985)
  ({"mt"}, {"2017"}, {35}, 0.98, 1.018)
  ({"mt"}, {"2018"}, {35}, 0.99, 1.009)
  ({"mt"}, {"2016"}, {36}, 0.962, 0.983)
  ({"mt"}, {"2017"}, {36}, 0.982, 1.019)
  ({"mt"}, {"2018"}, {36}, 0.99, 1.009)
  ({"mt"}, {"2016"}, {32}, 1.032, 0.923)
  ({"mt"}, {"2017"}, {32}, 1.138, 0.876)
  ({"mt"}, {"2018"}, {32}, 1.081, 0.928)
  ({"mt"}, {"2016"}, {33}, 1.034, 0.925)
  ({"mt"}, {"2017"}, {33}, 1.117, 0.874)
  ({"mt"}, {"2018"}, {33}, 1.08, 0.927)
  ({"tt"}, {"2016"}, {32}, 0.968, 0.982)
  ({"tt"}, {"2017"}, {32}, 0.991, 1.015)
  ({"tt"}, {"2018"}, {32}, 0.989, 1.009)
  ({"tt"}, {"2016"}, {35}, 1.02, 0.95)
  ({"tt"}, {"2017"}, {35}, 1.061, 0.892)
  ({"tt"}, {"2018"}, {35}, 1.084, 0.936)
  );
  
  cb.cp().process({"VVL"}).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id>::init
  ({"em"}, {"2016"}, {35}, 0.986, 1.01)
  ({"em"}, {"2017"}, {35}, 0.976, 1.019)
  ({"em"}, {"2018"}, {35}, 0.987, 1.013)
  ({"em"}, {"2016"}, {36}, 0.984, 1.012)
  ({"em"}, {"2017"}, {36}, 0.975, 1.024)
  ({"em"}, {"2018"}, {36}, 0.986, 1.014)
  ({"em"}, {"2016"}, {37}, 0.986, 1.013)
  ({"em"}, {"2017"}, {37}, 0.974, 1.023)
  ({"em"}, {"2018"}, {37}, 0.986, 1.014)
  ({"em"}, {"2016"}, {32}, 1.01, 0.991)
  ({"em"}, {"2017"}, {32}, 1.019, 0.981)
  ({"em"}, {"2018"}, {32}, 1.011, 0.986)
  ({"em"}, {"2016"}, {33}, 1.006, 0.994)
  ({"em"}, {"2017"}, {33}, 1.011, 0.987)
  ({"em"}, {"2018"}, {33}, 1.007, 0.992)
  ({"em"}, {"2016"}, {34}, 1.004, 0.995)
  ({"em"}, {"2017"}, {34}, 1.009, 0.991)
  ({"em"}, {"2018"}, {34}, 1.005, 0.994)
  ({"em"}, {"2016"}, {2}, 1.0, 1.0)
  ({"em"}, {"2017"}, {2}, 1.0, 1.0)
  ({"em"}, {"2018"}, {2}, 1.0, 1.0)
  ({"et"}, {"2016"}, {35}, 0.986, 1.015)
  ({"et"}, {"2017"}, {35}, 0.972, 1.033)
  ({"et"}, {"2018"}, {35}, 0.985, 1.016)
  ({"et"}, {"2016"}, {36}, 0.982, 1.018)
  ({"et"}, {"2017"}, {36}, 0.975, 1.027)
  ({"et"}, {"2018"}, {36}, 0.985, 1.013)
  ({"et"}, {"2016"}, {32}, 1.005, 0.994)
  ({"et"}, {"2017"}, {32}, 1.019, 0.978)
  ({"et"}, {"2018"}, {32}, 1.01, 0.99)
  ({"et"}, {"2016"}, {33}, 1.006, 0.994)
  ({"et"}, {"2017"}, {33}, 1.015, 0.984)
  ({"et"}, {"2018"}, {33}, 1.009, 0.993)
  ({"mt"}, {"2016"}, {35}, 0.984, 1.016)
  ({"mt"}, {"2017"}, {35}, 0.974, 1.024)
  ({"mt"}, {"2018"}, {35}, 0.987, 1.017)
  ({"mt"}, {"2016"}, {36}, 0.984, 1.013)
  ({"mt"}, {"2017"}, {36}, 0.969, 1.028)
  ({"mt"}, {"2018"}, {36}, 0.983, 1.016)
  ({"mt"}, {"2016"}, {32}, 1.005, 0.995)
  ({"mt"}, {"2017"}, {32}, 1.017, 0.984)
  ({"mt"}, {"2018"}, {32}, 1.008, 0.989)
  ({"mt"}, {"2016"}, {33}, 1.005, 0.996)
  ({"mt"}, {"2017"}, {33}, 1.019, 0.983)
  ({"mt"}, {"2018"}, {33}, 1.01, 0.991)
  ({"tt"}, {"2016"}, {32}, 0.976, 1.006)
  ({"tt"}, {"2017"}, {32}, 0.987, 1.033)
  ({"tt"}, {"2018"}, {32}, 0.984, 1.0)
  ({"tt"}, {"2016"}, {35}, 1.01, 0.998)
  ({"tt"}, {"2017"}, {35}, 1.01, 0.975)
  ({"tt"}, {"2018"}, {35}, 1.011, 1.0)
  );
  
  cb.cp().process(mssm_bbH_signals).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"60"}, 0.981, 1.002)
  ({"em"}, {"2017"}, btag_catagories, {"60"}, 0.957, 1.031)
  ({"em"}, {"2018"}, btag_catagories, {"60"}, 0.974, 0.977)
  ({"et"}, {"2016"}, btag_catagories, {"60"}, 1.0, 1.0)
  ({"et"}, {"2017"}, btag_catagories, {"60"}, 0.979, 0.977)
  ({"et"}, {"2018"}, btag_catagories, {"60"}, 0.936, 1.063)
  ({"mt"}, {"2016"}, btag_catagories, {"60"}, 0.959, 1.023)
  ({"mt"}, {"2017"}, btag_catagories, {"60"}, 0.961, 1.053)
  ({"mt"}, {"2018"}, btag_catagories, {"60"}, 0.965, 1.009)
  ({"tt"}, {"2016"}, btag_catagories, {"60"}, 0.856, 0.921)
  ({"tt"}, {"2017"}, btag_catagories, {"60"}, 0.964, 1.021)
  ({"tt"}, {"2018"}, btag_catagories, {"60"}, 1.0, 0.919)
  ({"em"}, {"2016"}, nobtag_catagories, {"60"}, 1.011, 0.999)
  ({"em"}, {"2017"}, nobtag_catagories, {"60"}, 1.019, 0.986)
  ({"em"}, {"2018"}, nobtag_catagories, {"60"}, 1.019, 1.017)
  ({"et"}, {"2016"}, nobtag_catagories, {"60"}, 1.0, 1.0)
  ({"et"}, {"2017"}, nobtag_catagories, {"60"}, 1.095, 1.105)
  ({"et"}, {"2018"}, nobtag_catagories, {"60"}, 1.132, 0.871)
  ({"mt"}, {"2016"}, nobtag_catagories, {"60"}, 1.08, 0.955)
  ({"mt"}, {"2017"}, nobtag_catagories, {"60"}, 1.065, 0.937)
  ({"mt"}, {"2018"}, nobtag_catagories, {"60"}, 1.038, 0.99)
  ({"tt"}, {"2016"}, nobtag_catagories, {"60"}, 1.496, 1.272)
  ({"tt"}, {"2017"}, nobtag_catagories, {"60"}, 1.107, 0.936)
  ({"tt"}, {"2018"}, nobtag_catagories, {"60"}, 1.0, 1.668)
  );
  
  cb.cp().process(mssm_bbH_signals).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"80"}, 0.969, 1.066)
  ({"em"}, {"2017"}, btag_catagories, {"80"}, 1.005, 1.014)
  ({"em"}, {"2018"}, btag_catagories, {"80"}, 0.989, 1.025)
  ({"et"}, {"2016"}, btag_catagories, {"80"}, 0.994, 1.019)
  ({"et"}, {"2017"}, btag_catagories, {"80"}, 0.961, 1.026)
  ({"et"}, {"2018"}, btag_catagories, {"80"}, 0.969, 1.048)
  ({"mt"}, {"2016"}, btag_catagories, {"80"}, 0.969, 1.023)
  ({"mt"}, {"2017"}, btag_catagories, {"80"}, 0.976, 1.033)
  ({"mt"}, {"2018"}, btag_catagories, {"80"}, 0.967, 1.03)
  ({"tt"}, {"2016"}, btag_catagories, {"80"}, 0.942, 1.023)
  ({"tt"}, {"2017"}, btag_catagories, {"80"}, 0.969, 0.992)
  ({"tt"}, {"2018"}, btag_catagories, {"80"}, 0.959, 1.012)
  ({"em"}, {"2016"}, nobtag_catagories, {"80"}, 1.007, 0.985)
  ({"em"}, {"2017"}, nobtag_catagories, {"80"}, 0.998, 0.994)
  ({"em"}, {"2018"}, nobtag_catagories, {"80"}, 1.004, 0.994)
  ({"et"}, {"2016"}, nobtag_catagories, {"80"}, 1.004, 0.989)
  ({"et"}, {"2017"}, nobtag_catagories, {"80"}, 1.025, 0.983)
  ({"et"}, {"2018"}, nobtag_catagories, {"80"}, 1.023, 0.965)
  ({"mt"}, {"2016"}, nobtag_catagories, {"80"}, 1.009, 0.993)
  ({"mt"}, {"2017"}, nobtag_catagories, {"80"}, 1.009, 0.987)
  ({"mt"}, {"2018"}, nobtag_catagories, {"80"}, 1.015, 0.986)
  ({"tt"}, {"2016"}, nobtag_catagories, {"80"}, 1.247, 0.9)
  ({"tt"}, {"2017"}, nobtag_catagories, {"80"}, 1.075, 1.019)
  ({"tt"}, {"2018"}, nobtag_catagories, {"80"}, 1.267, 0.92)
  );
  
  cb.cp().process(mssm_bbH_signals).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"100"}, 0.979, 1.015)
  ({"em"}, {"2017"}, btag_catagories, {"100"}, 0.969, 1.03)
  ({"em"}, {"2018"}, btag_catagories, {"100"}, 0.966, 1.036)
  ({"et"}, {"2016"}, btag_catagories, {"100"}, 0.984, 1.028)
  ({"et"}, {"2017"}, btag_catagories, {"100"}, 0.972, 1.043)
  ({"et"}, {"2018"}, btag_catagories, {"100"}, 0.979, 1.032)
  ({"mt"}, {"2016"}, btag_catagories, {"100"}, 0.973, 1.012)
  ({"mt"}, {"2017"}, btag_catagories, {"100"}, 0.939, 1.026)
  ({"mt"}, {"2018"}, btag_catagories, {"100"}, 0.974, 1.055)
  ({"tt"}, {"2016"}, btag_catagories, {"100"}, 0.998, 1.006)
  ({"tt"}, {"2017"}, btag_catagories, {"100"}, 0.979, 1.05)
  ({"tt"}, {"2018"}, btag_catagories, {"100"}, 0.957, 1.013)
  ({"em"}, {"2016"}, nobtag_catagories, {"100"}, 1.005, 0.996)
  ({"em"}, {"2017"}, nobtag_catagories, {"100"}, 1.01, 0.991)
  ({"em"}, {"2018"}, nobtag_catagories, {"100"}, 1.011, 0.988)
  ({"et"}, {"2016"}, nobtag_catagories, {"100"}, 1.005, 0.991)
  ({"et"}, {"2017"}, nobtag_catagories, {"100"}, 1.01, 0.988)
  ({"et"}, {"2018"}, nobtag_catagories, {"100"}, 1.006, 0.988)
  ({"mt"}, {"2016"}, nobtag_catagories, {"100"}, 1.007, 0.997)
  ({"mt"}, {"2017"}, nobtag_catagories, {"100"}, 1.019, 0.992)
  ({"mt"}, {"2018"}, nobtag_catagories, {"100"}, 1.008, 0.981)
  ({"tt"}, {"2016"}, nobtag_catagories, {"100"}, 1.002, 0.995)
  ({"tt"}, {"2017"}, nobtag_catagories, {"100"}, 1.015, 0.964)
  ({"tt"}, {"2018"}, nobtag_catagories, {"100"}, 1.038, 0.988)
  );
  
  cb.cp().process(mssm_bbH_signals).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"120"}, 0.974, 1.021)
  ({"em"}, {"2017"}, btag_catagories, {"120"}, 0.966, 1.067)
  ({"em"}, {"2018"}, btag_catagories, {"120"}, 0.979, 1.032)
  ({"et"}, {"2016"}, btag_catagories, {"120"}, 0.976, 1.008)
  ({"et"}, {"2017"}, btag_catagories, {"120"}, 0.972, 1.043)
  ({"et"}, {"2018"}, btag_catagories, {"120"}, 0.965, 1.033)
  ({"mt"}, {"2016"}, btag_catagories, {"120"}, 0.972, 1.027)
  ({"mt"}, {"2017"}, btag_catagories, {"120"}, 0.957, 1.043)
  ({"mt"}, {"2018"}, btag_catagories, {"120"}, 0.971, 1.036)
  ({"tt"}, {"2016"}, btag_catagories, {"120"}, 0.982, 1.019)
  ({"tt"}, {"2017"}, btag_catagories, {"120"}, 0.95, 1.032)
  ({"tt"}, {"2018"}, btag_catagories, {"120"}, 0.98, 1.034)
  ({"em"}, {"2016"}, nobtag_catagories, {"120"}, 1.006, 0.995)
  ({"em"}, {"2017"}, nobtag_catagories, {"120"}, 1.01, 0.98)
  ({"em"}, {"2018"}, nobtag_catagories, {"120"}, 1.007, 0.99)
  ({"et"}, {"2016"}, nobtag_catagories, {"120"}, 1.006, 0.998)
  ({"et"}, {"2017"}, nobtag_catagories, {"120"}, 1.009, 0.985)
  ({"et"}, {"2018"}, nobtag_catagories, {"120"}, 1.014, 0.988)
  ({"mt"}, {"2016"}, nobtag_catagories, {"120"}, 1.007, 0.993)
  ({"mt"}, {"2017"}, nobtag_catagories, {"120"}, 1.014, 0.986)
  ({"mt"}, {"2018"}, nobtag_catagories, {"120"}, 1.011, 0.986)
  ({"tt"}, {"2016"}, nobtag_catagories, {"120"}, 1.006, 0.993)
  ({"tt"}, {"2017"}, nobtag_catagories, {"120"}, 1.019, 0.988)
  ({"tt"}, {"2018"}, nobtag_catagories, {"120"}, 1.009, 0.986)
  );
  
  cb.cp().process(mssm_bbH_signals).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"125"}, 0.965, 1.029)
  ({"em"}, {"2017"}, btag_catagories, {"125"}, 0.96, 1.038)
  ({"em"}, {"2018"}, btag_catagories, {"125"}, 0.974, 1.021)
  ({"et"}, {"2016"}, btag_catagories, {"125"}, 0.973, 1.025)
  ({"et"}, {"2017"}, btag_catagories, {"125"}, 0.961, 1.024)
  ({"et"}, {"2018"}, btag_catagories, {"125"}, 0.966, 1.024)
  ({"mt"}, {"2016"}, btag_catagories, {"125"}, 0.965, 1.022)
  ({"mt"}, {"2017"}, btag_catagories, {"125"}, 0.957, 1.055)
  ({"mt"}, {"2018"}, btag_catagories, {"125"}, 0.966, 1.027)
  ({"tt"}, {"2016"}, btag_catagories, {"125"}, 0.971, 1.023)
  ({"tt"}, {"2017"}, btag_catagories, {"125"}, 0.973, 1.051)
  ({"tt"}, {"2018"}, btag_catagories, {"125"}, 0.983, 1.029)
  ({"em"}, {"2016"}, nobtag_catagories, {"125"}, 1.008, 0.994)
  ({"em"}, {"2017"}, nobtag_catagories, {"125"}, 1.012, 0.988)
  ({"em"}, {"2018"}, nobtag_catagories, {"125"}, 1.009, 0.993)
  ({"et"}, {"2016"}, nobtag_catagories, {"125"}, 1.008, 0.992)
  ({"et"}, {"2017"}, nobtag_catagories, {"125"}, 1.014, 0.993)
  ({"et"}, {"2018"}, nobtag_catagories, {"125"}, 1.014, 0.99)
  ({"mt"}, {"2016"}, nobtag_catagories, {"125"}, 1.009, 0.994)
  ({"mt"}, {"2017"}, nobtag_catagories, {"125"}, 1.015, 0.982)
  ({"mt"}, {"2018"}, nobtag_catagories, {"125"}, 1.012, 0.99)
  ({"tt"}, {"2016"}, nobtag_catagories, {"125"}, 1.01, 0.992)
  ({"tt"}, {"2017"}, nobtag_catagories, {"125"}, 1.011, 0.979)
  ({"tt"}, {"2018"}, nobtag_catagories, {"125"}, 1.007, 0.987)
  );
  
  cb.cp().process(mssm_bbH_signals).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"130"}, 0.97, 1.024)
  ({"em"}, {"2017"}, btag_catagories, {"130"}, 0.943, 1.033)
  ({"em"}, {"2018"}, btag_catagories, {"130"}, 0.975, 1.024)
  ({"et"}, {"2016"}, btag_catagories, {"130"}, 0.974, 1.029)
  ({"et"}, {"2017"}, btag_catagories, {"130"}, 0.943, 1.063)
  ({"et"}, {"2018"}, btag_catagories, {"130"}, 0.968, 1.031)
  ({"mt"}, {"2016"}, btag_catagories, {"130"}, 0.976, 1.021)
  ({"mt"}, {"2017"}, btag_catagories, {"130"}, 0.963, 1.039)
  ({"mt"}, {"2018"}, btag_catagories, {"130"}, 0.961, 1.028)
  ({"tt"}, {"2016"}, btag_catagories, {"130"}, 0.976, 1.03)
  ({"tt"}, {"2017"}, btag_catagories, {"130"}, 0.962, 1.044)
  ({"tt"}, {"2018"}, btag_catagories, {"130"}, 0.97, 1.037)
  ({"em"}, {"2016"}, nobtag_catagories, {"130"}, 1.008, 0.993)
  ({"em"}, {"2017"}, nobtag_catagories, {"130"}, 1.018, 0.99)
  ({"em"}, {"2018"}, nobtag_catagories, {"130"}, 1.008, 0.992)
  ({"et"}, {"2016"}, nobtag_catagories, {"130"}, 1.007, 0.992)
  ({"et"}, {"2017"}, nobtag_catagories, {"130"}, 1.018, 0.979)
  ({"et"}, {"2018"}, nobtag_catagories, {"130"}, 1.012, 0.988)
  ({"mt"}, {"2016"}, nobtag_catagories, {"130"}, 1.006, 0.994)
  ({"mt"}, {"2017"}, nobtag_catagories, {"130"}, 1.013, 0.988)
  ({"mt"}, {"2018"}, nobtag_catagories, {"130"}, 1.015, 0.989)
  ({"tt"}, {"2016"}, nobtag_catagories, {"130"}, 1.008, 0.99)
  ({"tt"}, {"2017"}, nobtag_catagories, {"130"}, 1.016, 0.981)
  ({"tt"}, {"2018"}, nobtag_catagories, {"130"}, 1.014, 0.983)
  );
  
  cb.cp().process(mssm_bbH_signals).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"140"}, 0.975, 1.022)
  ({"em"}, {"2017"}, btag_catagories, {"140"}, 0.964, 1.033)
  ({"em"}, {"2018"}, btag_catagories, {"140"}, 0.97, 1.041)
  ({"et"}, {"2016"}, btag_catagories, {"140"}, 0.972, 1.019)
  ({"et"}, {"2017"}, btag_catagories, {"140"}, 0.957, 1.035)
  ({"et"}, {"2018"}, btag_catagories, {"140"}, 0.97, 1.023)
  ({"mt"}, {"2016"}, btag_catagories, {"140"}, 0.97, 1.026)
  ({"mt"}, {"2017"}, btag_catagories, {"140"}, 0.952, 1.05)
  ({"mt"}, {"2018"}, btag_catagories, {"140"}, 0.964, 1.041)
  ({"tt"}, {"2016"}, btag_catagories, {"140"}, 0.978, 1.028)
  ({"tt"}, {"2017"}, btag_catagories, {"140"}, 0.972, 1.039)
  ({"tt"}, {"2018"}, btag_catagories, {"140"}, 0.975, 1.027)
  ({"em"}, {"2016"}, nobtag_catagories, {"140"}, 1.007, 0.994)
  ({"em"}, {"2017"}, nobtag_catagories, {"140"}, 1.012, 0.989)
  ({"em"}, {"2018"}, nobtag_catagories, {"140"}, 1.011, 0.985)
  ({"et"}, {"2016"}, nobtag_catagories, {"140"}, 1.007, 0.995)
  ({"et"}, {"2017"}, nobtag_catagories, {"140"}, 1.016, 0.987)
  ({"et"}, {"2018"}, nobtag_catagories, {"140"}, 1.011, 0.991)
  ({"mt"}, {"2016"}, nobtag_catagories, {"140"}, 1.008, 0.992)
  ({"mt"}, {"2017"}, nobtag_catagories, {"140"}, 1.017, 0.983)
  ({"mt"}, {"2018"}, nobtag_catagories, {"140"}, 1.014, 0.984)
  ({"tt"}, {"2016"}, nobtag_catagories, {"140"}, 1.007, 0.991)
  ({"tt"}, {"2017"}, nobtag_catagories, {"140"}, 1.012, 0.983)
  ({"tt"}, {"2018"}, nobtag_catagories, {"140"}, 1.01, 0.989)
  );
  
  cb.cp().process(mssm_bbH_signals).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"160"}, 0.972, 1.022)
  ({"em"}, {"2017"}, btag_catagories, {"160"}, 0.951, 1.033)
  ({"em"}, {"2018"}, btag_catagories, {"160"}, 0.985, 1.026)
  ({"et"}, {"2016"}, btag_catagories, {"160"}, 0.974, 1.033)
  ({"et"}, {"2017"}, btag_catagories, {"160"}, 0.966, 1.046)
  ({"et"}, {"2018"}, btag_catagories, {"160"}, 0.97, 1.037)
  ({"mt"}, {"2016"}, btag_catagories, {"160"}, 0.98, 1.017)
  ({"mt"}, {"2017"}, btag_catagories, {"160"}, 0.959, 1.046)
  ({"mt"}, {"2018"}, btag_catagories, {"160"}, 0.974, 1.03)
  ({"tt"}, {"2016"}, btag_catagories, {"160"}, 0.967, 1.022)
  ({"tt"}, {"2017"}, btag_catagories, {"160"}, 0.956, 1.05)
  ({"tt"}, {"2018"}, btag_catagories, {"160"}, 0.968, 1.041)
  ({"em"}, {"2016"}, nobtag_catagories, {"160"}, 1.008, 0.994)
  ({"em"}, {"2017"}, nobtag_catagories, {"160"}, 1.018, 0.988)
  ({"em"}, {"2018"}, nobtag_catagories, {"160"}, 1.006, 0.99)
  ({"et"}, {"2016"}, nobtag_catagories, {"160"}, 1.008, 0.99)
  ({"et"}, {"2017"}, nobtag_catagories, {"160"}, 1.013, 0.984)
  ({"et"}, {"2018"}, nobtag_catagories, {"160"}, 1.011, 0.985)
  ({"mt"}, {"2016"}, nobtag_catagories, {"160"}, 1.005, 0.995)
  ({"mt"}, {"2017"}, nobtag_catagories, {"160"}, 1.015, 0.983)
  ({"mt"}, {"2018"}, nobtag_catagories, {"160"}, 1.011, 0.988)
  ({"tt"}, {"2016"}, nobtag_catagories, {"160"}, 1.01, 0.993)
  ({"tt"}, {"2017"}, nobtag_catagories, {"160"}, 1.017, 0.981)
  ({"tt"}, {"2018"}, nobtag_catagories, {"160"}, 1.014, 0.982)
  );
  
  cb.cp().process(mssm_bbH_signals).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"180"}, 0.972, 1.026)
  ({"em"}, {"2017"}, btag_catagories, {"180"}, 0.965, 1.032)
  ({"em"}, {"2018"}, btag_catagories, {"180"}, 0.98, 1.033)
  ({"et"}, {"2016"}, btag_catagories, {"180"}, 0.977, 1.028)
  ({"et"}, {"2017"}, btag_catagories, {"180"}, 0.959, 1.038)
  ({"et"}, {"2018"}, btag_catagories, {"180"}, 0.979, 1.032)
  ({"mt"}, {"2016"}, btag_catagories, {"180"}, 0.973, 1.02)
  ({"mt"}, {"2017"}, btag_catagories, {"180"}, 0.959, 1.036)
  ({"mt"}, {"2018"}, btag_catagories, {"180"}, 0.97, 1.03)
  ({"tt"}, {"2016"}, btag_catagories, {"180"}, 0.981, 1.02)
  ({"tt"}, {"2017"}, btag_catagories, {"180"}, 0.955, 1.045)
  ({"tt"}, {"2018"}, btag_catagories, {"180"}, 0.979, 1.038)
  ({"em"}, {"2016"}, nobtag_catagories, {"180"}, 1.009, 0.992)
  ({"em"}, {"2017"}, nobtag_catagories, {"180"}, 1.014, 0.987)
  ({"em"}, {"2018"}, nobtag_catagories, {"180"}, 1.008, 0.987)
  ({"et"}, {"2016"}, nobtag_catagories, {"180"}, 1.007, 0.991)
  ({"et"}, {"2017"}, nobtag_catagories, {"180"}, 1.018, 0.984)
  ({"et"}, {"2018"}, nobtag_catagories, {"180"}, 1.009, 0.986)
  ({"mt"}, {"2016"}, nobtag_catagories, {"180"}, 1.009, 0.994)
  ({"mt"}, {"2017"}, nobtag_catagories, {"180"}, 1.016, 0.986)
  ({"mt"}, {"2018"}, nobtag_catagories, {"180"}, 1.012, 0.987)
  ({"tt"}, {"2016"}, nobtag_catagories, {"180"}, 1.006, 0.993)
  ({"tt"}, {"2017"}, nobtag_catagories, {"180"}, 1.02, 0.98)
  ({"tt"}, {"2018"}, nobtag_catagories, {"180"}, 1.009, 0.983)
  );
  
  cb.cp().process(mssm_bbH_signals).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"200"}, 0.981, 1.013)
  ({"em"}, {"2017"}, btag_catagories, {"200"}, 0.957, 1.047)
  ({"em"}, {"2018"}, btag_catagories, {"200"}, 0.971, 1.038)
  ({"et"}, {"2016"}, btag_catagories, {"200"}, 0.742, 0.775)
  ({"et"}, {"2017"}, btag_catagories, {"200"}, 0.744, 0.812)
  ({"et"}, {"2018"}, btag_catagories, {"200"}, 0.972, 1.035)
  ({"mt"}, {"2016"}, btag_catagories, {"200"}, 0.772, 0.805)
  ({"mt"}, {"2017"}, btag_catagories, {"200"}, 0.764, 0.823)
  ({"mt"}, {"2018"}, btag_catagories, {"200"}, 0.972, 1.03)
  ({"tt"}, {"2016"}, btag_catagories, {"200"}, 0.619, 0.65)
  ({"tt"}, {"2017"}, btag_catagories, {"200"}, 0.684, 0.742)
  ({"tt"}, {"2018"}, btag_catagories, {"200"}, 0.968, 1.034)
  ({"em"}, {"2016"}, nobtag_catagories, {"200"}, 1.006, 0.996)
  ({"em"}, {"2017"}, nobtag_catagories, {"200"}, 1.015, 0.984)
  ({"em"}, {"2018"}, nobtag_catagories, {"200"}, 1.012, 0.984)
  ({"et"}, {"2016"}, nobtag_catagories, {"200"}, 0.816, 0.805)
  ({"et"}, {"2017"}, nobtag_catagories, {"200"}, 0.875, 0.845)
  ({"et"}, {"2018"}, nobtag_catagories, {"200"}, 1.012, 0.985)
  ({"mt"}, {"2016"}, nobtag_catagories, {"200"}, 0.871, 0.859)
  ({"mt"}, {"2017"}, nobtag_catagories, {"200"}, 0.882, 0.857)
  ({"mt"}, {"2018"}, nobtag_catagories, {"200"}, 1.012, 0.987)
  ({"tt"}, {"2016"}, nobtag_catagories, {"200"}, 0.764, 0.751)
  ({"tt"}, {"2017"}, nobtag_catagories, {"200"}, 0.783, 0.756)
  ({"tt"}, {"2018"}, nobtag_catagories, {"200"}, 1.015, 0.984)
  );
  
  cb.cp().process(mssm_bbH_signals).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"250"}, 0.976, 1.028)
  ({"em"}, {"2017"}, btag_catagories, {"250"}, 0.966, 1.047)
  ({"em"}, {"2018"}, btag_catagories, {"250"}, 0.977, 1.032)
  ({"et"}, {"2016"}, btag_catagories, {"250"}, 0.98, 1.025)
  ({"et"}, {"2017"}, btag_catagories, {"250"}, 0.968, 1.04)
  ({"et"}, {"2018"}, btag_catagories, {"250"}, 0.97, 1.034)
  ({"mt"}, {"2016"}, btag_catagories, {"250"}, 0.979, 1.023)
  ({"mt"}, {"2017"}, btag_catagories, {"250"}, 0.95, 1.041)
  ({"mt"}, {"2018"}, btag_catagories, {"250"}, 0.975, 1.028)
  ({"tt"}, {"2016"}, btag_catagories, {"250"}, 0.973, 1.02)
  ({"tt"}, {"2017"}, btag_catagories, {"250"}, 0.957, 1.041)
  ({"tt"}, {"2018"}, btag_catagories, {"250"}, 0.975, 1.032)
  ({"em"}, {"2016"}, nobtag_catagories, {"250"}, 1.008, 0.99)
  ({"em"}, {"2017"}, nobtag_catagories, {"250"}, 1.015, 0.979)
  ({"em"}, {"2018"}, nobtag_catagories, {"250"}, 1.011, 0.985)
  ({"et"}, {"2016"}, nobtag_catagories, {"250"}, 1.008, 0.99)
  ({"et"}, {"2017"}, nobtag_catagories, {"250"}, 1.016, 0.981)
  ({"et"}, {"2018"}, nobtag_catagories, {"250"}, 1.015, 0.983)
  ({"mt"}, {"2016"}, nobtag_catagories, {"250"}, 1.008, 0.991)
  ({"mt"}, {"2017"}, nobtag_catagories, {"250"}, 1.022, 0.982)
  ({"mt"}, {"2018"}, nobtag_catagories, {"250"}, 1.013, 0.986)
  ({"tt"}, {"2016"}, nobtag_catagories, {"250"}, 1.011, 0.992)
  ({"tt"}, {"2017"}, nobtag_catagories, {"250"}, 1.022, 0.979)
  ({"tt"}, {"2018"}, nobtag_catagories, {"250"}, 1.013, 0.983)
  );
  
  cb.cp().process(mssm_bbH_signals).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"300"}, 1.0, 1.0)
  ({"em"}, {"2017"}, btag_catagories, {"300"}, 0.957, 1.037)
  ({"em"}, {"2018"}, btag_catagories, {"300"}, 0.98, 1.026)
  ({"et"}, {"2016"}, btag_catagories, {"300"}, 1.0, 1.0)
  ({"et"}, {"2017"}, btag_catagories, {"300"}, 0.965, 1.043)
  ({"et"}, {"2018"}, btag_catagories, {"300"}, 0.978, 1.025)
  ({"mt"}, {"2016"}, btag_catagories, {"300"}, 1.0, 1.0)
  ({"mt"}, {"2017"}, btag_catagories, {"300"}, 0.966, 1.036)
  ({"mt"}, {"2018"}, btag_catagories, {"300"}, 0.977, 1.026)
  ({"tt"}, {"2016"}, btag_catagories, {"300"}, 1.0, 1.0)
  ({"tt"}, {"2017"}, btag_catagories, {"300"}, 0.961, 1.039)
  ({"tt"}, {"2018"}, btag_catagories, {"300"}, 0.976, 1.024)
  ({"em"}, {"2016"}, nobtag_catagories, {"300"}, 1.0, 1.0)
  ({"em"}, {"2017"}, nobtag_catagories, {"300"}, 1.022, 0.981)
  ({"em"}, {"2018"}, nobtag_catagories, {"300"}, 1.01, 0.988)
  ({"et"}, {"2016"}, nobtag_catagories, {"300"}, 1.0, 1.0)
  ({"et"}, {"2017"}, nobtag_catagories, {"300"}, 1.018, 0.979)
  ({"et"}, {"2018"}, nobtag_catagories, {"300"}, 1.012, 0.986)
  ({"mt"}, {"2016"}, nobtag_catagories, {"300"}, 1.0, 1.0)
  ({"mt"}, {"2017"}, nobtag_catagories, {"300"}, 1.017, 0.983)
  ({"mt"}, {"2018"}, nobtag_catagories, {"300"}, 1.013, 0.985)
  ({"tt"}, {"2016"}, nobtag_catagories, {"300"}, 1.0, 1.0)
  ({"tt"}, {"2017"}, nobtag_catagories, {"300"}, 1.021, 0.979)
  ({"tt"}, {"2018"}, nobtag_catagories, {"300"}, 1.014, 0.986)
  );
  
  cb.cp().process(mssm_bbH_signals).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"350"}, 0.98, 1.025)
  ({"em"}, {"2017"}, btag_catagories, {"350"}, 0.967, 1.045)
  ({"em"}, {"2018"}, btag_catagories, {"350"}, 0.986, 1.028)
  ({"et"}, {"2016"}, btag_catagories, {"350"}, 0.98, 1.021)
  ({"et"}, {"2017"}, btag_catagories, {"350"}, 0.957, 1.041)
  ({"et"}, {"2018"}, btag_catagories, {"350"}, 0.977, 1.025)
  ({"mt"}, {"2016"}, btag_catagories, {"350"}, 0.974, 1.019)
  ({"mt"}, {"2017"}, btag_catagories, {"350"}, 0.96, 1.039)
  ({"mt"}, {"2018"}, btag_catagories, {"350"}, 0.976, 1.027)
  ({"tt"}, {"2016"}, btag_catagories, {"350"}, 0.981, 1.022)
  ({"tt"}, {"2017"}, btag_catagories, {"350"}, 0.962, 1.046)
  ({"tt"}, {"2018"}, btag_catagories, {"350"}, 0.976, 1.027)
  ({"em"}, {"2016"}, nobtag_catagories, {"350"}, 1.008, 0.99)
  ({"em"}, {"2017"}, nobtag_catagories, {"350"}, 1.017, 0.977)
  ({"em"}, {"2018"}, nobtag_catagories, {"350"}, 1.008, 0.984)
  ({"et"}, {"2016"}, nobtag_catagories, {"350"}, 1.009, 0.99)
  ({"et"}, {"2017"}, nobtag_catagories, {"350"}, 1.025, 0.977)
  ({"et"}, {"2018"}, nobtag_catagories, {"350"}, 1.014, 0.985)
  ({"mt"}, {"2016"}, nobtag_catagories, {"350"}, 1.011, 0.991)
  ({"mt"}, {"2017"}, nobtag_catagories, {"350"}, 1.022, 0.979)
  ({"mt"}, {"2018"}, nobtag_catagories, {"350"}, 1.014, 0.984)
  ({"tt"}, {"2016"}, nobtag_catagories, {"350"}, 1.01, 0.989)
  ({"tt"}, {"2017"}, nobtag_catagories, {"350"}, 1.022, 0.973)
  ({"tt"}, {"2018"}, nobtag_catagories, {"350"}, 1.016, 0.982)
  );
  
  cb.cp().process(mssm_bbH_signals).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"400"}, 0.979, 1.021)
  ({"em"}, {"2017"}, btag_catagories, {"400"}, 0.967, 1.034)
  ({"em"}, {"2018"}, btag_catagories, {"400"}, 0.973, 1.027)
  ({"et"}, {"2016"}, btag_catagories, {"400"}, 0.984, 1.022)
  ({"et"}, {"2017"}, btag_catagories, {"400"}, 0.959, 1.034)
  ({"et"}, {"2018"}, btag_catagories, {"400"}, 0.98, 1.025)
  ({"mt"}, {"2016"}, btag_catagories, {"400"}, 0.982, 1.018)
  ({"mt"}, {"2017"}, btag_catagories, {"400"}, 0.961, 1.028)
  ({"mt"}, {"2018"}, btag_catagories, {"400"}, 0.972, 1.028)
  ({"tt"}, {"2016"}, btag_catagories, {"400"}, 0.975, 1.019)
  ({"tt"}, {"2017"}, btag_catagories, {"400"}, 0.965, 1.036)
  ({"tt"}, {"2018"}, btag_catagories, {"400"}, 0.978, 1.026)
  ({"em"}, {"2016"}, nobtag_catagories, {"400"}, 1.009, 0.99)
  ({"em"}, {"2017"}, nobtag_catagories, {"400"}, 1.018, 0.982)
  ({"em"}, {"2018"}, nobtag_catagories, {"400"}, 1.016, 0.984)
  ({"et"}, {"2016"}, nobtag_catagories, {"400"}, 1.008, 0.989)
  ({"et"}, {"2017"}, nobtag_catagories, {"400"}, 1.024, 0.98)
  ({"et"}, {"2018"}, nobtag_catagories, {"400"}, 1.013, 0.985)
  ({"mt"}, {"2016"}, nobtag_catagories, {"400"}, 1.009, 0.991)
  ({"mt"}, {"2017"}, nobtag_catagories, {"400"}, 1.023, 0.983)
  ({"mt"}, {"2018"}, nobtag_catagories, {"400"}, 1.018, 0.982)
  ({"tt"}, {"2016"}, nobtag_catagories, {"400"}, 1.013, 0.99)
  ({"tt"}, {"2017"}, nobtag_catagories, {"400"}, 1.023, 0.976)
  ({"tt"}, {"2018"}, nobtag_catagories, {"400"}, 1.015, 0.982)
  );
  
  cb.cp().process(mssm_bbH_signals).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"450"}, 0.98, 1.021)
  ({"em"}, {"2017"}, btag_catagories, {"450"}, 0.961, 1.034)
  ({"em"}, {"2018"}, btag_catagories, {"450"}, 0.978, 1.028)
  ({"et"}, {"2016"}, btag_catagories, {"450"}, 0.974, 1.025)
  ({"et"}, {"2017"}, btag_catagories, {"450"}, 0.967, 1.038)
  ({"et"}, {"2018"}, btag_catagories, {"450"}, 0.976, 1.02)
  ({"mt"}, {"2016"}, btag_catagories, {"450"}, 0.978, 1.022)
  ({"mt"}, {"2017"}, btag_catagories, {"450"}, 0.969, 1.037)
  ({"mt"}, {"2018"}, btag_catagories, {"450"}, 0.972, 1.029)
  ({"tt"}, {"2016"}, btag_catagories, {"450"}, 0.979, 1.02)
  ({"tt"}, {"2017"}, btag_catagories, {"450"}, 0.966, 1.034)
  ({"tt"}, {"2018"}, btag_catagories, {"450"}, 0.977, 1.025)
  ({"em"}, {"2016"}, nobtag_catagories, {"450"}, 1.01, 0.99)
  ({"em"}, {"2017"}, nobtag_catagories, {"450"}, 1.023, 0.98)
  ({"em"}, {"2018"}, nobtag_catagories, {"450"}, 1.014, 0.982)
  ({"et"}, {"2016"}, nobtag_catagories, {"450"}, 1.014, 0.986)
  ({"et"}, {"2017"}, nobtag_catagories, {"450"}, 1.021, 0.977)
  ({"et"}, {"2018"}, nobtag_catagories, {"450"}, 1.015, 0.987)
  ({"mt"}, {"2016"}, nobtag_catagories, {"450"}, 1.011, 0.988)
  ({"mt"}, {"2017"}, nobtag_catagories, {"450"}, 1.02, 0.977)
  ({"mt"}, {"2018"}, nobtag_catagories, {"450"}, 1.019, 0.98)
  ({"tt"}, {"2016"}, nobtag_catagories, {"450"}, 1.011, 0.989)
  ({"tt"}, {"2017"}, nobtag_catagories, {"450"}, 1.023, 0.977)
  ({"tt"}, {"2018"}, nobtag_catagories, {"450"}, 1.017, 0.982)
  );
  
  cb.cp().process(mssm_bbH_signals).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"500"}, 0.978, 1.02)
  ({"em"}, {"2017"}, btag_catagories, {"500"}, 0.973, 1.034)
  ({"em"}, {"2018"}, btag_catagories, {"500"}, 0.982, 1.023)
  ({"et"}, {"2016"}, btag_catagories, {"500"}, 0.984, 1.018)
  ({"et"}, {"2017"}, btag_catagories, {"500"}, 0.969, 1.03)
  ({"et"}, {"2018"}, btag_catagories, {"500"}, 0.978, 1.023)
  ({"mt"}, {"2016"}, btag_catagories, {"500"}, 0.983, 1.027)
  ({"mt"}, {"2017"}, btag_catagories, {"500"}, 0.964, 1.035)
  ({"mt"}, {"2018"}, btag_catagories, {"500"}, 0.974, 1.026)
  ({"tt"}, {"2016"}, btag_catagories, {"500"}, 0.985, 1.021)
  ({"tt"}, {"2017"}, btag_catagories, {"500"}, 0.964, 1.033)
  ({"tt"}, {"2018"}, btag_catagories, {"500"}, 0.977, 1.022)
  ({"em"}, {"2016"}, nobtag_catagories, {"500"}, 1.011, 0.99)
  ({"em"}, {"2017"}, nobtag_catagories, {"500"}, 1.016, 0.981)
  ({"em"}, {"2018"}, nobtag_catagories, {"500"}, 1.012, 0.985)
  ({"et"}, {"2016"}, nobtag_catagories, {"500"}, 1.009, 0.991)
  ({"et"}, {"2017"}, nobtag_catagories, {"500"}, 1.02, 0.981)
  ({"et"}, {"2018"}, nobtag_catagories, {"500"}, 1.015, 0.984)
  ({"mt"}, {"2016"}, nobtag_catagories, {"500"}, 1.009, 0.985)
  ({"mt"}, {"2017"}, nobtag_catagories, {"500"}, 1.024, 0.978)
  ({"mt"}, {"2018"}, nobtag_catagories, {"500"}, 1.019, 0.982)
  ({"tt"}, {"2016"}, nobtag_catagories, {"500"}, 1.009, 0.988)
  ({"tt"}, {"2017"}, nobtag_catagories, {"500"}, 1.025, 0.977)
  ({"tt"}, {"2018"}, nobtag_catagories, {"500"}, 1.018, 0.983)
  );
  
  cb.cp().process(mssm_bbH_signals).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"600"}, 0.988, 1.017)
  ({"em"}, {"2017"}, btag_catagories, {"600"}, 0.972, 1.036)
  ({"em"}, {"2018"}, btag_catagories, {"600"}, 0.977, 1.017)
  ({"et"}, {"2016"}, btag_catagories, {"600"}, 0.975, 1.021)
  ({"et"}, {"2017"}, btag_catagories, {"600"}, 0.973, 1.031)
  ({"et"}, {"2018"}, btag_catagories, {"600"}, 0.98, 1.024)
  ({"mt"}, {"2016"}, btag_catagories, {"600"}, 0.981, 1.021)
  ({"mt"}, {"2017"}, btag_catagories, {"600"}, 0.965, 1.034)
  ({"mt"}, {"2018"}, btag_catagories, {"600"}, 0.981, 1.019)
  ({"tt"}, {"2016"}, btag_catagories, {"600"}, 0.984, 1.024)
  ({"tt"}, {"2017"}, btag_catagories, {"600"}, 0.967, 1.037)
  ({"tt"}, {"2018"}, btag_catagories, {"600"}, 0.978, 1.018)
  ({"em"}, {"2016"}, nobtag_catagories, {"600"}, 1.006, 0.992)
  ({"em"}, {"2017"}, nobtag_catagories, {"600"}, 1.018, 0.977)
  ({"em"}, {"2018"}, nobtag_catagories, {"600"}, 1.016, 0.988)
  ({"et"}, {"2016"}, nobtag_catagories, {"600"}, 1.015, 0.987)
  ({"et"}, {"2017"}, nobtag_catagories, {"600"}, 1.019, 0.978)
  ({"et"}, {"2018"}, nobtag_catagories, {"600"}, 1.015, 0.982)
  ({"mt"}, {"2016"}, nobtag_catagories, {"600"}, 1.011, 0.988)
  ({"mt"}, {"2017"}, nobtag_catagories, {"600"}, 1.023, 0.978)
  ({"mt"}, {"2018"}, nobtag_catagories, {"600"}, 1.015, 0.985)
  ({"tt"}, {"2016"}, nobtag_catagories, {"600"}, 1.01, 0.985)
  ({"tt"}, {"2017"}, nobtag_catagories, {"600"}, 1.024, 0.973)
  ({"tt"}, {"2018"}, nobtag_catagories, {"600"}, 1.018, 0.985)
  );
  
  cb.cp().process(mssm_bbH_signals).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"700"}, 1.0, 1.0)
  ({"em"}, {"2017"}, btag_catagories, {"700"}, 0.964, 1.031)
  ({"em"}, {"2018"}, btag_catagories, {"700"}, 0.98, 1.033)
  ({"et"}, {"2016"}, btag_catagories, {"700"}, 1.0, 1.0)
  ({"et"}, {"2017"}, btag_catagories, {"700"}, 0.966, 1.031)
  ({"et"}, {"2018"}, btag_catagories, {"700"}, 0.978, 1.031)
  ({"mt"}, {"2016"}, btag_catagories, {"700"}, 1.0, 1.0)
  ({"mt"}, {"2017"}, btag_catagories, {"700"}, 0.962, 1.029)
  ({"mt"}, {"2018"}, btag_catagories, {"700"}, 0.98, 1.015)
  ({"tt"}, {"2016"}, btag_catagories, {"700"}, 1.0, 1.0)
  ({"tt"}, {"2017"}, btag_catagories, {"700"}, 0.964, 1.035)
  ({"tt"}, {"2018"}, btag_catagories, {"700"}, 0.979, 1.018)
  ({"em"}, {"2016"}, nobtag_catagories, {"700"}, 1.0, 1.0)
  ({"em"}, {"2017"}, nobtag_catagories, {"700"}, 1.02, 0.981)
  ({"em"}, {"2018"}, nobtag_catagories, {"700"}, 1.015, 0.976)
  ({"et"}, {"2016"}, nobtag_catagories, {"700"}, 1.0, 1.0)
  ({"et"}, {"2017"}, nobtag_catagories, {"700"}, 1.025, 0.978)
  ({"et"}, {"2018"}, nobtag_catagories, {"700"}, 1.017, 0.975)
  ({"mt"}, {"2016"}, nobtag_catagories, {"700"}, 1.0, 1.0)
  ({"mt"}, {"2017"}, nobtag_catagories, {"700"}, 1.03, 0.979)
  ({"mt"}, {"2018"}, nobtag_catagories, {"700"}, 1.016, 0.988)
  ({"tt"}, {"2016"}, nobtag_catagories, {"700"}, 1.0, 1.0)
  ({"tt"}, {"2017"}, nobtag_catagories, {"700"}, 1.028, 0.974)
  ({"tt"}, {"2018"}, nobtag_catagories, {"700"}, 1.017, 0.985)
  );
  
  cb.cp().process(mssm_bbH_signals).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"800"}, 0.989, 1.017)
  ({"em"}, {"2017"}, btag_catagories, {"800"}, 0.961, 1.037)
  ({"em"}, {"2018"}, btag_catagories, {"800"}, 0.981, 1.018)
  ({"et"}, {"2016"}, btag_catagories, {"800"}, 0.976, 1.015)
  ({"et"}, {"2017"}, btag_catagories, {"800"}, 0.97, 1.034)
  ({"et"}, {"2018"}, btag_catagories, {"800"}, 0.983, 1.03)
  ({"mt"}, {"2016"}, btag_catagories, {"800"}, 0.985, 1.024)
  ({"mt"}, {"2017"}, btag_catagories, {"800"}, 0.969, 1.037)
  ({"mt"}, {"2018"}, btag_catagories, {"800"}, 0.981, 1.02)
  ({"tt"}, {"2016"}, btag_catagories, {"800"}, 0.982, 1.018)
  ({"tt"}, {"2017"}, btag_catagories, {"800"}, 0.975, 1.033)
  ({"tt"}, {"2018"}, btag_catagories, {"800"}, 0.979, 1.016)
  ({"em"}, {"2016"}, nobtag_catagories, {"800"}, 1.006, 0.989)
  ({"em"}, {"2017"}, nobtag_catagories, {"800"}, 1.024, 0.977)
  ({"em"}, {"2018"}, nobtag_catagories, {"800"}, 1.014, 0.986)
  ({"et"}, {"2016"}, nobtag_catagories, {"800"}, 1.014, 0.991)
  ({"et"}, {"2017"}, nobtag_catagories, {"800"}, 1.022, 0.976)
  ({"et"}, {"2018"}, nobtag_catagories, {"800"}, 1.013, 0.976)
  ({"mt"}, {"2016"}, nobtag_catagories, {"800"}, 1.01, 0.985)
  ({"mt"}, {"2017"}, nobtag_catagories, {"800"}, 1.022, 0.974)
  ({"mt"}, {"2018"}, nobtag_catagories, {"800"}, 1.015, 0.983)
  ({"tt"}, {"2016"}, nobtag_catagories, {"800"}, 1.013, 0.988)
  ({"tt"}, {"2017"}, nobtag_catagories, {"800"}, 1.02, 0.974)
  ({"tt"}, {"2018"}, nobtag_catagories, {"800"}, 1.019, 0.985)
  );
  
  cb.cp().process(mssm_bbH_signals).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"900"}, 0.984, 1.022)
  ({"em"}, {"2017"}, btag_catagories, {"900"}, 0.973, 1.027)
  ({"em"}, {"2018"}, btag_catagories, {"900"}, 0.977, 1.025)
  ({"et"}, {"2016"}, btag_catagories, {"900"}, 0.98, 1.017)
  ({"et"}, {"2017"}, btag_catagories, {"900"}, 0.967, 1.031)
  ({"et"}, {"2018"}, btag_catagories, {"900"}, 0.975, 1.021)
  ({"mt"}, {"2016"}, btag_catagories, {"900"}, 0.984, 1.016)
  ({"mt"}, {"2017"}, btag_catagories, {"900"}, 0.964, 1.04)
  ({"mt"}, {"2018"}, btag_catagories, {"900"}, 0.981, 1.025)
  ({"tt"}, {"2016"}, btag_catagories, {"900"}, 0.981, 1.016)
  ({"tt"}, {"2017"}, btag_catagories, {"900"}, 0.969, 1.032)
  ({"tt"}, {"2018"}, btag_catagories, {"900"}, 0.979, 1.023)
  ({"em"}, {"2016"}, nobtag_catagories, {"900"}, 1.009, 0.988)
  ({"em"}, {"2017"}, nobtag_catagories, {"900"}, 1.016, 0.982)
  ({"em"}, {"2018"}, nobtag_catagories, {"900"}, 1.016, 0.982)
  ({"et"}, {"2016"}, nobtag_catagories, {"900"}, 1.012, 0.99)
  ({"et"}, {"2017"}, nobtag_catagories, {"900"}, 1.024, 0.979)
  ({"et"}, {"2018"}, nobtag_catagories, {"900"}, 1.021, 0.982)
  ({"mt"}, {"2016"}, nobtag_catagories, {"900"}, 1.01, 0.99)
  ({"mt"}, {"2017"}, nobtag_catagories, {"900"}, 1.028, 0.97)
  ({"mt"}, {"2018"}, nobtag_catagories, {"900"}, 1.016, 0.979)
  ({"tt"}, {"2016"}, nobtag_catagories, {"900"}, 1.013, 0.989)
  ({"tt"}, {"2017"}, nobtag_catagories, {"900"}, 1.024, 0.975)
  ({"tt"}, {"2018"}, nobtag_catagories, {"900"}, 1.019, 0.979)
  );
  
  cb.cp().process(mssm_bbH_signals).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"1000"}, 1.0, 1.0)
  ({"em"}, {"2017"}, btag_catagories, {"1000"}, 0.97, 1.018)
  ({"em"}, {"2018"}, btag_catagories, {"1000"}, 1.0, 1.0)
  ({"et"}, {"2016"}, btag_catagories, {"1000"}, 1.0, 1.0)
  ({"et"}, {"2017"}, btag_catagories, {"1000"}, 0.966, 1.035)
  ({"et"}, {"2018"}, btag_catagories, {"1000"}, 0.987, 1.031)
  ({"mt"}, {"2016"}, btag_catagories, {"1000"}, 1.0, 1.0)
  ({"mt"}, {"2017"}, btag_catagories, {"1000"}, 0.972, 1.041)
  ({"mt"}, {"2018"}, btag_catagories, {"1000"}, 0.981, 1.018)
  ({"tt"}, {"2016"}, btag_catagories, {"1000"}, 1.0, 1.0)
  ({"tt"}, {"2017"}, btag_catagories, {"1000"}, 0.968, 1.03)
  ({"tt"}, {"2018"}, btag_catagories, {"1000"}, 0.982, 1.02)
  ({"em"}, {"2016"}, nobtag_catagories, {"1000"}, 1.0, 1.0)
  ({"em"}, {"2017"}, nobtag_catagories, {"1000"}, 1.021, 0.987)
  ({"em"}, {"2018"}, nobtag_catagories, {"1000"}, 1.0, 1.0)
  ({"et"}, {"2016"}, nobtag_catagories, {"1000"}, 1.0, 1.0)
  ({"et"}, {"2017"}, nobtag_catagories, {"1000"}, 1.026, 0.974)
  ({"et"}, {"2018"}, nobtag_catagories, {"1000"}, 1.009, 0.978)
  ({"mt"}, {"2016"}, nobtag_catagories, {"1000"}, 1.0, 1.0)
  ({"mt"}, {"2017"}, nobtag_catagories, {"1000"}, 1.02, 0.972)
  ({"mt"}, {"2018"}, nobtag_catagories, {"1000"}, 1.016, 0.985)
  ({"tt"}, {"2016"}, nobtag_catagories, {"1000"}, 1.0, 1.0)
  ({"tt"}, {"2017"}, nobtag_catagories, {"1000"}, 1.028, 0.975)
  ({"tt"}, {"2018"}, nobtag_catagories, {"1000"}, 1.016, 0.982)
  );
  
  cb.cp().process(mssm_bbH_signals).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"1200"}, 0.983, 1.016)
  ({"em"}, {"2017"}, btag_catagories, {"1200"}, 0.972, 1.03)
  ({"em"}, {"2018"}, btag_catagories, {"1200"}, 0.984, 1.02)
  ({"et"}, {"2016"}, btag_catagories, {"1200"}, 0.968, 1.02)
  ({"et"}, {"2017"}, btag_catagories, {"1200"}, 0.973, 1.026)
  ({"et"}, {"2018"}, btag_catagories, {"1200"}, 0.977, 1.023)
  ({"mt"}, {"2016"}, btag_catagories, {"1200"}, 0.982, 1.016)
  ({"mt"}, {"2017"}, btag_catagories, {"1200"}, 0.964, 1.03)
  ({"mt"}, {"2018"}, btag_catagories, {"1200"}, 0.98, 1.022)
  ({"tt"}, {"2016"}, btag_catagories, {"1200"}, 0.979, 1.02)
  ({"tt"}, {"2017"}, btag_catagories, {"1200"}, 0.972, 1.028)
  ({"tt"}, {"2018"}, btag_catagories, {"1200"}, 0.981, 1.016)
  ({"em"}, {"2016"}, nobtag_catagories, {"1200"}, 1.011, 0.989)
  ({"em"}, {"2017"}, nobtag_catagories, {"1200"}, 1.017, 0.978)
  ({"em"}, {"2018"}, nobtag_catagories, {"1200"}, 1.013, 0.984)
  ({"et"}, {"2016"}, nobtag_catagories, {"1200"}, 1.023, 0.986)
  ({"et"}, {"2017"}, nobtag_catagories, {"1200"}, 1.019, 0.981)
  ({"et"}, {"2018"}, nobtag_catagories, {"1200"}, 1.018, 0.982)
  ({"mt"}, {"2016"}, nobtag_catagories, {"1200"}, 1.012, 0.99)
  ({"mt"}, {"2017"}, nobtag_catagories, {"1200"}, 1.026, 0.977)
  ({"mt"}, {"2018"}, nobtag_catagories, {"1200"}, 1.016, 0.981)
  ({"tt"}, {"2016"}, nobtag_catagories, {"1200"}, 1.016, 0.985)
  ({"tt"}, {"2017"}, nobtag_catagories, {"1200"}, 1.023, 0.977)
  ({"tt"}, {"2018"}, nobtag_catagories, {"1200"}, 1.018, 0.984)
  );
  
  cb.cp().process(mssm_bbH_signals).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"1400"}, 0.985, 1.01)
  ({"em"}, {"2017"}, btag_catagories, {"1400"}, 0.984, 1.031)
  ({"em"}, {"2018"}, btag_catagories, {"1400"}, 0.985, 1.006)
  ({"et"}, {"2016"}, btag_catagories, {"1400"}, 0.993, 1.013)
  ({"et"}, {"2017"}, btag_catagories, {"1400"}, 0.966, 1.015)
  ({"et"}, {"2018"}, btag_catagories, {"1400"}, 0.982, 1.021)
  ({"mt"}, {"2016"}, btag_catagories, {"1400"}, 0.974, 1.013)
  ({"mt"}, {"2017"}, btag_catagories, {"1400"}, 0.97, 1.034)
  ({"mt"}, {"2018"}, btag_catagories, {"1400"}, 0.978, 1.019)
  ({"tt"}, {"2016"}, btag_catagories, {"1400"}, 0.98, 1.016)
  ({"tt"}, {"2017"}, btag_catagories, {"1400"}, 0.97, 1.023)
  ({"tt"}, {"2018"}, btag_catagories, {"1400"}, 0.982, 1.013)
  ({"em"}, {"2016"}, nobtag_catagories, {"1400"}, 1.011, 0.992)
  ({"em"}, {"2017"}, nobtag_catagories, {"1400"}, 1.011, 0.978)
  ({"em"}, {"2018"}, nobtag_catagories, {"1400"}, 1.012, 0.994)
  ({"et"}, {"2016"}, nobtag_catagories, {"1400"}, 1.005, 0.991)
  ({"et"}, {"2017"}, nobtag_catagories, {"1400"}, 1.028, 0.988)
  ({"et"}, {"2018"}, nobtag_catagories, {"1400"}, 1.015, 0.983)
  ({"mt"}, {"2016"}, nobtag_catagories, {"1400"}, 1.019, 0.991)
  ({"mt"}, {"2017"}, nobtag_catagories, {"1400"}, 1.025, 0.973)
  ({"mt"}, {"2018"}, nobtag_catagories, {"1400"}, 1.019, 0.984)
  ({"tt"}, {"2016"}, nobtag_catagories, {"1400"}, 1.015, 0.988)
  ({"tt"}, {"2017"}, nobtag_catagories, {"1400"}, 1.026, 0.981)
  ({"tt"}, {"2018"}, nobtag_catagories, {"1400"}, 1.018, 0.987)
  );
  
  cb.cp().process(mssm_bbH_signals).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"1600"}, 0.988, 1.008)
  ({"em"}, {"2017"}, btag_catagories, {"1600"}, 1.0, 1.0)
  ({"em"}, {"2018"}, btag_catagories, {"1600"}, 0.989, 1.016)
  ({"et"}, {"2016"}, btag_catagories, {"1600"}, 0.982, 1.029)
  ({"et"}, {"2017"}, btag_catagories, {"1600"}, 0.984, 1.035)
  ({"et"}, {"2018"}, btag_catagories, {"1600"}, 0.982, 1.02)
  ({"mt"}, {"2016"}, btag_catagories, {"1600"}, 0.984, 1.02)
  ({"mt"}, {"2017"}, btag_catagories, {"1600"}, 0.971, 1.034)
  ({"mt"}, {"2018"}, btag_catagories, {"1600"}, 0.988, 1.018)
  ({"tt"}, {"2016"}, btag_catagories, {"1600"}, 0.987, 1.014)
  ({"tt"}, {"2017"}, btag_catagories, {"1600"}, 0.973, 1.029)
  ({"tt"}, {"2018"}, btag_catagories, {"1600"}, 0.983, 1.014)
  ({"em"}, {"2016"}, nobtag_catagories, {"1600"}, 1.008, 0.996)
  ({"em"}, {"2017"}, nobtag_catagories, {"1600"}, 1.0, 1.0)
  ({"em"}, {"2018"}, nobtag_catagories, {"1600"}, 1.01, 0.985)
  ({"et"}, {"2016"}, nobtag_catagories, {"1600"}, 1.012, 0.98)
  ({"et"}, {"2017"}, nobtag_catagories, {"1600"}, 1.013, 0.973)
  ({"et"}, {"2018"}, nobtag_catagories, {"1600"}, 1.016, 0.982)
  ({"mt"}, {"2016"}, nobtag_catagories, {"1600"}, 1.011, 0.988)
  ({"mt"}, {"2017"}, nobtag_catagories, {"1600"}, 1.023, 0.973)
  ({"mt"}, {"2018"}, nobtag_catagories, {"1600"}, 1.01, 0.986)
  ({"tt"}, {"2016"}, nobtag_catagories, {"1600"}, 1.011, 0.989)
  ({"tt"}, {"2017"}, nobtag_catagories, {"1600"}, 1.025, 0.974)
  ({"tt"}, {"2018"}, nobtag_catagories, {"1600"}, 1.016, 0.986)
  );
  
  cb.cp().process(mssm_bbH_signals).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"1800"}, 0.975, 1.02)
  ({"em"}, {"2017"}, btag_catagories, {"1800"}, 1.0, 1.0)
  ({"em"}, {"2018"}, btag_catagories, {"1800"}, 0.994, 1.026)
  ({"et"}, {"2016"}, btag_catagories, {"1800"}, 0.982, 1.017)
  ({"et"}, {"2017"}, btag_catagories, {"1800"}, 0.973, 1.035)
  ({"et"}, {"2018"}, btag_catagories, {"1800"}, 0.984, 1.021)
  ({"mt"}, {"2016"}, btag_catagories, {"1800"}, 0.988, 1.02)
  ({"mt"}, {"2017"}, btag_catagories, {"1800"}, 0.966, 1.028)
  ({"mt"}, {"2018"}, btag_catagories, {"1800"}, 0.98, 1.012)
  ({"tt"}, {"2016"}, btag_catagories, {"1800"}, 0.99, 1.022)
  ({"tt"}, {"2017"}, btag_catagories, {"1800"}, 0.968, 1.026)
  ({"tt"}, {"2018"}, btag_catagories, {"1800"}, 0.984, 1.018)
  ({"em"}, {"2016"}, nobtag_catagories, {"1800"}, 1.016, 0.987)
  ({"em"}, {"2017"}, nobtag_catagories, {"1800"}, 1.0, 1.0)
  ({"em"}, {"2018"}, nobtag_catagories, {"1800"}, 1.007, 0.98)
  ({"et"}, {"2016"}, nobtag_catagories, {"1800"}, 1.014, 0.987)
  ({"et"}, {"2017"}, nobtag_catagories, {"1800"}, 1.021, 0.975)
  ({"et"}, {"2018"}, nobtag_catagories, {"1800"}, 1.014, 0.981)
  ({"mt"}, {"2016"}, nobtag_catagories, {"1800"}, 1.008, 0.987)
  ({"mt"}, {"2017"}, nobtag_catagories, {"1800"}, 1.025, 0.98)
  ({"mt"}, {"2018"}, nobtag_catagories, {"1800"}, 1.018, 0.99)
  ({"tt"}, {"2016"}, nobtag_catagories, {"1800"}, 1.007, 0.984)
  ({"tt"}, {"2017"}, nobtag_catagories, {"1800"}, 1.028, 0.977)
  ({"tt"}, {"2018"}, nobtag_catagories, {"1800"}, 1.017, 0.981)
  );
  
  cb.cp().process(mssm_bbH_signals).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"2000"}, 0.974, 0.998)
  ({"em"}, {"2017"}, btag_catagories, {"2000"}, 0.982, 1.021)
  ({"em"}, {"2018"}, btag_catagories, {"2000"}, 0.972, 1.008)
  ({"et"}, {"2016"}, btag_catagories, {"2000"}, 0.985, 1.021)
  ({"et"}, {"2017"}, btag_catagories, {"2000"}, 0.982, 1.04)
  ({"et"}, {"2018"}, btag_catagories, {"2000"}, 0.978, 1.012)
  ({"mt"}, {"2016"}, btag_catagories, {"2000"}, 0.976, 1.01)
  ({"mt"}, {"2017"}, btag_catagories, {"2000"}, 0.965, 1.024)
  ({"mt"}, {"2018"}, btag_catagories, {"2000"}, 0.98, 1.018)
  ({"tt"}, {"2016"}, btag_catagories, {"2000"}, 0.986, 1.016)
  ({"tt"}, {"2017"}, btag_catagories, {"2000"}, 0.974, 1.023)
  ({"tt"}, {"2018"}, btag_catagories, {"2000"}, 0.983, 1.019)
  ({"em"}, {"2016"}, nobtag_catagories, {"2000"}, 1.021, 1.002)
  ({"em"}, {"2017"}, nobtag_catagories, {"2000"}, 1.016, 0.987)
  ({"em"}, {"2018"}, nobtag_catagories, {"2000"}, 1.021, 0.998)
  ({"et"}, {"2016"}, nobtag_catagories, {"2000"}, 1.011, 0.986)
  ({"et"}, {"2017"}, nobtag_catagories, {"2000"}, 1.015, 0.968)
  ({"et"}, {"2018"}, nobtag_catagories, {"2000"}, 1.019, 0.99)
  ({"mt"}, {"2016"}, nobtag_catagories, {"2000"}, 1.016, 0.993)
  ({"mt"}, {"2017"}, nobtag_catagories, {"2000"}, 1.027, 0.98)
  ({"mt"}, {"2018"}, nobtag_catagories, {"2000"}, 1.018, 0.983)
  ({"tt"}, {"2016"}, nobtag_catagories, {"2000"}, 1.011, 0.987)
  ({"tt"}, {"2017"}, nobtag_catagories, {"2000"}, 1.023, 0.98)
  ({"tt"}, {"2018"}, nobtag_catagories, {"2000"}, 1.017, 0.981)
  );
  
  cb.cp().process(mssm_bbH_signals).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"2300"}, 0.984, 1.016)
  ({"em"}, {"2017"}, btag_catagories, {"2300"}, 0.99, 1.026)
  ({"em"}, {"2018"}, btag_catagories, {"2300"}, 0.976, 1.011)
  ({"et"}, {"2016"}, btag_catagories, {"2300"}, 0.974, 1.016)
  ({"et"}, {"2017"}, btag_catagories, {"2300"}, 0.969, 1.03)
  ({"et"}, {"2018"}, btag_catagories, {"2300"}, 0.984, 1.01)
  ({"mt"}, {"2016"}, btag_catagories, {"2300"}, 0.981, 1.014)
  ({"mt"}, {"2017"}, btag_catagories, {"2300"}, 0.966, 1.035)
  ({"mt"}, {"2018"}, btag_catagories, {"2300"}, 0.987, 1.021)
  ({"tt"}, {"2016"}, btag_catagories, {"2300"}, 0.983, 1.018)
  ({"tt"}, {"2017"}, btag_catagories, {"2300"}, 0.968, 1.03)
  ({"tt"}, {"2018"}, btag_catagories, {"2300"}, 0.985, 1.018)
  ({"em"}, {"2016"}, nobtag_catagories, {"2300"}, 1.009, 0.99)
  ({"em"}, {"2017"}, nobtag_catagories, {"2300"}, 1.007, 0.978)
  ({"em"}, {"2018"}, nobtag_catagories, {"2300"}, 1.022, 0.993)
  ({"et"}, {"2016"}, nobtag_catagories, {"2300"}, 1.018, 0.989)
  ({"et"}, {"2017"}, nobtag_catagories, {"2300"}, 1.025, 0.981)
  ({"et"}, {"2018"}, nobtag_catagories, {"2300"}, 1.014, 0.995)
  ({"mt"}, {"2016"}, nobtag_catagories, {"2300"}, 1.013, 0.989)
  ({"mt"}, {"2017"}, nobtag_catagories, {"2300"}, 1.028, 0.973)
  ({"mt"}, {"2018"}, nobtag_catagories, {"2300"}, 1.01, 0.982)
  ({"tt"}, {"2016"}, nobtag_catagories, {"2300"}, 1.014, 0.986)
  ({"tt"}, {"2017"}, nobtag_catagories, {"2300"}, 1.028, 0.973)
  ({"tt"}, {"2018"}, nobtag_catagories, {"2300"}, 1.016, 0.982)
  );
  
  cb.cp().process(mssm_bbH_signals).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"2600"}, 0.987, 1.01)
  ({"em"}, {"2017"}, btag_catagories, {"2600"}, 0.965, 1.016)
  ({"em"}, {"2018"}, btag_catagories, {"2600"}, 0.98, 1.014)
  ({"et"}, {"2016"}, btag_catagories, {"2600"}, 0.993, 1.013)
  ({"et"}, {"2017"}, btag_catagories, {"2600"}, 1.0, 1.0)
  ({"et"}, {"2018"}, btag_catagories, {"2600"}, 0.983, 1.012)
  ({"mt"}, {"2016"}, btag_catagories, {"2600"}, 0.985, 1.018)
  ({"mt"}, {"2017"}, btag_catagories, {"2600"}, 1.0, 1.0)
  ({"mt"}, {"2018"}, btag_catagories, {"2600"}, 0.982, 1.025)
  ({"tt"}, {"2016"}, btag_catagories, {"2600"}, 0.985, 1.02)
  ({"tt"}, {"2017"}, btag_catagories, {"2600"}, 1.0, 1.0)
  ({"tt"}, {"2018"}, btag_catagories, {"2600"}, 0.976, 1.015)
  ({"em"}, {"2016"}, nobtag_catagories, {"2600"}, 1.013, 0.996)
  ({"em"}, {"2017"}, nobtag_catagories, {"2600"}, 1.028, 0.985)
  ({"em"}, {"2018"}, nobtag_catagories, {"2600"}, 1.016, 0.986)
  ({"et"}, {"2016"}, nobtag_catagories, {"2600"}, 1.004, 0.991)
  ({"et"}, {"2017"}, nobtag_catagories, {"2600"}, 1.0, 1.0)
  ({"et"}, {"2018"}, nobtag_catagories, {"2600"}, 1.013, 0.99)
  ({"mt"}, {"2016"}, nobtag_catagories, {"2600"}, 1.011, 0.986)
  ({"mt"}, {"2017"}, nobtag_catagories, {"2600"}, 1.0, 1.0)
  ({"mt"}, {"2018"}, nobtag_catagories, {"2600"}, 1.016, 0.979)
  ({"tt"}, {"2016"}, nobtag_catagories, {"2600"}, 1.012, 0.984)
  ({"tt"}, {"2017"}, nobtag_catagories, {"2600"}, 1.0, 1.0)
  ({"tt"}, {"2018"}, nobtag_catagories, {"2600"}, 1.024, 0.985)
  );
  
  cb.cp().process(mssm_bbH_signals).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"2900"}, 0.985, 1.019)
  ({"em"}, {"2017"}, btag_catagories, {"2900"}, 1.0, 1.0)
  ({"em"}, {"2018"}, btag_catagories, {"2900"}, 0.989, 1.015)
  ({"et"}, {"2016"}, btag_catagories, {"2900"}, 0.977, 1.012)
  ({"et"}, {"2017"}, btag_catagories, {"2900"}, 0.968, 1.019)
  ({"et"}, {"2018"}, btag_catagories, {"2900"}, 0.986, 1.011)
  ({"mt"}, {"2016"}, btag_catagories, {"2900"}, 0.983, 1.021)
  ({"mt"}, {"2017"}, btag_catagories, {"2900"}, 0.961, 1.028)
  ({"mt"}, {"2018"}, btag_catagories, {"2900"}, 0.986, 1.018)
  ({"tt"}, {"2016"}, btag_catagories, {"2900"}, 0.982, 1.016)
  ({"tt"}, {"2017"}, btag_catagories, {"2900"}, 0.969, 1.033)
  ({"tt"}, {"2018"}, btag_catagories, {"2900"}, 0.979, 1.01)
  ({"em"}, {"2016"}, nobtag_catagories, {"2900"}, 1.01, 0.986)
  ({"em"}, {"2017"}, nobtag_catagories, {"2900"}, 1.0, 1.0)
  ({"em"}, {"2018"}, nobtag_catagories, {"2900"}, 1.013, 0.99)
  ({"et"}, {"2016"}, nobtag_catagories, {"2900"}, 1.015, 0.993)
  ({"et"}, {"2017"}, nobtag_catagories, {"2900"}, 1.028, 0.982)
  ({"et"}, {"2018"}, nobtag_catagories, {"2900"}, 1.012, 0.992)
  ({"mt"}, {"2016"}, nobtag_catagories, {"2900"}, 1.011, 0.985)
  ({"mt"}, {"2017"}, nobtag_catagories, {"2900"}, 1.031, 0.98)
  ({"mt"}, {"2018"}, nobtag_catagories, {"2900"}, 1.012, 0.984)
  ({"tt"}, {"2016"}, nobtag_catagories, {"2900"}, 1.015, 0.987)
  ({"tt"}, {"2017"}, nobtag_catagories, {"2900"}, 1.029, 0.969)
  ({"tt"}, {"2018"}, nobtag_catagories, {"2900"}, 1.021, 0.99)
  );
  
  cb.cp().process(mssm_bbH_signals).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"3200"}, 0.985, 1.01)
  ({"em"}, {"2017"}, btag_catagories, {"3200"}, 0.96, 1.01)
  ({"em"}, {"2018"}, btag_catagories, {"3200"}, 1.0, 1.0)
  ({"et"}, {"2016"}, btag_catagories, {"3200"}, 0.984, 1.024)
  ({"et"}, {"2017"}, btag_catagories, {"3200"}, 0.972, 1.016)
  ({"et"}, {"2018"}, btag_catagories, {"3200"}, 0.982, 1.019)
  ({"mt"}, {"2016"}, btag_catagories, {"3200"}, 0.991, 1.015)
  ({"mt"}, {"2017"}, btag_catagories, {"3200"}, 0.987, 1.03)
  ({"mt"}, {"2018"}, btag_catagories, {"3200"}, 0.983, 1.016)
  ({"tt"}, {"2016"}, btag_catagories, {"3200"}, 0.979, 1.016)
  ({"tt"}, {"2017"}, btag_catagories, {"3200"}, 0.98, 1.041)
  ({"tt"}, {"2018"}, btag_catagories, {"3200"}, 0.981, 1.017)
  ({"em"}, {"2016"}, nobtag_catagories, {"3200"}, 1.013, 0.995)
  ({"em"}, {"2017"}, nobtag_catagories, {"3200"}, 1.029, 0.995)
  ({"em"}, {"2018"}, nobtag_catagories, {"3200"}, 1.0, 1.0)
  ({"et"}, {"2016"}, nobtag_catagories, {"3200"}, 1.009, 0.985)
  ({"et"}, {"2017"}, nobtag_catagories, {"3200"}, 1.024, 0.986)
  ({"et"}, {"2018"}, nobtag_catagories, {"3200"}, 1.019, 0.981)
  ({"mt"}, {"2016"}, nobtag_catagories, {"3200"}, 1.007, 0.991)
  ({"mt"}, {"2017"}, nobtag_catagories, {"3200"}, 1.011, 0.981)
  ({"mt"}, {"2018"}, nobtag_catagories, {"3200"}, 1.016, 0.987)
  ({"tt"}, {"2016"}, nobtag_catagories, {"3200"}, 1.018, 0.986)
  ({"tt"}, {"2017"}, nobtag_catagories, {"3200"}, 1.017, 0.964)
  ({"tt"}, {"2018"}, nobtag_catagories, {"3200"}, 1.02, 0.983)
  );
  
  cb.cp().process(mssm_bbH_signals).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"3500"}, 0.991, 1.02)
  ({"em"}, {"2017"}, btag_catagories, {"3500"}, 0.969, 1.036)
  ({"em"}, {"2018"}, btag_catagories, {"3500"}, 0.973, 1.005)
  ({"et"}, {"2016"}, btag_catagories, {"3500"}, 0.982, 1.021)
  ({"et"}, {"2017"}, btag_catagories, {"3500"}, 0.967, 1.03)
  ({"et"}, {"2018"}, btag_catagories, {"3500"}, 0.974, 1.013)
  ({"mt"}, {"2016"}, btag_catagories, {"3500"}, 0.988, 1.013)
  ({"mt"}, {"2017"}, btag_catagories, {"3500"}, 0.979, 1.028)
  ({"mt"}, {"2018"}, btag_catagories, {"3500"}, 0.976, 1.025)
  ({"tt"}, {"2016"}, btag_catagories, {"3500"}, 0.983, 1.026)
  ({"tt"}, {"2017"}, btag_catagories, {"3500"}, 0.967, 1.03)
  ({"tt"}, {"2018"}, btag_catagories, {"3500"}, 0.978, 1.022)
  ({"em"}, {"2016"}, nobtag_catagories, {"3500"}, 1.008, 0.986)
  ({"em"}, {"2017"}, nobtag_catagories, {"3500"}, 1.025, 0.98)
  ({"em"}, {"2018"}, nobtag_catagories, {"3500"}, 1.022, 0.994)
  ({"et"}, {"2016"}, nobtag_catagories, {"3500"}, 1.013, 0.987)
  ({"et"}, {"2017"}, nobtag_catagories, {"3500"}, 1.028, 0.971)
  ({"et"}, {"2018"}, nobtag_catagories, {"3500"}, 1.023, 0.99)
  ({"mt"}, {"2016"}, nobtag_catagories, {"3500"}, 1.008, 0.991)
  ({"mt"}, {"2017"}, nobtag_catagories, {"3500"}, 1.018, 0.973)
  ({"mt"}, {"2018"}, nobtag_catagories, {"3500"}, 1.021, 0.978)
  ({"tt"}, {"2016"}, nobtag_catagories, {"3500"}, 1.014, 0.979)
  ({"tt"}, {"2017"}, nobtag_catagories, {"3500"}, 1.029, 0.974)
  ({"tt"}, {"2018"}, nobtag_catagories, {"3500"}, 1.023, 0.977)
  );
  
  cb.cp().process(mssm_ggH_signals).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"60"}, 1.0, 1.0)
  ({"em"}, {"2017"}, btag_catagories, {"60"}, 1.0, 1.0)
  ({"em"}, {"2018"}, btag_catagories, {"60"}, 1.0, 1.0)
  ({"et"}, {"2016"}, btag_catagories, {"60"}, 1.0, 1.0)
  ({"et"}, {"2017"}, btag_catagories, {"60"}, 0.972, 1.025)
  ({"et"}, {"2018"}, btag_catagories, {"60"}, 1.0, 1.06)
  ({"mt"}, {"2016"}, btag_catagories, {"60"}, 1.0, 1.025)
  ({"mt"}, {"2017"}, btag_catagories, {"60"}, 0.963, 1.075)
  ({"mt"}, {"2018"}, btag_catagories, {"60"}, 1.0, 1.016)
  ({"tt"}, {"2016"}, btag_catagories, {"60"}, 1.0, 1.0)
  ({"tt"}, {"2017"}, btag_catagories, {"60"}, 0.953, 1.033)
  ({"tt"}, {"2018"}, btag_catagories, {"60"}, 0.956, 1.034)
  ({"em"}, {"2016"}, nobtag_catagories, {"60"}, 1.0, 1.0)
  ({"em"}, {"2017"}, nobtag_catagories, {"60"}, 1.0, 1.0)
  ({"em"}, {"2018"}, nobtag_catagories, {"60"}, 1.0, 1.0)
  ({"et"}, {"2016"}, nobtag_catagories, {"60"}, 1.0, 1.0)
  ({"et"}, {"2017"}, nobtag_catagories, {"60"}, 1.011, 0.996)
  ({"et"}, {"2018"}, nobtag_catagories, {"60"}, 1.0, 0.995)
  ({"mt"}, {"2016"}, nobtag_catagories, {"60"}, 1.0, 0.998)
  ({"mt"}, {"2017"}, nobtag_catagories, {"60"}, 1.003, 0.997)
  ({"mt"}, {"2018"}, nobtag_catagories, {"60"}, 1.0, 0.998)
  ({"tt"}, {"2016"}, nobtag_catagories, {"60"}, 1.0, 1.0)
  ({"tt"}, {"2017"}, nobtag_catagories, {"60"}, 1.006, 0.996)
  ({"tt"}, {"2018"}, nobtag_catagories, {"60"}, 1.007, 0.994)
  );
  
  cb.cp().process(mssm_ggH_signals).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"80"}, 1.0, 1.0)
  ({"em"}, {"2017"}, btag_catagories, {"80"}, 1.0, 1.0)
  ({"em"}, {"2018"}, btag_catagories, {"80"}, 1.0, 1.0)
  ({"et"}, {"2016"}, btag_catagories, {"80"}, 0.975, 1.0)
  ({"et"}, {"2017"}, btag_catagories, {"80"}, 1.0, 1.035)
  ({"et"}, {"2018"}, btag_catagories, {"80"}, 0.935, 1.0)
  ({"mt"}, {"2016"}, btag_catagories, {"80"}, 0.978, 1.0)
  ({"mt"}, {"2017"}, btag_catagories, {"80"}, 0.954, 1.076)
  ({"mt"}, {"2018"}, btag_catagories, {"80"}, 0.983, 1.027)
  ({"tt"}, {"2016"}, btag_catagories, {"80"}, 0.982, 1.026)
  ({"tt"}, {"2017"}, btag_catagories, {"80"}, 0.936, 1.021)
  ({"tt"}, {"2018"}, btag_catagories, {"80"}, 1.0, 1.041)
  ({"em"}, {"2016"}, nobtag_catagories, {"80"}, 1.0, 1.0)
  ({"em"}, {"2017"}, nobtag_catagories, {"80"}, 1.0, 1.0)
  ({"em"}, {"2018"}, nobtag_catagories, {"80"}, 1.0, 1.0)
  ({"et"}, {"2016"}, nobtag_catagories, {"80"}, 1.001, 1.0)
  ({"et"}, {"2017"}, nobtag_catagories, {"80"}, 1.0, 0.998)
  ({"et"}, {"2018"}, nobtag_catagories, {"80"}, 1.005, 1.0)
  ({"mt"}, {"2016"}, nobtag_catagories, {"80"}, 1.001, 1.0)
  ({"mt"}, {"2017"}, nobtag_catagories, {"80"}, 1.002, 0.996)
  ({"mt"}, {"2018"}, nobtag_catagories, {"80"}, 1.001, 0.998)
  ({"tt"}, {"2016"}, nobtag_catagories, {"80"}, 1.002, 0.996)
  ({"tt"}, {"2017"}, nobtag_catagories, {"80"}, 1.01, 0.997)
  ({"tt"}, {"2018"}, nobtag_catagories, {"80"}, 1.0, 0.995)
  );
  
  cb.cp().process(mssm_ggH_signals).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"100"}, 0.987, 1.0)
  ({"em"}, {"2017"}, btag_catagories, {"100"}, 0.988, 1.032)
  ({"em"}, {"2018"}, btag_catagories, {"100"}, 1.0, 1.0)
  ({"et"}, {"2016"}, btag_catagories, {"100"}, 0.96, 1.043)
  ({"et"}, {"2017"}, btag_catagories, {"100"}, 0.942, 1.07)
  ({"et"}, {"2018"}, btag_catagories, {"100"}, 0.994, 1.051)
  ({"mt"}, {"2016"}, btag_catagories, {"100"}, 0.98, 1.014)
  ({"mt"}, {"2017"}, btag_catagories, {"100"}, 0.968, 1.034)
  ({"mt"}, {"2018"}, btag_catagories, {"100"}, 0.983, 1.026)
  ({"tt"}, {"2016"}, btag_catagories, {"100"}, 0.948, 1.025)
  ({"tt"}, {"2017"}, btag_catagories, {"100"}, 0.948, 1.027)
  ({"tt"}, {"2018"}, btag_catagories, {"100"}, 0.98, 1.02)
  ({"em"}, {"2016"}, nobtag_catagories, {"100"}, 1.0, 1.0)
  ({"em"}, {"2017"}, nobtag_catagories, {"100"}, 1.0, 0.999)
  ({"em"}, {"2018"}, nobtag_catagories, {"100"}, 1.0, 1.0)
  ({"et"}, {"2016"}, nobtag_catagories, {"100"}, 1.002, 0.998)
  ({"et"}, {"2017"}, nobtag_catagories, {"100"}, 1.002, 0.997)
  ({"et"}, {"2018"}, nobtag_catagories, {"100"}, 1.0, 0.997)
  ({"mt"}, {"2016"}, nobtag_catagories, {"100"}, 1.001, 0.999)
  ({"mt"}, {"2017"}, nobtag_catagories, {"100"}, 1.002, 0.999)
  ({"mt"}, {"2018"}, nobtag_catagories, {"100"}, 1.001, 0.999)
  ({"tt"}, {"2016"}, nobtag_catagories, {"100"}, 1.004, 0.998)
  ({"tt"}, {"2017"}, nobtag_catagories, {"100"}, 1.004, 0.998)
  ({"tt"}, {"2018"}, nobtag_catagories, {"100"}, 1.002, 0.998)
  );
  
  cb.cp().process(mssm_ggH_signals).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"120"}, 1.0, 1.0)
  ({"em"}, {"2017"}, btag_catagories, {"120"}, 0.995, 1.012)
  ({"em"}, {"2018"}, btag_catagories, {"120"}, 1.0, 1.02)
  ({"et"}, {"2016"}, btag_catagories, {"120"}, 0.982, 1.009)
  ({"et"}, {"2017"}, btag_catagories, {"120"}, 0.952, 1.047)
  ({"et"}, {"2018"}, btag_catagories, {"120"}, 0.983, 1.029)
  ({"mt"}, {"2016"}, btag_catagories, {"120"}, 0.981, 1.018)
  ({"mt"}, {"2017"}, btag_catagories, {"120"}, 0.98, 1.019)
  ({"mt"}, {"2018"}, btag_catagories, {"120"}, 0.99, 1.038)
  ({"tt"}, {"2016"}, btag_catagories, {"120"}, 0.981, 1.008)
  ({"tt"}, {"2017"}, btag_catagories, {"120"}, 0.983, 1.046)
  ({"tt"}, {"2018"}, btag_catagories, {"120"}, 0.977, 1.039)
  ({"em"}, {"2016"}, nobtag_catagories, {"120"}, 1.0, 1.0)
  ({"em"}, {"2017"}, nobtag_catagories, {"120"}, 1.0, 1.0)
  ({"em"}, {"2018"}, nobtag_catagories, {"120"}, 1.0, 0.999)
  ({"et"}, {"2016"}, nobtag_catagories, {"120"}, 1.001, 1.0)
  ({"et"}, {"2017"}, nobtag_catagories, {"120"}, 1.002, 0.999)
  ({"et"}, {"2018"}, nobtag_catagories, {"120"}, 1.0, 0.998)
  ({"mt"}, {"2016"}, nobtag_catagories, {"120"}, 1.001, 0.999)
  ({"mt"}, {"2017"}, nobtag_catagories, {"120"}, 1.001, 1.0)
  ({"mt"}, {"2018"}, nobtag_catagories, {"120"}, 1.0, 0.998)
  ({"tt"}, {"2016"}, nobtag_catagories, {"120"}, 1.001, 1.0)
  ({"tt"}, {"2017"}, nobtag_catagories, {"120"}, 1.001, 0.998)
  ({"tt"}, {"2018"}, nobtag_catagories, {"120"}, 1.002, 0.997)
  );
  
  cb.cp().process(mssm_ggH_signals).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"125"}, 0.985, 1.0)
  ({"em"}, {"2017"}, btag_catagories, {"125"}, 0.988, 1.004)
  ({"em"}, {"2018"}, btag_catagories, {"125"}, 0.991, 1.004)
  ({"et"}, {"2016"}, btag_catagories, {"125"}, 0.984, 1.047)
  ({"et"}, {"2017"}, btag_catagories, {"125"}, 0.958, 1.025)
  ({"et"}, {"2018"}, btag_catagories, {"125"}, 0.985, 1.039)
  ({"mt"}, {"2016"}, btag_catagories, {"125"}, 0.963, 1.023)
  ({"mt"}, {"2017"}, btag_catagories, {"125"}, 0.966, 1.031)
  ({"mt"}, {"2018"}, btag_catagories, {"125"}, 0.986, 1.015)
  ({"tt"}, {"2016"}, btag_catagories, {"125"}, 0.982, 1.017)
  ({"tt"}, {"2017"}, btag_catagories, {"125"}, 0.958, 1.033)
  ({"tt"}, {"2018"}, btag_catagories, {"125"}, 0.986, 1.016)
  ({"em"}, {"2016"}, nobtag_catagories, {"125"}, 1.0, 1.0)
  ({"em"}, {"2017"}, nobtag_catagories, {"125"}, 1.0, 1.0)
  ({"em"}, {"2018"}, nobtag_catagories, {"125"}, 1.0, 1.0)
  ({"et"}, {"2016"}, nobtag_catagories, {"125"}, 1.001, 0.998)
  ({"et"}, {"2017"}, nobtag_catagories, {"125"}, 1.002, 0.999)
  ({"et"}, {"2018"}, nobtag_catagories, {"125"}, 1.0, 0.998)
  ({"mt"}, {"2016"}, nobtag_catagories, {"125"}, 1.001, 0.999)
  ({"mt"}, {"2017"}, nobtag_catagories, {"125"}, 1.001, 0.999)
  ({"mt"}, {"2018"}, nobtag_catagories, {"125"}, 1.0, 0.999)
  ({"tt"}, {"2016"}, nobtag_catagories, {"125"}, 1.001, 0.999)
  ({"tt"}, {"2017"}, nobtag_catagories, {"125"}, 1.002, 0.998)
  ({"tt"}, {"2018"}, nobtag_catagories, {"125"}, 1.001, 0.999)
  );
  
  cb.cp().process(mssm_ggH_signals).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"130"}, 1.0, 1.011)
  ({"em"}, {"2017"}, btag_catagories, {"130"}, 0.992, 1.003)
  ({"em"}, {"2018"}, btag_catagories, {"130"}, 0.979, 1.0)
  ({"et"}, {"2016"}, btag_catagories, {"130"}, 0.974, 1.031)
  ({"et"}, {"2017"}, btag_catagories, {"130"}, 0.971, 1.041)
  ({"et"}, {"2018"}, btag_catagories, {"130"}, 0.988, 1.031)
  ({"mt"}, {"2016"}, btag_catagories, {"130"}, 0.978, 1.025)
  ({"mt"}, {"2017"}, btag_catagories, {"130"}, 0.966, 1.029)
  ({"mt"}, {"2018"}, btag_catagories, {"130"}, 0.976, 1.014)
  ({"tt"}, {"2016"}, btag_catagories, {"130"}, 0.97, 1.021)
  ({"tt"}, {"2017"}, btag_catagories, {"130"}, 0.948, 1.03)
  ({"tt"}, {"2018"}, btag_catagories, {"130"}, 0.965, 1.022)
  ({"em"}, {"2016"}, nobtag_catagories, {"130"}, 1.0, 1.0)
  ({"em"}, {"2017"}, nobtag_catagories, {"130"}, 1.0, 1.0)
  ({"em"}, {"2018"}, nobtag_catagories, {"130"}, 1.001, 1.0)
  ({"et"}, {"2016"}, nobtag_catagories, {"130"}, 1.0, 0.999)
  ({"et"}, {"2017"}, nobtag_catagories, {"130"}, 1.001, 0.999)
  ({"et"}, {"2018"}, nobtag_catagories, {"130"}, 1.0, 0.998)
  ({"mt"}, {"2016"}, nobtag_catagories, {"130"}, 1.001, 0.999)
  ({"mt"}, {"2017"}, nobtag_catagories, {"130"}, 1.001, 0.999)
  ({"mt"}, {"2018"}, nobtag_catagories, {"130"}, 1.001, 0.999)
  ({"tt"}, {"2016"}, nobtag_catagories, {"130"}, 1.001, 0.999)
  ({"tt"}, {"2017"}, nobtag_catagories, {"130"}, 1.003, 0.999)
  ({"tt"}, {"2018"}, nobtag_catagories, {"130"}, 1.002, 0.999)
  );
  
  cb.cp().process(mssm_ggH_signals).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"140"}, 0.997, 1.003)
  ({"em"}, {"2017"}, btag_catagories, {"140"}, 0.994, 1.008)
  ({"em"}, {"2018"}, btag_catagories, {"140"}, 1.0, 1.008)
  ({"et"}, {"2016"}, btag_catagories, {"140"}, 0.963, 1.018)
  ({"et"}, {"2017"}, btag_catagories, {"140"}, 0.961, 1.048)
  ({"et"}, {"2018"}, btag_catagories, {"140"}, 0.971, 1.021)
  ({"mt"}, {"2016"}, btag_catagories, {"140"}, 0.965, 1.03)
  ({"mt"}, {"2017"}, btag_catagories, {"140"}, 0.971, 1.034)
  ({"mt"}, {"2018"}, btag_catagories, {"140"}, 0.984, 1.012)
  ({"tt"}, {"2016"}, btag_catagories, {"140"}, 0.981, 1.04)
  ({"tt"}, {"2017"}, btag_catagories, {"140"}, 0.966, 1.043)
  ({"tt"}, {"2018"}, btag_catagories, {"140"}, 0.984, 1.03)
  ({"em"}, {"2016"}, nobtag_catagories, {"140"}, 1.0, 1.0)
  ({"em"}, {"2017"}, nobtag_catagories, {"140"}, 1.0, 1.0)
  ({"em"}, {"2018"}, nobtag_catagories, {"140"}, 1.0, 1.0)
  ({"et"}, {"2016"}, nobtag_catagories, {"140"}, 1.001, 1.0)
  ({"et"}, {"2017"}, nobtag_catagories, {"140"}, 1.002, 0.998)
  ({"et"}, {"2018"}, nobtag_catagories, {"140"}, 1.001, 0.999)
  ({"mt"}, {"2016"}, nobtag_catagories, {"140"}, 1.001, 0.999)
  ({"mt"}, {"2017"}, nobtag_catagories, {"140"}, 1.001, 0.999)
  ({"mt"}, {"2018"}, nobtag_catagories, {"140"}, 1.001, 0.999)
  ({"tt"}, {"2016"}, nobtag_catagories, {"140"}, 1.001, 0.999)
  ({"tt"}, {"2017"}, nobtag_catagories, {"140"}, 1.002, 0.998)
  ({"tt"}, {"2018"}, nobtag_catagories, {"140"}, 1.001, 0.998)
  );
  
  cb.cp().process(mssm_ggH_signals).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"160"}, 0.993, 1.011)
  ({"em"}, {"2017"}, btag_catagories, {"160"}, 0.99, 1.012)
  ({"em"}, {"2018"}, btag_catagories, {"160"}, 0.993, 1.003)
  ({"et"}, {"2016"}, btag_catagories, {"160"}, 0.972, 1.026)
  ({"et"}, {"2017"}, btag_catagories, {"160"}, 0.961, 1.019)
  ({"et"}, {"2018"}, btag_catagories, {"160"}, 0.984, 1.025)
  ({"mt"}, {"2016"}, btag_catagories, {"160"}, 0.98, 1.02)
  ({"mt"}, {"2017"}, btag_catagories, {"160"}, 0.959, 1.022)
  ({"mt"}, {"2018"}, btag_catagories, {"160"}, 0.975, 1.023)
  ({"tt"}, {"2016"}, btag_catagories, {"160"}, 0.979, 1.033)
  ({"tt"}, {"2017"}, btag_catagories, {"160"}, 0.967, 1.035)
  ({"tt"}, {"2018"}, btag_catagories, {"160"}, 0.985, 1.023)
  ({"em"}, {"2016"}, nobtag_catagories, {"160"}, 1.0, 1.0)
  ({"em"}, {"2017"}, nobtag_catagories, {"160"}, 1.0, 1.0)
  ({"em"}, {"2018"}, nobtag_catagories, {"160"}, 1.0, 1.0)
  ({"et"}, {"2016"}, nobtag_catagories, {"160"}, 1.001, 0.999)
  ({"et"}, {"2017"}, nobtag_catagories, {"160"}, 1.002, 0.999)
  ({"et"}, {"2018"}, nobtag_catagories, {"160"}, 1.001, 0.999)
  ({"mt"}, {"2016"}, nobtag_catagories, {"160"}, 1.001, 0.999)
  ({"mt"}, {"2017"}, nobtag_catagories, {"160"}, 1.002, 0.999)
  ({"mt"}, {"2018"}, nobtag_catagories, {"160"}, 1.001, 0.998)
  ({"tt"}, {"2016"}, nobtag_catagories, {"160"}, 1.001, 0.999)
  ({"tt"}, {"2017"}, nobtag_catagories, {"160"}, 1.002, 0.998)
  ({"tt"}, {"2018"}, nobtag_catagories, {"160"}, 1.001, 0.999)
  );
  
  cb.cp().process(mssm_ggH_signals).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"180"}, 0.993, 1.007)
  ({"em"}, {"2017"}, btag_catagories, {"180"}, 0.995, 1.0)
  ({"em"}, {"2018"}, btag_catagories, {"180"}, 0.995, 1.015)
  ({"et"}, {"2016"}, btag_catagories, {"180"}, 0.977, 1.031)
  ({"et"}, {"2017"}, btag_catagories, {"180"}, 0.972, 1.041)
  ({"et"}, {"2018"}, btag_catagories, {"180"}, 0.993, 1.022)
  ({"mt"}, {"2016"}, btag_catagories, {"180"}, 0.989, 1.025)
  ({"mt"}, {"2017"}, btag_catagories, {"180"}, 0.964, 1.033)
  ({"mt"}, {"2018"}, btag_catagories, {"180"}, 0.979, 1.025)
  ({"tt"}, {"2016"}, btag_catagories, {"180"}, 0.975, 1.021)
  ({"tt"}, {"2017"}, btag_catagories, {"180"}, 0.965, 1.028)
  ({"tt"}, {"2018"}, btag_catagories, {"180"}, 0.976, 1.024)
  ({"em"}, {"2016"}, nobtag_catagories, {"180"}, 1.0, 1.0)
  ({"em"}, {"2017"}, nobtag_catagories, {"180"}, 1.0, 1.0)
  ({"em"}, {"2018"}, nobtag_catagories, {"180"}, 1.0, 0.999)
  ({"et"}, {"2016"}, nobtag_catagories, {"180"}, 1.001, 0.999)
  ({"et"}, {"2017"}, nobtag_catagories, {"180"}, 1.002, 0.998)
  ({"et"}, {"2018"}, nobtag_catagories, {"180"}, 1.0, 0.999)
  ({"mt"}, {"2016"}, nobtag_catagories, {"180"}, 1.0, 0.999)
  ({"mt"}, {"2017"}, nobtag_catagories, {"180"}, 1.002, 0.999)
  ({"mt"}, {"2018"}, nobtag_catagories, {"180"}, 1.001, 0.999)
  ({"tt"}, {"2016"}, nobtag_catagories, {"180"}, 1.001, 0.999)
  ({"tt"}, {"2017"}, nobtag_catagories, {"180"}, 1.002, 0.999)
  ({"tt"}, {"2018"}, nobtag_catagories, {"180"}, 1.001, 0.999)
  );
  
  cb.cp().process(mssm_ggH_signals).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"200"}, 0.991, 1.003)
  ({"em"}, {"2017"}, btag_catagories, {"200"}, 0.996, 1.003)
  ({"em"}, {"2018"}, btag_catagories, {"200"}, 0.997, 1.003)
  ({"et"}, {"2016"}, btag_catagories, {"200"}, 0.974, 1.015)
  ({"et"}, {"2017"}, btag_catagories, {"200"}, 0.959, 1.043)
  ({"et"}, {"2018"}, btag_catagories, {"200"}, 0.968, 1.014)
  ({"mt"}, {"2016"}, btag_catagories, {"200"}, 0.979, 1.018)
  ({"mt"}, {"2017"}, btag_catagories, {"200"}, 0.964, 1.049)
  ({"mt"}, {"2018"}, btag_catagories, {"200"}, 0.981, 1.025)
  ({"tt"}, {"2016"}, btag_catagories, {"200"}, 0.97, 1.028)
  ({"tt"}, {"2017"}, btag_catagories, {"200"}, 0.96, 1.031)
  ({"tt"}, {"2018"}, btag_catagories, {"200"}, 0.976, 1.027)
  ({"em"}, {"2016"}, nobtag_catagories, {"200"}, 1.0, 1.0)
  ({"em"}, {"2017"}, nobtag_catagories, {"200"}, 1.0, 1.0)
  ({"em"}, {"2018"}, nobtag_catagories, {"200"}, 1.0, 1.0)
  ({"et"}, {"2016"}, nobtag_catagories, {"200"}, 1.001, 1.0)
  ({"et"}, {"2017"}, nobtag_catagories, {"200"}, 1.002, 0.998)
  ({"et"}, {"2018"}, nobtag_catagories, {"200"}, 1.001, 0.999)
  ({"mt"}, {"2016"}, nobtag_catagories, {"200"}, 1.001, 0.999)
  ({"mt"}, {"2017"}, nobtag_catagories, {"200"}, 1.002, 0.998)
  ({"mt"}, {"2018"}, nobtag_catagories, {"200"}, 1.001, 0.998)
  ({"tt"}, {"2016"}, nobtag_catagories, {"200"}, 1.001, 0.999)
  ({"tt"}, {"2017"}, nobtag_catagories, {"200"}, 1.002, 0.999)
  ({"tt"}, {"2018"}, nobtag_catagories, {"200"}, 1.001, 0.998)
  );
  
  cb.cp().process(mssm_ggH_signals).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"250"}, 0.992, 1.012)
  ({"em"}, {"2017"}, btag_catagories, {"250"}, 0.987, 1.005)
  ({"em"}, {"2018"}, btag_catagories, {"250"}, 1.0, 1.015)
  ({"et"}, {"2016"}, btag_catagories, {"250"}, 0.984, 1.019)
  ({"et"}, {"2017"}, btag_catagories, {"250"}, 0.982, 1.021)
  ({"et"}, {"2018"}, btag_catagories, {"250"}, 0.99, 1.02)
  ({"mt"}, {"2016"}, btag_catagories, {"250"}, 0.972, 1.02)
  ({"mt"}, {"2017"}, btag_catagories, {"250"}, 0.975, 1.045)
  ({"mt"}, {"2018"}, btag_catagories, {"250"}, 0.979, 1.026)
  ({"tt"}, {"2016"}, btag_catagories, {"250"}, 0.981, 1.021)
  ({"tt"}, {"2017"}, btag_catagories, {"250"}, 0.968, 1.034)
  ({"tt"}, {"2018"}, btag_catagories, {"250"}, 0.983, 1.024)
  ({"em"}, {"2016"}, nobtag_catagories, {"250"}, 1.0, 1.0)
  ({"em"}, {"2017"}, nobtag_catagories, {"250"}, 1.001, 1.0)
  ({"em"}, {"2018"}, nobtag_catagories, {"250"}, 1.0, 0.999)
  ({"et"}, {"2016"}, nobtag_catagories, {"250"}, 1.001, 0.999)
  ({"et"}, {"2017"}, nobtag_catagories, {"250"}, 1.001, 0.999)
  ({"et"}, {"2018"}, nobtag_catagories, {"250"}, 1.001, 0.999)
  ({"mt"}, {"2016"}, nobtag_catagories, {"250"}, 1.001, 0.999)
  ({"mt"}, {"2017"}, nobtag_catagories, {"250"}, 1.001, 0.998)
  ({"mt"}, {"2018"}, nobtag_catagories, {"250"}, 1.001, 0.999)
  ({"tt"}, {"2016"}, nobtag_catagories, {"250"}, 1.001, 0.999)
  ({"tt"}, {"2017"}, nobtag_catagories, {"250"}, 1.002, 0.998)
  ({"tt"}, {"2018"}, nobtag_catagories, {"250"}, 1.001, 0.999)
  );
  
  cb.cp().process(mssm_ggH_signals).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"300"}, 0.997, 1.006)
  ({"em"}, {"2017"}, btag_catagories, {"300"}, 0.995, 1.013)
  ({"em"}, {"2018"}, btag_catagories, {"300"}, 0.999, 1.002)
  ({"et"}, {"2016"}, btag_catagories, {"300"}, 0.983, 1.021)
  ({"et"}, {"2017"}, btag_catagories, {"300"}, 0.961, 1.046)
  ({"et"}, {"2018"}, btag_catagories, {"300"}, 0.972, 1.021)
  ({"mt"}, {"2016"}, btag_catagories, {"300"}, 0.974, 1.025)
  ({"mt"}, {"2017"}, btag_catagories, {"300"}, 0.969, 1.046)
  ({"mt"}, {"2018"}, btag_catagories, {"300"}, 0.98, 1.022)
  ({"tt"}, {"2016"}, btag_catagories, {"300"}, 0.974, 1.029)
  ({"tt"}, {"2017"}, btag_catagories, {"300"}, 0.968, 1.029)
  ({"tt"}, {"2018"}, btag_catagories, {"300"}, 0.976, 1.02)
  ({"em"}, {"2016"}, nobtag_catagories, {"300"}, 1.0, 1.0)
  ({"em"}, {"2017"}, nobtag_catagories, {"300"}, 1.0, 0.999)
  ({"em"}, {"2018"}, nobtag_catagories, {"300"}, 1.0, 1.0)
  ({"et"}, {"2016"}, nobtag_catagories, {"300"}, 1.001, 0.999)
  ({"et"}, {"2017"}, nobtag_catagories, {"300"}, 1.002, 0.998)
  ({"et"}, {"2018"}, nobtag_catagories, {"300"}, 1.001, 0.999)
  ({"mt"}, {"2016"}, nobtag_catagories, {"300"}, 1.001, 0.999)
  ({"mt"}, {"2017"}, nobtag_catagories, {"300"}, 1.001, 0.998)
  ({"mt"}, {"2018"}, nobtag_catagories, {"300"}, 1.001, 0.998)
  ({"tt"}, {"2016"}, nobtag_catagories, {"300"}, 1.001, 0.999)
  ({"tt"}, {"2017"}, nobtag_catagories, {"300"}, 1.002, 0.999)
  ({"tt"}, {"2018"}, nobtag_catagories, {"300"}, 1.002, 0.999)
  );
  
  cb.cp().process(mssm_ggH_signals).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"350"}, 0.994, 1.014)
  ({"em"}, {"2017"}, btag_catagories, {"350"}, 0.998, 1.009)
  ({"em"}, {"2018"}, btag_catagories, {"350"}, 0.999, 1.012)
  ({"et"}, {"2016"}, btag_catagories, {"350"}, 0.977, 1.022)
  ({"et"}, {"2017"}, btag_catagories, {"350"}, 0.969, 1.049)
  ({"et"}, {"2018"}, btag_catagories, {"350"}, 1.0, 1.0)
  ({"mt"}, {"2016"}, btag_catagories, {"350"}, 0.976, 1.019)
  ({"mt"}, {"2017"}, btag_catagories, {"350"}, 0.959, 1.031)
  ({"mt"}, {"2018"}, btag_catagories, {"350"}, 1.0, 1.0)
  ({"tt"}, {"2016"}, btag_catagories, {"350"}, 0.975, 1.027)
  ({"tt"}, {"2017"}, btag_catagories, {"350"}, 0.97, 1.039)
  ({"tt"}, {"2018"}, btag_catagories, {"350"}, 1.0, 1.0)
  ({"em"}, {"2016"}, nobtag_catagories, {"350"}, 1.0, 1.0)
  ({"em"}, {"2017"}, nobtag_catagories, {"350"}, 1.0, 1.0)
  ({"em"}, {"2018"}, nobtag_catagories, {"350"}, 1.0, 0.999)
  ({"et"}, {"2016"}, nobtag_catagories, {"350"}, 1.001, 0.999)
  ({"et"}, {"2017"}, nobtag_catagories, {"350"}, 1.002, 0.998)
  ({"et"}, {"2018"}, nobtag_catagories, {"350"}, 1.0, 1.0)
  ({"mt"}, {"2016"}, nobtag_catagories, {"350"}, 1.001, 0.999)
  ({"mt"}, {"2017"}, nobtag_catagories, {"350"}, 1.002, 0.999)
  ({"mt"}, {"2018"}, nobtag_catagories, {"350"}, 1.0, 1.0)
  ({"tt"}, {"2016"}, nobtag_catagories, {"350"}, 1.001, 0.999)
  ({"tt"}, {"2017"}, nobtag_catagories, {"350"}, 1.002, 0.998)
  ({"tt"}, {"2018"}, nobtag_catagories, {"350"}, 1.0, 1.0)
  );
  
  cb.cp().process(mssm_ggH_signals).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"400"}, 0.992, 1.0)
  ({"em"}, {"2017"}, btag_catagories, {"400"}, 0.998, 1.003)
  ({"em"}, {"2018"}, btag_catagories, {"400"}, 0.997, 1.009)
  ({"et"}, {"2016"}, btag_catagories, {"400"}, 0.969, 1.025)
  ({"et"}, {"2017"}, btag_catagories, {"400"}, 0.952, 1.041)
  ({"et"}, {"2018"}, btag_catagories, {"400"}, 1.0, 1.0)
  ({"mt"}, {"2016"}, btag_catagories, {"400"}, 0.968, 1.027)
  ({"mt"}, {"2017"}, btag_catagories, {"400"}, 0.964, 1.036)
  ({"mt"}, {"2018"}, btag_catagories, {"400"}, 1.0, 1.0)
  ({"tt"}, {"2016"}, btag_catagories, {"400"}, 0.978, 1.024)
  ({"tt"}, {"2017"}, btag_catagories, {"400"}, 0.969, 1.037)
  ({"tt"}, {"2018"}, btag_catagories, {"400"}, 1.0, 1.0)
  ({"em"}, {"2016"}, nobtag_catagories, {"400"}, 1.0, 1.0)
  ({"em"}, {"2017"}, nobtag_catagories, {"400"}, 1.0, 1.0)
  ({"em"}, {"2018"}, nobtag_catagories, {"400"}, 1.0, 0.999)
  ({"et"}, {"2016"}, nobtag_catagories, {"400"}, 1.001, 0.999)
  ({"et"}, {"2017"}, nobtag_catagories, {"400"}, 1.002, 0.998)
  ({"et"}, {"2018"}, nobtag_catagories, {"400"}, 1.0, 1.0)
  ({"mt"}, {"2016"}, nobtag_catagories, {"400"}, 1.001, 0.999)
  ({"mt"}, {"2017"}, nobtag_catagories, {"400"}, 1.002, 0.998)
  ({"mt"}, {"2018"}, nobtag_catagories, {"400"}, 1.0, 1.0)
  ({"tt"}, {"2016"}, nobtag_catagories, {"400"}, 1.001, 0.999)
  ({"tt"}, {"2017"}, nobtag_catagories, {"400"}, 1.002, 0.998)
  ({"tt"}, {"2018"}, nobtag_catagories, {"400"}, 1.0, 1.0)
  );
  
  cb.cp().process(mssm_ggH_signals).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"450"}, 0.991, 1.008)
  ({"em"}, {"2017"}, btag_catagories, {"450"}, 0.991, 1.01)
  ({"em"}, {"2018"}, btag_catagories, {"450"}, 0.997, 1.008)
  ({"et"}, {"2016"}, btag_catagories, {"450"}, 0.97, 1.022)
  ({"et"}, {"2017"}, btag_catagories, {"450"}, 0.964, 1.035)
  ({"et"}, {"2018"}, btag_catagories, {"450"}, 0.119, 0.12)
  ({"mt"}, {"2016"}, btag_catagories, {"450"}, 0.984, 1.019)
  ({"mt"}, {"2017"}, btag_catagories, {"450"}, 0.967, 1.041)
  ({"mt"}, {"2018"}, btag_catagories, {"450"}, 0.126, 0.132)
  ({"tt"}, {"2016"}, btag_catagories, {"450"}, 0.974, 1.026)
  ({"tt"}, {"2017"}, btag_catagories, {"450"}, 0.962, 1.033)
  ({"tt"}, {"2018"}, btag_catagories, {"450"}, 0.622, 0.646)
  ({"em"}, {"2016"}, nobtag_catagories, {"450"}, 1.0, 1.0)
  ({"em"}, {"2017"}, nobtag_catagories, {"450"}, 1.0, 1.0)
  ({"em"}, {"2018"}, nobtag_catagories, {"450"}, 1.0, 0.999)
  ({"et"}, {"2016"}, nobtag_catagories, {"450"}, 1.001, 0.999)
  ({"et"}, {"2017"}, nobtag_catagories, {"450"}, 1.002, 0.998)
  ({"et"}, {"2018"}, nobtag_catagories, {"450"}, 0.12, 0.12)
  ({"mt"}, {"2016"}, nobtag_catagories, {"450"}, 1.001, 0.999)
  ({"mt"}, {"2017"}, nobtag_catagories, {"450"}, 1.002, 0.998)
  ({"mt"}, {"2018"}, nobtag_catagories, {"450"}, 0.117, 0.117)
  ({"tt"}, {"2016"}, nobtag_catagories, {"450"}, 1.001, 0.999)
  ({"tt"}, {"2017"}, nobtag_catagories, {"450"}, 1.002, 0.998)
  ({"tt"}, {"2018"}, nobtag_catagories, {"450"}, 0.636, 0.635)
  );
  
  cb.cp().process(mssm_ggH_signals).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"500"}, 0.997, 1.016)
  ({"em"}, {"2017"}, btag_catagories, {"500"}, 0.997, 1.016)
  ({"em"}, {"2018"}, btag_catagories, {"500"}, 0.992, 1.009)
  ({"et"}, {"2016"}, btag_catagories, {"500"}, 0.981, 1.033)
  ({"et"}, {"2017"}, btag_catagories, {"500"}, 0.944, 1.047)
  ({"et"}, {"2018"}, btag_catagories, {"500"}, 1.0, 1.0)
  ({"mt"}, {"2016"}, btag_catagories, {"500"}, 0.973, 1.022)
  ({"mt"}, {"2017"}, btag_catagories, {"500"}, 0.976, 1.017)
  ({"mt"}, {"2018"}, btag_catagories, {"500"}, 1.0, 1.0)
  ({"tt"}, {"2016"}, btag_catagories, {"500"}, 0.977, 1.027)
  ({"tt"}, {"2017"}, btag_catagories, {"500"}, 0.969, 1.045)
  ({"tt"}, {"2018"}, btag_catagories, {"500"}, 1.0, 1.0)
  ({"em"}, {"2016"}, nobtag_catagories, {"500"}, 1.0, 0.999)
  ({"em"}, {"2017"}, nobtag_catagories, {"500"}, 1.0, 0.999)
  ({"em"}, {"2018"}, nobtag_catagories, {"500"}, 1.0, 0.999)
  ({"et"}, {"2016"}, nobtag_catagories, {"500"}, 1.001, 0.999)
  ({"et"}, {"2017"}, nobtag_catagories, {"500"}, 1.003, 0.997)
  ({"et"}, {"2018"}, nobtag_catagories, {"500"}, 1.0, 1.0)
  ({"mt"}, {"2016"}, nobtag_catagories, {"500"}, 1.001, 0.999)
  ({"mt"}, {"2017"}, nobtag_catagories, {"500"}, 1.001, 0.999)
  ({"mt"}, {"2018"}, nobtag_catagories, {"500"}, 1.0, 1.0)
  ({"tt"}, {"2016"}, nobtag_catagories, {"500"}, 1.001, 0.999)
  ({"tt"}, {"2017"}, nobtag_catagories, {"500"}, 1.002, 0.998)
  ({"tt"}, {"2018"}, nobtag_catagories, {"500"}, 1.0, 1.0)
  );
  
  cb.cp().process(mssm_ggH_signals).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"600"}, 1.0, 1.018)
  ({"em"}, {"2017"}, btag_catagories, {"600"}, 1.0, 1.0)
  ({"em"}, {"2018"}, btag_catagories, {"600"}, 1.0, 1.007)
  ({"et"}, {"2016"}, btag_catagories, {"600"}, 0.985, 1.022)
  ({"et"}, {"2017"}, btag_catagories, {"600"}, 0.962, 1.029)
  ({"et"}, {"2018"}, btag_catagories, {"600"}, 0.991, 1.097)
  ({"mt"}, {"2016"}, btag_catagories, {"600"}, 0.992, 1.025)
  ({"mt"}, {"2017"}, btag_catagories, {"600"}, 0.966, 1.013)
  ({"mt"}, {"2018"}, btag_catagories, {"600"}, 0.981, 1.031)
  ({"tt"}, {"2016"}, btag_catagories, {"600"}, 0.973, 1.015)
  ({"tt"}, {"2017"}, btag_catagories, {"600"}, 0.983, 1.038)
  ({"tt"}, {"2018"}, btag_catagories, {"600"}, 0.983, 1.013)
  ({"em"}, {"2016"}, nobtag_catagories, {"600"}, 1.0, 0.999)
  ({"em"}, {"2017"}, nobtag_catagories, {"600"}, 1.0, 1.0)
  ({"em"}, {"2018"}, nobtag_catagories, {"600"}, 1.0, 1.0)
  ({"et"}, {"2016"}, nobtag_catagories, {"600"}, 1.001, 0.999)
  ({"et"}, {"2017"}, nobtag_catagories, {"600"}, 1.002, 0.997)
  ({"et"}, {"2018"}, nobtag_catagories, {"600"}, 1.001, 0.994)
  ({"mt"}, {"2016"}, nobtag_catagories, {"600"}, 1.0, 0.999)
  ({"mt"}, {"2017"}, nobtag_catagories, {"600"}, 1.002, 1.0)
  ({"mt"}, {"2018"}, nobtag_catagories, {"600"}, 1.002, 0.998)
  ({"tt"}, {"2016"}, nobtag_catagories, {"600"}, 1.001, 0.999)
  ({"tt"}, {"2017"}, nobtag_catagories, {"600"}, 1.001, 0.998)
  ({"tt"}, {"2018"}, nobtag_catagories, {"600"}, 1.001, 0.999)
  );
  
  cb.cp().process(mssm_ggH_signals).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"700"}, 0.975, 1.02)
  ({"em"}, {"2017"}, btag_catagories, {"700"}, 0.971, 1.008)
  ({"em"}, {"2018"}, btag_catagories, {"700"}, 0.988, 1.014)
  ({"et"}, {"2016"}, btag_catagories, {"700"}, 0.976, 1.016)
  ({"et"}, {"2017"}, btag_catagories, {"700"}, 0.968, 1.025)
  ({"et"}, {"2018"}, btag_catagories, {"700"}, 0.983, 1.025)
  ({"mt"}, {"2016"}, btag_catagories, {"700"}, 0.964, 1.028)
  ({"mt"}, {"2017"}, btag_catagories, {"700"}, 0.955, 1.052)
  ({"mt"}, {"2018"}, btag_catagories, {"700"}, 0.972, 1.017)
  ({"tt"}, {"2016"}, btag_catagories, {"700"}, 0.983, 1.028)
  ({"tt"}, {"2017"}, btag_catagories, {"700"}, 0.981, 1.036)
  ({"tt"}, {"2018"}, btag_catagories, {"700"}, 0.976, 1.017)
  ({"em"}, {"2016"}, nobtag_catagories, {"700"}, 1.001, 0.999)
  ({"em"}, {"2017"}, nobtag_catagories, {"700"}, 1.002, 1.0)
  ({"em"}, {"2018"}, nobtag_catagories, {"700"}, 1.001, 0.999)
  ({"et"}, {"2016"}, nobtag_catagories, {"700"}, 1.001, 0.999)
  ({"et"}, {"2017"}, nobtag_catagories, {"700"}, 1.002, 0.998)
  ({"et"}, {"2018"}, nobtag_catagories, {"700"}, 1.001, 0.999)
  ({"mt"}, {"2016"}, nobtag_catagories, {"700"}, 1.002, 0.998)
  ({"mt"}, {"2017"}, nobtag_catagories, {"700"}, 1.002, 0.997)
  ({"mt"}, {"2018"}, nobtag_catagories, {"700"}, 1.002, 0.999)
  ({"tt"}, {"2016"}, nobtag_catagories, {"700"}, 1.001, 0.999)
  ({"tt"}, {"2017"}, nobtag_catagories, {"700"}, 1.001, 0.998)
  ({"tt"}, {"2018"}, nobtag_catagories, {"700"}, 1.002, 0.999)
  );
  
  cb.cp().process(mssm_ggH_signals).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"800"}, 0.992, 1.0)
  ({"em"}, {"2017"}, btag_catagories, {"800"}, 1.0, 1.007)
  ({"em"}, {"2018"}, btag_catagories, {"800"}, 1.0, 1.0)
  ({"et"}, {"2016"}, btag_catagories, {"800"}, 0.976, 1.04)
  ({"et"}, {"2017"}, btag_catagories, {"800"}, 0.941, 1.045)
  ({"et"}, {"2018"}, btag_catagories, {"800"}, 0.976, 1.026)
  ({"mt"}, {"2016"}, btag_catagories, {"800"}, 0.974, 1.03)
  ({"mt"}, {"2017"}, btag_catagories, {"800"}, 0.969, 1.051)
  ({"mt"}, {"2018"}, btag_catagories, {"800"}, 0.986, 1.026)
  ({"tt"}, {"2016"}, btag_catagories, {"800"}, 0.977, 1.021)
  ({"tt"}, {"2017"}, btag_catagories, {"800"}, 0.971, 1.027)
  ({"tt"}, {"2018"}, btag_catagories, {"800"}, 0.983, 1.023)
  ({"em"}, {"2016"}, nobtag_catagories, {"800"}, 1.0, 1.0)
  ({"em"}, {"2017"}, nobtag_catagories, {"800"}, 1.0, 1.0)
  ({"em"}, {"2018"}, nobtag_catagories, {"800"}, 1.0, 1.0)
  ({"et"}, {"2016"}, nobtag_catagories, {"800"}, 1.001, 0.998)
  ({"et"}, {"2017"}, nobtag_catagories, {"800"}, 1.004, 0.998)
  ({"et"}, {"2018"}, nobtag_catagories, {"800"}, 1.002, 0.998)
  ({"mt"}, {"2016"}, nobtag_catagories, {"800"}, 1.001, 0.998)
  ({"mt"}, {"2017"}, nobtag_catagories, {"800"}, 1.002, 0.997)
  ({"mt"}, {"2018"}, nobtag_catagories, {"800"}, 1.001, 0.998)
  ({"tt"}, {"2016"}, nobtag_catagories, {"800"}, 1.001, 0.999)
  ({"tt"}, {"2017"}, nobtag_catagories, {"800"}, 1.002, 0.998)
  ({"tt"}, {"2018"}, nobtag_catagories, {"800"}, 1.001, 0.998)
  );
  
  cb.cp().process(mssm_ggH_signals).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"900"}, 0.997, 1.0)
  ({"em"}, {"2017"}, btag_catagories, {"900"}, 0.999, 1.008)
  ({"em"}, {"2018"}, btag_catagories, {"900"}, 0.994, 1.005)
  ({"et"}, {"2016"}, btag_catagories, {"900"}, 0.979, 1.024)
  ({"et"}, {"2017"}, btag_catagories, {"900"}, 0.981, 1.047)
  ({"et"}, {"2018"}, btag_catagories, {"900"}, 0.989, 1.015)
  ({"mt"}, {"2016"}, btag_catagories, {"900"}, 0.978, 1.024)
  ({"mt"}, {"2017"}, btag_catagories, {"900"}, 0.954, 1.042)
  ({"mt"}, {"2018"}, btag_catagories, {"900"}, 0.978, 1.042)
  ({"tt"}, {"2016"}, btag_catagories, {"900"}, 0.968, 1.022)
  ({"tt"}, {"2017"}, btag_catagories, {"900"}, 0.972, 1.042)
  ({"tt"}, {"2018"}, btag_catagories, {"900"}, 0.983, 1.018)
  ({"em"}, {"2016"}, nobtag_catagories, {"900"}, 1.0, 1.0)
  ({"em"}, {"2017"}, nobtag_catagories, {"900"}, 1.0, 1.0)
  ({"em"}, {"2018"}, nobtag_catagories, {"900"}, 1.0, 1.0)
  ({"et"}, {"2016"}, nobtag_catagories, {"900"}, 1.001, 0.999)
  ({"et"}, {"2017"}, nobtag_catagories, {"900"}, 1.001, 0.997)
  ({"et"}, {"2018"}, nobtag_catagories, {"900"}, 1.001, 0.999)
  ({"mt"}, {"2016"}, nobtag_catagories, {"900"}, 1.001, 0.999)
  ({"mt"}, {"2017"}, nobtag_catagories, {"900"}, 1.003, 0.997)
  ({"mt"}, {"2018"}, nobtag_catagories, {"900"}, 1.002, 0.997)
  ({"tt"}, {"2016"}, nobtag_catagories, {"900"}, 1.002, 0.999)
  ({"tt"}, {"2017"}, nobtag_catagories, {"900"}, 1.002, 0.997)
  ({"tt"}, {"2018"}, nobtag_catagories, {"900"}, 1.002, 0.998)
  );
  
  cb.cp().process(mssm_ggH_signals).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"1000"}, 0.998, 1.013)
  ({"em"}, {"2017"}, btag_catagories, {"1000"}, 0.992, 1.008)
  ({"em"}, {"2018"}, btag_catagories, {"1000"}, 0.988, 1.006)
  ({"et"}, {"2016"}, btag_catagories, {"1000"}, 0.986, 1.031)
  ({"et"}, {"2017"}, btag_catagories, {"1000"}, 0.98, 1.047)
  ({"et"}, {"2018"}, btag_catagories, {"1000"}, 0.977, 1.02)
  ({"mt"}, {"2016"}, btag_catagories, {"1000"}, 0.981, 1.019)
  ({"mt"}, {"2017"}, btag_catagories, {"1000"}, 0.966, 1.046)
  ({"mt"}, {"2018"}, btag_catagories, {"1000"}, 0.981, 1.036)
  ({"tt"}, {"2016"}, btag_catagories, {"1000"}, 0.977, 1.018)
  ({"tt"}, {"2017"}, btag_catagories, {"1000"}, 0.955, 1.043)
  ({"tt"}, {"2018"}, btag_catagories, {"1000"}, 0.98, 1.024)
  ({"em"}, {"2016"}, nobtag_catagories, {"1000"}, 1.0, 0.999)
  ({"em"}, {"2017"}, nobtag_catagories, {"1000"}, 1.001, 0.999)
  ({"em"}, {"2018"}, nobtag_catagories, {"1000"}, 1.001, 1.0)
  ({"et"}, {"2016"}, nobtag_catagories, {"1000"}, 1.001, 0.998)
  ({"et"}, {"2017"}, nobtag_catagories, {"1000"}, 1.002, 0.997)
  ({"et"}, {"2018"}, nobtag_catagories, {"1000"}, 1.002, 0.998)
  ({"mt"}, {"2016"}, nobtag_catagories, {"1000"}, 1.001, 0.999)
  ({"mt"}, {"2017"}, nobtag_catagories, {"1000"}, 1.002, 0.997)
  ({"mt"}, {"2018"}, nobtag_catagories, {"1000"}, 1.001, 0.997)
  ({"tt"}, {"2016"}, nobtag_catagories, {"1000"}, 1.001, 0.999)
  ({"tt"}, {"2017"}, nobtag_catagories, {"1000"}, 1.003, 0.997)
  ({"tt"}, {"2018"}, nobtag_catagories, {"1000"}, 1.002, 0.998)
  );
  
  cb.cp().process(mssm_ggH_signals).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"1200"}, 0.99, 1.01)
  ({"em"}, {"2017"}, btag_catagories, {"1200"}, 1.0, 1.0)
  ({"em"}, {"2018"}, btag_catagories, {"1200"}, 0.995, 1.007)
  ({"et"}, {"2016"}, btag_catagories, {"1200"}, 0.963, 1.016)
  ({"et"}, {"2017"}, btag_catagories, {"1200"}, 0.986, 1.039)
  ({"et"}, {"2018"}, btag_catagories, {"1200"}, 0.986, 1.031)
  ({"mt"}, {"2016"}, btag_catagories, {"1200"}, 0.976, 1.03)
  ({"mt"}, {"2017"}, btag_catagories, {"1200"}, 0.96, 1.014)
  ({"mt"}, {"2018"}, btag_catagories, {"1200"}, 0.973, 1.019)
  ({"tt"}, {"2016"}, btag_catagories, {"1200"}, 0.974, 1.017)
  ({"tt"}, {"2017"}, btag_catagories, {"1200"}, 0.964, 1.035)
  ({"tt"}, {"2018"}, btag_catagories, {"1200"}, 0.981, 1.028)
  ({"em"}, {"2016"}, nobtag_catagories, {"1200"}, 1.001, 0.999)
  ({"em"}, {"2017"}, nobtag_catagories, {"1200"}, 1.0, 1.0)
  ({"em"}, {"2018"}, nobtag_catagories, {"1200"}, 1.0, 0.999)
  ({"et"}, {"2016"}, nobtag_catagories, {"1200"}, 1.002, 0.999)
  ({"et"}, {"2017"}, nobtag_catagories, {"1200"}, 1.001, 0.997)
  ({"et"}, {"2018"}, nobtag_catagories, {"1200"}, 1.001, 0.997)
  ({"mt"}, {"2016"}, nobtag_catagories, {"1200"}, 1.001, 0.998)
  ({"mt"}, {"2017"}, nobtag_catagories, {"1200"}, 1.003, 0.999)
  ({"mt"}, {"2018"}, nobtag_catagories, {"1200"}, 1.002, 0.998)
  ({"tt"}, {"2016"}, nobtag_catagories, {"1200"}, 1.002, 0.999)
  ({"tt"}, {"2017"}, nobtag_catagories, {"1200"}, 1.003, 0.997)
  ({"tt"}, {"2018"}, nobtag_catagories, {"1200"}, 1.002, 0.998)
  );
  
  cb.cp().process(mssm_ggH_signals).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"1400"}, 0.969, 0.991)
  ({"em"}, {"2017"}, btag_catagories, {"1400"}, 0.987, 1.013)
  ({"em"}, {"2018"}, btag_catagories, {"1400"}, 1.0, 1.003)
  ({"et"}, {"2016"}, btag_catagories, {"1400"}, 0.978, 1.003)
  ({"et"}, {"2017"}, btag_catagories, {"1400"}, 0.962, 0.995)
  ({"et"}, {"2018"}, btag_catagories, {"1400"}, 0.945, 1.02)
  ({"mt"}, {"2016"}, btag_catagories, {"1400"}, 0.981, 1.015)
  ({"mt"}, {"2017"}, btag_catagories, {"1400"}, 0.956, 1.036)
  ({"mt"}, {"2018"}, btag_catagories, {"1400"}, 0.977, 1.017)
  ({"tt"}, {"2016"}, btag_catagories, {"1400"}, 0.983, 1.038)
  ({"tt"}, {"2017"}, btag_catagories, {"1400"}, 0.966, 1.041)
  ({"tt"}, {"2018"}, btag_catagories, {"1400"}, 0.992, 1.014)
  ({"em"}, {"2016"}, nobtag_catagories, {"1400"}, 1.002, 1.001)
  ({"em"}, {"2017"}, nobtag_catagories, {"1400"}, 1.001, 0.999)
  ({"em"}, {"2018"}, nobtag_catagories, {"1400"}, 1.0, 1.0)
  ({"et"}, {"2016"}, nobtag_catagories, {"1400"}, 1.001, 1.0)
  ({"et"}, {"2017"}, nobtag_catagories, {"1400"}, 1.004, 0.999)
  ({"et"}, {"2018"}, nobtag_catagories, {"1400"}, 1.005, 0.998)
  ({"mt"}, {"2016"}, nobtag_catagories, {"1400"}, 1.001, 0.999)
  ({"mt"}, {"2017"}, nobtag_catagories, {"1400"}, 1.003, 0.997)
  ({"mt"}, {"2018"}, nobtag_catagories, {"1400"}, 1.002, 0.998)
  ({"tt"}, {"2016"}, nobtag_catagories, {"1400"}, 1.001, 0.997)
  ({"tt"}, {"2017"}, nobtag_catagories, {"1400"}, 1.003, 0.997)
  ({"tt"}, {"2018"}, nobtag_catagories, {"1400"}, 1.001, 0.999)
  );
  
  cb.cp().process(mssm_ggH_signals).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"1600"}, 1.0, 1.001)
  ({"em"}, {"2017"}, btag_catagories, {"1600"}, 0.996, 1.012)
  ({"em"}, {"2018"}, btag_catagories, {"1600"}, 0.988, 0.996)
  ({"et"}, {"2016"}, btag_catagories, {"1600"}, 0.97, 1.026)
  ({"et"}, {"2017"}, btag_catagories, {"1600"}, 0.961, 1.04)
  ({"et"}, {"2018"}, btag_catagories, {"1600"}, 0.974, 1.021)
  ({"mt"}, {"2016"}, btag_catagories, {"1600"}, 0.973, 1.04)
  ({"mt"}, {"2017"}, btag_catagories, {"1600"}, 0.972, 1.048)
  ({"mt"}, {"2018"}, btag_catagories, {"1600"}, 0.971, 1.026)
  ({"tt"}, {"2016"}, btag_catagories, {"1600"}, 0.984, 1.013)
  ({"tt"}, {"2017"}, btag_catagories, {"1600"}, 0.948, 1.045)
  ({"tt"}, {"2018"}, btag_catagories, {"1600"}, 0.979, 1.034)
  ({"em"}, {"2016"}, nobtag_catagories, {"1600"}, 1.0, 1.0)
  ({"em"}, {"2017"}, nobtag_catagories, {"1600"}, 1.0, 0.999)
  ({"em"}, {"2018"}, nobtag_catagories, {"1600"}, 1.001, 1.0)
  ({"et"}, {"2016"}, nobtag_catagories, {"1600"}, 1.002, 0.998)
  ({"et"}, {"2017"}, nobtag_catagories, {"1600"}, 1.002, 0.997)
  ({"et"}, {"2018"}, nobtag_catagories, {"1600"}, 1.002, 0.998)
  ({"mt"}, {"2016"}, nobtag_catagories, {"1600"}, 1.002, 0.997)
  ({"mt"}, {"2017"}, nobtag_catagories, {"1600"}, 1.002, 0.997)
  ({"mt"}, {"2018"}, nobtag_catagories, {"1600"}, 1.003, 0.997)
  ({"tt"}, {"2016"}, nobtag_catagories, {"1600"}, 1.001, 0.999)
  ({"tt"}, {"2017"}, nobtag_catagories, {"1600"}, 1.004, 0.996)
  ({"tt"}, {"2018"}, nobtag_catagories, {"1600"}, 1.002, 0.997)
  );
  
  cb.cp().process(mssm_ggH_signals).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"1800"}, 0.993, 1.0)
  ({"em"}, {"2017"}, btag_catagories, {"1800"}, 0.981, 1.003)
  ({"em"}, {"2018"}, btag_catagories, {"1800"}, 0.988, 1.004)
  ({"et"}, {"2016"}, btag_catagories, {"1800"}, 0.978, 1.018)
  ({"et"}, {"2017"}, btag_catagories, {"1800"}, 0.967, 1.025)
  ({"et"}, {"2018"}, btag_catagories, {"1800"}, 0.977, 1.018)
  ({"mt"}, {"2016"}, btag_catagories, {"1800"}, 0.98, 1.033)
  ({"mt"}, {"2017"}, btag_catagories, {"1800"}, 0.952, 1.028)
  ({"mt"}, {"2018"}, btag_catagories, {"1800"}, 0.974, 1.026)
  ({"tt"}, {"2016"}, btag_catagories, {"1800"}, 0.975, 1.012)
  ({"tt"}, {"2017"}, btag_catagories, {"1800"}, 0.942, 1.031)
  ({"tt"}, {"2018"}, btag_catagories, {"1800"}, 0.983, 1.025)
  ({"em"}, {"2016"}, nobtag_catagories, {"1800"}, 1.001, 1.0)
  ({"em"}, {"2017"}, nobtag_catagories, {"1800"}, 1.001, 1.0)
  ({"em"}, {"2018"}, nobtag_catagories, {"1800"}, 1.001, 1.0)
  ({"et"}, {"2016"}, nobtag_catagories, {"1800"}, 1.001, 0.999)
  ({"et"}, {"2017"}, nobtag_catagories, {"1800"}, 1.003, 0.999)
  ({"et"}, {"2018"}, nobtag_catagories, {"1800"}, 1.002, 0.998)
  ({"mt"}, {"2016"}, nobtag_catagories, {"1800"}, 1.001, 0.998)
  ({"mt"}, {"2017"}, nobtag_catagories, {"1800"}, 1.004, 0.998)
  ({"mt"}, {"2018"}, nobtag_catagories, {"1800"}, 1.002, 0.997)
  ({"tt"}, {"2016"}, nobtag_catagories, {"1800"}, 1.002, 0.999)
  ({"tt"}, {"2017"}, nobtag_catagories, {"1800"}, 1.006, 0.997)
  ({"tt"}, {"2018"}, nobtag_catagories, {"1800"}, 1.002, 0.997)
  );
  
  cb.cp().process(mssm_ggH_signals).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"2000"}, 1.0, 1.008)
  ({"em"}, {"2017"}, btag_catagories, {"2000"}, 0.968, 1.016)
  ({"em"}, {"2018"}, btag_catagories, {"2000"}, 0.996, 1.006)
  ({"et"}, {"2016"}, btag_catagories, {"2000"}, 0.981, 1.015)
  ({"et"}, {"2017"}, btag_catagories, {"2000"}, 0.97, 1.019)
  ({"et"}, {"2018"}, btag_catagories, {"2000"}, 0.98, 1.024)
  ({"mt"}, {"2016"}, btag_catagories, {"2000"}, 0.988, 1.024)
  ({"mt"}, {"2017"}, btag_catagories, {"2000"}, 0.972, 1.047)
  ({"mt"}, {"2018"}, btag_catagories, {"2000"}, 0.982, 1.015)
  ({"tt"}, {"2016"}, btag_catagories, {"2000"}, 0.964, 1.02)
  ({"tt"}, {"2017"}, btag_catagories, {"2000"}, 0.963, 1.035)
  ({"tt"}, {"2018"}, btag_catagories, {"2000"}, 0.974, 1.024)
  ({"em"}, {"2016"}, nobtag_catagories, {"2000"}, 1.0, 0.999)
  ({"em"}, {"2017"}, nobtag_catagories, {"2000"}, 1.003, 0.999)
  ({"em"}, {"2018"}, nobtag_catagories, {"2000"}, 1.0, 0.999)
  ({"et"}, {"2016"}, nobtag_catagories, {"2000"}, 1.001, 0.999)
  ({"et"}, {"2017"}, nobtag_catagories, {"2000"}, 1.002, 0.998)
  ({"et"}, {"2018"}, nobtag_catagories, {"2000"}, 1.002, 0.998)
  ({"mt"}, {"2016"}, nobtag_catagories, {"2000"}, 1.001, 0.998)
  ({"mt"}, {"2017"}, nobtag_catagories, {"2000"}, 1.002, 0.996)
  ({"mt"}, {"2018"}, nobtag_catagories, {"2000"}, 1.002, 0.998)
  ({"tt"}, {"2016"}, nobtag_catagories, {"2000"}, 1.003, 0.998)
  ({"tt"}, {"2017"}, nobtag_catagories, {"2000"}, 1.004, 0.997)
  ({"tt"}, {"2018"}, nobtag_catagories, {"2000"}, 1.003, 0.997)
  );
  
  cb.cp().process(mssm_ggH_signals).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"2300"}, 0.991, 1.004)
  ({"em"}, {"2017"}, btag_catagories, {"2300"}, 0.979, 1.021)
  ({"em"}, {"2018"}, btag_catagories, {"2300"}, 0.992, 1.003)
  ({"et"}, {"2016"}, btag_catagories, {"2300"}, 0.976, 1.035)
  ({"et"}, {"2017"}, btag_catagories, {"2300"}, 0.963, 1.04)
  ({"et"}, {"2018"}, btag_catagories, {"2300"}, 0.979, 1.009)
  ({"mt"}, {"2016"}, btag_catagories, {"2300"}, 0.986, 1.017)
  ({"mt"}, {"2017"}, btag_catagories, {"2300"}, 0.963, 1.021)
  ({"mt"}, {"2018"}, btag_catagories, {"2300"}, 0.979, 1.014)
  ({"tt"}, {"2016"}, btag_catagories, {"2300"}, 0.976, 1.021)
  ({"tt"}, {"2017"}, btag_catagories, {"2300"}, 0.951, 1.044)
  ({"tt"}, {"2018"}, btag_catagories, {"2300"}, 0.982, 1.017)
  ({"em"}, {"2016"}, nobtag_catagories, {"2300"}, 1.001, 1.0)
  ({"em"}, {"2017"}, nobtag_catagories, {"2300"}, 1.002, 0.998)
  ({"em"}, {"2018"}, nobtag_catagories, {"2300"}, 1.001, 1.0)
  ({"et"}, {"2016"}, nobtag_catagories, {"2300"}, 1.002, 0.997)
  ({"et"}, {"2017"}, nobtag_catagories, {"2300"}, 1.003, 0.997)
  ({"et"}, {"2018"}, nobtag_catagories, {"2300"}, 1.002, 0.999)
  ({"mt"}, {"2016"}, nobtag_catagories, {"2300"}, 1.001, 0.999)
  ({"mt"}, {"2017"}, nobtag_catagories, {"2300"}, 1.003, 0.998)
  ({"mt"}, {"2018"}, nobtag_catagories, {"2300"}, 1.002, 0.998)
  ({"tt"}, {"2016"}, nobtag_catagories, {"2300"}, 1.002, 0.998)
  ({"tt"}, {"2017"}, nobtag_catagories, {"2300"}, 1.004, 0.996)
  ({"tt"}, {"2018"}, nobtag_catagories, {"2300"}, 1.002, 0.998)
  );
  
  cb.cp().process(mssm_ggH_signals).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"2600"}, 0.984, 1.0)
  ({"em"}, {"2017"}, btag_catagories, {"2600"}, 0.993, 1.029)
  ({"em"}, {"2018"}, btag_catagories, {"2600"}, 0.996, 1.011)
  ({"et"}, {"2016"}, btag_catagories, {"2600"}, 0.985, 1.016)
  ({"et"}, {"2017"}, btag_catagories, {"2600"}, 0.969, 1.029)
  ({"et"}, {"2018"}, btag_catagories, {"2600"}, 0.997, 1.022)
  ({"mt"}, {"2016"}, btag_catagories, {"2600"}, 0.988, 1.024)
  ({"mt"}, {"2017"}, btag_catagories, {"2600"}, 0.97, 1.037)
  ({"mt"}, {"2018"}, btag_catagories, {"2600"}, 0.984, 1.02)
  ({"tt"}, {"2016"}, btag_catagories, {"2600"}, 0.968, 1.029)
  ({"tt"}, {"2017"}, btag_catagories, {"2600"}, 0.948, 1.04)
  ({"tt"}, {"2018"}, btag_catagories, {"2600"}, 0.976, 1.025)
  ({"em"}, {"2016"}, nobtag_catagories, {"2600"}, 1.001, 1.0)
  ({"em"}, {"2017"}, nobtag_catagories, {"2600"}, 1.001, 0.997)
  ({"em"}, {"2018"}, nobtag_catagories, {"2600"}, 1.001, 0.999)
  ({"et"}, {"2016"}, nobtag_catagories, {"2600"}, 1.001, 0.999)
  ({"et"}, {"2017"}, nobtag_catagories, {"2600"}, 1.002, 0.998)
  ({"et"}, {"2018"}, nobtag_catagories, {"2600"}, 1.0, 0.998)
  ({"mt"}, {"2016"}, nobtag_catagories, {"2600"}, 1.0, 0.998)
  ({"mt"}, {"2017"}, nobtag_catagories, {"2600"}, 1.003, 0.996)
  ({"mt"}, {"2018"}, nobtag_catagories, {"2600"}, 1.002, 0.998)
  ({"tt"}, {"2016"}, nobtag_catagories, {"2600"}, 1.003, 0.998)
  ({"tt"}, {"2017"}, nobtag_catagories, {"2600"}, 1.005, 0.996)
  ({"tt"}, {"2018"}, nobtag_catagories, {"2600"}, 1.003, 0.997)
  );
  
  cb.cp().process(mssm_ggH_signals).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"2900"}, 0.998, 1.01)
  ({"em"}, {"2017"}, btag_catagories, {"2900"}, 0.992, 1.006)
  ({"em"}, {"2018"}, btag_catagories, {"2900"}, 0.987, 1.013)
  ({"et"}, {"2016"}, btag_catagories, {"2900"}, 0.976, 1.047)
  ({"et"}, {"2017"}, btag_catagories, {"2900"}, 0.963, 1.037)
  ({"et"}, {"2018"}, btag_catagories, {"2900"}, 0.982, 1.013)
  ({"mt"}, {"2016"}, btag_catagories, {"2900"}, 0.974, 1.024)
  ({"mt"}, {"2017"}, btag_catagories, {"2900"}, 0.969, 1.062)
  ({"mt"}, {"2018"}, btag_catagories, {"2900"}, 0.977, 1.018)
  ({"tt"}, {"2016"}, btag_catagories, {"2900"}, 0.977, 1.019)
  ({"tt"}, {"2017"}, btag_catagories, {"2900"}, 0.966, 1.047)
  ({"tt"}, {"2018"}, btag_catagories, {"2900"}, 0.982, 1.016)
  ({"em"}, {"2016"}, nobtag_catagories, {"2900"}, 1.0, 0.999)
  ({"em"}, {"2017"}, nobtag_catagories, {"2900"}, 1.001, 0.999)
  ({"em"}, {"2018"}, nobtag_catagories, {"2900"}, 1.002, 0.998)
  ({"et"}, {"2016"}, nobtag_catagories, {"2900"}, 1.002, 0.996)
  ({"et"}, {"2017"}, nobtag_catagories, {"2900"}, 1.003, 0.997)
  ({"et"}, {"2018"}, nobtag_catagories, {"2900"}, 1.002, 0.999)
  ({"mt"}, {"2016"}, nobtag_catagories, {"2900"}, 1.003, 0.998)
  ({"mt"}, {"2017"}, nobtag_catagories, {"2900"}, 1.003, 0.994)
  ({"mt"}, {"2018"}, nobtag_catagories, {"2900"}, 1.003, 0.998)
  ({"tt"}, {"2016"}, nobtag_catagories, {"2900"}, 1.002, 0.998)
  ({"tt"}, {"2017"}, nobtag_catagories, {"2900"}, 1.003, 0.995)
  ({"tt"}, {"2018"}, nobtag_catagories, {"2900"}, 1.002, 0.998)
  );
  
  cb.cp().process(mssm_ggH_signals).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"3200"}, 0.985, 1.006)
  ({"em"}, {"2017"}, btag_catagories, {"3200"}, 0.991, 1.014)
  ({"em"}, {"2018"}, btag_catagories, {"3200"}, 0.998, 1.001)
  ({"et"}, {"2016"}, btag_catagories, {"3200"}, 0.987, 1.032)
  ({"et"}, {"2017"}, btag_catagories, {"3200"}, 0.957, 1.038)
  ({"et"}, {"2018"}, btag_catagories, {"3200"}, 0.967, 1.019)
  ({"mt"}, {"2016"}, btag_catagories, {"3200"}, 0.976, 1.022)
  ({"mt"}, {"2017"}, btag_catagories, {"3200"}, 0.972, 1.026)
  ({"mt"}, {"2018"}, btag_catagories, {"3200"}, 0.986, 1.025)
  ({"tt"}, {"2016"}, btag_catagories, {"3200"}, 0.979, 1.02)
  ({"tt"}, {"2017"}, btag_catagories, {"3200"}, 0.958, 1.037)
  ({"tt"}, {"2018"}, btag_catagories, {"3200"}, 0.976, 1.02)
  ({"em"}, {"2016"}, nobtag_catagories, {"3200"}, 1.001, 0.999)
  ({"em"}, {"2017"}, nobtag_catagories, {"3200"}, 1.001, 0.998)
  ({"em"}, {"2018"}, nobtag_catagories, {"3200"}, 1.0, 1.0)
  ({"et"}, {"2016"}, nobtag_catagories, {"3200"}, 1.001, 0.997)
  ({"et"}, {"2017"}, nobtag_catagories, {"3200"}, 1.004, 0.997)
  ({"et"}, {"2018"}, nobtag_catagories, {"3200"}, 1.004, 0.998)
  ({"mt"}, {"2016"}, nobtag_catagories, {"3200"}, 1.002, 0.998)
  ({"mt"}, {"2017"}, nobtag_catagories, {"3200"}, 1.003, 0.996)
  ({"mt"}, {"2018"}, nobtag_catagories, {"3200"}, 1.002, 0.997)
  ({"tt"}, {"2016"}, nobtag_catagories, {"3200"}, 1.002, 0.998)
  ({"tt"}, {"2017"}, nobtag_catagories, {"3200"}, 1.005, 0.996)
  ({"tt"}, {"2018"}, nobtag_catagories, {"3200"}, 1.003, 0.998)
  );
  
  cb.cp().process(mssm_ggH_signals).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"3500"}, 1.0, 1.01)
  ({"em"}, {"2017"}, btag_catagories, {"3500"}, 0.985, 1.0)
  ({"em"}, {"2018"}, btag_catagories, {"3500"}, 0.995, 1.002)
  ({"et"}, {"2016"}, btag_catagories, {"3500"}, 0.969, 1.027)
  ({"et"}, {"2017"}, btag_catagories, {"3500"}, 0.964, 1.041)
  ({"et"}, {"2018"}, btag_catagories, {"3500"}, 0.983, 1.016)
  ({"mt"}, {"2016"}, btag_catagories, {"3500"}, 0.975, 1.01)
  ({"mt"}, {"2017"}, btag_catagories, {"3500"}, 0.97, 1.044)
  ({"mt"}, {"2018"}, btag_catagories, {"3500"}, 0.968, 1.032)
  ({"tt"}, {"2016"}, btag_catagories, {"3500"}, 0.975, 1.021)
  ({"tt"}, {"2017"}, btag_catagories, {"3500"}, 0.962, 1.03)
  ({"tt"}, {"2018"}, btag_catagories, {"3500"}, 0.972, 1.028)
  ({"em"}, {"2016"}, nobtag_catagories, {"3500"}, 1.0, 0.999)
  ({"em"}, {"2017"}, nobtag_catagories, {"3500"}, 1.001, 1.0)
  ({"em"}, {"2018"}, nobtag_catagories, {"3500"}, 1.001, 1.0)
  ({"et"}, {"2016"}, nobtag_catagories, {"3500"}, 1.003, 0.998)
  ({"et"}, {"2017"}, nobtag_catagories, {"3500"}, 1.003, 0.995)
  ({"et"}, {"2018"}, nobtag_catagories, {"3500"}, 1.002, 0.998)
  ({"mt"}, {"2016"}, nobtag_catagories, {"3500"}, 1.002, 0.999)
  ({"mt"}, {"2017"}, nobtag_catagories, {"3500"}, 1.003, 0.996)
  ({"mt"}, {"2018"}, nobtag_catagories, {"3500"}, 1.004, 0.997)
  ({"tt"}, {"2016"}, nobtag_catagories, {"3500"}, 1.002, 0.998)
  ({"tt"}, {"2017"}, nobtag_catagories, {"3500"}, 1.004, 0.997)
  ({"tt"}, {"2018"}, nobtag_catagories, {"3500"}, 1.004, 0.997)
  );
  
  cb.cp().process({"W"}).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id>::init
  ({"em"}, {"2016"}, {35}, 1.0, 1.076)
  ({"em"}, {"2017"}, {35}, 1.0, 1.039)
  ({"em"}, {"2018"}, {35}, 0.945, 1.189)
  ({"em"}, {"2016"}, {36}, 1.0, 1.125)
  ({"em"}, {"2017"}, {36}, 0.956, 1.057)
  ({"em"}, {"2018"}, {36}, 0.86, 1.055)
  ({"em"}, {"2016"}, {37}, 0.753, 1.065)
  ({"em"}, {"2017"}, {37}, 0.929, 1.054)
  ({"em"}, {"2018"}, {37}, 0.932, 1.03)
  ({"em"}, {"2016"}, {32}, 1.0, 0.991)
  ({"em"}, {"2017"}, {32}, 1.0, 0.994)
  ({"em"}, {"2018"}, {32}, 1.008, 0.972)
  ({"em"}, {"2016"}, {33}, 1.0, 0.992)
  ({"em"}, {"2017"}, {33}, 1.004, 0.995)
  ({"em"}, {"2018"}, {33}, 1.01, 0.996)
  ({"em"}, {"2016"}, {34}, 1.013, 0.994)
  ({"em"}, {"2017"}, {34}, 1.004, 0.997)
  ({"em"}, {"2018"}, {34}, 1.004, 0.998)
  ({"em"}, {"2016"}, {2}, 1.001, 1.001)
  ({"em"}, {"2017"}, {2}, 1.0, 1.0)
  ({"em"}, {"2018"}, {2}, 1.0, 1.0)
  );
  
  cb.cp().process({"ZL"}).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id>::init
  ({"em"}, {"2016"}, {35}, 0.855, 1.0)
  ({"em"}, {"2017"}, {35}, 1.0, 1.023)
  ({"em"}, {"2018"}, {35}, 0.983, 1.139)
  ({"em"}, {"2016"}, {36}, 0.949, 1.019)
  ({"em"}, {"2017"}, {36}, 0.965, 1.135)
  ({"em"}, {"2018"}, {36}, 0.915, 1.086)
  ({"em"}, {"2016"}, {37}, 0.879, 1.025)
  ({"em"}, {"2017"}, {37}, 0.971, 1.046)
  ({"em"}, {"2018"}, {37}, 0.935, 1.064)
  ({"em"}, {"2016"}, {32}, 1.014, 1.0)
  ({"em"}, {"2017"}, {32}, 1.0, 0.997)
  ({"em"}, {"2018"}, {32}, 1.002, 0.983)
  ({"em"}, {"2016"}, {33}, 1.001, 0.998)
  ({"em"}, {"2017"}, {33}, 1.001, 0.995)
  ({"em"}, {"2018"}, {33}, 1.004, 0.996)
  ({"em"}, {"2016"}, {34}, 1.005, 1.0)
  ({"em"}, {"2017"}, {34}, 1.003, 0.998)
  ({"em"}, {"2018"}, {34}, 1.005, 0.994)
  ({"em"}, {"2016"}, {2}, 1.002, 1.002)
  ({"em"}, {"2017"}, {2}, 0.999, 0.999)
  ({"em"}, {"2018"}, {2}, 1.001, 1.001)
  ({"et"}, {"2016"}, {35}, 0.944, 1.074)
  ({"et"}, {"2017"}, {35}, 0.938, 1.093)
  ({"et"}, {"2018"}, {35}, 0.874, 1.103)
  ({"et"}, {"2016"}, {36}, 0.959, 1.055)
  ({"et"}, {"2017"}, {36}, 0.97, 1.063)
  ({"et"}, {"2018"}, {36}, 0.847, 1.05)
  ({"et"}, {"2016"}, {32}, 1.001, 0.998)
  ({"et"}, {"2017"}, {32}, 1.004, 0.999)
  ({"et"}, {"2018"}, {32}, 1.005, 0.996)
  ({"et"}, {"2016"}, {33}, 1.001, 0.998)
  ({"et"}, {"2017"}, {33}, 1.005, 1.001)
  ({"et"}, {"2018"}, {33}, 1.01, 0.997)
  ({"mt"}, {"2016"}, {35}, 0.926, 1.086)
  ({"mt"}, {"2017"}, {35}, 0.96, 1.098)
  ({"mt"}, {"2018"}, {35}, 0.942, 1.191)
  ({"mt"}, {"2016"}, {36}, 0.891, 1.114)
  ({"mt"}, {"2017"}, {36}, 0.931, 1.023)
  ({"mt"}, {"2018"}, {36}, 0.946, 1.093)
  ({"mt"}, {"2016"}, {32}, 1.002, 0.998)
  ({"mt"}, {"2017"}, {32}, 1.004, 0.998)
  ({"mt"}, {"2018"}, {32}, 1.002, 0.993)
  ({"mt"}, {"2016"}, {33}, 1.003, 0.996)
  ({"mt"}, {"2017"}, {33}, 1.006, 1.001)
  ({"mt"}, {"2018"}, {33}, 1.004, 0.994)
  ({"tt"}, {"2016"}, {32}, 0.949, 1.018)
  ({"tt"}, {"2017"}, {32}, 0.924, 1.098)
  ({"tt"}, {"2018"}, {32}, 0.867, 1.054)
  ({"tt"}, {"2016"}, {35}, 1.003, 0.999)
  ({"tt"}, {"2017"}, {35}, 1.007, 0.999)
  ({"tt"}, {"2018"}, {35}, 1.01, 0.996)
  );
  
  cb.cp().process({"TTL"}).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id>::init
  ({"em"}, {"2016"}, {35}, 0.999, 1.001)
  ({"em"}, {"2017"}, {35}, 0.998, 1.001)
  ({"em"}, {"2018"}, {35}, 0.998, 1.001)
  ({"em"}, {"2016"}, {36}, 0.999, 1.001)
  ({"em"}, {"2017"}, {36}, 0.998, 1.001)
  ({"em"}, {"2018"}, {36}, 0.998, 1.001)
  ({"em"}, {"2016"}, {37}, 0.999, 1.0)
  ({"em"}, {"2017"}, {37}, 0.998, 1.001)
  ({"em"}, {"2018"}, {37}, 0.998, 1.001)
  ({"em"}, {"2016"}, {32}, 1.005, 0.993)
  ({"em"}, {"2017"}, {32}, 1.008, 0.991)
  ({"em"}, {"2018"}, {32}, 1.013, 0.988)
  ({"em"}, {"2016"}, {33}, 1.004, 0.995)
  ({"em"}, {"2017"}, {33}, 1.006, 0.99)
  ({"em"}, {"2018"}, {33}, 1.01, 0.989)
  ({"em"}, {"2016"}, {34}, 1.004, 0.995)
  ({"em"}, {"2017"}, {34}, 1.007, 0.992)
  ({"em"}, {"2018"}, {34}, 1.01, 0.989)
  ({"em"}, {"2016"}, {2}, 1.0, 1.0)
  ({"em"}, {"2017"}, {2}, 1.0, 1.0)
  ({"em"}, {"2018"}, {2}, 1.0, 1.0)
  ({"et"}, {"2016"}, {35}, 0.973, 0.975)
  ({"et"}, {"2017"}, {35}, 0.999, 1.001)
  ({"et"}, {"2018"}, {35}, 0.999, 1.001)
  ({"et"}, {"2016"}, {36}, 0.97, 0.971)
  ({"et"}, {"2017"}, {36}, 0.999, 1.0)
  ({"et"}, {"2018"}, {36}, 0.999, 1.001)
  ({"et"}, {"2016"}, {32}, 0.976, 0.968)
  ({"et"}, {"2017"}, {32}, 1.006, 0.99)
  ({"et"}, {"2018"}, {32}, 1.009, 0.992)
  ({"et"}, {"2016"}, {33}, 0.977, 0.971)
  ({"et"}, {"2017"}, {33}, 1.007, 0.997)
  ({"et"}, {"2018"}, {33}, 1.009, 0.991)
  ({"mt"}, {"2016"}, {35}, 0.974, 0.975)
  ({"mt"}, {"2017"}, {35}, 0.999, 1.001)
  ({"mt"}, {"2018"}, {35}, 0.999, 1.001)
  ({"mt"}, {"2016"}, {36}, 0.973, 0.974)
  ({"mt"}, {"2017"}, {36}, 0.999, 1.001)
  ({"mt"}, {"2018"}, {36}, 0.999, 1.001)
  ({"mt"}, {"2016"}, {32}, 0.981, 0.973)
  ({"mt"}, {"2017"}, {32}, 1.005, 0.995)
  ({"mt"}, {"2018"}, {32}, 1.01, 0.99)
  ({"mt"}, {"2016"}, {33}, 0.978, 0.97)
  ({"mt"}, {"2017"}, {33}, 1.006, 0.995)
  ({"mt"}, {"2018"}, {33}, 1.01, 0.991)
  ({"tt"}, {"2016"}, {32}, 0.976, 0.976)
  ({"tt"}, {"2017"}, {32}, 1.0, 1.0)
  ({"tt"}, {"2018"}, {32}, 0.999, 1.001)
  ({"tt"}, {"2016"}, {35}, 0.982, 0.982)
  ({"tt"}, {"2017"}, {35}, 1.0, 1.0)
  ({"tt"}, {"2018"}, {35}, 1.011, 0.99)
  );
  
  cb.cp().process({"VVL"}).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id>::init
  ({"em"}, {"2016"}, {35}, 0.993, 1.008)
  ({"em"}, {"2017"}, {35}, 0.99, 1.009)
  ({"em"}, {"2018"}, {35}, 0.987, 1.013)
  ({"em"}, {"2016"}, {36}, 0.992, 1.006)
  ({"em"}, {"2017"}, {36}, 0.989, 1.01)
  ({"em"}, {"2018"}, {36}, 0.986, 1.014)
  ({"em"}, {"2016"}, {37}, 0.992, 1.008)
  ({"em"}, {"2017"}, {37}, 0.989, 1.012)
  ({"em"}, {"2018"}, {37}, 0.984, 1.015)
  ({"em"}, {"2016"}, {32}, 1.005, 0.993)
  ({"em"}, {"2017"}, {32}, 1.006, 0.99)
  ({"em"}, {"2018"}, {32}, 1.011, 0.986)
  ({"em"}, {"2016"}, {33}, 1.003, 0.996)
  ({"em"}, {"2017"}, {33}, 1.004, 0.994)
  ({"em"}, {"2018"}, {33}, 1.007, 0.992)
  ({"em"}, {"2016"}, {34}, 1.002, 0.997)
  ({"em"}, {"2017"}, {34}, 1.004, 0.995)
  ({"em"}, {"2018"}, {34}, 1.006, 0.994)
  ({"em"}, {"2016"}, {2}, 1.0, 1.0)
  ({"em"}, {"2017"}, {2}, 1.0, 1.0)
  ({"em"}, {"2018"}, {2}, 1.0, 1.0)
  ({"et"}, {"2016"}, {35}, 0.99, 1.009)
  ({"et"}, {"2017"}, {35}, 0.993, 1.006)
  ({"et"}, {"2018"}, {35}, 0.991, 1.012)
  ({"et"}, {"2016"}, {36}, 0.99, 1.007)
  ({"et"}, {"2017"}, {36}, 0.994, 1.005)
  ({"et"}, {"2018"}, {36}, 0.988, 1.01)
  ({"et"}, {"2016"}, {32}, 1.004, 0.997)
  ({"et"}, {"2017"}, {32}, 1.005, 0.996)
  ({"et"}, {"2018"}, {32}, 1.006, 0.992)
  ({"et"}, {"2016"}, {33}, 1.003, 0.998)
  ({"et"}, {"2017"}, {33}, 1.004, 0.997)
  ({"et"}, {"2018"}, {33}, 1.007, 0.994)
  ({"mt"}, {"2016"}, {35}, 0.992, 1.011)
  ({"mt"}, {"2017"}, {35}, 0.993, 1.005)
  ({"mt"}, {"2018"}, {35}, 0.991, 1.011)
  ({"mt"}, {"2016"}, {36}, 0.987, 1.006)
  ({"mt"}, {"2017"}, {36}, 0.995, 1.005)
  ({"mt"}, {"2018"}, {36}, 0.99, 1.009)
  ({"mt"}, {"2016"}, {32}, 1.003, 0.996)
  ({"mt"}, {"2017"}, {32}, 1.005, 0.997)
  ({"mt"}, {"2018"}, {32}, 1.006, 0.993)
  ({"mt"}, {"2016"}, {33}, 1.004, 0.998)
  ({"mt"}, {"2017"}, {33}, 1.003, 0.997)
  ({"mt"}, {"2018"}, {33}, 1.006, 0.995)
  ({"tt"}, {"2016"}, {32}, 1.001, 1.01)
  ({"tt"}, {"2017"}, {32}, 1.0, 0.998)
  ({"tt"}, {"2018"}, {32}, 0.992, 1.009)
  ({"tt"}, {"2016"}, {35}, 1.0, 0.996)
  ({"tt"}, {"2017"}, {35}, 1.0, 1.001)
  ({"tt"}, {"2018"}, {35}, 1.005, 0.994)
  );
  
  cb.cp().process(mssm_bbH_signals).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"60"}, 0.981, 1.0)
  ({"em"}, {"2017"}, btag_catagories, {"60"}, 1.0, 1.0)
  ({"em"}, {"2018"}, btag_catagories, {"60"}, 0.976, 1.0)
  ({"et"}, {"2016"}, btag_catagories, {"60"}, 1.0, 1.0)
  ({"et"}, {"2017"}, btag_catagories, {"60"}, 1.0, 1.0)
  ({"et"}, {"2018"}, btag_catagories, {"60"}, 0.931, 1.0)
  ({"mt"}, {"2016"}, btag_catagories, {"60"}, 0.988, 1.001)
  ({"mt"}, {"2017"}, btag_catagories, {"60"}, 0.983, 1.0)
  ({"mt"}, {"2018"}, btag_catagories, {"60"}, 1.0, 1.0)
  ({"tt"}, {"2016"}, btag_catagories, {"60"}, 1.0, 1.0)
  ({"tt"}, {"2017"}, btag_catagories, {"60"}, 1.0, 1.0)
  ({"tt"}, {"2018"}, btag_catagories, {"60"}, 1.0, 1.0)
  ({"em"}, {"2016"}, nobtag_catagories, {"60"}, 1.011, 1.0)
  ({"em"}, {"2017"}, nobtag_catagories, {"60"}, 1.0, 1.0)
  ({"em"}, {"2018"}, nobtag_catagories, {"60"}, 1.017, 1.0)
  ({"et"}, {"2016"}, nobtag_catagories, {"60"}, 1.0, 1.0)
  ({"et"}, {"2017"}, nobtag_catagories, {"60"}, 1.0, 1.0)
  ({"et"}, {"2018"}, nobtag_catagories, {"60"}, 1.001, 1.0)
  ({"mt"}, {"2016"}, nobtag_catagories, {"60"}, 1.023, 0.998)
  ({"mt"}, {"2017"}, nobtag_catagories, {"60"}, 1.033, 1.0)
  ({"mt"}, {"2018"}, nobtag_catagories, {"60"}, 1.0, 1.0)
  ({"tt"}, {"2016"}, nobtag_catagories, {"60"}, 1.0, 1.0)
  ({"tt"}, {"2017"}, nobtag_catagories, {"60"}, 1.0, 1.0)
  ({"tt"}, {"2018"}, nobtag_catagories, {"60"}, 1.0, 1.0)
  );
  
  cb.cp().process(mssm_bbH_signals).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"80"}, 1.0, 1.006)
  ({"em"}, {"2017"}, btag_catagories, {"80"}, 0.997, 1.01)
  ({"em"}, {"2018"}, btag_catagories, {"80"}, 0.997, 1.011)
  ({"et"}, {"2016"}, btag_catagories, {"80"}, 0.977, 1.003)
  ({"et"}, {"2017"}, btag_catagories, {"80"}, 1.0, 0.993)
  ({"et"}, {"2018"}, btag_catagories, {"80"}, 1.002, 1.01)
  ({"mt"}, {"2016"}, btag_catagories, {"80"}, 0.996, 1.01)
  ({"mt"}, {"2017"}, btag_catagories, {"80"}, 0.996, 1.0)
  ({"mt"}, {"2018"}, btag_catagories, {"80"}, 0.993, 1.012)
  ({"tt"}, {"2016"}, btag_catagories, {"80"}, 1.0, 1.0)
  ({"tt"}, {"2017"}, btag_catagories, {"80"}, 1.0, 1.0)
  ({"tt"}, {"2018"}, btag_catagories, {"80"}, 1.0, 1.0)
  ({"em"}, {"2016"}, nobtag_catagories, {"80"}, 1.0, 0.999)
  ({"em"}, {"2017"}, nobtag_catagories, {"80"}, 1.001, 0.997)
  ({"em"}, {"2018"}, nobtag_catagories, {"80"}, 1.001, 0.997)
  ({"et"}, {"2016"}, nobtag_catagories, {"80"}, 1.013, 0.998)
  ({"et"}, {"2017"}, nobtag_catagories, {"80"}, 1.0, 1.005)
  ({"et"}, {"2018"}, nobtag_catagories, {"80"}, 0.999, 0.993)
  ({"mt"}, {"2016"}, nobtag_catagories, {"80"}, 1.001, 0.997)
  ({"mt"}, {"2017"}, nobtag_catagories, {"80"}, 1.002, 1.0)
  ({"mt"}, {"2018"}, nobtag_catagories, {"80"}, 1.003, 0.994)
  ({"tt"}, {"2016"}, nobtag_catagories, {"80"}, 1.0, 1.0)
  ({"tt"}, {"2017"}, nobtag_catagories, {"80"}, 1.0, 1.0)
  ({"tt"}, {"2018"}, nobtag_catagories, {"80"}, 1.0, 1.0)
  );
  
  cb.cp().process(mssm_bbH_signals).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"100"}, 0.999, 1.005)
  ({"em"}, {"2017"}, btag_catagories, {"100"}, 0.997, 1.008)
  ({"em"}, {"2018"}, btag_catagories, {"100"}, 0.984, 1.008)
  ({"et"}, {"2016"}, btag_catagories, {"100"}, 0.997, 1.007)
  ({"et"}, {"2017"}, btag_catagories, {"100"}, 0.988, 1.004)
  ({"et"}, {"2018"}, btag_catagories, {"100"}, 0.998, 1.012)
  ({"mt"}, {"2016"}, btag_catagories, {"100"}, 0.999, 1.003)
  ({"mt"}, {"2017"}, btag_catagories, {"100"}, 0.996, 1.002)
  ({"mt"}, {"2018"}, btag_catagories, {"100"}, 0.989, 1.01)
  ({"tt"}, {"2016"}, btag_catagories, {"100"}, 0.993, 1.005)
  ({"tt"}, {"2017"}, btag_catagories, {"100"}, 1.0, 0.999)
  ({"tt"}, {"2018"}, btag_catagories, {"100"}, 0.989, 1.005)
  ({"em"}, {"2016"}, nobtag_catagories, {"100"}, 1.0, 0.999)
  ({"em"}, {"2017"}, nobtag_catagories, {"100"}, 1.001, 0.997)
  ({"em"}, {"2018"}, nobtag_catagories, {"100"}, 1.005, 0.998)
  ({"et"}, {"2016"}, nobtag_catagories, {"100"}, 1.001, 0.998)
  ({"et"}, {"2017"}, nobtag_catagories, {"100"}, 1.004, 1.0)
  ({"et"}, {"2018"}, nobtag_catagories, {"100"}, 1.0, 0.995)
  ({"mt"}, {"2016"}, nobtag_catagories, {"100"}, 1.001, 0.999)
  ({"mt"}, {"2017"}, nobtag_catagories, {"100"}, 1.001, 0.999)
  ({"mt"}, {"2018"}, nobtag_catagories, {"100"}, 1.003, 0.997)
  ({"tt"}, {"2016"}, nobtag_catagories, {"100"}, 1.005, 0.996)
  ({"tt"}, {"2017"}, nobtag_catagories, {"100"}, 1.0, 1.001)
  ({"tt"}, {"2018"}, nobtag_catagories, {"100"}, 1.01, 0.995)
  );
  
  cb.cp().process(mssm_bbH_signals).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"120"}, 0.997, 1.005)
  ({"em"}, {"2017"}, btag_catagories, {"120"}, 0.988, 1.013)
  ({"em"}, {"2018"}, btag_catagories, {"120"}, 0.988, 1.011)
  ({"et"}, {"2016"}, btag_catagories, {"120"}, 0.995, 1.009)
  ({"et"}, {"2017"}, btag_catagories, {"120"}, 0.99, 1.01)
  ({"et"}, {"2018"}, btag_catagories, {"120"}, 0.991, 1.01)
  ({"mt"}, {"2016"}, btag_catagories, {"120"}, 0.996, 1.004)
  ({"mt"}, {"2017"}, btag_catagories, {"120"}, 0.996, 1.006)
  ({"mt"}, {"2018"}, btag_catagories, {"120"}, 0.991, 1.004)
  ({"tt"}, {"2016"}, btag_catagories, {"120"}, 0.998, 1.004)
  ({"tt"}, {"2017"}, btag_catagories, {"120"}, 0.998, 1.003)
  ({"tt"}, {"2018"}, btag_catagories, {"120"}, 0.993, 1.015)
  ({"em"}, {"2016"}, nobtag_catagories, {"120"}, 1.001, 0.999)
  ({"em"}, {"2017"}, nobtag_catagories, {"120"}, 1.004, 0.996)
  ({"em"}, {"2018"}, nobtag_catagories, {"120"}, 1.004, 0.996)
  ({"et"}, {"2016"}, nobtag_catagories, {"120"}, 1.001, 0.998)
  ({"et"}, {"2017"}, nobtag_catagories, {"120"}, 1.003, 0.997)
  ({"et"}, {"2018"}, nobtag_catagories, {"120"}, 1.004, 0.997)
  ({"mt"}, {"2016"}, nobtag_catagories, {"120"}, 1.001, 0.999)
  ({"mt"}, {"2017"}, nobtag_catagories, {"120"}, 1.002, 0.998)
  ({"mt"}, {"2018"}, nobtag_catagories, {"120"}, 1.003, 0.998)
  ({"tt"}, {"2016"}, nobtag_catagories, {"120"}, 1.001, 0.998)
  ({"tt"}, {"2017"}, nobtag_catagories, {"120"}, 1.001, 0.999)
  ({"tt"}, {"2018"}, nobtag_catagories, {"120"}, 1.003, 0.994)
  );
  
  cb.cp().process(mssm_bbH_signals).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"125"}, 1.0, 1.011)
  ({"em"}, {"2017"}, btag_catagories, {"125"}, 0.991, 1.007)
  ({"em"}, {"2018"}, btag_catagories, {"125"}, 0.991, 1.011)
  ({"et"}, {"2016"}, btag_catagories, {"125"}, 0.999, 1.009)
  ({"et"}, {"2017"}, btag_catagories, {"125"}, 0.996, 1.006)
  ({"et"}, {"2018"}, btag_catagories, {"125"}, 0.992, 1.004)
  ({"mt"}, {"2016"}, btag_catagories, {"125"}, 0.994, 1.006)
  ({"mt"}, {"2017"}, btag_catagories, {"125"}, 0.996, 1.007)
  ({"mt"}, {"2018"}, btag_catagories, {"125"}, 0.996, 1.012)
  ({"tt"}, {"2016"}, btag_catagories, {"125"}, 0.995, 1.006)
  ({"tt"}, {"2017"}, btag_catagories, {"125"}, 0.996, 1.001)
  ({"tt"}, {"2018"}, btag_catagories, {"125"}, 0.99, 1.009)
  ({"em"}, {"2016"}, nobtag_catagories, {"125"}, 1.0, 0.998)
  ({"em"}, {"2017"}, nobtag_catagories, {"125"}, 1.003, 0.998)
  ({"em"}, {"2018"}, nobtag_catagories, {"125"}, 1.003, 0.996)
  ({"et"}, {"2016"}, nobtag_catagories, {"125"}, 1.001, 0.998)
  ({"et"}, {"2017"}, nobtag_catagories, {"125"}, 1.002, 0.998)
  ({"et"}, {"2018"}, nobtag_catagories, {"125"}, 1.004, 0.998)
  ({"mt"}, {"2016"}, nobtag_catagories, {"125"}, 1.001, 0.999)
  ({"mt"}, {"2017"}, nobtag_catagories, {"125"}, 1.001, 0.998)
  ({"mt"}, {"2018"}, nobtag_catagories, {"125"}, 1.003, 0.995)
  ({"tt"}, {"2016"}, nobtag_catagories, {"125"}, 1.002, 0.998)
  ({"tt"}, {"2017"}, nobtag_catagories, {"125"}, 1.002, 1.0)
  ({"tt"}, {"2018"}, nobtag_catagories, {"125"}, 1.004, 0.996)
  );
  
  cb.cp().process(mssm_bbH_signals).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"130"}, 0.996, 1.01)
  ({"em"}, {"2017"}, btag_catagories, {"130"}, 0.996, 1.001)
  ({"em"}, {"2018"}, btag_catagories, {"130"}, 0.993, 1.006)
  ({"et"}, {"2016"}, btag_catagories, {"130"}, 0.994, 1.004)
  ({"et"}, {"2017"}, btag_catagories, {"130"}, 0.995, 1.005)
  ({"et"}, {"2018"}, btag_catagories, {"130"}, 0.989, 1.008)
  ({"mt"}, {"2016"}, btag_catagories, {"130"}, 0.999, 1.004)
  ({"mt"}, {"2017"}, btag_catagories, {"130"}, 0.997, 1.006)
  ({"mt"}, {"2018"}, btag_catagories, {"130"}, 0.991, 1.008)
  ({"tt"}, {"2016"}, btag_catagories, {"130"}, 0.994, 1.003)
  ({"tt"}, {"2017"}, btag_catagories, {"130"}, 0.997, 1.003)
  ({"tt"}, {"2018"}, btag_catagories, {"130"}, 0.99, 1.001)
  ({"em"}, {"2016"}, nobtag_catagories, {"130"}, 1.001, 0.997)
  ({"em"}, {"2017"}, nobtag_catagories, {"130"}, 1.001, 1.0)
  ({"em"}, {"2018"}, nobtag_catagories, {"130"}, 1.002, 0.998)
  ({"et"}, {"2016"}, nobtag_catagories, {"130"}, 1.002, 0.999)
  ({"et"}, {"2017"}, nobtag_catagories, {"130"}, 1.002, 0.999)
  ({"et"}, {"2018"}, nobtag_catagories, {"130"}, 1.004, 0.997)
  ({"mt"}, {"2016"}, nobtag_catagories, {"130"}, 1.0, 0.999)
  ({"mt"}, {"2017"}, nobtag_catagories, {"130"}, 1.001, 0.998)
  ({"mt"}, {"2018"}, nobtag_catagories, {"130"}, 1.002, 0.996)
  ({"tt"}, {"2016"}, nobtag_catagories, {"130"}, 1.002, 0.999)
  ({"tt"}, {"2017"}, nobtag_catagories, {"130"}, 1.001, 0.999)
  ({"tt"}, {"2018"}, nobtag_catagories, {"130"}, 1.005, 1.0)
  );
  
  cb.cp().process(mssm_bbH_signals).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"140"}, 0.994, 1.002)
  ({"em"}, {"2017"}, btag_catagories, {"140"}, 0.995, 1.004)
  ({"em"}, {"2018"}, btag_catagories, {"140"}, 0.987, 1.007)
  ({"et"}, {"2016"}, btag_catagories, {"140"}, 0.999, 1.004)
  ({"et"}, {"2017"}, btag_catagories, {"140"}, 0.992, 1.005)
  ({"et"}, {"2018"}, btag_catagories, {"140"}, 0.999, 1.009)
  ({"mt"}, {"2016"}, btag_catagories, {"140"}, 0.999, 1.005)
  ({"mt"}, {"2017"}, btag_catagories, {"140"}, 0.994, 1.006)
  ({"mt"}, {"2018"}, btag_catagories, {"140"}, 0.989, 1.005)
  ({"tt"}, {"2016"}, btag_catagories, {"140"}, 0.994, 0.999)
  ({"tt"}, {"2017"}, btag_catagories, {"140"}, 0.996, 1.009)
  ({"tt"}, {"2018"}, btag_catagories, {"140"}, 0.996, 1.006)
  ({"em"}, {"2016"}, nobtag_catagories, {"140"}, 1.002, 1.0)
  ({"em"}, {"2017"}, nobtag_catagories, {"140"}, 1.002, 0.999)
  ({"em"}, {"2018"}, nobtag_catagories, {"140"}, 1.005, 0.997)
  ({"et"}, {"2016"}, nobtag_catagories, {"140"}, 1.0, 0.999)
  ({"et"}, {"2017"}, nobtag_catagories, {"140"}, 1.003, 0.999)
  ({"et"}, {"2018"}, nobtag_catagories, {"140"}, 0.999, 0.996)
  ({"mt"}, {"2016"}, nobtag_catagories, {"140"}, 1.0, 0.998)
  ({"mt"}, {"2017"}, nobtag_catagories, {"140"}, 1.002, 0.998)
  ({"mt"}, {"2018"}, nobtag_catagories, {"140"}, 1.003, 0.998)
  ({"tt"}, {"2016"}, nobtag_catagories, {"140"}, 1.002, 1.0)
  ({"tt"}, {"2017"}, nobtag_catagories, {"140"}, 1.002, 0.996)
  ({"tt"}, {"2018"}, nobtag_catagories, {"140"}, 1.002, 0.997)
  );
  
  cb.cp().process(mssm_bbH_signals).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"160"}, 0.995, 1.004)
  ({"em"}, {"2017"}, btag_catagories, {"160"}, 0.987, 1.008)
  ({"em"}, {"2018"}, btag_catagories, {"160"}, 0.992, 1.014)
  ({"et"}, {"2016"}, btag_catagories, {"160"}, 0.993, 1.006)
  ({"et"}, {"2017"}, btag_catagories, {"160"}, 0.989, 1.01)
  ({"et"}, {"2018"}, btag_catagories, {"160"}, 0.988, 1.006)
  ({"mt"}, {"2016"}, btag_catagories, {"160"}, 0.996, 1.004)
  ({"mt"}, {"2017"}, btag_catagories, {"160"}, 0.994, 1.007)
  ({"mt"}, {"2018"}, btag_catagories, {"160"}, 0.993, 1.008)
  ({"tt"}, {"2016"}, btag_catagories, {"160"}, 0.994, 1.005)
  ({"tt"}, {"2017"}, btag_catagories, {"160"}, 0.996, 1.005)
  ({"tt"}, {"2018"}, btag_catagories, {"160"}, 0.994, 1.008)
  ({"em"}, {"2016"}, nobtag_catagories, {"160"}, 1.001, 0.999)
  ({"em"}, {"2017"}, nobtag_catagories, {"160"}, 1.005, 0.997)
  ({"em"}, {"2018"}, nobtag_catagories, {"160"}, 1.003, 0.994)
  ({"et"}, {"2016"}, nobtag_catagories, {"160"}, 1.002, 0.998)
  ({"et"}, {"2017"}, nobtag_catagories, {"160"}, 1.004, 0.996)
  ({"et"}, {"2018"}, nobtag_catagories, {"160"}, 1.005, 0.997)
  ({"mt"}, {"2016"}, nobtag_catagories, {"160"}, 1.001, 0.999)
  ({"mt"}, {"2017"}, nobtag_catagories, {"160"}, 1.002, 0.998)
  ({"mt"}, {"2018"}, nobtag_catagories, {"160"}, 1.003, 0.996)
  ({"tt"}, {"2016"}, nobtag_catagories, {"160"}, 1.002, 0.998)
  ({"tt"}, {"2017"}, nobtag_catagories, {"160"}, 1.002, 0.998)
  ({"tt"}, {"2018"}, nobtag_catagories, {"160"}, 1.003, 0.996)
  );
  
  cb.cp().process(mssm_bbH_signals).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"180"}, 0.994, 1.004)
  ({"em"}, {"2017"}, btag_catagories, {"180"}, 0.991, 1.008)
  ({"em"}, {"2018"}, btag_catagories, {"180"}, 0.989, 1.006)
  ({"et"}, {"2016"}, btag_catagories, {"180"}, 0.995, 1.003)
  ({"et"}, {"2017"}, btag_catagories, {"180"}, 0.993, 1.005)
  ({"et"}, {"2018"}, btag_catagories, {"180"}, 0.988, 1.006)
  ({"mt"}, {"2016"}, btag_catagories, {"180"}, 0.996, 1.003)
  ({"mt"}, {"2017"}, btag_catagories, {"180"}, 0.993, 1.003)
  ({"mt"}, {"2018"}, btag_catagories, {"180"}, 0.992, 1.008)
  ({"tt"}, {"2016"}, btag_catagories, {"180"}, 0.998, 1.004)
  ({"tt"}, {"2017"}, btag_catagories, {"180"}, 0.998, 1.005)
  ({"tt"}, {"2018"}, btag_catagories, {"180"}, 0.991, 1.006)
  ({"em"}, {"2016"}, nobtag_catagories, {"180"}, 1.002, 0.999)
  ({"em"}, {"2017"}, nobtag_catagories, {"180"}, 1.003, 0.997)
  ({"em"}, {"2018"}, nobtag_catagories, {"180"}, 1.004, 0.998)
  ({"et"}, {"2016"}, nobtag_catagories, {"180"}, 1.001, 0.999)
  ({"et"}, {"2017"}, nobtag_catagories, {"180"}, 1.003, 0.998)
  ({"et"}, {"2018"}, nobtag_catagories, {"180"}, 1.005, 0.997)
  ({"mt"}, {"2016"}, nobtag_catagories, {"180"}, 1.001, 0.999)
  ({"mt"}, {"2017"}, nobtag_catagories, {"180"}, 1.003, 0.999)
  ({"mt"}, {"2018"}, nobtag_catagories, {"180"}, 1.003, 0.996)
  ({"tt"}, {"2016"}, nobtag_catagories, {"180"}, 1.001, 0.999)
  ({"tt"}, {"2017"}, nobtag_catagories, {"180"}, 1.001, 0.998)
  ({"tt"}, {"2018"}, nobtag_catagories, {"180"}, 1.004, 0.997)
  );
  
  cb.cp().process(mssm_bbH_signals).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"200"}, 0.997, 1.009)
  ({"em"}, {"2017"}, btag_catagories, {"200"}, 0.994, 1.008)
  ({"em"}, {"2018"}, btag_catagories, {"200"}, 0.992, 1.01)
  ({"et"}, {"2016"}, btag_catagories, {"200"}, 0.758, 0.765)
  ({"et"}, {"2017"}, btag_catagories, {"200"}, 0.772, 0.777)
  ({"et"}, {"2018"}, btag_catagories, {"200"}, 0.993, 1.008)
  ({"mt"}, {"2016"}, btag_catagories, {"200"}, 0.783, 0.79)
  ({"mt"}, {"2017"}, btag_catagories, {"200"}, 0.791, 0.797)
  ({"mt"}, {"2018"}, btag_catagories, {"200"}, 0.99, 1.006)
  ({"tt"}, {"2016"}, btag_catagories, {"200"}, 0.632, 0.638)
  ({"tt"}, {"2017"}, btag_catagories, {"200"}, 0.707, 0.714)
  ({"tt"}, {"2018"}, btag_catagories, {"200"}, 0.988, 1.008)
  ({"em"}, {"2016"}, nobtag_catagories, {"200"}, 1.001, 0.997)
  ({"em"}, {"2017"}, nobtag_catagories, {"200"}, 1.002, 0.997)
  ({"em"}, {"2018"}, nobtag_catagories, {"200"}, 1.003, 0.996)
  ({"et"}, {"2016"}, nobtag_catagories, {"200"}, 0.81, 0.808)
  ({"et"}, {"2017"}, nobtag_catagories, {"200"}, 0.862, 0.86)
  ({"et"}, {"2018"}, nobtag_catagories, {"200"}, 1.003, 0.996)
  ({"mt"}, {"2016"}, nobtag_catagories, {"200"}, 0.867, 0.865)
  ({"mt"}, {"2017"}, nobtag_catagories, {"200"}, 0.871, 0.868)
  ({"mt"}, {"2018"}, nobtag_catagories, {"200"}, 1.004, 0.997)
  ({"tt"}, {"2016"}, nobtag_catagories, {"200"}, 0.758, 0.756)
  ({"tt"}, {"2017"}, nobtag_catagories, {"200"}, 0.772, 0.769)
  ({"tt"}, {"2018"}, nobtag_catagories, {"200"}, 1.006, 0.996)
  );
  
  cb.cp().process(mssm_bbH_signals).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"250"}, 0.998, 1.005)
  ({"em"}, {"2017"}, btag_catagories, {"250"}, 0.996, 1.006)
  ({"em"}, {"2018"}, btag_catagories, {"250"}, 0.99, 1.015)
  ({"et"}, {"2016"}, btag_catagories, {"250"}, 1.003, 1.005)
  ({"et"}, {"2017"}, btag_catagories, {"250"}, 0.997, 1.007)
  ({"et"}, {"2018"}, btag_catagories, {"250"}, 0.99, 1.012)
  ({"mt"}, {"2016"}, btag_catagories, {"250"}, 0.997, 1.004)
  ({"mt"}, {"2017"}, btag_catagories, {"250"}, 0.996, 1.004)
  ({"mt"}, {"2018"}, btag_catagories, {"250"}, 0.988, 1.011)
  ({"tt"}, {"2016"}, btag_catagories, {"250"}, 0.995, 1.003)
  ({"tt"}, {"2017"}, btag_catagories, {"250"}, 0.996, 1.004)
  ({"tt"}, {"2018"}, btag_catagories, {"250"}, 0.992, 1.008)
  ({"em"}, {"2016"}, nobtag_catagories, {"250"}, 1.001, 0.998)
  ({"em"}, {"2017"}, nobtag_catagories, {"250"}, 1.002, 0.998)
  ({"em"}, {"2018"}, nobtag_catagories, {"250"}, 1.005, 0.993)
  ({"et"}, {"2016"}, nobtag_catagories, {"250"}, 0.999, 0.998)
  ({"et"}, {"2017"}, nobtag_catagories, {"250"}, 1.002, 0.997)
  ({"et"}, {"2018"}, nobtag_catagories, {"250"}, 1.004, 0.993)
  ({"mt"}, {"2016"}, nobtag_catagories, {"250"}, 1.001, 0.999)
  ({"mt"}, {"2017"}, nobtag_catagories, {"250"}, 1.001, 0.998)
  ({"mt"}, {"2018"}, nobtag_catagories, {"250"}, 1.006, 0.994)
  ({"tt"}, {"2016"}, nobtag_catagories, {"250"}, 1.002, 0.999)
  ({"tt"}, {"2017"}, nobtag_catagories, {"250"}, 1.002, 0.998)
  ({"tt"}, {"2018"}, nobtag_catagories, {"250"}, 1.004, 0.996)
  );
  
  cb.cp().process(mssm_bbH_signals).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"300"}, 1.0, 1.0)
  ({"em"}, {"2017"}, btag_catagories, {"300"}, 0.993, 1.006)
  ({"em"}, {"2018"}, btag_catagories, {"300"}, 0.993, 1.006)
  ({"et"}, {"2016"}, btag_catagories, {"300"}, 1.0, 1.0)
  ({"et"}, {"2017"}, btag_catagories, {"300"}, 0.992, 1.005)
  ({"et"}, {"2018"}, btag_catagories, {"300"}, 0.987, 1.008)
  ({"mt"}, {"2016"}, btag_catagories, {"300"}, 1.0, 1.0)
  ({"mt"}, {"2017"}, btag_catagories, {"300"}, 0.997, 1.004)
  ({"mt"}, {"2018"}, btag_catagories, {"300"}, 0.992, 1.008)
  ({"tt"}, {"2016"}, btag_catagories, {"300"}, 1.0, 1.0)
  ({"tt"}, {"2017"}, btag_catagories, {"300"}, 0.995, 1.007)
  ({"tt"}, {"2018"}, btag_catagories, {"300"}, 0.992, 1.008)
  ({"em"}, {"2016"}, nobtag_catagories, {"300"}, 1.0, 1.0)
  ({"em"}, {"2017"}, nobtag_catagories, {"300"}, 1.003, 0.997)
  ({"em"}, {"2018"}, nobtag_catagories, {"300"}, 1.004, 0.997)
  ({"et"}, {"2016"}, nobtag_catagories, {"300"}, 1.0, 1.0)
  ({"et"}, {"2017"}, nobtag_catagories, {"300"}, 1.004, 0.998)
  ({"et"}, {"2018"}, nobtag_catagories, {"300"}, 1.007, 0.996)
  ({"mt"}, {"2016"}, nobtag_catagories, {"300"}, 1.0, 1.0)
  ({"mt"}, {"2017"}, nobtag_catagories, {"300"}, 1.002, 0.998)
  ({"mt"}, {"2018"}, nobtag_catagories, {"300"}, 1.004, 0.996)
  ({"tt"}, {"2016"}, nobtag_catagories, {"300"}, 1.0, 1.0)
  ({"tt"}, {"2017"}, nobtag_catagories, {"300"}, 1.003, 0.996)
  ({"tt"}, {"2018"}, nobtag_catagories, {"300"}, 1.005, 0.995)
  );
  
  cb.cp().process(mssm_bbH_signals).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"350"}, 0.992, 1.006)
  ({"em"}, {"2017"}, btag_catagories, {"350"}, 0.991, 1.002)
  ({"em"}, {"2018"}, btag_catagories, {"350"}, 0.992, 1.011)
  ({"et"}, {"2016"}, btag_catagories, {"350"}, 0.998, 1.003)
  ({"et"}, {"2017"}, btag_catagories, {"350"}, 0.996, 1.006)
  ({"et"}, {"2018"}, btag_catagories, {"350"}, 0.992, 1.006)
  ({"mt"}, {"2016"}, btag_catagories, {"350"}, 0.996, 1.005)
  ({"mt"}, {"2017"}, btag_catagories, {"350"}, 0.996, 1.007)
  ({"mt"}, {"2018"}, btag_catagories, {"350"}, 0.994, 1.007)
  ({"tt"}, {"2016"}, btag_catagories, {"350"}, 0.996, 1.004)
  ({"tt"}, {"2017"}, btag_catagories, {"350"}, 0.996, 1.006)
  ({"tt"}, {"2018"}, btag_catagories, {"350"}, 0.997, 1.007)
  ({"em"}, {"2016"}, nobtag_catagories, {"350"}, 1.003, 0.997)
  ({"em"}, {"2017"}, nobtag_catagories, {"350"}, 1.005, 0.999)
  ({"em"}, {"2018"}, nobtag_catagories, {"350"}, 1.005, 0.994)
  ({"et"}, {"2016"}, nobtag_catagories, {"350"}, 1.001, 0.999)
  ({"et"}, {"2017"}, nobtag_catagories, {"350"}, 1.003, 0.997)
  ({"et"}, {"2018"}, nobtag_catagories, {"350"}, 1.005, 0.996)
  ({"mt"}, {"2016"}, nobtag_catagories, {"350"}, 1.002, 0.998)
  ({"mt"}, {"2017"}, nobtag_catagories, {"350"}, 1.002, 0.996)
  ({"mt"}, {"2018"}, nobtag_catagories, {"350"}, 1.004, 0.996)
  ({"tt"}, {"2016"}, nobtag_catagories, {"350"}, 1.002, 0.998)
  ({"tt"}, {"2017"}, nobtag_catagories, {"350"}, 1.002, 0.997)
  ({"tt"}, {"2018"}, nobtag_catagories, {"350"}, 1.002, 0.995)
  );
  
  cb.cp().process(mssm_bbH_signals).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"400"}, 0.993, 1.007)
  ({"em"}, {"2017"}, btag_catagories, {"400"}, 0.995, 1.008)
  ({"em"}, {"2018"}, btag_catagories, {"400"}, 0.991, 1.007)
  ({"et"}, {"2016"}, btag_catagories, {"400"}, 0.996, 1.002)
  ({"et"}, {"2017"}, btag_catagories, {"400"}, 0.997, 1.006)
  ({"et"}, {"2018"}, btag_catagories, {"400"}, 0.992, 1.009)
  ({"mt"}, {"2016"}, btag_catagories, {"400"}, 0.996, 1.003)
  ({"mt"}, {"2017"}, btag_catagories, {"400"}, 0.996, 1.005)
  ({"mt"}, {"2018"}, btag_catagories, {"400"}, 0.993, 1.006)
  ({"tt"}, {"2016"}, btag_catagories, {"400"}, 0.996, 1.006)
  ({"tt"}, {"2017"}, btag_catagories, {"400"}, 0.995, 1.003)
  ({"tt"}, {"2018"}, btag_catagories, {"400"}, 0.993, 1.006)
  ({"em"}, {"2016"}, nobtag_catagories, {"400"}, 1.003, 0.997)
  ({"em"}, {"2017"}, nobtag_catagories, {"400"}, 1.003, 0.996)
  ({"em"}, {"2018"}, nobtag_catagories, {"400"}, 1.005, 0.996)
  ({"et"}, {"2016"}, nobtag_catagories, {"400"}, 1.002, 0.999)
  ({"et"}, {"2017"}, nobtag_catagories, {"400"}, 1.002, 0.997)
  ({"et"}, {"2018"}, nobtag_catagories, {"400"}, 1.005, 0.994)
  ({"mt"}, {"2016"}, nobtag_catagories, {"400"}, 1.002, 0.999)
  ({"mt"}, {"2017"}, nobtag_catagories, {"400"}, 1.002, 0.997)
  ({"mt"}, {"2018"}, nobtag_catagories, {"400"}, 1.005, 0.996)
  ({"tt"}, {"2016"}, nobtag_catagories, {"400"}, 1.002, 0.997)
  ({"tt"}, {"2017"}, nobtag_catagories, {"400"}, 1.003, 0.998)
  ({"tt"}, {"2018"}, nobtag_catagories, {"400"}, 1.005, 0.995)
  );
  
  cb.cp().process(mssm_bbH_signals).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"450"}, 0.992, 1.003)
  ({"em"}, {"2017"}, btag_catagories, {"450"}, 0.993, 1.002)
  ({"em"}, {"2018"}, btag_catagories, {"450"}, 0.991, 1.006)
  ({"et"}, {"2016"}, btag_catagories, {"450"}, 0.997, 1.004)
  ({"et"}, {"2017"}, btag_catagories, {"450"}, 0.996, 1.005)
  ({"et"}, {"2018"}, btag_catagories, {"450"}, 0.99, 1.008)
  ({"mt"}, {"2016"}, btag_catagories, {"450"}, 0.996, 1.003)
  ({"mt"}, {"2017"}, btag_catagories, {"450"}, 0.995, 1.003)
  ({"mt"}, {"2018"}, btag_catagories, {"450"}, 0.993, 1.008)
  ({"tt"}, {"2016"}, btag_catagories, {"450"}, 0.996, 1.004)
  ({"tt"}, {"2017"}, btag_catagories, {"450"}, 0.996, 1.004)
  ({"tt"}, {"2018"}, btag_catagories, {"450"}, 0.994, 1.008)
  ({"em"}, {"2016"}, nobtag_catagories, {"450"}, 1.004, 0.999)
  ({"em"}, {"2017"}, nobtag_catagories, {"450"}, 1.004, 0.999)
  ({"em"}, {"2018"}, nobtag_catagories, {"450"}, 1.006, 0.996)
  ({"et"}, {"2016"}, nobtag_catagories, {"450"}, 1.001, 0.998)
  ({"et"}, {"2017"}, nobtag_catagories, {"450"}, 1.003, 0.998)
  ({"et"}, {"2018"}, nobtag_catagories, {"450"}, 1.006, 0.994)
  ({"mt"}, {"2016"}, nobtag_catagories, {"450"}, 1.002, 0.998)
  ({"mt"}, {"2017"}, nobtag_catagories, {"450"}, 1.003, 0.998)
  ({"mt"}, {"2018"}, nobtag_catagories, {"450"}, 1.005, 0.995)
  ({"tt"}, {"2016"}, nobtag_catagories, {"450"}, 1.002, 0.998)
  ({"tt"}, {"2017"}, nobtag_catagories, {"450"}, 1.002, 0.997)
  ({"tt"}, {"2018"}, nobtag_catagories, {"450"}, 1.005, 0.994)
  );
  
  cb.cp().process(mssm_bbH_signals).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"500"}, 0.998, 1.003)
  ({"em"}, {"2017"}, btag_catagories, {"500"}, 0.996, 1.005)
  ({"em"}, {"2018"}, btag_catagories, {"500"}, 0.99, 1.008)
  ({"et"}, {"2016"}, btag_catagories, {"500"}, 0.995, 1.003)
  ({"et"}, {"2017"}, btag_catagories, {"500"}, 0.996, 1.003)
  ({"et"}, {"2018"}, btag_catagories, {"500"}, 0.992, 1.008)
  ({"mt"}, {"2016"}, btag_catagories, {"500"}, 0.996, 1.005)
  ({"mt"}, {"2017"}, btag_catagories, {"500"}, 0.995, 1.008)
  ({"mt"}, {"2018"}, btag_catagories, {"500"}, 0.994, 1.007)
  ({"tt"}, {"2016"}, btag_catagories, {"500"}, 0.995, 1.005)
  ({"tt"}, {"2017"}, btag_catagories, {"500"}, 0.996, 1.004)
  ({"tt"}, {"2018"}, btag_catagories, {"500"}, 0.993, 1.007)
  ({"em"}, {"2016"}, nobtag_catagories, {"500"}, 1.001, 0.998)
  ({"em"}, {"2017"}, nobtag_catagories, {"500"}, 1.002, 0.997)
  ({"em"}, {"2018"}, nobtag_catagories, {"500"}, 1.006, 0.995)
  ({"et"}, {"2016"}, nobtag_catagories, {"500"}, 1.003, 0.998)
  ({"et"}, {"2017"}, nobtag_catagories, {"500"}, 1.003, 0.998)
  ({"et"}, {"2018"}, nobtag_catagories, {"500"}, 1.005, 0.993)
  ({"mt"}, {"2016"}, nobtag_catagories, {"500"}, 1.002, 0.997)
  ({"mt"}, {"2017"}, nobtag_catagories, {"500"}, 1.003, 0.996)
  ({"mt"}, {"2018"}, nobtag_catagories, {"500"}, 1.005, 0.995)
  ({"tt"}, {"2016"}, nobtag_catagories, {"500"}, 1.003, 0.997)
  ({"tt"}, {"2017"}, nobtag_catagories, {"500"}, 1.003, 0.997)
  ({"tt"}, {"2018"}, nobtag_catagories, {"500"}, 1.006, 0.994)
  );
  
  cb.cp().process(mssm_bbH_signals).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"600"}, 0.99, 1.007)
  ({"em"}, {"2017"}, btag_catagories, {"600"}, 0.993, 1.007)
  ({"em"}, {"2018"}, btag_catagories, {"600"}, 0.988, 1.007)
  ({"et"}, {"2016"}, btag_catagories, {"600"}, 0.995, 1.0)
  ({"et"}, {"2017"}, btag_catagories, {"600"}, 0.997, 1.0)
  ({"et"}, {"2018"}, btag_catagories, {"600"}, 0.989, 1.009)
  ({"mt"}, {"2016"}, btag_catagories, {"600"}, 0.994, 1.005)
  ({"mt"}, {"2017"}, btag_catagories, {"600"}, 0.995, 1.003)
  ({"mt"}, {"2018"}, btag_catagories, {"600"}, 0.993, 1.008)
  ({"tt"}, {"2016"}, btag_catagories, {"600"}, 0.998, 1.005)
  ({"tt"}, {"2017"}, btag_catagories, {"600"}, 0.996, 1.005)
  ({"tt"}, {"2018"}, btag_catagories, {"600"}, 0.995, 1.009)
  ({"em"}, {"2016"}, nobtag_catagories, {"600"}, 1.005, 0.996)
  ({"em"}, {"2017"}, nobtag_catagories, {"600"}, 1.004, 0.996)
  ({"em"}, {"2018"}, nobtag_catagories, {"600"}, 1.008, 0.995)
  ({"et"}, {"2016"}, nobtag_catagories, {"600"}, 1.003, 1.0)
  ({"et"}, {"2017"}, nobtag_catagories, {"600"}, 1.003, 1.0)
  ({"et"}, {"2018"}, nobtag_catagories, {"600"}, 1.008, 0.993)
  ({"mt"}, {"2016"}, nobtag_catagories, {"600"}, 1.003, 0.997)
  ({"mt"}, {"2017"}, nobtag_catagories, {"600"}, 1.004, 0.998)
  ({"mt"}, {"2018"}, nobtag_catagories, {"600"}, 1.006, 0.994)
  ({"tt"}, {"2016"}, nobtag_catagories, {"600"}, 1.001, 0.997)
  ({"tt"}, {"2017"}, nobtag_catagories, {"600"}, 1.003, 0.996)
  ({"tt"}, {"2018"}, nobtag_catagories, {"600"}, 1.004, 0.993)
  );
  
  cb.cp().process(mssm_bbH_signals).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"700"}, 1.0, 1.0)
  ({"em"}, {"2017"}, btag_catagories, {"700"}, 0.991, 1.003)
  ({"em"}, {"2018"}, btag_catagories, {"700"}, 0.994, 1.008)
  ({"et"}, {"2016"}, btag_catagories, {"700"}, 1.0, 1.0)
  ({"et"}, {"2017"}, btag_catagories, {"700"}, 0.997, 1.001)
  ({"et"}, {"2018"}, btag_catagories, {"700"}, 0.99, 1.01)
  ({"mt"}, {"2016"}, btag_catagories, {"700"}, 1.0, 1.0)
  ({"mt"}, {"2017"}, btag_catagories, {"700"}, 0.995, 1.002)
  ({"mt"}, {"2018"}, btag_catagories, {"700"}, 0.993, 1.005)
  ({"tt"}, {"2016"}, btag_catagories, {"700"}, 1.0, 1.0)
  ({"tt"}, {"2017"}, btag_catagories, {"700"}, 0.996, 1.005)
  ({"tt"}, {"2018"}, btag_catagories, {"700"}, 0.991, 1.007)
  ({"em"}, {"2016"}, nobtag_catagories, {"700"}, 1.0, 1.0)
  ({"em"}, {"2017"}, nobtag_catagories, {"700"}, 1.005, 0.998)
  ({"em"}, {"2018"}, nobtag_catagories, {"700"}, 1.004, 0.994)
  ({"et"}, {"2016"}, nobtag_catagories, {"700"}, 1.0, 1.0)
  ({"et"}, {"2017"}, nobtag_catagories, {"700"}, 1.001, 0.999)
  ({"et"}, {"2018"}, nobtag_catagories, {"700"}, 1.008, 0.992)
  ({"mt"}, {"2016"}, nobtag_catagories, {"700"}, 1.0, 1.0)
  ({"mt"}, {"2017"}, nobtag_catagories, {"700"}, 1.005, 0.999)
  ({"mt"}, {"2018"}, nobtag_catagories, {"700"}, 1.006, 0.996)
  ({"tt"}, {"2016"}, nobtag_catagories, {"700"}, 1.0, 1.0)
  ({"tt"}, {"2017"}, nobtag_catagories, {"700"}, 1.003, 0.996)
  ({"tt"}, {"2018"}, nobtag_catagories, {"700"}, 1.007, 0.995)
  );
  
  cb.cp().process(mssm_bbH_signals).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"800"}, 0.989, 1.004)
  ({"em"}, {"2017"}, btag_catagories, {"800"}, 0.989, 1.006)
  ({"em"}, {"2018"}, btag_catagories, {"800"}, 0.985, 1.01)
  ({"et"}, {"2016"}, btag_catagories, {"800"}, 0.996, 1.002)
  ({"et"}, {"2017"}, btag_catagories, {"800"}, 0.996, 1.006)
  ({"et"}, {"2018"}, btag_catagories, {"800"}, 0.992, 1.005)
  ({"mt"}, {"2016"}, btag_catagories, {"800"}, 0.997, 1.006)
  ({"mt"}, {"2017"}, btag_catagories, {"800"}, 0.995, 1.005)
  ({"mt"}, {"2018"}, btag_catagories, {"800"}, 0.992, 1.009)
  ({"tt"}, {"2016"}, btag_catagories, {"800"}, 0.996, 1.006)
  ({"tt"}, {"2017"}, btag_catagories, {"800"}, 0.996, 1.008)
  ({"tt"}, {"2018"}, btag_catagories, {"800"}, 0.991, 1.007)
  ({"em"}, {"2016"}, nobtag_catagories, {"800"}, 1.007, 0.997)
  ({"em"}, {"2017"}, nobtag_catagories, {"800"}, 1.007, 0.996)
  ({"em"}, {"2018"}, nobtag_catagories, {"800"}, 1.011, 0.993)
  ({"et"}, {"2016"}, nobtag_catagories, {"800"}, 1.003, 0.999)
  ({"et"}, {"2017"}, nobtag_catagories, {"800"}, 1.002, 0.996)
  ({"et"}, {"2018"}, nobtag_catagories, {"800"}, 1.007, 0.996)
  ({"mt"}, {"2016"}, nobtag_catagories, {"800"}, 1.002, 0.996)
  ({"mt"}, {"2017"}, nobtag_catagories, {"800"}, 1.004, 0.996)
  ({"mt"}, {"2018"}, nobtag_catagories, {"800"}, 1.006, 0.992)
  ({"tt"}, {"2016"}, nobtag_catagories, {"800"}, 1.002, 0.996)
  ({"tt"}, {"2017"}, nobtag_catagories, {"800"}, 1.003, 0.994)
  ({"tt"}, {"2018"}, nobtag_catagories, {"800"}, 1.008, 0.994)
  );
  
  cb.cp().process(mssm_bbH_signals).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"900"}, 0.994, 1.006)
  ({"em"}, {"2017"}, btag_catagories, {"900"}, 0.993, 1.009)
  ({"em"}, {"2018"}, btag_catagories, {"900"}, 0.989, 1.006)
  ({"et"}, {"2016"}, btag_catagories, {"900"}, 0.997, 1.007)
  ({"et"}, {"2017"}, btag_catagories, {"900"}, 0.995, 1.005)
  ({"et"}, {"2018"}, btag_catagories, {"900"}, 0.995, 1.005)
  ({"mt"}, {"2016"}, btag_catagories, {"900"}, 0.997, 1.003)
  ({"mt"}, {"2017"}, btag_catagories, {"900"}, 0.996, 1.007)
  ({"mt"}, {"2018"}, btag_catagories, {"900"}, 0.995, 1.007)
  ({"tt"}, {"2016"}, btag_catagories, {"900"}, 0.997, 1.003)
  ({"tt"}, {"2017"}, btag_catagories, {"900"}, 0.996, 1.005)
  ({"tt"}, {"2018"}, btag_catagories, {"900"}, 0.99, 1.008)
  ({"em"}, {"2016"}, nobtag_catagories, {"900"}, 1.003, 0.997)
  ({"em"}, {"2017"}, nobtag_catagories, {"900"}, 1.004, 0.993)
  ({"em"}, {"2018"}, nobtag_catagories, {"900"}, 1.008, 0.996)
  ({"et"}, {"2016"}, nobtag_catagories, {"900"}, 1.002, 0.996)
  ({"et"}, {"2017"}, nobtag_catagories, {"900"}, 1.004, 0.997)
  ({"et"}, {"2018"}, nobtag_catagories, {"900"}, 1.003, 0.995)
  ({"mt"}, {"2016"}, nobtag_catagories, {"900"}, 1.002, 0.998)
  ({"mt"}, {"2017"}, nobtag_catagories, {"900"}, 1.003, 0.995)
  ({"mt"}, {"2018"}, nobtag_catagories, {"900"}, 1.004, 0.994)
  ({"tt"}, {"2016"}, nobtag_catagories, {"900"}, 1.002, 0.998)
  ({"tt"}, {"2017"}, nobtag_catagories, {"900"}, 1.003, 0.996)
  ({"tt"}, {"2018"}, nobtag_catagories, {"900"}, 1.009, 0.993)
  );
  
  cb.cp().process(mssm_bbH_signals).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"1000"}, 1.0, 1.0)
  ({"em"}, {"2017"}, btag_catagories, {"1000"}, 0.994, 1.007)
  ({"em"}, {"2018"}, btag_catagories, {"1000"}, 1.0, 1.0)
  ({"et"}, {"2016"}, btag_catagories, {"1000"}, 1.0, 1.0)
  ({"et"}, {"2017"}, btag_catagories, {"1000"}, 0.996, 1.007)
  ({"et"}, {"2018"}, btag_catagories, {"1000"}, 0.997, 1.008)
  ({"mt"}, {"2016"}, btag_catagories, {"1000"}, 1.0, 1.0)
  ({"mt"}, {"2017"}, btag_catagories, {"1000"}, 0.995, 1.005)
  ({"mt"}, {"2018"}, btag_catagories, {"1000"}, 0.99, 1.007)
  ({"tt"}, {"2016"}, btag_catagories, {"1000"}, 1.0, 1.0)
  ({"tt"}, {"2017"}, btag_catagories, {"1000"}, 0.997, 1.005)
  ({"tt"}, {"2018"}, btag_catagories, {"1000"}, 0.995, 1.008)
  ({"em"}, {"2016"}, nobtag_catagories, {"1000"}, 1.0, 1.0)
  ({"em"}, {"2017"}, nobtag_catagories, {"1000"}, 1.005, 0.996)
  ({"em"}, {"2018"}, nobtag_catagories, {"1000"}, 1.0, 1.0)
  ({"et"}, {"2016"}, nobtag_catagories, {"1000"}, 1.0, 1.0)
  ({"et"}, {"2017"}, nobtag_catagories, {"1000"}, 1.004, 0.995)
  ({"et"}, {"2018"}, nobtag_catagories, {"1000"}, 1.001, 0.995)
  ({"mt"}, {"2016"}, nobtag_catagories, {"1000"}, 1.0, 1.0)
  ({"mt"}, {"2017"}, nobtag_catagories, {"1000"}, 1.004, 0.996)
  ({"mt"}, {"2018"}, nobtag_catagories, {"1000"}, 1.008, 0.995)
  ({"tt"}, {"2016"}, nobtag_catagories, {"1000"}, 1.0, 1.0)
  ({"tt"}, {"2017"}, nobtag_catagories, {"1000"}, 1.003, 0.996)
  ({"tt"}, {"2018"}, nobtag_catagories, {"1000"}, 1.005, 0.993)
  );
  
  cb.cp().process(mssm_bbH_signals).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"1200"}, 0.996, 1.006)
  ({"em"}, {"2017"}, btag_catagories, {"1200"}, 0.991, 1.004)
  ({"em"}, {"2018"}, btag_catagories, {"1200"}, 0.99, 1.007)
  ({"et"}, {"2016"}, btag_catagories, {"1200"}, 0.994, 1.005)
  ({"et"}, {"2017"}, btag_catagories, {"1200"}, 0.996, 1.003)
  ({"et"}, {"2018"}, btag_catagories, {"1200"}, 0.996, 1.007)
  ({"mt"}, {"2016"}, btag_catagories, {"1200"}, 0.996, 1.005)
  ({"mt"}, {"2017"}, btag_catagories, {"1200"}, 0.997, 1.004)
  ({"mt"}, {"2018"}, btag_catagories, {"1200"}, 0.986, 1.01)
  ({"tt"}, {"2016"}, btag_catagories, {"1200"}, 0.997, 1.001)
  ({"tt"}, {"2017"}, btag_catagories, {"1200"}, 0.998, 1.004)
  ({"tt"}, {"2018"}, btag_catagories, {"1200"}, 0.993, 1.008)
  ({"em"}, {"2016"}, nobtag_catagories, {"1200"}, 1.002, 0.996)
  ({"em"}, {"2017"}, nobtag_catagories, {"1200"}, 1.005, 0.996)
  ({"em"}, {"2018"}, nobtag_catagories, {"1200"}, 1.008, 0.993)
  ({"et"}, {"2016"}, nobtag_catagories, {"1200"}, 1.005, 0.996)
  ({"et"}, {"2017"}, nobtag_catagories, {"1200"}, 1.002, 0.997)
  ({"et"}, {"2018"}, nobtag_catagories, {"1200"}, 1.003, 0.995)
  ({"mt"}, {"2016"}, nobtag_catagories, {"1200"}, 1.003, 0.997)
  ({"mt"}, {"2017"}, nobtag_catagories, {"1200"}, 1.002, 0.996)
  ({"mt"}, {"2018"}, nobtag_catagories, {"1200"}, 1.012, 0.991)
  ({"tt"}, {"2016"}, nobtag_catagories, {"1200"}, 1.002, 0.999)
  ({"tt"}, {"2017"}, nobtag_catagories, {"1200"}, 1.002, 0.997)
  ({"tt"}, {"2018"}, nobtag_catagories, {"1200"}, 1.006, 0.992)
  );
  
  cb.cp().process(mssm_bbH_signals).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"1400"}, 0.988, 1.001)
  ({"em"}, {"2017"}, btag_catagories, {"1400"}, 0.993, 1.009)
  ({"em"}, {"2018"}, btag_catagories, {"1400"}, 0.989, 1.006)
  ({"et"}, {"2016"}, btag_catagories, {"1400"}, 0.999, 1.001)
  ({"et"}, {"2017"}, btag_catagories, {"1400"}, 0.998, 1.004)
  ({"et"}, {"2018"}, btag_catagories, {"1400"}, 0.992, 1.009)
  ({"mt"}, {"2016"}, btag_catagories, {"1400"}, 0.995, 1.002)
  ({"mt"}, {"2017"}, btag_catagories, {"1400"}, 0.996, 1.005)
  ({"mt"}, {"2018"}, btag_catagories, {"1400"}, 0.989, 1.006)
  ({"tt"}, {"2016"}, btag_catagories, {"1400"}, 0.998, 1.004)
  ({"tt"}, {"2017"}, btag_catagories, {"1400"}, 0.993, 1.004)
  ({"tt"}, {"2018"}, btag_catagories, {"1400"}, 0.991, 1.011)
  ({"em"}, {"2016"}, nobtag_catagories, {"1400"}, 1.008, 0.999)
  ({"em"}, {"2017"}, nobtag_catagories, {"1400"}, 1.006, 0.994)
  ({"em"}, {"2018"}, nobtag_catagories, {"1400"}, 1.01, 0.994)
  ({"et"}, {"2016"}, nobtag_catagories, {"1400"}, 1.001, 0.999)
  ({"et"}, {"2017"}, nobtag_catagories, {"1400"}, 1.003, 0.996)
  ({"et"}, {"2018"}, nobtag_catagories, {"1400"}, 1.007, 0.992)
  ({"mt"}, {"2016"}, nobtag_catagories, {"1400"}, 1.004, 0.999)
  ({"mt"}, {"2017"}, nobtag_catagories, {"1400"}, 1.004, 0.996)
  ({"mt"}, {"2018"}, nobtag_catagories, {"1400"}, 1.009, 0.995)
  ({"tt"}, {"2016"}, nobtag_catagories, {"1400"}, 1.002, 0.997)
  ({"tt"}, {"2017"}, nobtag_catagories, {"1400"}, 1.006, 0.997)
  ({"tt"}, {"2018"}, nobtag_catagories, {"1400"}, 1.009, 0.989)
  );
  
  cb.cp().process(mssm_bbH_signals).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"1600"}, 0.992, 1.002)
  ({"em"}, {"2017"}, btag_catagories, {"1600"}, 1.0, 1.0)
  ({"em"}, {"2018"}, btag_catagories, {"1600"}, 0.981, 1.018)
  ({"et"}, {"2016"}, btag_catagories, {"1600"}, 0.998, 1.006)
  ({"et"}, {"2017"}, btag_catagories, {"1600"}, 0.994, 1.006)
  ({"et"}, {"2018"}, btag_catagories, {"1600"}, 0.988, 1.008)
  ({"mt"}, {"2016"}, btag_catagories, {"1600"}, 0.997, 1.004)
  ({"mt"}, {"2017"}, btag_catagories, {"1600"}, 0.998, 1.004)
  ({"mt"}, {"2018"}, btag_catagories, {"1600"}, 0.991, 1.011)
  ({"tt"}, {"2016"}, btag_catagories, {"1600"}, 0.996, 1.007)
  ({"tt"}, {"2017"}, btag_catagories, {"1600"}, 0.994, 1.005)
  ({"tt"}, {"2018"}, btag_catagories, {"1600"}, 0.992, 1.01)
  ({"em"}, {"2016"}, nobtag_catagories, {"1600"}, 1.005, 0.999)
  ({"em"}, {"2017"}, nobtag_catagories, {"1600"}, 1.0, 1.0)
  ({"em"}, {"2018"}, nobtag_catagories, {"1600"}, 1.015, 0.984)
  ({"et"}, {"2016"}, nobtag_catagories, {"1600"}, 1.002, 0.996)
  ({"et"}, {"2017"}, nobtag_catagories, {"1600"}, 1.004, 0.997)
  ({"et"}, {"2018"}, nobtag_catagories, {"1600"}, 1.011, 0.992)
  ({"mt"}, {"2016"}, nobtag_catagories, {"1600"}, 1.002, 0.997)
  ({"mt"}, {"2017"}, nobtag_catagories, {"1600"}, 1.002, 0.997)
  ({"mt"}, {"2018"}, nobtag_catagories, {"1600"}, 1.007, 0.991)
  ({"tt"}, {"2016"}, nobtag_catagories, {"1600"}, 1.003, 0.995)
  ({"tt"}, {"2017"}, nobtag_catagories, {"1600"}, 1.005, 0.996)
  ({"tt"}, {"2018"}, nobtag_catagories, {"1600"}, 1.007, 0.99)
  );
  
  cb.cp().process(mssm_bbH_signals).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"1800"}, 0.994, 1.004)
  ({"em"}, {"2017"}, btag_catagories, {"1800"}, 1.0, 1.0)
  ({"em"}, {"2018"}, btag_catagories, {"1800"}, 0.993, 1.003)
  ({"et"}, {"2016"}, btag_catagories, {"1800"}, 0.999, 1.002)
  ({"et"}, {"2017"}, btag_catagories, {"1800"}, 0.996, 1.004)
  ({"et"}, {"2018"}, btag_catagories, {"1800"}, 0.991, 1.011)
  ({"mt"}, {"2016"}, btag_catagories, {"1800"}, 0.997, 1.005)
  ({"mt"}, {"2017"}, btag_catagories, {"1800"}, 0.995, 1.005)
  ({"mt"}, {"2018"}, btag_catagories, {"1800"}, 0.99, 1.008)
  ({"tt"}, {"2016"}, btag_catagories, {"1800"}, 0.996, 1.007)
  ({"tt"}, {"2017"}, btag_catagories, {"1800"}, 0.996, 1.004)
  ({"tt"}, {"2018"}, btag_catagories, {"1800"}, 0.991, 1.007)
  ({"em"}, {"2016"}, nobtag_catagories, {"1800"}, 1.002, 0.996)
  ({"em"}, {"2017"}, nobtag_catagories, {"1800"}, 1.0, 1.0)
  ({"em"}, {"2018"}, nobtag_catagories, {"1800"}, 1.008, 0.999)
  ({"et"}, {"2016"}, nobtag_catagories, {"1800"}, 1.0, 0.998)
  ({"et"}, {"2017"}, nobtag_catagories, {"1800"}, 1.004, 0.996)
  ({"et"}, {"2018"}, nobtag_catagories, {"1800"}, 1.009, 0.99)
  ({"mt"}, {"2016"}, nobtag_catagories, {"1800"}, 1.002, 0.997)
  ({"mt"}, {"2017"}, nobtag_catagories, {"1800"}, 1.003, 0.994)
  ({"mt"}, {"2018"}, nobtag_catagories, {"1800"}, 1.008, 0.993)
  ({"tt"}, {"2016"}, nobtag_catagories, {"1800"}, 1.003, 0.994)
  ({"tt"}, {"2017"}, nobtag_catagories, {"1800"}, 1.004, 0.997)
  ({"tt"}, {"2018"}, nobtag_catagories, {"1800"}, 1.01, 0.992)
  );
  
  cb.cp().process(mssm_bbH_signals).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"2000"}, 0.988, 1.001)
  ({"em"}, {"2017"}, btag_catagories, {"2000"}, 0.989, 1.004)
  ({"em"}, {"2018"}, btag_catagories, {"2000"}, 0.995, 1.004)
  ({"et"}, {"2016"}, btag_catagories, {"2000"}, 0.992, 1.008)
  ({"et"}, {"2017"}, btag_catagories, {"2000"}, 0.991, 1.005)
  ({"et"}, {"2018"}, btag_catagories, {"2000"}, 0.993, 1.011)
  ({"mt"}, {"2016"}, btag_catagories, {"2000"}, 0.996, 1.003)
  ({"mt"}, {"2017"}, btag_catagories, {"2000"}, 0.996, 1.009)
  ({"mt"}, {"2018"}, btag_catagories, {"2000"}, 0.99, 1.01)
  ({"tt"}, {"2016"}, btag_catagories, {"2000"}, 0.993, 1.005)
  ({"tt"}, {"2017"}, btag_catagories, {"2000"}, 0.994, 1.006)
  ({"tt"}, {"2018"}, btag_catagories, {"2000"}, 0.989, 1.009)
  ({"em"}, {"2016"}, nobtag_catagories, {"2000"}, 1.009, 1.0)
  ({"em"}, {"2017"}, nobtag_catagories, {"2000"}, 1.009, 0.999)
  ({"em"}, {"2018"}, nobtag_catagories, {"2000"}, 1.005, 0.997)
  ({"et"}, {"2016"}, nobtag_catagories, {"2000"}, 1.005, 0.994)
  ({"et"}, {"2017"}, nobtag_catagories, {"2000"}, 1.006, 0.997)
  ({"et"}, {"2018"}, nobtag_catagories, {"2000"}, 1.006, 0.99)
  ({"mt"}, {"2016"}, nobtag_catagories, {"2000"}, 1.003, 0.997)
  ({"mt"}, {"2017"}, nobtag_catagories, {"2000"}, 1.004, 0.993)
  ({"mt"}, {"2018"}, nobtag_catagories, {"2000"}, 1.009, 0.99)
  ({"tt"}, {"2016"}, nobtag_catagories, {"2000"}, 1.005, 0.996)
  ({"tt"}, {"2017"}, nobtag_catagories, {"2000"}, 1.005, 0.995)
  ({"tt"}, {"2018"}, nobtag_catagories, {"2000"}, 1.011, 0.991)
  );
  
  cb.cp().process(mssm_bbH_signals).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"2300"}, 1.001, 1.007)
  ({"em"}, {"2017"}, btag_catagories, {"2300"}, 0.994, 1.008)
  ({"em"}, {"2018"}, btag_catagories, {"2300"}, 0.986, 1.011)
  ({"et"}, {"2016"}, btag_catagories, {"2300"}, 0.993, 1.005)
  ({"et"}, {"2017"}, btag_catagories, {"2300"}, 0.992, 1.009)
  ({"et"}, {"2018"}, btag_catagories, {"2300"}, 0.993, 1.008)
  ({"mt"}, {"2016"}, btag_catagories, {"2300"}, 0.993, 1.002)
  ({"mt"}, {"2017"}, btag_catagories, {"2300"}, 0.994, 1.003)
  ({"mt"}, {"2018"}, btag_catagories, {"2300"}, 0.99, 1.003)
  ({"tt"}, {"2016"}, btag_catagories, {"2300"}, 0.997, 1.005)
  ({"tt"}, {"2017"}, btag_catagories, {"2300"}, 0.995, 1.005)
  ({"tt"}, {"2018"}, btag_catagories, {"2300"}, 0.992, 1.005)
  ({"em"}, {"2016"}, nobtag_catagories, {"2300"}, 0.999, 0.995)
  ({"em"}, {"2017"}, nobtag_catagories, {"2300"}, 1.005, 0.993)
  ({"em"}, {"2018"}, nobtag_catagories, {"2300"}, 1.013, 0.991)
  ({"et"}, {"2016"}, nobtag_catagories, {"2300"}, 1.006, 0.997)
  ({"et"}, {"2017"}, nobtag_catagories, {"2300"}, 1.007, 0.993)
  ({"et"}, {"2018"}, nobtag_catagories, {"2300"}, 1.006, 0.993)
  ({"mt"}, {"2016"}, nobtag_catagories, {"2300"}, 1.006, 0.999)
  ({"mt"}, {"2017"}, nobtag_catagories, {"2300"}, 1.005, 0.997)
  ({"mt"}, {"2018"}, nobtag_catagories, {"2300"}, 1.008, 0.995)
  ({"tt"}, {"2016"}, nobtag_catagories, {"2300"}, 1.003, 0.996)
  ({"tt"}, {"2017"}, nobtag_catagories, {"2300"}, 1.005, 0.995)
  ({"tt"}, {"2018"}, nobtag_catagories, {"2300"}, 1.008, 0.994)
  );
  
  cb.cp().process(mssm_bbH_signals).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"2600"}, 0.996, 1.008)
  ({"em"}, {"2017"}, btag_catagories, {"2600"}, 0.986, 0.998)
  ({"em"}, {"2018"}, btag_catagories, {"2600"}, 0.992, 1.006)
  ({"et"}, {"2016"}, btag_catagories, {"2600"}, 0.993, 1.007)
  ({"et"}, {"2017"}, btag_catagories, {"2600"}, 1.0, 1.0)
  ({"et"}, {"2018"}, btag_catagories, {"2600"}, 0.994, 1.012)
  ({"mt"}, {"2016"}, btag_catagories, {"2600"}, 0.994, 1.005)
  ({"mt"}, {"2017"}, btag_catagories, {"2600"}, 1.0, 1.0)
  ({"mt"}, {"2018"}, btag_catagories, {"2600"}, 0.99, 1.009)
  ({"tt"}, {"2016"}, btag_catagories, {"2600"}, 0.996, 1.007)
  ({"tt"}, {"2017"}, btag_catagories, {"2600"}, 1.0, 1.0)
  ({"tt"}, {"2018"}, btag_catagories, {"2600"}, 0.996, 1.007)
  ({"em"}, {"2016"}, nobtag_catagories, {"2600"}, 1.007, 0.997)
  ({"em"}, {"2017"}, nobtag_catagories, {"2600"}, 1.012, 1.001)
  ({"em"}, {"2018"}, nobtag_catagories, {"2600"}, 1.008, 0.996)
  ({"et"}, {"2016"}, nobtag_catagories, {"2600"}, 1.005, 0.995)
  ({"et"}, {"2017"}, nobtag_catagories, {"2600"}, 1.0, 1.0)
  ({"et"}, {"2018"}, nobtag_catagories, {"2600"}, 1.004, 0.99)
  ({"mt"}, {"2016"}, nobtag_catagories, {"2600"}, 1.005, 0.996)
  ({"mt"}, {"2017"}, nobtag_catagories, {"2600"}, 1.0, 1.0)
  ({"mt"}, {"2018"}, nobtag_catagories, {"2600"}, 1.008, 0.992)
  ({"tt"}, {"2016"}, nobtag_catagories, {"2600"}, 1.003, 0.995)
  ({"tt"}, {"2017"}, nobtag_catagories, {"2600"}, 1.0, 1.0)
  ({"tt"}, {"2018"}, nobtag_catagories, {"2600"}, 1.004, 0.993)
  );
  
  cb.cp().process(mssm_bbH_signals).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"2900"}, 0.99, 1.003)
  ({"em"}, {"2017"}, btag_catagories, {"2900"}, 1.0, 1.0)
  ({"em"}, {"2018"}, btag_catagories, {"2900"}, 0.986, 1.002)
  ({"et"}, {"2016"}, btag_catagories, {"2900"}, 0.991, 1.002)
  ({"et"}, {"2017"}, btag_catagories, {"2900"}, 0.997, 1.01)
  ({"et"}, {"2018"}, btag_catagories, {"2900"}, 0.989, 1.01)
  ({"mt"}, {"2016"}, btag_catagories, {"2900"}, 0.995, 1.003)
  ({"mt"}, {"2017"}, btag_catagories, {"2900"}, 0.998, 1.007)
  ({"mt"}, {"2018"}, btag_catagories, {"2900"}, 0.993, 1.014)
  ({"tt"}, {"2016"}, btag_catagories, {"2900"}, 0.997, 1.004)
  ({"tt"}, {"2017"}, btag_catagories, {"2900"}, 0.994, 1.004)
  ({"tt"}, {"2018"}, btag_catagories, {"2900"}, 0.99, 1.011)
  ({"em"}, {"2016"}, nobtag_catagories, {"2900"}, 1.005, 0.996)
  ({"em"}, {"2017"}, nobtag_catagories, {"2900"}, 1.0, 1.0)
  ({"em"}, {"2018"}, nobtag_catagories, {"2900"}, 1.015, 1.001)
  ({"et"}, {"2016"}, nobtag_catagories, {"2900"}, 1.006, 0.999)
  ({"et"}, {"2017"}, nobtag_catagories, {"2900"}, 1.006, 0.993)
  ({"et"}, {"2018"}, nobtag_catagories, {"2900"}, 1.01, 0.993)
  ({"mt"}, {"2016"}, nobtag_catagories, {"2900"}, 1.004, 0.997)
  ({"mt"}, {"2017"}, nobtag_catagories, {"2900"}, 1.005, 0.996)
  ({"mt"}, {"2018"}, nobtag_catagories, {"2900"}, 1.007, 0.988)
  ({"tt"}, {"2016"}, nobtag_catagories, {"2900"}, 1.002, 0.997)
  ({"tt"}, {"2017"}, nobtag_catagories, {"2900"}, 1.005, 0.997)
  ({"tt"}, {"2018"}, nobtag_catagories, {"2900"}, 1.01, 0.989)
  );
  
  cb.cp().process(mssm_bbH_signals).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"3200"}, 0.995, 1.009)
  ({"em"}, {"2017"}, btag_catagories, {"3200"}, 0.989, 1.002)
  ({"em"}, {"2018"}, btag_catagories, {"3200"}, 1.0, 1.0)
  ({"et"}, {"2016"}, btag_catagories, {"3200"}, 0.997, 1.005)
  ({"et"}, {"2017"}, btag_catagories, {"3200"}, 0.996, 0.998)
  ({"et"}, {"2018"}, btag_catagories, {"3200"}, 0.986, 1.011)
  ({"mt"}, {"2016"}, btag_catagories, {"3200"}, 0.998, 1.006)
  ({"mt"}, {"2017"}, btag_catagories, {"3200"}, 0.99, 1.003)
  ({"mt"}, {"2018"}, btag_catagories, {"3200"}, 0.991, 1.01)
  ({"tt"}, {"2016"}, btag_catagories, {"3200"}, 0.996, 1.007)
  ({"tt"}, {"2017"}, btag_catagories, {"3200"}, 0.991, 1.009)
  ({"tt"}, {"2018"}, btag_catagories, {"3200"}, 0.996, 1.008)
  ({"em"}, {"2016"}, nobtag_catagories, {"3200"}, 1.006, 0.995)
  ({"em"}, {"2017"}, nobtag_catagories, {"3200"}, 1.007, 0.998)
  ({"em"}, {"2018"}, nobtag_catagories, {"3200"}, 1.0, 1.0)
  ({"et"}, {"2016"}, nobtag_catagories, {"3200"}, 1.002, 0.999)
  ({"et"}, {"2017"}, nobtag_catagories, {"3200"}, 1.003, 0.999)
  ({"et"}, {"2018"}, nobtag_catagories, {"3200"}, 1.014, 0.988)
  ({"mt"}, {"2016"}, nobtag_catagories, {"3200"}, 1.003, 0.998)
  ({"mt"}, {"2017"}, nobtag_catagories, {"3200"}, 1.008, 0.998)
  ({"mt"}, {"2018"}, nobtag_catagories, {"3200"}, 1.008, 0.993)
  ({"tt"}, {"2016"}, nobtag_catagories, {"3200"}, 1.003, 0.995)
  ({"tt"}, {"2017"}, nobtag_catagories, {"3200"}, 1.008, 0.993)
  ({"tt"}, {"2018"}, nobtag_catagories, {"3200"}, 1.004, 0.992)
  );
  
  cb.cp().process(mssm_bbH_signals).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"3500"}, 1.001, 1.001)
  ({"em"}, {"2017"}, btag_catagories, {"3500"}, 0.997, 1.009)
  ({"em"}, {"2018"}, btag_catagories, {"3500"}, 0.983, 1.005)
  ({"et"}, {"2016"}, btag_catagories, {"3500"}, 0.992, 1.005)
  ({"et"}, {"2017"}, btag_catagories, {"3500"}, 0.995, 1.012)
  ({"et"}, {"2018"}, btag_catagories, {"3500"}, 0.985, 1.008)
  ({"mt"}, {"2016"}, btag_catagories, {"3500"}, 0.993, 1.004)
  ({"mt"}, {"2017"}, btag_catagories, {"3500"}, 0.99, 1.009)
  ({"mt"}, {"2018"}, btag_catagories, {"3500"}, 0.984, 1.007)
  ({"tt"}, {"2016"}, btag_catagories, {"3500"}, 0.997, 1.008)
  ({"tt"}, {"2017"}, btag_catagories, {"3500"}, 0.994, 1.003)
  ({"tt"}, {"2018"}, btag_catagories, {"3500"}, 0.991, 1.008)
  ({"em"}, {"2016"}, nobtag_catagories, {"3500"}, 1.0, 0.998)
  ({"em"}, {"2017"}, nobtag_catagories, {"3500"}, 1.007, 0.998)
  ({"em"}, {"2018"}, nobtag_catagories, {"3500"}, 1.015, 0.991)
  ({"et"}, {"2016"}, nobtag_catagories, {"3500"}, 1.007, 0.996)
  ({"et"}, {"2017"}, nobtag_catagories, {"3500"}, 1.002, 0.988)
  ({"et"}, {"2018"}, nobtag_catagories, {"3500"}, 1.015, 0.995)
  ({"mt"}, {"2016"}, nobtag_catagories, {"3500"}, 1.005, 0.997)
  ({"mt"}, {"2017"}, nobtag_catagories, {"3500"}, 1.007, 0.995)
  ({"mt"}, {"2018"}, nobtag_catagories, {"3500"}, 1.013, 0.995)
  ({"tt"}, {"2016"}, nobtag_catagories, {"3500"}, 1.002, 0.994)
  ({"tt"}, {"2017"}, nobtag_catagories, {"3500"}, 1.005, 0.998)
  ({"tt"}, {"2018"}, nobtag_catagories, {"3500"}, 1.009, 0.992)
  );
  
  cb.cp().process(mssm_ggH_signals).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"60"}, 1.0, 1.0)
  ({"em"}, {"2017"}, btag_catagories, {"60"}, 1.0, 1.0)
  ({"em"}, {"2018"}, btag_catagories, {"60"}, 1.0, 1.0)
  ({"et"}, {"2016"}, btag_catagories, {"60"}, 0.941, 1.071)
  ({"et"}, {"2017"}, btag_catagories, {"60"}, 0.899, 1.131)
  ({"et"}, {"2018"}, btag_catagories, {"60"}, 0.81, 1.149)
  ({"mt"}, {"2016"}, btag_catagories, {"60"}, 0.95, 1.087)
  ({"mt"}, {"2017"}, btag_catagories, {"60"}, 0.982, 1.056)
  ({"mt"}, {"2018"}, btag_catagories, {"60"}, 0.883, 1.074)
  ({"tt"}, {"2016"}, btag_catagories, {"60"}, 1.0, 1.0)
  ({"tt"}, {"2017"}, btag_catagories, {"60"}, 1.0, 1.0)
  ({"tt"}, {"2018"}, btag_catagories, {"60"}, 0.972, 1.034)
  ({"em"}, {"2016"}, nobtag_catagories, {"60"}, 1.0, 1.0)
  ({"em"}, {"2017"}, nobtag_catagories, {"60"}, 1.0, 1.0)
  ({"em"}, {"2018"}, nobtag_catagories, {"60"}, 1.0, 1.0)
  ({"et"}, {"2016"}, nobtag_catagories, {"60"}, 1.006, 0.993)
  ({"et"}, {"2017"}, nobtag_catagories, {"60"}, 1.015, 0.981)
  ({"et"}, {"2018"}, nobtag_catagories, {"60"}, 1.015, 0.988)
  ({"mt"}, {"2016"}, nobtag_catagories, {"60"}, 1.003, 0.994)
  ({"mt"}, {"2017"}, nobtag_catagories, {"60"}, 1.004, 0.996)
  ({"mt"}, {"2018"}, nobtag_catagories, {"60"}, 1.012, 0.992)
  ({"tt"}, {"2016"}, nobtag_catagories, {"60"}, 1.0, 1.0)
  ({"tt"}, {"2017"}, nobtag_catagories, {"60"}, 1.0, 1.0)
  ({"tt"}, {"2018"}, nobtag_catagories, {"60"}, 1.005, 0.994)
  );
  
  cb.cp().process(mssm_ggH_signals).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"80"}, 0.879, 1.081)
  ({"em"}, {"2017"}, btag_catagories, {"80"}, 0.88, 1.091)
  ({"em"}, {"2018"}, btag_catagories, {"80"}, 0.917, 1.142)
  ({"et"}, {"2016"}, btag_catagories, {"80"}, 0.939, 1.032)
  ({"et"}, {"2017"}, btag_catagories, {"80"}, 0.916, 1.102)
  ({"et"}, {"2018"}, btag_catagories, {"80"}, 1.004, 1.056)
  ({"mt"}, {"2016"}, btag_catagories, {"80"}, 0.922, 1.073)
  ({"mt"}, {"2017"}, btag_catagories, {"80"}, 0.926, 1.061)
  ({"mt"}, {"2018"}, btag_catagories, {"80"}, 0.938, 1.094)
  ({"tt"}, {"2016"}, btag_catagories, {"80"}, 0.944, 1.054)
  ({"tt"}, {"2017"}, btag_catagories, {"80"}, 0.956, 1.012)
  ({"tt"}, {"2018"}, btag_catagories, {"80"}, 0.877, 1.195)
  ({"em"}, {"2016"}, nobtag_catagories, {"80"}, 1.005, 0.997)
  ({"em"}, {"2017"}, nobtag_catagories, {"80"}, 1.006, 0.996)
  ({"em"}, {"2018"}, nobtag_catagories, {"80"}, 1.004, 0.993)
  ({"et"}, {"2016"}, nobtag_catagories, {"80"}, 1.004, 0.998)
  ({"et"}, {"2017"}, nobtag_catagories, {"80"}, 1.006, 0.993)
  ({"et"}, {"2018"}, nobtag_catagories, {"80"}, 1.0, 0.995)
  ({"mt"}, {"2016"}, nobtag_catagories, {"80"}, 1.004, 0.997)
  ({"mt"}, {"2017"}, nobtag_catagories, {"80"}, 1.004, 0.997)
  ({"mt"}, {"2018"}, nobtag_catagories, {"80"}, 1.003, 0.994)
  ({"tt"}, {"2016"}, nobtag_catagories, {"80"}, 1.008, 0.993)
  ({"tt"}, {"2017"}, nobtag_catagories, {"80"}, 1.007, 0.998)
  ({"tt"}, {"2018"}, nobtag_catagories, {"80"}, 1.014, 0.978)
  );
  
  cb.cp().process(mssm_ggH_signals).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"100"}, 0.913, 1.132)
  ({"em"}, {"2017"}, btag_catagories, {"100"}, 0.901, 1.121)
  ({"em"}, {"2018"}, btag_catagories, {"100"}, 0.869, 1.172)
  ({"et"}, {"2016"}, btag_catagories, {"100"}, 0.959, 1.081)
  ({"et"}, {"2017"}, btag_catagories, {"100"}, 0.925, 1.072)
  ({"et"}, {"2018"}, btag_catagories, {"100"}, 0.876, 1.053)
  ({"mt"}, {"2016"}, btag_catagories, {"100"}, 0.969, 1.102)
  ({"mt"}, {"2017"}, btag_catagories, {"100"}, 0.956, 1.033)
  ({"mt"}, {"2018"}, btag_catagories, {"100"}, 0.891, 1.114)
  ({"tt"}, {"2016"}, btag_catagories, {"100"}, 0.937, 1.063)
  ({"tt"}, {"2017"}, btag_catagories, {"100"}, 0.942, 1.089)
  ({"tt"}, {"2018"}, btag_catagories, {"100"}, 0.915, 1.107)
  ({"em"}, {"2016"}, nobtag_catagories, {"100"}, 1.002, 0.996)
  ({"em"}, {"2017"}, nobtag_catagories, {"100"}, 1.004, 0.996)
  ({"em"}, {"2018"}, nobtag_catagories, {"100"}, 1.005, 0.993)
  ({"et"}, {"2016"}, nobtag_catagories, {"100"}, 1.002, 0.997)
  ({"et"}, {"2017"}, nobtag_catagories, {"100"}, 1.003, 0.997)
  ({"et"}, {"2018"}, nobtag_catagories, {"100"}, 1.007, 0.997)
  ({"mt"}, {"2016"}, nobtag_catagories, {"100"}, 1.001, 0.996)
  ({"mt"}, {"2017"}, nobtag_catagories, {"100"}, 1.002, 0.999)
  ({"mt"}, {"2018"}, nobtag_catagories, {"100"}, 1.006, 0.994)
  ({"tt"}, {"2016"}, nobtag_catagories, {"100"}, 1.005, 0.995)
  ({"tt"}, {"2017"}, nobtag_catagories, {"100"}, 1.004, 0.993)
  ({"tt"}, {"2018"}, nobtag_catagories, {"100"}, 1.008, 0.99)
  );
  
  cb.cp().process(mssm_ggH_signals).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"120"}, 0.883, 1.093)
  ({"em"}, {"2017"}, btag_catagories, {"120"}, 0.906, 1.116)
  ({"em"}, {"2018"}, btag_catagories, {"120"}, 0.829, 1.105)
  ({"et"}, {"2016"}, btag_catagories, {"120"}, 0.907, 1.099)
  ({"et"}, {"2017"}, btag_catagories, {"120"}, 0.929, 1.07)
  ({"et"}, {"2018"}, btag_catagories, {"120"}, 0.845, 1.117)
  ({"mt"}, {"2016"}, btag_catagories, {"120"}, 0.942, 1.105)
  ({"mt"}, {"2017"}, btag_catagories, {"120"}, 0.927, 1.089)
  ({"mt"}, {"2018"}, btag_catagories, {"120"}, 0.887, 1.127)
  ({"tt"}, {"2016"}, btag_catagories, {"120"}, 0.957, 1.063)
  ({"tt"}, {"2017"}, btag_catagories, {"120"}, 0.947, 1.052)
  ({"tt"}, {"2018"}, btag_catagories, {"120"}, 0.885, 1.099)
  ({"em"}, {"2016"}, nobtag_catagories, {"120"}, 1.004, 0.997)
  ({"em"}, {"2017"}, nobtag_catagories, {"120"}, 1.004, 0.996)
  ({"em"}, {"2018"}, nobtag_catagories, {"120"}, 1.007, 0.996)
  ({"et"}, {"2016"}, nobtag_catagories, {"120"}, 1.003, 0.997)
  ({"et"}, {"2017"}, nobtag_catagories, {"120"}, 1.003, 0.998)
  ({"et"}, {"2018"}, nobtag_catagories, {"120"}, 1.008, 0.993)
  ({"mt"}, {"2016"}, nobtag_catagories, {"120"}, 1.002, 0.997)
  ({"mt"}, {"2017"}, nobtag_catagories, {"120"}, 1.003, 0.996)
  ({"mt"}, {"2018"}, nobtag_catagories, {"120"}, 1.005, 0.994)
  ({"tt"}, {"2016"}, nobtag_catagories, {"120"}, 1.002, 0.997)
  ({"tt"}, {"2017"}, nobtag_catagories, {"120"}, 1.003, 0.997)
  ({"tt"}, {"2018"}, nobtag_catagories, {"120"}, 1.009, 0.992)
  );
  
  cb.cp().process(mssm_ggH_signals).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"125"}, 0.936, 1.087)
  ({"em"}, {"2017"}, btag_catagories, {"125"}, 0.927, 1.124)
  ({"em"}, {"2018"}, btag_catagories, {"125"}, 0.877, 1.21)
  ({"et"}, {"2016"}, btag_catagories, {"125"}, 0.961, 1.079)
  ({"et"}, {"2017"}, btag_catagories, {"125"}, 0.899, 1.062)
  ({"et"}, {"2018"}, btag_catagories, {"125"}, 0.876, 1.101)
  ({"mt"}, {"2016"}, btag_catagories, {"125"}, 0.946, 1.093)
  ({"mt"}, {"2017"}, btag_catagories, {"125"}, 0.925, 1.073)
  ({"mt"}, {"2018"}, btag_catagories, {"125"}, 0.883, 1.125)
  ({"tt"}, {"2016"}, btag_catagories, {"125"}, 0.943, 1.036)
  ({"tt"}, {"2017"}, btag_catagories, {"125"}, 0.916, 1.082)
  ({"tt"}, {"2018"}, btag_catagories, {"125"}, 0.901, 1.097)
  ({"em"}, {"2016"}, nobtag_catagories, {"125"}, 1.002, 0.997)
  ({"em"}, {"2017"}, nobtag_catagories, {"125"}, 1.003, 0.995)
  ({"em"}, {"2018"}, nobtag_catagories, {"125"}, 1.006, 0.99)
  ({"et"}, {"2016"}, nobtag_catagories, {"125"}, 1.001, 0.997)
  ({"et"}, {"2017"}, nobtag_catagories, {"125"}, 1.005, 0.997)
  ({"et"}, {"2018"}, nobtag_catagories, {"125"}, 1.006, 0.995)
  ({"mt"}, {"2016"}, nobtag_catagories, {"125"}, 1.002, 0.997)
  ({"mt"}, {"2017"}, nobtag_catagories, {"125"}, 1.003, 0.997)
  ({"mt"}, {"2018"}, nobtag_catagories, {"125"}, 1.005, 0.993)
  ({"tt"}, {"2016"}, nobtag_catagories, {"125"}, 1.003, 0.998)
  ({"tt"}, {"2017"}, nobtag_catagories, {"125"}, 1.004, 0.996)
  ({"tt"}, {"2018"}, nobtag_catagories, {"125"}, 1.006, 0.994)
  );
  
  cb.cp().process(mssm_ggH_signals).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"130"}, 0.88, 1.093)
  ({"em"}, {"2017"}, btag_catagories, {"130"}, 0.875, 1.143)
  ({"em"}, {"2018"}, btag_catagories, {"130"}, 0.875, 1.121)
  ({"et"}, {"2016"}, btag_catagories, {"130"}, 0.945, 1.069)
  ({"et"}, {"2017"}, btag_catagories, {"130"}, 0.942, 1.081)
  ({"et"}, {"2018"}, btag_catagories, {"130"}, 0.899, 1.17)
  ({"mt"}, {"2016"}, btag_catagories, {"130"}, 0.916, 1.13)
  ({"mt"}, {"2017"}, btag_catagories, {"130"}, 0.923, 1.058)
  ({"mt"}, {"2018"}, btag_catagories, {"130"}, 0.88, 1.12)
  ({"tt"}, {"2016"}, btag_catagories, {"130"}, 0.967, 1.081)
  ({"tt"}, {"2017"}, btag_catagories, {"130"}, 0.916, 1.079)
  ({"tt"}, {"2018"}, btag_catagories, {"130"}, 0.889, 1.114)
  ({"em"}, {"2016"}, nobtag_catagories, {"130"}, 1.003, 0.998)
  ({"em"}, {"2017"}, nobtag_catagories, {"130"}, 1.005, 0.994)
  ({"em"}, {"2018"}, nobtag_catagories, {"130"}, 1.006, 0.995)
  ({"et"}, {"2016"}, nobtag_catagories, {"130"}, 1.001, 0.998)
  ({"et"}, {"2017"}, nobtag_catagories, {"130"}, 1.003, 0.996)
  ({"et"}, {"2018"}, nobtag_catagories, {"130"}, 1.004, 0.992)
  ({"mt"}, {"2016"}, nobtag_catagories, {"130"}, 1.002, 0.996)
  ({"mt"}, {"2017"}, nobtag_catagories, {"130"}, 1.003, 0.998)
  ({"mt"}, {"2018"}, nobtag_catagories, {"130"}, 1.006, 0.994)
  ({"tt"}, {"2016"}, nobtag_catagories, {"130"}, 1.002, 0.996)
  ({"tt"}, {"2017"}, nobtag_catagories, {"130"}, 1.004, 0.996)
  ({"tt"}, {"2018"}, nobtag_catagories, {"130"}, 1.007, 0.993)
  );
  
  cb.cp().process(mssm_ggH_signals).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"140"}, 0.937, 1.104)
  ({"em"}, {"2017"}, btag_catagories, {"140"}, 0.852, 1.096)
  ({"em"}, {"2018"}, btag_catagories, {"140"}, 0.863, 1.149)
  ({"et"}, {"2016"}, btag_catagories, {"140"}, 0.943, 1.066)
  ({"et"}, {"2017"}, btag_catagories, {"140"}, 0.914, 1.076)
  ({"et"}, {"2018"}, btag_catagories, {"140"}, 0.862, 1.142)
  ({"mt"}, {"2016"}, btag_catagories, {"140"}, 0.937, 1.069)
  ({"mt"}, {"2017"}, btag_catagories, {"140"}, 0.901, 1.076)
  ({"mt"}, {"2018"}, btag_catagories, {"140"}, 0.893, 1.139)
  ({"tt"}, {"2016"}, btag_catagories, {"140"}, 0.959, 1.106)
  ({"tt"}, {"2017"}, btag_catagories, {"140"}, 0.926, 1.061)
  ({"tt"}, {"2018"}, btag_catagories, {"140"}, 0.909, 1.128)
  ({"em"}, {"2016"}, nobtag_catagories, {"140"}, 1.002, 0.997)
  ({"em"}, {"2017"}, nobtag_catagories, {"140"}, 1.006, 0.996)
  ({"em"}, {"2018"}, nobtag_catagories, {"140"}, 1.006, 0.994)
  ({"et"}, {"2016"}, nobtag_catagories, {"140"}, 1.002, 0.998)
  ({"et"}, {"2017"}, nobtag_catagories, {"140"}, 1.004, 0.997)
  ({"et"}, {"2018"}, nobtag_catagories, {"140"}, 1.007, 0.992)
  ({"mt"}, {"2016"}, nobtag_catagories, {"140"}, 1.002, 0.998)
  ({"mt"}, {"2017"}, nobtag_catagories, {"140"}, 1.004, 0.997)
  ({"mt"}, {"2018"}, nobtag_catagories, {"140"}, 1.005, 0.993)
  ({"tt"}, {"2016"}, nobtag_catagories, {"140"}, 1.002, 0.996)
  ({"tt"}, {"2017"}, nobtag_catagories, {"140"}, 1.004, 0.997)
  ({"tt"}, {"2018"}, nobtag_catagories, {"140"}, 1.005, 0.993)
  );
  
  cb.cp().process(mssm_ggH_signals).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"160"}, 0.885, 1.117)
  ({"em"}, {"2017"}, btag_catagories, {"160"}, 0.886, 1.115)
  ({"em"}, {"2018"}, btag_catagories, {"160"}, 0.86, 1.137)
  ({"et"}, {"2016"}, btag_catagories, {"160"}, 0.927, 1.109)
  ({"et"}, {"2017"}, btag_catagories, {"160"}, 0.897, 1.059)
  ({"et"}, {"2018"}, btag_catagories, {"160"}, 0.867, 1.127)
  ({"mt"}, {"2016"}, btag_catagories, {"160"}, 0.93, 1.079)
  ({"mt"}, {"2017"}, btag_catagories, {"160"}, 0.909, 1.059)
  ({"mt"}, {"2018"}, btag_catagories, {"160"}, 0.874, 1.11)
  ({"tt"}, {"2016"}, btag_catagories, {"160"}, 0.924, 1.086)
  ({"tt"}, {"2017"}, btag_catagories, {"160"}, 0.916, 1.089)
  ({"tt"}, {"2018"}, btag_catagories, {"160"}, 0.88, 1.103)
  ({"em"}, {"2016"}, nobtag_catagories, {"160"}, 1.004, 0.996)
  ({"em"}, {"2017"}, nobtag_catagories, {"160"}, 1.005, 0.995)
  ({"em"}, {"2018"}, nobtag_catagories, {"160"}, 1.007, 0.993)
  ({"et"}, {"2016"}, nobtag_catagories, {"160"}, 1.002, 0.996)
  ({"et"}, {"2017"}, nobtag_catagories, {"160"}, 1.005, 0.998)
  ({"et"}, {"2018"}, nobtag_catagories, {"160"}, 1.007, 0.993)
  ({"mt"}, {"2016"}, nobtag_catagories, {"160"}, 1.002, 0.997)
  ({"mt"}, {"2017"}, nobtag_catagories, {"160"}, 1.004, 0.998)
  ({"mt"}, {"2018"}, nobtag_catagories, {"160"}, 1.006, 0.994)
  ({"tt"}, {"2016"}, nobtag_catagories, {"160"}, 1.003, 0.997)
  ({"tt"}, {"2017"}, nobtag_catagories, {"160"}, 1.004, 0.996)
  ({"tt"}, {"2018"}, nobtag_catagories, {"160"}, 1.007, 0.994)
  );
  
  cb.cp().process(mssm_ggH_signals).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"180"}, 0.912, 1.15)
  ({"em"}, {"2017"}, btag_catagories, {"180"}, 0.913, 1.09)
  ({"em"}, {"2018"}, btag_catagories, {"180"}, 0.862, 1.166)
  ({"et"}, {"2016"}, btag_catagories, {"180"}, 0.927, 1.063)
  ({"et"}, {"2017"}, btag_catagories, {"180"}, 0.926, 1.077)
  ({"et"}, {"2018"}, btag_catagories, {"180"}, 0.862, 1.098)
  ({"mt"}, {"2016"}, btag_catagories, {"180"}, 0.946, 1.089)
  ({"mt"}, {"2017"}, btag_catagories, {"180"}, 0.934, 1.076)
  ({"mt"}, {"2018"}, btag_catagories, {"180"}, 0.868, 1.118)
  ({"tt"}, {"2016"}, btag_catagories, {"180"}, 0.923, 1.09)
  ({"tt"}, {"2017"}, btag_catagories, {"180"}, 0.922, 1.07)
  ({"tt"}, {"2018"}, btag_catagories, {"180"}, 0.902, 1.122)
  ({"em"}, {"2016"}, nobtag_catagories, {"180"}, 1.003, 0.996)
  ({"em"}, {"2017"}, nobtag_catagories, {"180"}, 1.004, 0.996)
  ({"em"}, {"2018"}, nobtag_catagories, {"180"}, 1.007, 0.991)
  ({"et"}, {"2016"}, nobtag_catagories, {"180"}, 1.003, 0.998)
  ({"et"}, {"2017"}, nobtag_catagories, {"180"}, 1.004, 0.997)
  ({"et"}, {"2018"}, nobtag_catagories, {"180"}, 1.007, 0.995)
  ({"mt"}, {"2016"}, nobtag_catagories, {"180"}, 1.002, 0.997)
  ({"mt"}, {"2017"}, nobtag_catagories, {"180"}, 1.003, 0.997)
  ({"mt"}, {"2018"}, nobtag_catagories, {"180"}, 1.006, 0.994)
  ({"tt"}, {"2016"}, nobtag_catagories, {"180"}, 1.003, 0.997)
  ({"tt"}, {"2017"}, nobtag_catagories, {"180"}, 1.004, 0.997)
  ({"tt"}, {"2018"}, nobtag_catagories, {"180"}, 1.006, 0.993)
  );
  
  cb.cp().process(mssm_ggH_signals).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"200"}, 0.902, 1.093)
  ({"em"}, {"2017"}, btag_catagories, {"200"}, 0.891, 1.107)
  ({"em"}, {"2018"}, btag_catagories, {"200"}, 0.865, 1.152)
  ({"et"}, {"2016"}, btag_catagories, {"200"}, 0.927, 1.085)
  ({"et"}, {"2017"}, btag_catagories, {"200"}, 0.929, 1.067)
  ({"et"}, {"2018"}, btag_catagories, {"200"}, 0.889, 1.111)
  ({"mt"}, {"2016"}, btag_catagories, {"200"}, 0.926, 1.082)
  ({"mt"}, {"2017"}, btag_catagories, {"200"}, 0.917, 1.073)
  ({"mt"}, {"2018"}, btag_catagories, {"200"}, 0.88, 1.108)
  ({"tt"}, {"2016"}, btag_catagories, {"200"}, 0.94, 1.087)
  ({"tt"}, {"2017"}, btag_catagories, {"200"}, 0.915, 1.086)
  ({"tt"}, {"2018"}, btag_catagories, {"200"}, 0.87, 1.113)
  ({"em"}, {"2016"}, nobtag_catagories, {"200"}, 1.003, 0.997)
  ({"em"}, {"2017"}, nobtag_catagories, {"200"}, 1.005, 0.995)
  ({"em"}, {"2018"}, nobtag_catagories, {"200"}, 1.007, 0.992)
  ({"et"}, {"2016"}, nobtag_catagories, {"200"}, 1.003, 0.997)
  ({"et"}, {"2017"}, nobtag_catagories, {"200"}, 1.003, 0.997)
  ({"et"}, {"2018"}, nobtag_catagories, {"200"}, 1.006, 0.993)
  ({"mt"}, {"2016"}, nobtag_catagories, {"200"}, 1.003, 0.997)
  ({"mt"}, {"2017"}, nobtag_catagories, {"200"}, 1.004, 0.997)
  ({"mt"}, {"2018"}, nobtag_catagories, {"200"}, 1.007, 0.993)
  ({"tt"}, {"2016"}, nobtag_catagories, {"200"}, 1.002, 0.997)
  ({"tt"}, {"2017"}, nobtag_catagories, {"200"}, 1.004, 0.996)
  ({"tt"}, {"2018"}, nobtag_catagories, {"200"}, 1.007, 0.994)
  );
  
  cb.cp().process(mssm_ggH_signals).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"250"}, 0.907, 1.102)
  ({"em"}, {"2017"}, btag_catagories, {"250"}, 0.901, 1.111)
  ({"em"}, {"2018"}, btag_catagories, {"250"}, 0.844, 1.145)
  ({"et"}, {"2016"}, btag_catagories, {"250"}, 0.92, 1.053)
  ({"et"}, {"2017"}, btag_catagories, {"250"}, 0.921, 1.088)
  ({"et"}, {"2018"}, btag_catagories, {"250"}, 0.871, 1.077)
  ({"mt"}, {"2016"}, btag_catagories, {"250"}, 0.938, 1.098)
  ({"mt"}, {"2017"}, btag_catagories, {"250"}, 0.908, 1.109)
  ({"mt"}, {"2018"}, btag_catagories, {"250"}, 0.888, 1.145)
  ({"tt"}, {"2016"}, btag_catagories, {"250"}, 0.916, 1.064)
  ({"tt"}, {"2017"}, btag_catagories, {"250"}, 0.925, 1.067)
  ({"tt"}, {"2018"}, btag_catagories, {"250"}, 0.887, 1.12)
  ({"em"}, {"2016"}, nobtag_catagories, {"250"}, 1.003, 0.996)
  ({"em"}, {"2017"}, nobtag_catagories, {"250"}, 1.005, 0.994)
  ({"em"}, {"2018"}, nobtag_catagories, {"250"}, 1.009, 0.991)
  ({"et"}, {"2016"}, nobtag_catagories, {"250"}, 1.003, 0.998)
  ({"et"}, {"2017"}, nobtag_catagories, {"250"}, 1.004, 0.996)
  ({"et"}, {"2018"}, nobtag_catagories, {"250"}, 1.008, 0.995)
  ({"mt"}, {"2016"}, nobtag_catagories, {"250"}, 1.002, 0.996)
  ({"mt"}, {"2017"}, nobtag_catagories, {"250"}, 1.004, 0.995)
  ({"mt"}, {"2018"}, nobtag_catagories, {"250"}, 1.007, 0.991)
  ({"tt"}, {"2016"}, nobtag_catagories, {"250"}, 1.003, 0.997)
  ({"tt"}, {"2017"}, nobtag_catagories, {"250"}, 1.004, 0.997)
  ({"tt"}, {"2018"}, nobtag_catagories, {"250"}, 1.007, 0.993)
  );
  
  cb.cp().process(mssm_ggH_signals).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"300"}, 0.913, 1.085)
  ({"em"}, {"2017"}, btag_catagories, {"300"}, 0.87, 1.087)
  ({"em"}, {"2018"}, btag_catagories, {"300"}, 0.876, 1.158)
  ({"et"}, {"2016"}, btag_catagories, {"300"}, 0.943, 1.091)
  ({"et"}, {"2017"}, btag_catagories, {"300"}, 0.912, 1.077)
  ({"et"}, {"2018"}, btag_catagories, {"300"}, 0.895, 1.143)
  ({"mt"}, {"2016"}, btag_catagories, {"300"}, 0.91, 1.086)
  ({"mt"}, {"2017"}, btag_catagories, {"300"}, 0.918, 1.073)
  ({"mt"}, {"2018"}, btag_catagories, {"300"}, 0.881, 1.122)
  ({"tt"}, {"2016"}, btag_catagories, {"300"}, 0.945, 1.069)
  ({"tt"}, {"2017"}, btag_catagories, {"300"}, 0.918, 1.066)
  ({"tt"}, {"2018"}, btag_catagories, {"300"}, 0.871, 1.13)
  ({"em"}, {"2016"}, nobtag_catagories, {"300"}, 1.003, 0.997)
  ({"em"}, {"2017"}, nobtag_catagories, {"300"}, 1.006, 0.996)
  ({"em"}, {"2018"}, nobtag_catagories, {"300"}, 1.008, 0.99)
  ({"et"}, {"2016"}, nobtag_catagories, {"300"}, 1.002, 0.997)
  ({"et"}, {"2017"}, nobtag_catagories, {"300"}, 1.005, 0.996)
  ({"et"}, {"2018"}, nobtag_catagories, {"300"}, 1.006, 0.991)
  ({"mt"}, {"2016"}, nobtag_catagories, {"300"}, 1.003, 0.997)
  ({"mt"}, {"2017"}, nobtag_catagories, {"300"}, 1.004, 0.996)
  ({"mt"}, {"2018"}, nobtag_catagories, {"300"}, 1.007, 0.992)
  ({"tt"}, {"2016"}, nobtag_catagories, {"300"}, 1.002, 0.997)
  ({"tt"}, {"2017"}, nobtag_catagories, {"300"}, 1.004, 0.997)
  ({"tt"}, {"2018"}, nobtag_catagories, {"300"}, 1.008, 0.992)
  );
  
  cb.cp().process(mssm_ggH_signals).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"350"}, 0.904, 1.098)
  ({"em"}, {"2017"}, btag_catagories, {"350"}, 0.898, 1.14)
  ({"em"}, {"2018"}, btag_catagories, {"350"}, 0.878, 1.154)
  ({"et"}, {"2016"}, btag_catagories, {"350"}, 0.924, 1.095)
  ({"et"}, {"2017"}, btag_catagories, {"350"}, 0.926, 1.083)
  ({"et"}, {"2018"}, btag_catagories, {"350"}, 1.0, 1.0)
  ({"mt"}, {"2016"}, btag_catagories, {"350"}, 0.94, 1.093)
  ({"mt"}, {"2017"}, btag_catagories, {"350"}, 0.918, 1.089)
  ({"mt"}, {"2018"}, btag_catagories, {"350"}, 1.0, 1.0)
  ({"tt"}, {"2016"}, btag_catagories, {"350"}, 0.933, 1.091)
  ({"tt"}, {"2017"}, btag_catagories, {"350"}, 0.917, 1.094)
  ({"tt"}, {"2018"}, btag_catagories, {"350"}, 1.0, 1.0)
  ({"em"}, {"2016"}, nobtag_catagories, {"350"}, 1.003, 0.997)
  ({"em"}, {"2017"}, nobtag_catagories, {"350"}, 1.004, 0.994)
  ({"em"}, {"2018"}, nobtag_catagories, {"350"}, 1.007, 0.991)
  ({"et"}, {"2016"}, nobtag_catagories, {"350"}, 1.003, 0.996)
  ({"et"}, {"2017"}, nobtag_catagories, {"350"}, 1.004, 0.996)
  ({"et"}, {"2018"}, nobtag_catagories, {"350"}, 1.0, 1.0)
  ({"mt"}, {"2016"}, nobtag_catagories, {"350"}, 1.002, 0.997)
  ({"mt"}, {"2017"}, nobtag_catagories, {"350"}, 1.004, 0.996)
  ({"mt"}, {"2018"}, nobtag_catagories, {"350"}, 1.0, 1.0)
  ({"tt"}, {"2016"}, nobtag_catagories, {"350"}, 1.003, 0.997)
  ({"tt"}, {"2017"}, nobtag_catagories, {"350"}, 1.004, 0.995)
  ({"tt"}, {"2018"}, nobtag_catagories, {"350"}, 1.0, 1.0)
  );
  
  cb.cp().process(mssm_ggH_signals).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"400"}, 0.892, 1.113)
  ({"em"}, {"2017"}, btag_catagories, {"400"}, 0.914, 1.116)
  ({"em"}, {"2018"}, btag_catagories, {"400"}, 0.853, 1.144)
  ({"et"}, {"2016"}, btag_catagories, {"400"}, 0.926, 1.063)
  ({"et"}, {"2017"}, btag_catagories, {"400"}, 0.912, 1.072)
  ({"et"}, {"2018"}, btag_catagories, {"400"}, 1.0, 1.0)
  ({"mt"}, {"2016"}, btag_catagories, {"400"}, 0.929, 1.086)
  ({"mt"}, {"2017"}, btag_catagories, {"400"}, 0.916, 1.081)
  ({"mt"}, {"2018"}, btag_catagories, {"400"}, 1.0, 1.0)
  ({"tt"}, {"2016"}, btag_catagories, {"400"}, 0.939, 1.095)
  ({"tt"}, {"2017"}, btag_catagories, {"400"}, 0.916, 1.073)
  ({"tt"}, {"2018"}, btag_catagories, {"400"}, 1.0, 1.0)
  ({"em"}, {"2016"}, nobtag_catagories, {"400"}, 1.004, 0.996)
  ({"em"}, {"2017"}, nobtag_catagories, {"400"}, 1.004, 0.995)
  ({"em"}, {"2018"}, nobtag_catagories, {"400"}, 1.008, 0.992)
  ({"et"}, {"2016"}, nobtag_catagories, {"400"}, 1.003, 0.997)
  ({"et"}, {"2017"}, nobtag_catagories, {"400"}, 1.004, 0.996)
  ({"et"}, {"2018"}, nobtag_catagories, {"400"}, 1.0, 1.0)
  ({"mt"}, {"2016"}, nobtag_catagories, {"400"}, 1.003, 0.997)
  ({"mt"}, {"2017"}, nobtag_catagories, {"400"}, 1.004, 0.996)
  ({"mt"}, {"2018"}, nobtag_catagories, {"400"}, 1.0, 1.0)
  ({"tt"}, {"2016"}, nobtag_catagories, {"400"}, 1.002, 0.996)
  ({"tt"}, {"2017"}, nobtag_catagories, {"400"}, 1.004, 0.996)
  ({"tt"}, {"2018"}, nobtag_catagories, {"400"}, 1.0, 1.0)
  );
  
  cb.cp().process(mssm_ggH_signals).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"450"}, 0.929, 1.103)
  ({"em"}, {"2017"}, btag_catagories, {"450"}, 0.89, 1.121)
  ({"em"}, {"2018"}, btag_catagories, {"450"}, 0.845, 1.126)
  ({"et"}, {"2016"}, btag_catagories, {"450"}, 0.926, 1.072)
  ({"et"}, {"2017"}, btag_catagories, {"450"}, 0.912, 1.099)
  ({"et"}, {"2018"}, btag_catagories, {"450"}, 0.106, 0.13)
  ({"mt"}, {"2016"}, btag_catagories, {"450"}, 0.924, 1.073)
  ({"mt"}, {"2017"}, btag_catagories, {"450"}, 0.913, 1.092)
  ({"mt"}, {"2018"}, btag_catagories, {"450"}, 0.112, 0.145)
  ({"tt"}, {"2016"}, btag_catagories, {"450"}, 0.943, 1.084)
  ({"tt"}, {"2017"}, btag_catagories, {"450"}, 0.917, 1.092)
  ({"tt"}, {"2018"}, btag_catagories, {"450"}, 0.558, 0.696)
  ({"em"}, {"2016"}, nobtag_catagories, {"450"}, 1.003, 0.996)
  ({"em"}, {"2017"}, nobtag_catagories, {"450"}, 1.005, 0.994)
  ({"em"}, {"2018"}, nobtag_catagories, {"450"}, 1.009, 0.992)
  ({"et"}, {"2016"}, nobtag_catagories, {"450"}, 1.003, 0.997)
  ({"et"}, {"2017"}, nobtag_catagories, {"450"}, 1.005, 0.995)
  ({"et"}, {"2018"}, nobtag_catagories, {"450"}, 0.121, 0.12)
  ({"mt"}, {"2016"}, nobtag_catagories, {"450"}, 1.003, 0.997)
  ({"mt"}, {"2017"}, nobtag_catagories, {"450"}, 1.004, 0.995)
  ({"mt"}, {"2018"}, nobtag_catagories, {"450"}, 0.118, 0.116)
  ({"tt"}, {"2016"}, nobtag_catagories, {"450"}, 1.002, 0.997)
  ({"tt"}, {"2017"}, nobtag_catagories, {"450"}, 1.004, 0.995)
  ({"tt"}, {"2018"}, nobtag_catagories, {"450"}, 0.641, 0.631)
  );
  
  cb.cp().process(mssm_ggH_signals).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"500"}, 0.911, 1.107)
  ({"em"}, {"2017"}, btag_catagories, {"500"}, 0.917, 1.13)
  ({"em"}, {"2018"}, btag_catagories, {"500"}, 0.866, 1.159)
  ({"et"}, {"2016"}, btag_catagories, {"500"}, 0.914, 1.09)
  ({"et"}, {"2017"}, btag_catagories, {"500"}, 0.916, 1.07)
  ({"et"}, {"2018"}, btag_catagories, {"500"}, 1.0, 1.0)
  ({"mt"}, {"2016"}, btag_catagories, {"500"}, 0.921, 1.048)
  ({"mt"}, {"2017"}, btag_catagories, {"500"}, 0.917, 1.088)
  ({"mt"}, {"2018"}, btag_catagories, {"500"}, 1.0, 1.0)
  ({"tt"}, {"2016"}, btag_catagories, {"500"}, 0.942, 1.087)
  ({"tt"}, {"2017"}, btag_catagories, {"500"}, 0.926, 1.08)
  ({"tt"}, {"2018"}, btag_catagories, {"500"}, 1.0, 1.0)
  ({"em"}, {"2016"}, nobtag_catagories, {"500"}, 1.004, 0.996)
  ({"em"}, {"2017"}, nobtag_catagories, {"500"}, 1.004, 0.994)
  ({"em"}, {"2018"}, nobtag_catagories, {"500"}, 1.008, 0.991)
  ({"et"}, {"2016"}, nobtag_catagories, {"500"}, 1.003, 0.997)
  ({"et"}, {"2017"}, nobtag_catagories, {"500"}, 1.005, 0.996)
  ({"et"}, {"2018"}, nobtag_catagories, {"500"}, 1.0, 1.0)
  ({"mt"}, {"2016"}, nobtag_catagories, {"500"}, 1.003, 0.998)
  ({"mt"}, {"2017"}, nobtag_catagories, {"500"}, 1.005, 0.995)
  ({"mt"}, {"2018"}, nobtag_catagories, {"500"}, 1.0, 1.0)
  ({"tt"}, {"2016"}, nobtag_catagories, {"500"}, 1.002, 0.996)
  ({"tt"}, {"2017"}, nobtag_catagories, {"500"}, 1.004, 0.996)
  ({"tt"}, {"2018"}, nobtag_catagories, {"500"}, 1.0, 1.0)
  );
  
  cb.cp().process(mssm_ggH_signals).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"600"}, 0.867, 1.071)
  ({"em"}, {"2017"}, btag_catagories, {"600"}, 0.944, 1.078)
  ({"em"}, {"2018"}, btag_catagories, {"600"}, 0.849, 1.196)
  ({"et"}, {"2016"}, btag_catagories, {"600"}, 0.953, 1.083)
  ({"et"}, {"2017"}, btag_catagories, {"600"}, 0.905, 1.059)
  ({"et"}, {"2018"}, btag_catagories, {"600"}, 0.86, 1.119)
  ({"mt"}, {"2016"}, btag_catagories, {"600"}, 0.913, 1.107)
  ({"mt"}, {"2017"}, btag_catagories, {"600"}, 0.928, 1.088)
  ({"mt"}, {"2018"}, btag_catagories, {"600"}, 0.901, 1.138)
  ({"tt"}, {"2016"}, btag_catagories, {"600"}, 0.967, 1.067)
  ({"tt"}, {"2017"}, btag_catagories, {"600"}, 0.931, 1.032)
  ({"tt"}, {"2018"}, btag_catagories, {"600"}, 0.868, 1.107)
  ({"em"}, {"2016"}, nobtag_catagories, {"600"}, 1.006, 0.997)
  ({"em"}, {"2017"}, nobtag_catagories, {"600"}, 1.002, 0.997)
  ({"em"}, {"2018"}, nobtag_catagories, {"600"}, 1.009, 0.989)
  ({"et"}, {"2016"}, nobtag_catagories, {"600"}, 1.001, 0.996)
  ({"et"}, {"2017"}, nobtag_catagories, {"600"}, 1.006, 0.996)
  ({"et"}, {"2018"}, nobtag_catagories, {"600"}, 1.008, 0.992)
  ({"mt"}, {"2016"}, nobtag_catagories, {"600"}, 1.003, 0.995)
  ({"mt"}, {"2017"}, nobtag_catagories, {"600"}, 1.004, 0.996)
  ({"mt"}, {"2018"}, nobtag_catagories, {"600"}, 1.008, 0.989)
  ({"tt"}, {"2016"}, nobtag_catagories, {"600"}, 1.002, 0.997)
  ({"tt"}, {"2017"}, nobtag_catagories, {"600"}, 1.004, 0.998)
  ({"tt"}, {"2018"}, nobtag_catagories, {"600"}, 1.009, 0.993)
  );
  
  cb.cp().process(mssm_ggH_signals).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"700"}, 0.932, 1.091)
  ({"em"}, {"2017"}, btag_catagories, {"700"}, 0.944, 1.092)
  ({"em"}, {"2018"}, btag_catagories, {"700"}, 0.817, 1.16)
  ({"et"}, {"2016"}, btag_catagories, {"700"}, 0.9, 1.095)
  ({"et"}, {"2017"}, btag_catagories, {"700"}, 0.938, 1.064)
  ({"et"}, {"2018"}, btag_catagories, {"700"}, 0.874, 1.121)
  ({"mt"}, {"2016"}, btag_catagories, {"700"}, 0.946, 1.064)
  ({"mt"}, {"2017"}, btag_catagories, {"700"}, 0.945, 1.096)
  ({"mt"}, {"2018"}, btag_catagories, {"700"}, 0.881, 1.11)
  ({"tt"}, {"2016"}, btag_catagories, {"700"}, 0.911, 1.068)
  ({"tt"}, {"2017"}, btag_catagories, {"700"}, 0.899, 1.06)
  ({"tt"}, {"2018"}, btag_catagories, {"700"}, 0.88, 1.134)
  ({"em"}, {"2016"}, nobtag_catagories, {"700"}, 1.003, 0.996)
  ({"em"}, {"2017"}, nobtag_catagories, {"700"}, 1.003, 0.994)
  ({"em"}, {"2018"}, nobtag_catagories, {"700"}, 1.011, 0.989)
  ({"et"}, {"2016"}, nobtag_catagories, {"700"}, 1.005, 0.996)
  ({"et"}, {"2017"}, nobtag_catagories, {"700"}, 1.003, 0.997)
  ({"et"}, {"2018"}, nobtag_catagories, {"700"}, 1.008, 0.992)
  ({"mt"}, {"2016"}, nobtag_catagories, {"700"}, 1.003, 0.997)
  ({"mt"}, {"2017"}, nobtag_catagories, {"700"}, 1.003, 0.995)
  ({"mt"}, {"2018"}, nobtag_catagories, {"700"}, 1.008, 0.992)
  ({"tt"}, {"2016"}, nobtag_catagories, {"700"}, 1.005, 0.996)
  ({"tt"}, {"2017"}, nobtag_catagories, {"700"}, 1.006, 0.996)
  ({"tt"}, {"2018"}, nobtag_catagories, {"700"}, 1.009, 0.99)
  );
  
  cb.cp().process(mssm_ggH_signals).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"800"}, 0.906, 1.122)
  ({"em"}, {"2017"}, btag_catagories, {"800"}, 0.862, 1.08)
  ({"em"}, {"2018"}, btag_catagories, {"800"}, 0.893, 1.186)
  ({"et"}, {"2016"}, btag_catagories, {"800"}, 0.935, 1.083)
  ({"et"}, {"2017"}, btag_catagories, {"800"}, 0.907, 1.088)
  ({"et"}, {"2018"}, btag_catagories, {"800"}, 0.845, 1.115)
  ({"mt"}, {"2016"}, btag_catagories, {"800"}, 0.96, 1.121)
  ({"mt"}, {"2017"}, btag_catagories, {"800"}, 0.93, 1.084)
  ({"mt"}, {"2018"}, btag_catagories, {"800"}, 0.882, 1.145)
  ({"tt"}, {"2016"}, btag_catagories, {"800"}, 0.928, 1.071)
  ({"tt"}, {"2017"}, btag_catagories, {"800"}, 0.923, 1.083)
  ({"tt"}, {"2018"}, btag_catagories, {"800"}, 0.886, 1.105)
  ({"em"}, {"2016"}, nobtag_catagories, {"800"}, 1.005, 0.994)
  ({"em"}, {"2017"}, nobtag_catagories, {"800"}, 1.007, 0.996)
  ({"em"}, {"2018"}, nobtag_catagories, {"800"}, 1.008, 0.986)
  ({"et"}, {"2016"}, nobtag_catagories, {"800"}, 1.003, 0.996)
  ({"et"}, {"2017"}, nobtag_catagories, {"800"}, 1.006, 0.994)
  ({"et"}, {"2018"}, nobtag_catagories, {"800"}, 1.011, 0.992)
  ({"mt"}, {"2016"}, nobtag_catagories, {"800"}, 1.002, 0.994)
  ({"mt"}, {"2017"}, nobtag_catagories, {"800"}, 1.004, 0.995)
  ({"mt"}, {"2018"}, nobtag_catagories, {"800"}, 1.008, 0.99)
  ({"tt"}, {"2016"}, nobtag_catagories, {"800"}, 1.004, 0.996)
  ({"tt"}, {"2017"}, nobtag_catagories, {"800"}, 1.005, 0.994)
  ({"tt"}, {"2018"}, nobtag_catagories, {"800"}, 1.009, 0.992)
  );
  
  cb.cp().process(mssm_ggH_signals).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"900"}, 0.926, 1.119)
  ({"em"}, {"2017"}, btag_catagories, {"900"}, 0.901, 1.135)
  ({"em"}, {"2018"}, btag_catagories, {"900"}, 0.842, 1.13)
  ({"et"}, {"2016"}, btag_catagories, {"900"}, 0.943, 1.08)
  ({"et"}, {"2017"}, btag_catagories, {"900"}, 0.914, 1.082)
  ({"et"}, {"2018"}, btag_catagories, {"900"}, 0.853, 1.143)
  ({"mt"}, {"2016"}, btag_catagories, {"900"}, 0.923, 1.078)
  ({"mt"}, {"2017"}, btag_catagories, {"900"}, 0.922, 1.081)
  ({"mt"}, {"2018"}, btag_catagories, {"900"}, 0.891, 1.127)
  ({"tt"}, {"2016"}, btag_catagories, {"900"}, 0.94, 1.099)
  ({"tt"}, {"2017"}, btag_catagories, {"900"}, 0.925, 1.082)
  ({"tt"}, {"2018"}, btag_catagories, {"900"}, 0.878, 1.117)
  ({"em"}, {"2016"}, nobtag_catagories, {"900"}, 1.004, 0.993)
  ({"em"}, {"2017"}, nobtag_catagories, {"900"}, 1.006, 0.991)
  ({"em"}, {"2018"}, nobtag_catagories, {"900"}, 1.013, 0.989)
  ({"et"}, {"2016"}, nobtag_catagories, {"900"}, 1.003, 0.995)
  ({"et"}, {"2017"}, nobtag_catagories, {"900"}, 1.005, 0.995)
  ({"et"}, {"2018"}, nobtag_catagories, {"900"}, 1.012, 0.988)
  ({"mt"}, {"2016"}, nobtag_catagories, {"900"}, 1.004, 0.996)
  ({"mt"}, {"2017"}, nobtag_catagories, {"900"}, 1.005, 0.995)
  ({"mt"}, {"2018"}, nobtag_catagories, {"900"}, 1.009, 0.99)
  ({"tt"}, {"2016"}, nobtag_catagories, {"900"}, 1.003, 0.995)
  ({"tt"}, {"2017"}, nobtag_catagories, {"900"}, 1.005, 0.995)
  ({"tt"}, {"2018"}, nobtag_catagories, {"900"}, 1.011, 0.99)
  );
  
  cb.cp().process(mssm_ggH_signals).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"1000"}, 0.909, 1.078)
  ({"em"}, {"2017"}, btag_catagories, {"1000"}, 0.913, 1.124)
  ({"em"}, {"2018"}, btag_catagories, {"1000"}, 0.829, 1.127)
  ({"et"}, {"2016"}, btag_catagories, {"1000"}, 0.942, 1.068)
  ({"et"}, {"2017"}, btag_catagories, {"1000"}, 0.922, 1.079)
  ({"et"}, {"2018"}, btag_catagories, {"1000"}, 0.876, 1.096)
  ({"mt"}, {"2016"}, btag_catagories, {"1000"}, 0.933, 1.064)
  ({"mt"}, {"2017"}, btag_catagories, {"1000"}, 0.921, 1.086)
  ({"mt"}, {"2018"}, btag_catagories, {"1000"}, 0.868, 1.132)
  ({"tt"}, {"2016"}, btag_catagories, {"1000"}, 0.931, 1.086)
  ({"tt"}, {"2017"}, btag_catagories, {"1000"}, 0.93, 1.079)
  ({"tt"}, {"2018"}, btag_catagories, {"1000"}, 0.879, 1.113)
  ({"em"}, {"2016"}, nobtag_catagories, {"1000"}, 1.005, 0.996)
  ({"em"}, {"2017"}, nobtag_catagories, {"1000"}, 1.006, 0.992)
  ({"em"}, {"2018"}, nobtag_catagories, {"1000"}, 1.014, 0.989)
  ({"et"}, {"2016"}, nobtag_catagories, {"1000"}, 1.003, 0.996)
  ({"et"}, {"2017"}, nobtag_catagories, {"1000"}, 1.005, 0.995)
  ({"et"}, {"2018"}, nobtag_catagories, {"1000"}, 1.011, 0.992)
  ({"mt"}, {"2016"}, nobtag_catagories, {"1000"}, 1.004, 0.996)
  ({"mt"}, {"2017"}, nobtag_catagories, {"1000"}, 1.006, 0.994)
  ({"mt"}, {"2018"}, nobtag_catagories, {"1000"}, 1.011, 0.989)
  ({"tt"}, {"2016"}, nobtag_catagories, {"1000"}, 1.004, 0.995)
  ({"tt"}, {"2017"}, nobtag_catagories, {"1000"}, 1.005, 0.994)
  ({"tt"}, {"2018"}, nobtag_catagories, {"1000"}, 1.011, 0.99)
  );
  
  cb.cp().process(mssm_ggH_signals).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"1200"}, 0.886, 1.135)
  ({"em"}, {"2017"}, btag_catagories, {"1200"}, 0.91, 1.109)
  ({"em"}, {"2018"}, btag_catagories, {"1200"}, 0.875, 1.095)
  ({"et"}, {"2016"}, btag_catagories, {"1200"}, 0.934, 1.101)
  ({"et"}, {"2017"}, btag_catagories, {"1200"}, 0.935, 1.112)
  ({"et"}, {"2018"}, btag_catagories, {"1200"}, 0.866, 1.134)
  ({"mt"}, {"2016"}, btag_catagories, {"1200"}, 0.938, 1.084)
  ({"mt"}, {"2017"}, btag_catagories, {"1200"}, 0.957, 1.082)
  ({"mt"}, {"2018"}, btag_catagories, {"1200"}, 0.863, 1.101)
  ({"tt"}, {"2016"}, btag_catagories, {"1200"}, 0.926, 1.074)
  ({"tt"}, {"2017"}, btag_catagories, {"1200"}, 0.918, 1.072)
  ({"tt"}, {"2018"}, btag_catagories, {"1200"}, 0.877, 1.134)
  ({"em"}, {"2016"}, nobtag_catagories, {"1200"}, 1.006, 0.993)
  ({"em"}, {"2017"}, nobtag_catagories, {"1200"}, 1.006, 0.993)
  ({"em"}, {"2018"}, nobtag_catagories, {"1200"}, 1.011, 0.992)
  ({"et"}, {"2016"}, nobtag_catagories, {"1200"}, 1.004, 0.993)
  ({"et"}, {"2017"}, nobtag_catagories, {"1200"}, 1.005, 0.992)
  ({"et"}, {"2018"}, nobtag_catagories, {"1200"}, 1.011, 0.989)
  ({"mt"}, {"2016"}, nobtag_catagories, {"1200"}, 1.004, 0.995)
  ({"mt"}, {"2017"}, nobtag_catagories, {"1200"}, 1.004, 0.995)
  ({"mt"}, {"2018"}, nobtag_catagories, {"1200"}, 1.012, 0.991)
  ({"tt"}, {"2016"}, nobtag_catagories, {"1200"}, 1.005, 0.995)
  ({"tt"}, {"2017"}, nobtag_catagories, {"1200"}, 1.006, 0.995)
  ({"tt"}, {"2018"}, nobtag_catagories, {"1200"}, 1.011, 0.988)
  );
  
  cb.cp().process(mssm_ggH_signals).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"1400"}, 0.954, 1.134)
  ({"em"}, {"2017"}, btag_catagories, {"1400"}, 0.871, 1.172)
  ({"em"}, {"2018"}, btag_catagories, {"1400"}, 0.883, 1.181)
  ({"et"}, {"2016"}, btag_catagories, {"1400"}, 0.967, 1.06)
  ({"et"}, {"2017"}, btag_catagories, {"1400"}, 0.933, 1.065)
  ({"et"}, {"2018"}, btag_catagories, {"1400"}, 0.88, 1.092)
  ({"mt"}, {"2016"}, btag_catagories, {"1400"}, 0.961, 1.085)
  ({"mt"}, {"2017"}, btag_catagories, {"1400"}, 0.926, 1.084)
  ({"mt"}, {"2018"}, btag_catagories, {"1400"}, 0.906, 1.117)
  ({"tt"}, {"2016"}, btag_catagories, {"1400"}, 0.947, 1.086)
  ({"tt"}, {"2017"}, btag_catagories, {"1400"}, 0.926, 1.09)
  ({"tt"}, {"2018"}, btag_catagories, {"1400"}, 0.87, 1.116)
  ({"em"}, {"2016"}, nobtag_catagories, {"1400"}, 1.003, 0.991)
  ({"em"}, {"2017"}, nobtag_catagories, {"1400"}, 1.008, 0.988)
  ({"em"}, {"2018"}, nobtag_catagories, {"1400"}, 1.01, 0.985)
  ({"et"}, {"2016"}, nobtag_catagories, {"1400"}, 1.002, 0.996)
  ({"et"}, {"2017"}, nobtag_catagories, {"1400"}, 1.006, 0.993)
  ({"et"}, {"2018"}, nobtag_catagories, {"1400"}, 1.011, 0.991)
  ({"mt"}, {"2016"}, nobtag_catagories, {"1400"}, 1.003, 0.995)
  ({"mt"}, {"2017"}, nobtag_catagories, {"1400"}, 1.005, 0.994)
  ({"mt"}, {"2018"}, nobtag_catagories, {"1400"}, 1.009, 0.988)
  ({"tt"}, {"2016"}, nobtag_catagories, {"1400"}, 1.004, 0.994)
  ({"tt"}, {"2017"}, nobtag_catagories, {"1400"}, 1.006, 0.993)
  ({"tt"}, {"2018"}, nobtag_catagories, {"1400"}, 1.013, 0.988)
  );
  
  cb.cp().process(mssm_ggH_signals).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"1600"}, 0.915, 1.079)
  ({"em"}, {"2017"}, btag_catagories, {"1600"}, 0.931, 1.134)
  ({"em"}, {"2018"}, btag_catagories, {"1600"}, 0.833, 1.145)
  ({"et"}, {"2016"}, btag_catagories, {"1600"}, 0.943, 1.112)
  ({"et"}, {"2017"}, btag_catagories, {"1600"}, 0.936, 1.063)
  ({"et"}, {"2018"}, btag_catagories, {"1600"}, 0.906, 1.143)
  ({"mt"}, {"2016"}, btag_catagories, {"1600"}, 0.944, 1.086)
  ({"mt"}, {"2017"}, btag_catagories, {"1600"}, 0.924, 1.088)
  ({"mt"}, {"2018"}, btag_catagories, {"1600"}, 0.875, 1.107)
  ({"tt"}, {"2016"}, btag_catagories, {"1600"}, 0.922, 1.056)
  ({"tt"}, {"2017"}, btag_catagories, {"1600"}, 0.932, 1.095)
  ({"tt"}, {"2018"}, btag_catagories, {"1600"}, 0.874, 1.123)
  ({"em"}, {"2016"}, nobtag_catagories, {"1600"}, 1.006, 0.994)
  ({"em"}, {"2017"}, nobtag_catagories, {"1600"}, 1.005, 0.99)
  ({"em"}, {"2018"}, nobtag_catagories, {"1600"}, 1.016, 0.986)
  ({"et"}, {"2016"}, nobtag_catagories, {"1600"}, 1.004, 0.992)
  ({"et"}, {"2017"}, nobtag_catagories, {"1600"}, 1.005, 0.995)
  ({"et"}, {"2018"}, nobtag_catagories, {"1600"}, 1.009, 0.987)
  ({"mt"}, {"2016"}, nobtag_catagories, {"1600"}, 1.004, 0.994)
  ({"mt"}, {"2017"}, nobtag_catagories, {"1600"}, 1.006, 0.993)
  ({"mt"}, {"2018"}, nobtag_catagories, {"1600"}, 1.013, 0.989)
  ({"tt"}, {"2016"}, nobtag_catagories, {"1600"}, 1.006, 0.996)
  ({"tt"}, {"2017"}, nobtag_catagories, {"1600"}, 1.005, 0.992)
  ({"tt"}, {"2018"}, nobtag_catagories, {"1600"}, 1.013, 0.988)
  );
  
  cb.cp().process(mssm_ggH_signals).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"1800"}, 0.916, 1.049)
  ({"em"}, {"2017"}, btag_catagories, {"1800"}, 0.928, 1.07)
  ({"em"}, {"2018"}, btag_catagories, {"1800"}, 0.83, 1.145)
  ({"et"}, {"2016"}, btag_catagories, {"1800"}, 0.958, 1.097)
  ({"et"}, {"2017"}, btag_catagories, {"1800"}, 0.928, 1.064)
  ({"et"}, {"2018"}, btag_catagories, {"1800"}, 0.891, 1.121)
  ({"mt"}, {"2016"}, btag_catagories, {"1800"}, 0.932, 1.072)
  ({"mt"}, {"2017"}, btag_catagories, {"1800"}, 0.931, 1.112)
  ({"mt"}, {"2018"}, btag_catagories, {"1800"}, 0.896, 1.131)
  ({"tt"}, {"2016"}, btag_catagories, {"1800"}, 0.931, 1.091)
  ({"tt"}, {"2017"}, btag_catagories, {"1800"}, 0.941, 1.051)
  ({"tt"}, {"2018"}, btag_catagories, {"1800"}, 0.881, 1.119)
  ({"em"}, {"2016"}, nobtag_catagories, {"1800"}, 1.006, 0.996)
  ({"em"}, {"2017"}, nobtag_catagories, {"1800"}, 1.005, 0.995)
  ({"em"}, {"2018"}, nobtag_catagories, {"1800"}, 1.017, 0.986)
  ({"et"}, {"2016"}, nobtag_catagories, {"1800"}, 1.003, 0.994)
  ({"et"}, {"2017"}, nobtag_catagories, {"1800"}, 1.007, 0.993)
  ({"et"}, {"2018"}, nobtag_catagories, {"1800"}, 1.011, 0.986)
  ({"mt"}, {"2016"}, nobtag_catagories, {"1800"}, 1.005, 0.995)
  ({"mt"}, {"2017"}, nobtag_catagories, {"1800"}, 1.006, 0.991)
  ({"mt"}, {"2018"}, nobtag_catagories, {"1800"}, 1.01, 0.987)
  ({"tt"}, {"2016"}, nobtag_catagories, {"1800"}, 1.006, 0.992)
  ({"tt"}, {"2017"}, nobtag_catagories, {"1800"}, 1.006, 0.995)
  ({"tt"}, {"2018"}, nobtag_catagories, {"1800"}, 1.013, 0.987)
  );
  
  cb.cp().process(mssm_ggH_signals).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"2000"}, 0.937, 1.102)
  ({"em"}, {"2017"}, btag_catagories, {"2000"}, 0.89, 1.081)
  ({"em"}, {"2018"}, btag_catagories, {"2000"}, 0.821, 1.154)
  ({"et"}, {"2016"}, btag_catagories, {"2000"}, 0.96, 1.053)
  ({"et"}, {"2017"}, btag_catagories, {"2000"}, 0.902, 1.084)
  ({"et"}, {"2018"}, btag_catagories, {"2000"}, 0.886, 1.115)
  ({"mt"}, {"2016"}, btag_catagories, {"2000"}, 0.92, 1.075)
  ({"mt"}, {"2017"}, btag_catagories, {"2000"}, 0.952, 1.068)
  ({"mt"}, {"2018"}, btag_catagories, {"2000"}, 0.884, 1.121)
  ({"tt"}, {"2016"}, btag_catagories, {"2000"}, 0.939, 1.066)
  ({"tt"}, {"2017"}, btag_catagories, {"2000"}, 0.936, 1.057)
  ({"tt"}, {"2018"}, btag_catagories, {"2000"}, 0.884, 1.088)
  ({"em"}, {"2016"}, nobtag_catagories, {"2000"}, 1.005, 0.993)
  ({"em"}, {"2017"}, nobtag_catagories, {"2000"}, 1.01, 0.993)
  ({"em"}, {"2018"}, nobtag_catagories, {"2000"}, 1.018, 0.984)
  ({"et"}, {"2016"}, nobtag_catagories, {"2000"}, 1.003, 0.996)
  ({"et"}, {"2017"}, nobtag_catagories, {"2000"}, 1.007, 0.994)
  ({"et"}, {"2018"}, nobtag_catagories, {"2000"}, 1.011, 0.989)
  ({"mt"}, {"2016"}, nobtag_catagories, {"2000"}, 1.006, 0.995)
  ({"mt"}, {"2017"}, nobtag_catagories, {"2000"}, 1.004, 0.994)
  ({"mt"}, {"2018"}, nobtag_catagories, {"2000"}, 1.011, 0.988)
  ({"tt"}, {"2016"}, nobtag_catagories, {"2000"}, 1.005, 0.994)
  ({"tt"}, {"2017"}, nobtag_catagories, {"2000"}, 1.006, 0.994)
  ({"tt"}, {"2018"}, nobtag_catagories, {"2000"}, 1.013, 0.99)
  );
  
  cb.cp().process(mssm_ggH_signals).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"2300"}, 0.946, 1.136)
  ({"em"}, {"2017"}, btag_catagories, {"2300"}, 0.902, 1.102)
  ({"em"}, {"2018"}, btag_catagories, {"2300"}, 0.872, 1.133)
  ({"et"}, {"2016"}, btag_catagories, {"2300"}, 0.917, 1.081)
  ({"et"}, {"2017"}, btag_catagories, {"2300"}, 0.887, 1.097)
  ({"et"}, {"2018"}, btag_catagories, {"2300"}, 0.889, 1.101)
  ({"mt"}, {"2016"}, btag_catagories, {"2300"}, 0.949, 1.09)
  ({"mt"}, {"2017"}, btag_catagories, {"2300"}, 0.913, 1.066)
  ({"mt"}, {"2018"}, btag_catagories, {"2300"}, 0.875, 1.114)
  ({"tt"}, {"2016"}, btag_catagories, {"2300"}, 0.926, 1.081)
  ({"tt"}, {"2017"}, btag_catagories, {"2300"}, 0.924, 1.078)
  ({"tt"}, {"2018"}, btag_catagories, {"2300"}, 0.89, 1.105)
  ({"em"}, {"2016"}, nobtag_catagories, {"2300"}, 1.004, 0.991)
  ({"em"}, {"2017"}, nobtag_catagories, {"2300"}, 1.009, 0.99)
  ({"em"}, {"2018"}, nobtag_catagories, {"2300"}, 1.013, 0.987)
  ({"et"}, {"2016"}, nobtag_catagories, {"2300"}, 1.007, 0.993)
  ({"et"}, {"2017"}, nobtag_catagories, {"2300"}, 1.009, 0.993)
  ({"et"}, {"2018"}, nobtag_catagories, {"2300"}, 1.012, 0.989)
  ({"mt"}, {"2016"}, nobtag_catagories, {"2300"}, 1.004, 0.992)
  ({"mt"}, {"2017"}, nobtag_catagories, {"2300"}, 1.008, 0.995)
  ({"mt"}, {"2018"}, nobtag_catagories, {"2300"}, 1.014, 0.987)
  ({"tt"}, {"2016"}, nobtag_catagories, {"2300"}, 1.006, 0.993)
  ({"tt"}, {"2017"}, nobtag_catagories, {"2300"}, 1.007, 0.993)
  ({"tt"}, {"2018"}, nobtag_catagories, {"2300"}, 1.013, 0.988)
  );
  
  cb.cp().process(mssm_ggH_signals).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"2600"}, 0.933, 1.086)
  ({"em"}, {"2017"}, btag_catagories, {"2600"}, 0.911, 1.107)
  ({"em"}, {"2018"}, btag_catagories, {"2600"}, 0.843, 1.1)
  ({"et"}, {"2016"}, btag_catagories, {"2600"}, 0.969, 1.061)
  ({"et"}, {"2017"}, btag_catagories, {"2600"}, 0.92, 1.086)
  ({"et"}, {"2018"}, btag_catagories, {"2600"}, 0.897, 1.096)
  ({"mt"}, {"2016"}, btag_catagories, {"2600"}, 0.954, 1.064)
  ({"mt"}, {"2017"}, btag_catagories, {"2600"}, 0.927, 1.071)
  ({"mt"}, {"2018"}, btag_catagories, {"2600"}, 0.9, 1.119)
  ({"tt"}, {"2016"}, btag_catagories, {"2600"}, 0.955, 1.073)
  ({"tt"}, {"2017"}, btag_catagories, {"2600"}, 0.923, 1.064)
  ({"tt"}, {"2018"}, btag_catagories, {"2600"}, 0.891, 1.106)
  ({"em"}, {"2016"}, nobtag_catagories, {"2600"}, 1.004, 0.993)
  ({"em"}, {"2017"}, nobtag_catagories, {"2600"}, 1.009, 0.989)
  ({"em"}, {"2018"}, nobtag_catagories, {"2600"}, 1.019, 0.988)
  ({"et"}, {"2016"}, nobtag_catagories, {"2600"}, 1.002, 0.995)
  ({"et"}, {"2017"}, nobtag_catagories, {"2600"}, 1.008, 0.993)
  ({"et"}, {"2018"}, nobtag_catagories, {"2600"}, 1.011, 0.99)
  ({"mt"}, {"2016"}, nobtag_catagories, {"2600"}, 1.004, 0.993)
  ({"mt"}, {"2017"}, nobtag_catagories, {"2600"}, 1.007, 0.993)
  ({"mt"}, {"2018"}, nobtag_catagories, {"2600"}, 1.012, 0.985)
  ({"tt"}, {"2016"}, nobtag_catagories, {"2600"}, 1.004, 0.994)
  ({"tt"}, {"2017"}, nobtag_catagories, {"2600"}, 1.008, 0.994)
  ({"tt"}, {"2018"}, nobtag_catagories, {"2600"}, 1.013, 0.987)
  );
  
  cb.cp().process(mssm_ggH_signals).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"2900"}, 0.901, 1.063)
  ({"em"}, {"2017"}, btag_catagories, {"2900"}, 0.914, 1.073)
  ({"em"}, {"2018"}, btag_catagories, {"2900"}, 0.88, 1.128)
  ({"et"}, {"2016"}, btag_catagories, {"2900"}, 0.948, 1.102)
  ({"et"}, {"2017"}, btag_catagories, {"2900"}, 0.927, 1.096)
  ({"et"}, {"2018"}, btag_catagories, {"2900"}, 0.894, 1.089)
  ({"mt"}, {"2016"}, btag_catagories, {"2900"}, 0.938, 1.056)
  ({"mt"}, {"2017"}, btag_catagories, {"2900"}, 0.909, 1.081)
  ({"mt"}, {"2018"}, btag_catagories, {"2900"}, 0.887, 1.108)
  ({"tt"}, {"2016"}, btag_catagories, {"2900"}, 0.939, 1.067)
  ({"tt"}, {"2017"}, btag_catagories, {"2900"}, 0.927, 1.053)
  ({"tt"}, {"2018"}, btag_catagories, {"2900"}, 0.895, 1.122)
  ({"em"}, {"2016"}, nobtag_catagories, {"2900"}, 1.01, 0.994)
  ({"em"}, {"2017"}, nobtag_catagories, {"2900"}, 1.008, 0.993)
  ({"em"}, {"2018"}, nobtag_catagories, {"2900"}, 1.014, 0.985)
  ({"et"}, {"2016"}, nobtag_catagories, {"2900"}, 1.004, 0.992)
  ({"et"}, {"2017"}, nobtag_catagories, {"2900"}, 1.006, 0.992)
  ({"et"}, {"2018"}, nobtag_catagories, {"2900"}, 1.012, 0.99)
  ({"mt"}, {"2016"}, nobtag_catagories, {"2900"}, 1.006, 0.994)
  ({"mt"}, {"2017"}, nobtag_catagories, {"2900"}, 1.008, 0.993)
  ({"mt"}, {"2018"}, nobtag_catagories, {"2900"}, 1.013, 0.988)
  ({"tt"}, {"2016"}, nobtag_catagories, {"2900"}, 1.006, 0.994)
  ({"tt"}, {"2017"}, nobtag_catagories, {"2900"}, 1.007, 0.995)
  ({"tt"}, {"2018"}, nobtag_catagories, {"2900"}, 1.013, 0.985)
  );
  
  cb.cp().process(mssm_ggH_signals).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"3200"}, 0.902, 1.06)
  ({"em"}, {"2017"}, btag_catagories, {"3200"}, 0.877, 1.082)
  ({"em"}, {"2018"}, btag_catagories, {"3200"}, 0.908, 1.134)
  ({"et"}, {"2016"}, btag_catagories, {"3200"}, 0.932, 1.085)
  ({"et"}, {"2017"}, btag_catagories, {"3200"}, 0.929, 1.048)
  ({"et"}, {"2018"}, btag_catagories, {"3200"}, 0.877, 1.075)
  ({"mt"}, {"2016"}, btag_catagories, {"3200"}, 0.936, 1.058)
  ({"mt"}, {"2017"}, btag_catagories, {"3200"}, 0.929, 1.057)
  ({"mt"}, {"2018"}, btag_catagories, {"3200"}, 0.901, 1.081)
  ({"tt"}, {"2016"}, btag_catagories, {"3200"}, 0.946, 1.076)
  ({"tt"}, {"2017"}, btag_catagories, {"3200"}, 0.933, 1.075)
  ({"tt"}, {"2018"}, btag_catagories, {"3200"}, 0.883, 1.105)
  ({"em"}, {"2016"}, nobtag_catagories, {"3200"}, 1.009, 0.995)
  ({"em"}, {"2017"}, nobtag_catagories, {"3200"}, 1.013, 0.991)
  ({"em"}, {"2018"}, nobtag_catagories, {"3200"}, 1.01, 0.985)
  ({"et"}, {"2016"}, nobtag_catagories, {"3200"}, 1.005, 0.993)
  ({"et"}, {"2017"}, nobtag_catagories, {"3200"}, 1.006, 0.996)
  ({"et"}, {"2018"}, nobtag_catagories, {"3200"}, 1.014, 0.991)
  ({"mt"}, {"2016"}, nobtag_catagories, {"3200"}, 1.005, 0.995)
  ({"mt"}, {"2017"}, nobtag_catagories, {"3200"}, 1.006, 0.994)
  ({"mt"}, {"2018"}, nobtag_catagories, {"3200"}, 1.012, 0.991)
  ({"tt"}, {"2016"}, nobtag_catagories, {"3200"}, 1.005, 0.993)
  ({"tt"}, {"2017"}, nobtag_catagories, {"3200"}, 1.007, 0.992)
  ({"tt"}, {"2018"}, nobtag_catagories, {"3200"}, 1.014, 0.987)
  );
  
  cb.cp().process(mssm_ggH_signals).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
  ({"em"}, {"2016"}, btag_catagories, {"3500"}, 0.893, 1.086)
  ({"em"}, {"2017"}, btag_catagories, {"3500"}, 0.889, 1.077)
  ({"em"}, {"2018"}, btag_catagories, {"3500"}, 0.885, 1.114)
  ({"et"}, {"2016"}, btag_catagories, {"3500"}, 0.933, 1.079)
  ({"et"}, {"2017"}, btag_catagories, {"3500"}, 0.92, 1.071)
  ({"et"}, {"2018"}, btag_catagories, {"3500"}, 0.907, 1.116)
  ({"mt"}, {"2016"}, btag_catagories, {"3500"}, 0.931, 1.077)
  ({"mt"}, {"2017"}, btag_catagories, {"3500"}, 0.952, 1.089)
  ({"mt"}, {"2018"}, btag_catagories, {"3500"}, 0.902, 1.096)
  ({"tt"}, {"2016"}, btag_catagories, {"3500"}, 0.945, 1.084)
  ({"tt"}, {"2017"}, btag_catagories, {"3500"}, 0.924, 1.054)
  ({"tt"}, {"2018"}, btag_catagories, {"3500"}, 0.876, 1.108)
  ({"em"}, {"2016"}, nobtag_catagories, {"3500"}, 1.009, 0.993)
  ({"em"}, {"2017"}, nobtag_catagories, {"3500"}, 1.009, 0.994)
  ({"em"}, {"2018"}, nobtag_catagories, {"3500"}, 1.014, 0.986)
  ({"et"}, {"2016"}, nobtag_catagories, {"3500"}, 1.006, 0.993)
  ({"et"}, {"2017"}, nobtag_catagories, {"3500"}, 1.008, 0.993)
  ({"et"}, {"2018"}, nobtag_catagories, {"3500"}, 1.01, 0.987)
  ({"mt"}, {"2016"}, nobtag_catagories, {"3500"}, 1.006, 0.994)
  ({"mt"}, {"2017"}, nobtag_catagories, {"3500"}, 1.005, 0.991)
  ({"mt"}, {"2018"}, nobtag_catagories, {"3500"}, 1.01, 0.99)
  ({"tt"}, {"2016"}, nobtag_catagories, {"3500"}, 1.005, 0.992)
  ({"tt"}, {"2017"}, nobtag_catagories, {"3500"}, 1.008, 0.994)
  ({"tt"}, {"2018"}, nobtag_catagories, {"3500"}, 1.015, 0.987)
  );

  // ##########################################################################
  // Uncertainty: Electron energy scale
  // References:
  // - MC: https://twiki.cern.ch/twiki/bin/view/CMS/EgammaRunIIRecommendations#E_gamma_Energy_Corrections
  // - Embedding: ?
  // Notes:
  // - FIXME: References for embedding missing, need proper correlation accross years for mc, see here: https://twiki.cern.ch/twiki/bin/view/CMS/EgammaRunIIRecommendations#Recommendations_on_Combining_Sys
  // ##########################################################################

  // MC uncorrelated uncertainty

  cb.cp()
      .channel({"em", "et"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_scale_e", "shape", SystMap<>::init(1.00));


  // Only using electron resolution for SM categories
  cb.cp()
      .channel({"em", "et"})
      .process(mc_processes)
      .bin_id(mssm_categories, false)
      .AddSyst(cb, "CMS_res_e", "shape", SystMap<>::init(1.00));

  // Embedded uncorrelated uncertainty

  cb.cp()
      .channel({"em", "et"})
      .process({"EMB"})
      .AddSyst(cb, "CMS_scale_e_emb", "shape", SystMap<>::init(1.00));


  // ##########################################################################
  // Uncertainty: Tau energy scale
  // References:
  // Notes:
  // - Tau energy scale is split by decay mode.
  // - FIXME: References?
  // - FIXME: Need it for H->WW in mt, et, (and tt)?
  // ##########################################################################


  // Common component acting on MC
  std::vector<std::string> tau_es_processes = JoinStr({{"ZTT", "TTT", "TTL", "VVT", "VVL"}, signals, signals_HWW, mssm_signals});
  std::vector<std::string> tau_es_processes_emb = {"EMB"};
  if (sm){
      std::vector<std::string> tau_es_processes = JoinStr({{"ZTT", "TTT", "TTL", "VVT", "VVL"}, signals, signals_HWW, mssm_signals, jetFakes});
      std::vector<std::string> tau_es_processes_emb =  JoinStr({{"EMB"}, jetFakes});
  }
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process(tau_es_processes)
      .AddSyst(cb, "CMS_scale_t_1prong_$ERA","shape", SystMap<>::init(1.0));

  cb.cp()
      .channel({"et", "mt", "tt"})
      .process(tau_es_processes)
      .AddSyst(cb, "CMS_scale_t_1prong1pizero_$ERA","shape", SystMap<>::init(1.0));

  cb.cp()
      .channel({"et", "mt", "tt"})
      .process(tau_es_processes)
      .AddSyst(cb, "CMS_scale_t_3prong_$ERA", "shape", SystMap<>::init(1.0));

  cb.cp()
      .channel({"et", "mt", "tt"})
      .process(tau_es_processes)
      .AddSyst(cb, "CMS_scale_t_3prong1pizero_$ERA", "shape", SystMap<>::init(1.0));

  // Component for EMB only
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process(tau_es_processes_emb)
      .AddSyst(cb, "CMS_scale_t_emb_1prong_$ERA", "shape", SystMap<>::init(0.866));

  cb.cp()
      .channel({"et", "mt", "tt"})
      .process(tau_es_processes_emb)
      .AddSyst(cb, "CMS_scale_t_emb_1prong1pizero_$ERA", "shape", SystMap<>::init(0.866));

  cb.cp()
      .channel({"et", "mt", "tt"})
      .process(tau_es_processes_emb)
      .AddSyst(cb, "CMS_scale_t_emb_3prong_$ERA", "shape", SystMap<>::init(0.866));

  cb.cp()
      .channel({"et", "mt", "tt"})
      .process(tau_es_processes_emb)
      .AddSyst(cb, "CMS_scale_t_emb_3prong1pizero_$ERA", "shape", SystMap<>::init(0.866));

  // Common component acting on EMB
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"EMB"})
      .AddSyst(cb, "CMS_scale_t_1prong_$ERA", "shape", SystMap<>::init(0.5));

  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"EMB"})
      .AddSyst(cb, "CMS_scale_t_1prong1pizero_$ERA", "shape", SystMap<>::init(0.5));

  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"EMB"})
      .AddSyst(cb, "CMS_scale_t_3prong_$ERA", "shape", SystMap<>::init(0.5));

  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"EMB"})
      .AddSyst(cb, "CMS_scale_t_3prong1pizero_$ERA", "shape", SystMap<>::init(0.5));

  // ##########################################################################
  // Uncertainty: Jet energy scale
  // References:
  // - Talk in CMS Htt meeting by Daniel Winterbottom about regional JES splits:
  //   https://indico.cern.ch/event/740094/contributions/3055870/
  // Notes:
  // ##########################################################################

  // Regional JES
    // uncorrelated between eras
    cb.cp()
        .channel({"et", "mt", "tt", "em"})
        .process(mc_processes)
        .AddSyst(cb, "CMS_scale_j_Absolute_$ERA", "shape", SystMap<>::init(1.00));

    cb.cp()
        .channel({"et", "mt", "tt", "em"})
        .process(mc_processes)
        .AddSyst(cb, "CMS_scale_j_BBEC1_$ERA", "shape", SystMap<>::init(1.00));

    cb.cp()
        .channel({"et", "mt", "tt", "em"})
        .process(mc_processes)
        .AddSyst(cb, "CMS_scale_j_EC2_$ERA", "shape", SystMap<>::init(1.00));

    cb.cp()
        .channel({"et", "mt", "tt", "em"})
        .process(mc_processes)
        .AddSyst(cb, "CMS_scale_j_HF_$ERA", "shape", SystMap<>::init(1.00));

    cb.cp()
        .channel({"et", "mt", "tt", "em"})
        .process(mc_processes)
        .AddSyst(cb, "CMS_scale_j_RelativeSample_$ERA", "shape", SystMap<>::init(1.00));
    // correlated between eras
    cb.cp()
        .channel({"et", "mt", "tt", "em"})
        .process(mc_processes)
        .AddSyst(cb, "CMS_scale_j_Absolute", "shape", SystMap<>::init(1.00));

    cb.cp()
        .channel({"et", "mt", "tt", "em"})
        .process(mc_processes)
        .AddSyst(cb, "CMS_scale_j_BBEC1", "shape", SystMap<>::init(1.00));

    cb.cp()
        .channel({"et", "mt", "tt", "em"})
        .process(mc_processes)
        .AddSyst(cb, "CMS_scale_j_EC2", "shape", SystMap<>::init(1.00));

    cb.cp()
        .channel({"et", "mt", "tt", "em"})
        .process(mc_processes)
        .AddSyst(cb, "CMS_scale_j_HF", "shape", SystMap<>::init(1.00));

    cb.cp()
        .channel({"et", "mt", "tt", "em"})
        .process(mc_processes)
        .AddSyst(cb, "CMS_scale_j_FlavorQCD", "shape", SystMap<>::init(1.00));

    cb.cp()
        .channel({"et", "mt", "tt", "em"})
        .process(mc_processes)
        .AddSyst(cb, "CMS_scale_j_RelativeBal", "shape", SystMap<>::init(1.00));

  // JER
  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_res_j_$ERA", "shape", SystMap<>::init(1.00));

  // ##########################################################################
  // Uncertainty: MET energy scale and Recoil
  // References:
  // Notes:
  // - FIXME: Clustered vs unclustered MET? Inclusion of JES splitting?
  // - FIXME: References?
  // ##########################################################################
  if (sm){
    cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process({"TT", "TTT", "TTL", "TTJ", "VV", "VVT", "VVL", "VVJ", "ST"})
      .AddSyst(cb, "CMS_scale_met_unclustered", "shape", SystMap<>::init(1.00));
  }
  else{
    cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process({"TT", "TTT", "TTL", "TTJ", "VV", "VVT", "VVL", "VVJ", "ST"})
      .AddSyst(cb, "CMS_scale_met_unclustered_$ERA", "shape", SystMap<>::init(1.00));
  }

  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process(JoinStr({signals, signals_HWW, mssm_signals, {"ZTT", "ZL", "ZJ", "W"}}))
      .AddSyst(cb, "CMS_htt_boson_scale_met_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process(JoinStr({signals, signals_HWW, mssm_signals, {"ZTT", "ZL", "ZJ", "W"}}))
      .AddSyst(cb, "CMS_htt_boson_res_met_$ERA", "shape", SystMap<>::init(1.00));

  // met uncertainty templates are included from taking 100% variation in the correction
  // these are scaled here to take the correct 1-sigma ranges

  // small uncertainty decorrelated by channel to account for statistical uncertainties on corrections
  cb.cp()
      .process({"EMB"})
      .channel({"et", "mt", "tt"})
      .bin_id(mssm_categories)
      .AddSyst(cb, "scale_embed_met_$CHANNEL_$ERA", "shape", SystMap<>::init(0.1)); 

  // the other component of the uncertainty is systematic and correlated between channels (but decorrelated by era) 

  cb.cp()
      .process({"EMB"})
      .bin_id(mssm_categories)
      .channel({"tt"})
      .era({"2016"})
      .AddSyst(cb, "scale_embed_met_$ERA", "shape", SystMap<>::init(0.22));
  cb.cp()
      .process({"EMB"})
      .bin_id(mssm_categories)
      .channel({"tt"})
      .era({"2017"})
      .AddSyst(cb, "scale_embed_met_$ERA", "shape", SystMap<>::init(0.25));
  cb.cp()
      .process({"EMB"})
      .bin_id(mssm_categories)
      .channel({"tt"})
      .era({"2018"})
      .AddSyst(cb, "scale_embed_met_$ERA", "shape", SystMap<>::init(0.2));

  cb.cp()
      .process({"EMB"})
      .bin_id(mssm_categories)
      .channel({"mt"})
      .era({"2016"})
      .AddSyst(cb, "scale_embed_met_$ERA", "shape", SystMap<>::init(1.0));
  cb.cp()
      .process({"EMB"})
      .bin_id(mssm_categories)
      .channel({"mt"})
      .era({"2017"})
      .AddSyst(cb, "scale_embed_met_$ERA", "shape", SystMap<>::init(0.67));
  cb.cp()
      .process({"EMB"})
      .bin_id(mssm_categories)
      .channel({"mt"})
      .era({"2018"})
      .AddSyst(cb, "scale_embed_met_$ERA", "shape", SystMap<>::init(0.85));

  cb.cp()
      .process({"EMB"})
      .bin_id(mssm_categories)
      .channel({"et"})
      .era({"2016"})
      .AddSyst(cb, "scale_embed_met_$ERA", "shape", SystMap<>::init(0.84));
  cb.cp()
      .process({"EMB"})
      .bin_id(mssm_categories)
      .channel({"et"})
      .era({"2017"})
      .AddSyst(cb, "scale_embed_met_$ERA", "shape", SystMap<>::init(0.63));
  cb.cp()
      .process({"EMB"})
      .bin_id(mssm_categories)
      .channel({"et"})
      .era({"2018"})
      .AddSyst(cb, "scale_embed_met_$ERA", "shape", SystMap<>::init(0.73));

  // ##########################################################################
  // Uncertainty: Background normalizations
  // References:
  // Notes:
  // - FIXME: Remeasure QCD extrapolation factors for SS and ABCD methods?
  //          Current values are measured by KIT.
  // - FIXME: Adapt for fake factor and embedding
  // - FIXME: W uncertainties: Do we need lnN uncertainties based on the Ersatz
  //          study in Run1 (found in HIG-16043 uncertainty model)
  // - FIXME: References?
  // ##########################################################################

  // VV
  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process({"VVT", "VVJ", "VVL", "VV", "ST"})
      .AddSyst(cb, "CMS_htt_vvXsec", "lnN", SystMap<>::init(1.05));

  // TT
//  cb.cp()
//      .channel({"et", "mt", "tt", "em"})
//      .process({"TTT", "TTL", "TTJ", "TT"})
//      .AddSyst(cb, "CMS_htt_tjXsec", "lnN", SystMap<>::init(1.06));
//
  // use unconstrained rate parameter for ttbar yield
  // We don't need above uncertainty on cross section if using the rate parameter
  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process({"TTT", "TTL", "TTJ", "TT"})
      .AddSyst(cb, "rate_ttbar","rateParam",SystMap<>::init(1.0));
  cb.GetParameter("rate_ttbar")->set_range(0.5,1.5);

  // We can also remove the lumi and em trigger uncertainties for ttbar if using the rate parameter
  cb.FilterSysts([](ch::Systematic *syst) {
      return (syst->name().find("lumi") != string::npos || syst->name().find("CMS_eff_trigger_em") != string::npos) &&
        (syst->process() == "TT" || syst->process() == "TTT" || syst->process() == "TTL" || syst->process() == "TTJ");
  });

  // W
  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process({"W"})
      .AddSyst(cb, "CMS_htt_wjXsec", "lnN", SystMap<>::init(1.04));

  // Z
  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process({"ZTT", "ZL", "ZJ"})
      .AddSyst(cb, "CMS_htt_zjXsec", "lnN", SystMap<>::init(1.02));

  cb.cp()
      .channel({"em"})
      .process({"QCD"})
      .AddSyst(cb, "CMS_htt_qcd_0jet_rate_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"em"})
      .process({"QCD"})
      .AddSyst(cb, "CMS_htt_qcd_0jet_shape_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"em"})
      .process({"QCD"})
      .AddSyst(cb, "CMS_htt_qcd_0jet_shape2_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"em"})
      .process({"QCD"})
      .AddSyst(cb, "CMS_htt_qcd_1jet_rate_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"em"})
      .process({"QCD"})
      .AddSyst(cb, "CMS_htt_qcd_1jet_shape_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"em"})
      .process({"QCD"})
      .AddSyst(cb, "CMS_htt_qcd_1jet_shape2_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
     .channel({"em"})
     .process({"QCD"})
     .AddSyst(cb, "CMS_htt_qcd_2jet_rate_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
     .channel({"em"})
     .process({"QCD"})
     .AddSyst(cb, "CMS_htt_qcd_2jet_shape_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
     .channel({"em"})
     .process({"QCD"})
     .AddSyst(cb, "CMS_htt_qcd_2jet_shape2_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
     .channel({"em"})
     .process({"QCD"})
     .AddSyst(cb, "CMS_htt_qcd_iso", "shape", SystMap<>::init(1.00));
  if (!sm){
      // em closure correction uncertainties in btag categories.
      // from https://indico.cern.ch/event/999841/contributions/4199219/attachments/2176453/3675294/MSSM_HTT_20210122.pdf
    cb.cp()
      .channel({"em"})
      .process({"QCD"})
      .AddSyst(cb, "CMS_htt_qcd_nbtag_closure_stat_$ERA", "lnN", SystMap<bin_id>::init
              ({35,36,37}, 1.07));
    cb.cp()
      .channel({"em"})
      .process({"QCD"})
      .AddSyst(cb, "CMS_htt_qcd_nbtag_closure_syst", "lnN", SystMap<bin_id>::init
              ({35,36,37}, 1.05));
  }

  // ##########################################################################
  // Uncertainty: Drell-Yan LO->NLO reweighting
  // References:
  // Notes:
  // - FIXME: References?
  // ##########################################################################

  if (era == 2016) {
      cb.cp()
          .channel({"et", "mt", "tt", "em"})
          .process({"ZTT", "ZL", "ZJ"})
          .AddSyst(cb, "CMS_htt_dyShape_$ERA", "shape", SystMap<>::init(0.10));
  } else {
      cb.cp()
          .channel({"et", "mt", "tt", "em"})
          .process({"ZTT", "ZL", "ZJ"})
          .AddSyst(cb, "CMS_htt_dyShape", "shape", SystMap<>::init(0.10));
  }

  // ##########################################################################
  // Uncertainty: TT shape reweighting
  // References:
  // Notes:
  // - FIXME: References?
  // ##########################################################################

  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process({"TTT", "TTL", "TTJ", "TT"})
      .AddSyst(cb, "CMS_htt_ttbarShape", "shapeU", SystMap<>::init(1.00));
  cb.GetParameter("CMS_htt_ttbarShape")->set_range(-1.0,1.0);

  // ##########################################################################
  // Uncertainty: Electron/muon to tau fakes and ZL energy scale
  // References:
  // Notes:
  // - FIXME: References?
  // ##########################################################################

  // ZL energy scale split by decay mode
  cb.cp()
      .channel({"mt"})
      .process({"ZL"})
      .AddSyst(cb, "CMS_ZLShape_$CHANNEL_1prong_$ERA", "shape",
               SystMap<>::init(1.00));

  cb.cp()
      .channel({"mt"})
      .process({"ZL"})
      .AddSyst(cb, "CMS_ZLShape_$CHANNEL_1prong1pizero_$ERA", "shape",
               SystMap<>::init(1.00));

  // split by eta for SM categories to match HIG-19-010, for MSSM it is included as a single uncertainty so will be decorrelated (should be fine as the m_sv cut removes all the ZL anyway for the MSSM+SM combinations)
  cb.cp()
      .channel({"et"})
      .process({"ZL"})
      .bin_id(mssm_categories, false)
      .AddSyst(cb, "CMS_ZLShape_$CHANNEL_1prong_barrel_$ERA", "shape",
               SystMap<>::init(1.00));

  cb.cp()
      .channel({"et"})
      .process({"ZL"})
      .bin_id(mssm_categories, false)
      .AddSyst(cb, "CMS_ZLShape_$CHANNEL_1prong1pizero_barrel_$ERA", "shape",
               SystMap<>::init(1.00));
  cb.cp()
      .channel({"et"})
      .process({"ZL"})
      .bin_id(mssm_categories, false)
      .AddSyst(cb, "CMS_ZLShape_$CHANNEL_1prong_endcap_$ERA", "shape",
               SystMap<>::init(1.00));

  cb.cp()
      .channel({"et"})
      .process({"ZL"})
      .bin_id(mssm_categories, false)
      .AddSyst(cb, "CMS_ZLShape_$CHANNEL_1prong1pizero_endcap_$ERA", "shape",
               SystMap<>::init(1.00));

  //single eta bin for MSSM cats: 
  cb.cp()
      .channel({"et"})
      .process({"ZL"})
      .bin_id(mssm_categories)
      .AddSyst(cb, "CMS_ZLShape_$CHANNEL_1prong_$ERA", "shape",
               SystMap<>::init(1.00));

  cb.cp()
      .channel({"et"})
      .process({"ZL"})
      .bin_id(mssm_categories)
      .AddSyst(cb, "CMS_ZLShape_$CHANNEL_1prong1pizero_$ERA", "shape",
               SystMap<>::init(1.00));


  // Lepton fake rate uncertainties are kept as shape uncertainties for SM categories to match HIG-19-010 but converted to lnN uncertainties for MSSM categories

  // Electron fakes
  cb.cp()
      .channel({"et"})
      .process({"ZL"})
      .bin_id({mssm_categories}, false)
      .AddSyst(cb, "CMS_fake_e_BA_$ERA", "shape",
               SystMap<>::init(1.00));
  cb.cp()
      .channel({"et"})
      .process({"ZL"})
      .bin_id({mssm_categories}, false)
      .AddSyst(cb, "CMS_fake_e_EC_$ERA", "shape",
               SystMap<>::init(1.00));

  cb.cp().channel({"et"}).process({"ZL"}).bin_id(mssm_categories).AddSyst(cb, "CMS_fake_e_BA_$ERA", "lnN", SystMap<ch::syst::era>::init
    ({"2016"}, 1.234)
    ({"2017"}, 1.202)
    ({"2018"}, 1.149));
  
  cb.cp().channel({"et"}).process({"ZL"}).bin_id(mssm_categories).AddSyst(cb, "CMS_fake_e_EC_$ERA", "lnN", SystMap<ch::syst::era>::init
    ({"2016"}, 1.052)
    ({"2017"}, 1.084)
    ({"2018"}, 1.053));

  // Muon fakes
  cb.cp()
      .channel({"mt"})
      .process({"ZL"})
      .bin_id({mssm_categories}, false)
      .AddSyst(cb, "CMS_fake_m_WH1_$ERA", "shape",
               SystMap<>::init(1.00));
  cb.cp()
      .channel({"mt"})
      .process({"ZL"})
      .bin_id({mssm_categories}, false)
      .AddSyst(cb, "CMS_fake_m_WH2_$ERA", "shape",
               SystMap<>::init(1.00));
  cb.cp()
      .channel({"mt"})
      .process({"ZL"})
      .bin_id({mssm_categories}, false)
      .AddSyst(cb, "CMS_fake_m_WH3_$ERA", "shape",
               SystMap<>::init(1.00));
  cb.cp()
      .channel({"mt"})
      .process({"ZL"})
      .bin_id({mssm_categories}, false)
      .AddSyst(cb, "CMS_fake_m_WH4_$ERA", "shape",
               SystMap<>::init(1.00));
  cb.cp()
      .channel({"mt"})
      .process({"ZL"})
      .bin_id({mssm_categories}, false)
      .AddSyst(cb, "CMS_fake_m_WH5_$ERA", "shape",
               SystMap<>::init(1.00));

  cb.cp().channel({"mt"}).process({"ZL"}).bin_id(mssm_categories).AddSyst(cb, "CMS_fake_m_WH1_$ERA", "lnN", SystMap<ch::syst::era>::init
    ({"2016"}, 1.042)
    ({"2017"}, 1.058)
    ({"2018"}, 1.045));
  
  cb.cp().channel({"mt"}).process({"ZL"}).bin_id(mssm_categories).AddSyst(cb, "CMS_fake_m_WH2_$ERA", "lnN", SystMap<ch::syst::era>::init
    ({"2016"}, 1.033)
    ({"2017"}, 1.048)
    ({"2018"}, 1.039));
  
  cb.cp().channel({"mt"}).process({"ZL"}).bin_id(mssm_categories).AddSyst(cb, "CMS_fake_m_WH3_$ERA", "lnN", SystMap<ch::syst::era>::init
    ({"2016"}, 1.029)
    ({"2017"}, 1.054)
    ({"2018"}, 1.038));
  
  cb.cp().channel({"mt"}).process({"ZL"}).bin_id(mssm_categories).AddSyst(cb, "CMS_fake_m_WH4_$ERA", "lnN", SystMap<ch::syst::era>::init
    ({"2016"}, 1.054)
    ({"2017"}, 1.055)
    ({"2018"}, 1.047));
  
  cb.cp().channel({"mt"}).process({"ZL"}).bin_id(mssm_categories).AddSyst(cb, "CMS_fake_m_WH5_$ERA", "lnN", SystMap<ch::syst::era>::init
    ({"2016"}, 1.016)
    ({"2017"}, 1.046)
    ({"2018"}, 1.049));

  // ##########################################################################
  // Uncertainty: Theory uncertainties
  // References:
  // - Gluon-fusion WG1 uncertainty scheme:
  //   https://twiki.cern.ch/twiki/bin/view/CMS/HiggsWG/SignalModelingTools
  // Notes:
  // - FIXME: WG1 scheme currently NOT applied to ggHWW -> on purpose?
  // - FIXME: Add TopMassTreatment from HIG-16043 uncertainty model
  // - FIXME: Compare to HIG-16043 uncertainty model:
  //           - PDF uncertainties split by category?
  //           - QCDUnc uncertainties?
  //           - UEPS uncertainties?
  // - FIXME: Check VH QCD scale uncertainty
  // - FIXME: References?
  // ##########################################################################

  // Uncertainty on branching ratio for HTT at 125 GeV
  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process(signals)
      .AddSyst(cb, "BR_htt_THU", "lnN", SystMap<>::init(1.017));
  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process(signals)
      .AddSyst(cb, "BR_htt_PU_mq", "lnN", SystMap<>::init(1.0099));
  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process(signals)
      .AddSyst(cb, "BR_htt_PU_alphas", "lnN", SystMap<>::init(1.0062));
  // Uncertainty on branching ratio for HWW at 125 GeV
  cb.cp()
      .channel({"et", "mt", "tt", "em"})
     .process(signals_HWW)
     .AddSyst(cb, "BR_hww_THU", "lnN", SystMap<>::init(1.0099));
  cb.cp()
      .channel({"et", "mt", "tt", "em"})
     .process(signals_HWW)
     .AddSyst(cb, "BR_hww_PU_mq", "lnN", SystMap<>::init(1.0099));
  cb.cp()
      .channel({"et", "mt", "tt", "em"})
     .process(signals_HWW)
     .AddSyst(cb, "BR_hww_PU_alphas", "lnN", SystMap<>::init(1.0066));
  // QCD scale (non ggH & qqH signals)
  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process({"ZH125", "ZHWW125"})
      .AddSyst(cb, "QCDScale_VH", "lnN", SystMap<>::init(1.009));
  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process({"WH125", "WHWW125"})
      .AddSyst(cb, "QCDScale_VH", "lnN", SystMap<>::init(1.008));
  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process({"ttH125"})
      .AddSyst(cb, "QCDScale_ttH", "lnN", SystMap<>::init(1.08));

  // PDF
  cb.cp()
      .channel({"et", "mt", "tt", "em"})
     .process(JoinStr({signals_ggH,signals_ggHToWW, {"ggh"}}))
      .AddSyst(cb, "pdf_Higgs_gg", "lnN", SystMap<>::init(1.032));
  cb.cp()
     .channel({"et", "mt", "tt", "em"})
     .process(JoinStr({signals_qqH,signals_qqHToWW, {"qqh"}}))
     .AddSyst(cb, "pdf_Higgs_qqbar", "lnN", SystMap<>::init(1.021));
  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process({"ZH125", "ZHWW125"})
      .AddSyst(cb, "pdf_Higgs_VH", "lnN", SystMap<>::init(1.013));
  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process({"WH125", "WHWW125"})
      .AddSyst(cb, "pdf_Higgs_VH", "lnN", SystMap<>::init(1.018));
  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process({"ttH125"})
      .AddSyst(cb, "pdf_Higgs_ttH", "lnN", SystMap<>::init(1.036));
  if (sm)
  {
    // Gluon-fusion WG1 uncertainty scheme
    cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process(JoinStr({signals_ggH, {"ggh"}}))
      .AddSyst(cb, "THU_ggH_Mig01", "shape", SystMap<>::init(1.00));
    cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process(JoinStr({signals_ggH, {"ggh"}}))
      .AddSyst(cb, "THU_ggH_Mig12", "shape", SystMap<>::init(1.00));
    cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process(JoinStr({signals_ggH, {"ggh"}}))
      .AddSyst(cb, "THU_ggH_Mu", "shape", SystMap<>::init(1.00));
    cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process(JoinStr({signals_ggH, {"ggh"}}))
      .AddSyst(cb, "THU_ggH_PT120", "shape", SystMap<>::init(1.00));
    cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process(JoinStr({signals_ggH, {"ggh"}}))
      .AddSyst(cb, "THU_ggH_PT60", "shape", SystMap<>::init(1.00));
    cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process(JoinStr({signals_ggH, {"ggh"}}))
      .AddSyst(cb, "THU_ggH_Res", "shape", SystMap<>::init(1.00));
    cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process(JoinStr({signals_ggH, {"ggh"}}))
      .AddSyst(cb, "THU_ggH_VBF2j", "shape", SystMap<>::init(1.00));
    cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process(JoinStr({signals_ggH, {"ggh"}}))
      .AddSyst(cb, "THU_ggH_VBF3j", "shape", SystMap<>::init(1.00));
    cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process(JoinStr({signals_ggH, {"ggh"}}))
      .AddSyst(cb, "THU_ggH_qmtop", "shape", SystMap<>::init(1.00));
    cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process({signals_ggHToWW})
      .AddSyst(cb, "QCDScale_ggHWW", "lnN", SystMap<>::init(1.039));
    // VBF WG1 uncertainty scheme
   cb.cp()
     .channel({"et", "mt", "tt", "em"})
     .process(JoinStr({signals_qqH, {"qqh"}}))
     .AddSyst(cb, "THU_qqH_TOT", "shape", SystMap<>::init(1.00));
   cb.cp()
     .channel({"et", "mt", "tt", "em"})
     .process(JoinStr({signals_qqH, {"qqh"}}))
     .AddSyst(cb, "THU_qqH_PTH200", "shape", SystMap<>::init(1.00));
   cb.cp()
     .channel({"et", "mt", "tt", "em"})
     .process(JoinStr({signals_qqH, {"qqh"}}))
     .AddSyst(cb, "THU_qqH_Mjj60", "shape", SystMap<>::init(1.00));
   cb.cp()
     .channel({"et", "mt", "tt", "em"})
     .process(JoinStr({signals_qqH, {"qqh"}}))
     .AddSyst(cb, "THU_qqH_Mjj120", "shape", SystMap<>::init(1.00));
   cb.cp()
     .channel({"et", "mt", "tt", "em"})
     .process(JoinStr({signals_qqH, {"qqh"}}))
     .AddSyst(cb, "THU_qqH_Mjj350", "shape", SystMap<>::init(1.00));
   cb.cp()
     .channel({"et", "mt", "tt", "em"})
     .process(JoinStr({signals_qqH, {"qqh"}}))
     .AddSyst(cb, "THU_qqH_Mjj700", "shape", SystMap<>::init(1.00));
   cb.cp()
     .channel({"et", "mt", "tt", "em"})
     .process(JoinStr({signals_qqH, {"qqh"}}))
     .AddSyst(cb, "THU_qqH_Mjj1000", "shape", SystMap<>::init(1.00));
   cb.cp()
     .channel({"et", "mt", "tt", "em"})
     .process(JoinStr({signals_qqH, {"qqh"}}))
     .AddSyst(cb, "THU_qqH_Mjj1500", "shape", SystMap<>::init(1.00));
   cb.cp()
     .channel({"et", "mt", "tt", "em"})
     .process(JoinStr({signals_qqH, {"qqh"}}))
     .AddSyst(cb, "THU_qqH_25", "shape", SystMap<>::init(1.00));
   cb.cp()
     .channel({"et", "mt", "tt", "em"})
     .process(JoinStr({signals_qqH, {"qqh"}}))
     .AddSyst(cb, "THU_qqH_JET01", "shape", SystMap<>::init(1.00));
  }
   cb.cp()
     .channel({"et", "mt", "tt", "em"})
     .process({signals_qqHToWW})
     .AddSyst(cb, "QCDScale_qqH", "lnN", SystMap<>::init(1.005));
  // ##########################################################################
  // Uncertainty: Embedded events
  // References:
  // - https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauTauEmbeddingSamples2016
  // Notes:
  // ##########################################################################

  // Embedded Normalization: No Lumi, Zjxsec information used, instead derived from data using dimuon selection efficiency
  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process({"EMB"})
      .AddSyst(cb, "CMS_htt_doublemutrg_$ERA", "lnN", SystMap<>::init(1.04));

  // TTbar contamination in embedded events: 10% shape uncertainty of assumed ttbar->tautau event shape
  cb.cp()
    .channel({"et", "mt", "tt", "em"})
    .process({"EMB"})
    .AddSyst(cb, "CMS_htt_emb_ttbar_$ERA", "shape", SystMap<>::init(1.00));

  // Uncertainty of hadronic tau track efficiency correction
  // uncorrelated between eras
  cb.cp()
    .channel({"et", "mt", "tt"})
    .process({"EMB"})
    .AddSyst(cb, "CMS_3ProngEff_$ERA", "shape", SystMap<>::init(1.0));

  cb.cp()
    .channel({"et", "mt", "tt"})
    .process({"EMB"})
    .AddSyst(cb, "CMS_1ProngPi0Eff_$ERA", "shape", SystMap<>::init(1.0));

  // ##########################################################################
  // Uncertainty: Jet fakes MSSM part
  // References:
  // - https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauJet2TauFakes
  // Notes:
  // - FIXME: add 2017 norm uncertainties, and properly correlate across years
  // ##########################################################################
  if (!sm){
    std::string jet_bins[2] = {"njet0", "njet1"};
    std::string jet_bins_lt[2] = {"njets0", "njets1"};
    std::string unc_regions[3] = {"low", "med", "high"};
    std::string extra_uncs[2] = {"unc1", "unc2"};
    std::string qcd_tt_uncs[3] = {"unc1", "unc2", "unc3"};
    std::string wjets_uncs[4] = {"unc1", "unc2", "unc3", "unc4"};


    // QCD shape stat.
    for (auto njet: jet_bins) {
        //only add njets0 uncerts for nobtag categories
        std::vector<int> bins = mssm_categories;
        if(njet=="njet0") bins = mssm_nobtag_catagories;
        for (auto reg: unc_regions) {
            for (auto unc: qcd_tt_uncs) {
                cb.cp()
                    .bin_id(bins)
                    .channel({"et", "mt", "tt"})
                    .process({"jetFakes"})
                    .AddSyst(cb, "CMS_ff_total_qcd_stat_"+njet+"_jet_pt_"+reg+"_"+unc+"_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));
            }
        }
    }

    // W shape stat.
    for (auto njet: jet_bins) {
        //only add njets0 uncerts for nobtag categories
        std::vector<int> bins = mssm_categories;
        if(njet=="njet0") bins = mssm_nobtag_catagories;
        for (auto reg: unc_regions) {
            for (auto unc: wjets_uncs) {
                cb.cp()
                    .bin_id(bins)
                    .channel({"et", "mt"})
                    .process({"jetFakes"})
                    .AddSyst(cb, "CMS_ff_total_wjets_stat_"+njet+"_jet_pt_"+reg+"_"+unc+"_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));
            }
        }
    }

    // TT shape stat.
    for (auto reg: unc_regions) {
        for (auto unc: qcd_tt_uncs) {
            cb.cp()
                .channel({"et", "mt"})
                .process({"jetFakes"})
                .AddSyst(cb, "CMS_ff_total_ttbar_stat_jet_pt_"+reg+"_"+unc+"_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));
        }
    }

    for (auto unc: extra_uncs) {
        for (auto njet: jet_bins_lt) {
          //only add njets0 uncerts for nobtag categories
          std::vector<int> bins = mssm_categories;
          if(njet=="njet0") bins = mssm_nobtag_catagories;
            cb.cp()
                .bin_id(bins)
                .channel({"et", "mt"})
                .process({"jetFakes"})
                .AddSyst(cb, "CMS_ff_total_qcd_stat_ss_"+njet+"_"+unc+"_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));

            cb.cp()
                .channel({"et", "mt"})
                .process({"jetFakes"})
                .AddSyst(cb, "CMS_ff_total_wjets_stat_met_"+njet+"_"+unc+"_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));

            cb.cp()
                .channel({"et", "mt"})
                .process({"jetFakes"})
                .AddSyst(cb, "CMS_ff_total_wjets_stat_l_pt_"+njet+"_"+unc+"_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));
        }

        cb.cp()
            .channel({"et", "mt"})
            .process({"jetFakes"})
            .AddSyst(cb, "CMS_ff_total_qcd_stat_l_pt_"+unc+"_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));

        cb.cp()
            .channel({"et", "mt"})
            .process({"jetFakes"})
            .AddSyst(cb, "CMS_ff_total_qcd_stat_iso_"+unc+"_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));

        cb.cp()
            .channel({"et", "mt"})
            .process({"jetFakes"})
            .AddSyst(cb, "CMS_ff_total_qcd_stat_os_"+unc+"_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));

        cb.cp()
            .channel({"et", "mt"})
            .process({"jetFakes"})
            .AddSyst(cb, "CMS_ff_total_wjets_stat_extrap_"+unc+"_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));

        cb.cp()
            .channel({"et", "mt"})
            .process({"jetFakes"})
            .AddSyst(cb, "CMS_ff_total_ttbar_stat_met_"+unc+"_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));
        cb.cp()
            .channel({"et", "mt"})
            .process({"jetFakes"})
            .AddSyst(cb, "CMS_ff_total_ttbar_stat_l_pt_"+unc+"_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));

        cb.cp()
            .channel({"tt"})
            .process({"jetFakes"})
            .AddSyst(cb, "CMS_ff_total_qcd_stat_dR_"+unc+"_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));

        cb.cp()
            .channel({"tt"})
            .process({"jetFakes"})
            .AddSyst(cb, "CMS_ff_total_qcd_stat_pt_"+unc+"_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));
    }
    cb.cp()
        .channel({"et", "mt", "tt"})
        .process({"jetFakes"})
        .AddSyst(cb, "CMS_ff_total_qcd_syst_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));

    cb.cp()
        .channel({"et", "mt"})
        .process({"jetFakes"})
        .AddSyst(cb, "CMS_ff_total_qcd_syst_iso_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));

    cb.cp()
        .channel({"et", "mt", "tt"})
        .process({"jetFakes"})
        .AddSyst(cb, "CMS_ff_total_wjets_syst_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));

    cb.cp()
        .channel({"et", "mt"})
        .process({"jetFakes"})
        .AddSyst(cb, "CMS_ff_total_wjets_frac_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));

    cb.cp()
        .channel({"et", "mt", "tt"})
        .process({"jetFakes"})
        .AddSyst(cb, "CMS_ff_total_ttbar_syst_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));

    cb.cp()
        .channel({"et", "mt"})
        .process({"jetFakes"})
        .AddSyst(cb, "CMS_ff_total_ttbar_frac_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));

    cb.cp()
        .channel({"et", "mt"})
        .process({"jetFakes"})
        .AddSyst(cb, "CMS_ff_total_low_pt_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));

    cb.cp()
        .channel({"tt"})
        .process({"wFakes"})
        .AddSyst(cb, "CMS_ff_total_wFakesNorm_$ERA", "lnN", SystMap<>::init(1.2));

    cb.cp()
        .channel({"et", "mt", "tt"})
        .process({"jetFakes"})
        .AddSyst(cb, "CMS_ff_total_sub_syst_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));

    // new additional uncertainties

    cb.cp()
        .channel({"et", "mt", "tt"})
        .process({"jetFakes"})
        .AddSyst(cb, "CMS_ff_total_syst_alt_func_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));

    cb.cp()
        .channel({"tt"})
        .process({"jetFakes"})
        .AddSyst(cb, "CMS_ff_total_qcd_syst_pt_2_closure_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));

    cb.cp()
        .channel({"tt"})
        .process({"jetFakes"})
        .AddSyst(cb, "CMS_ff_total_qcd_syst_dr_closure_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));

    cb.cp()
        .channel({"et", "mt", "tt"})
        .process({"jetFakes"})
        .AddSyst(cb, "CMS_ff_total_qcd_syst_met_closure_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));

    cb.cp()
        .channel({"et", "mt"})
        .process({"jetFakes"})
        .AddSyst(cb, "CMS_ff_total_ttbar_syst_met_closure_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));

    cb.cp()
        .channel({"et", "mt"})
        .process({"jetFakes"})
        .AddSyst(cb, "CMS_ff_total_wjets_syst_met_closure_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));

    cb.cp()
        .channel({"et", "mt"})
        .process({"jetFakes"})
        .AddSyst(cb, "CMS_ff_total_ttbar_syst_l_pt_closure_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));

    cb.cp()
        .channel({"et", "mt"})
        .process({"jetFakes"})
        .AddSyst(cb, "CMS_ff_total_wjets_syst_l_pt_closure_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));

    cb.cp()
        .channel({"et", "mt"})
        .process({"jetFakes"})
        .AddSyst(cb, "CMS_ff_total_wjets_syst_bkg_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));

    cb.cp()
        .channel({"et", "mt"})
        .process({"jetFakes"})
        .AddSyst(cb, "CMS_ff_total_qcd_syst_bkg_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));

  }
  else {
    // QCD shape stat.
    cb.cp()
        .channel({"et", "mt", "tt"})
        .process({"jetFakesSM"})
        .AddSyst(cb, "CMS_ffSM_qcd_njet0_morphed_stat_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));
    cb.cp()
        .channel({"et", "mt", "tt"})
        .process({"jetFakesSM"})
        .AddSyst(cb, "CMS_ffSM_qcd_njet1_morphed_stat_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));
    cb.cp()
        .channel({"et", "mt", "tt"})
        .process({"jetFakesSM"})
        .AddSyst(cb, "CMS_ffSM_qcd_njet2_morphed_stat_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));


    // W shape stat.
    cb.cp()
        .channel({"et", "mt"})
        .process({"jetFakesSM"})
        .AddSyst(cb, "CMS_ffSM_w_lowdR_njet0_morphed_stat_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));
    cb.cp()
        .channel({"et", "mt"})
        .process({"jetFakesSM"})
        .AddSyst(cb, "CMS_ffSM_w_lowdR_njet1_morphed_stat_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));
    cb.cp()
        .channel({"et", "mt"})
        .process({"jetFakesSM"})
        .AddSyst(cb, "CMS_ffSM_w_lowdR_njet2_morphed_stat_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));
    cb.cp()
        .channel({"et", "mt"})
        .process({"jetFakesSM"})
        .AddSyst(cb, "CMS_ffSM_w_highdR_njet0_morphed_stat_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));
    cb.cp()
        .channel({"et", "mt"})
        .process({"jetFakesSM"})
        .AddSyst(cb, "CMS_ffSM_w_highdR_njet1_morphed_stat_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));
    cb.cp()
        .channel({"et", "mt"})
        .process({"jetFakesSM"})
        .AddSyst(cb, "CMS_ffSM_w_highdR_njet2_morphed_stat_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));



    // TT shape stat.
    cb.cp()
        .channel({"et", "mt"})
        .process({"jetFakesSM"})
        .AddSyst(cb, "CMS_ffSM_tt_njet0_morphed_stat_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));

    cb.cp()
        .channel({"et", "mt"})
        .process({"jetFakesSM"})
        .AddSyst(cb, "CMS_ffSM_tt_njet1_morphed_stat_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));
    
    // MC subtraction uncertainty
    // uncorrelated between eras
    cb.cp()
        .channel({"et", "mt"})
        .process({"jetFakesSM"})
        .AddSyst(cb, "CMS_ffSM_w_mc_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
    cb.cp()
        .channel({"et", "mt", "tt"})
        .process({"jetFakesSM"})
        .AddSyst(cb, "CMS_ffSM_qcd_mc_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
    cb.cp()
        .channel({"et", "mt"})
        .process({"jetFakesSM"})
        .AddSyst(cb, "CMS_ffSM_frac_w_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));

        
    // Shape syst. of different contributions (QCD/W/tt)
    // uncorrelated between eras
    cb.cp()
        .channel({"et", "mt", "tt"})
        .process({"jetFakesSM"})
        .AddSyst(cb, "CMS_ffSM_qcd_mvis_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
    cb.cp()
        .channel({"et", "mt", "tt"})
        .process({"jetFakesSM"})
        .AddSyst(cb, "CMS_ffSM_qcd_mvis_osss_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
    cb.cp()
        .channel({"et", "mt", "tt"})
        .process({"jetFakesSM"})
        .AddSyst(cb, "CMS_ffSM_corr_qcd_mvis_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
    cb.cp()
        .channel({"et", "mt", "tt"})
        .process({"jetFakesSM"})
        .AddSyst(cb, "CMS_ffSM_corr_qcd_mvis_osss_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
    cb.cp()
        .channel({"et", "mt"})
        .process({"jetFakesSM"})
        .AddSyst(cb, "CMS_ffSM_qcd_muiso_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
    cb.cp()
        .channel({"et", "mt"})
        .process({"jetFakesSM"})
        .AddSyst(cb, "CMS_ffSM_corr_qcd_muiso_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
    cb.cp()
        .channel({"tt"})
        .process({"jetFakesSM"})
        .AddSyst(cb, "CMS_ffSM_qcd_tau2_pt_0jet_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
    cb.cp()
        .channel({"tt"})
        .process({"jetFakesSM"})
        .AddSyst(cb, "CMS_ffSM_corr_qcd_tau2_pt_0jet_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
    cb.cp()
        .channel({"tt"})
        .process({"jetFakesSM"})
        .AddSyst(cb, "CMS_ffSM_qcd_tau2_pt_1jet_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
    cb.cp()
        .channel({"tt"})
        .process({"jetFakesSM"})
        .AddSyst(cb, "CMS_ffSM_corr_qcd_tau2_pt_1jet_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
    cb.cp()
        .channel({"tt"})
        .process({"jetFakesSM"})
        .AddSyst(cb, "CMS_ffSM_w_syst_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
    cb.cp()
        .channel({"tt"})
        .process({"jetFakesSM"})
        .AddSyst(cb, "CMS_ffSM_tt_syst_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
    cb.cp()
        .channel({"et", "mt"})
        .process({"jetFakesSM"})
        .AddSyst(cb, "CMS_ffSM_tt_morphed_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
    cb.cp()
        .channel({"et", "mt"})
        .process({"jetFakesSM"})
        .AddSyst(cb, "CMS_ffSM_corr_tt_syst_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
    cb.cp()
        .channel({"et", "mt"})
        .process({"jetFakesSM"})
        .AddSyst(cb, "CMS_ffSM_tt_sf_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
    cb.cp()
        .channel({"et", "mt"})
        .process({"jetFakesSM"})
        .AddSyst(cb, "CMS_ffSM_w_lepPt_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
    cb.cp()
        .channel({"et", "mt"})
        .process({"jetFakesSM"})
        .AddSyst(cb, "CMS_ffSM_corr_w_lepPt_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
    cb.cp()
        .channel({"et", "mt"})
        .process({"jetFakesSM"})
        .AddSyst(cb, "CMS_ffSM_w_mt_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
    cb.cp()
        .channel({"et", "mt"})
        .process({"jetFakesSM"})
        .AddSyst(cb, "CMS_ffSM_corr_w_mt_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));


    //below: jetFakesSM norm uncertainties. Current values are for 2016, which are probably a good approx. for 2017. To be updated.


    // Stat. norm (uncorrelated across years)
    cb.cp()
        .channel({"et", "mt", "tt"})
        .process({"jetFakesSM"})
        .AddSyst(cb, "CMS_ffSM_norm_stat_$CHANNEL_$BIN_$ERA", "lnN", SystMap<channel, bin_id>::init
            ({"mt"}, {11},  1.04) //w
            ({"mt"}, {12},  1.052) //ztt
            ({"mt"}, {13},  1.051) //tt
            ({"mt"}, {14},  1.047) //ss
            ({"mt"}, {15},  1.04) //zll
            ({"mt"}, {16},  1.059) //misc
            ({"mt"}, {20},  1.052) //emb
            ({"mt"}, {21},  1.047) //ff
            ({"mt"}, {300}, 1.037) //incl
            ({"et"}, {11},  1.066) //w
            ({"et"}, {12},  1.095) //ztt
            ({"et"}, {13},  1.083) //tt
            ({"et"}, {14},  1.054) //ss
            ({"et"}, {15},  1.095) //zll
            ({"et"}, {16},  1.107) //misc
            ({"et"}, {20},  1.095) //emb
            ({"et"}, {21},  1.066) //ff
            ({"et"}, {300}, 1.065) //incl
            ({"tt"}, {12},  1.049) //ztt
            ({"tt"}, {16},  1.028) //misc
            ({"tt"}, {17},  1.041) //noniso
            ({"tt"}, {20},  1.049) //emb
            ({"tt"}, {21},  1.041) //ff
            ({"tt"}, {300}, 1.041) //incl
            );
    // ggH and qqH categories
    cb.cp()
        .channel({"et", "mt", "tt"})
        .process({"jetFakesSM"})
        .AddSyst(cb, "CMS_ffSM_norm_stat_$CHANNEL_ggH_$ERA", "lnN", SystMap<channel, bin_id>::init
            ({"mt"}, {1, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110},  1.049)
            ({"et"}, {1, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110},  1.074)
            ({"tt"}, {1, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110},  1.041)
            );

    cb.cp()
        .channel({"et", "mt", "tt"})
        .process({"jetFakesSM"})
        .AddSyst(cb, "CMS_ffSM_norm_stat_$CHANNEL_qqH_$ERA", "lnN", SystMap<channel, bin_id>::init
            ({"mt"}, {2, 200, 201, 202, 203},  1.068)
            ({"et"}, {2, 200, 201, 202, 203},  1.112)
            ({"tt"}, {2, 200, 201, 202, 203},  1.052)
            );

    // Syst. norm: Bin-correlated
    // uncorrelated between eras

    cb.cp()
        .channel({"et", "mt", "tt"})
        .process({"jetFakesSM"})
        .AddSyst(cb, "CMS_ffSM_jetbinned_stat_0jet_norm_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
    cb.cp()
        .channel({"et", "mt", "tt"})
        .process({"jetFakesSM"})
        .AddSyst(cb, "CMS_ffSM_jetbinned_stat_1jet_norm_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
    cb.cp()
        .channel({"et", "mt", "tt"})
        .process({"jetFakesSM"})
        .AddSyst(cb, "CMS_ffSM_jetbinned_stat_2jet_norm_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));

    // Syst. norm: Bin-dependent, correlated across years
    // uncorrelated between eras
    cb.cp()
        .channel({"et", "mt", "tt"})
        .process({"jetFakesSM"})
        .AddSyst(cb, "CMS_ffSM_sub_syst_$CHANNEL_$BIN_$ERA", "lnN", SystMap<channel, bin_id>::init
            ({"mt"}, {11},  1.025) //w
            ({"mt"}, {12},  1.045) //ztt
            ({"mt"}, {13},  1.03) //tt
            ({"mt"}, {14},  1.02) //ss
            ({"mt"}, {15},  1.04) //zll
            ({"mt"}, {16},  1.035) //misc
            ({"mt"}, {20},  1.045) //emb
            ({"mt"}, {21},  1.024) //ss
            ({"mt"}, {300}, 1.035) //incl
            ({"et"}, {11},  1.02) //w
            ({"et"}, {12},  1.04) //ztt
            ({"et"}, {13},  1.03) //tt
            ({"et"}, {14},  1.02) //ss
            ({"et"}, {15},  1.04) //zll
            ({"et"}, {16},  1.035) //misc
            ({"et"}, {20},  1.04) //emb
            ({"et"}, {21},  1.023) //ff
            ({"et"}, {300}, 1.035) //incl
            ({"tt"}, {12},  1.035) //ztt
            ({"tt"}, {16},  1.03) //misc
            ({"tt"}, {17},  1.02) //noniso
            ({"tt"}, {20},  1.035) //emb
            ({"tt"}, {21},  1.02) //ff
            ({"tt"}, {300}, 1.03) //incl
            );

    // ggH and qqH categories
    cb.cp()
        .channel({"et", "mt", "tt"})
        .process({"jetFakesSM"})
        .AddSyst(cb, "CMS_ffSM_sub_syst_$CHANNEL_ggH_$ERA", "lnN", SystMap<channel, bin_id>::init
            ({"mt"}, {1, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110},  1.04)
            ({"et"}, {1, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110},  1.04)
            ({"tt"}, {1, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110},  1.03)
            );

    cb.cp()
        .channel({"et", "mt", "tt"})
        .process({"jetFakesSM"})
        .AddSyst(cb, "CMS_ffSM_sub_syst_$CHANNEL_qqH_$ERA", "lnN", SystMap<channel, bin_id>::init
            ({"mt"}, {2, 200, 201, 202, 203},  1.04)
            ({"et"}, {2, 200, 201, 202, 203},  1.035)
            ({"tt"}, {2, 200, 201, 202, 203},  1.03)
            );
        }
}
} // namespace ch
