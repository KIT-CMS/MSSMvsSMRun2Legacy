#include "CombineHarvester/MSSMvsSMRun2Legacy/interface/HttSystematics_MSSMvsSMRun2.h"
#include "CombineHarvester/CombineTools/interface/Process.h"
#include "CombineHarvester/CombineTools/interface/Systematics.h"
#include "CombineHarvester/CombineTools/interface/Utilities.h"
#include <string>
#include <vector>
#include <fstream>
#include "CombineHarvester/CombineTools/interface/JsonTools.h"

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
      "ggH125",
      "ggH",
      "ggh",
      "ggH1",
      "ggX",
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
      "qqH125", 
      "qqH",
      "qqh",
      "qqH1",
      "qqX",
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
      "WH125", "Wh", "WH", "WH1", "ZH125", "Zh", "ZH", "ZH1", "ttH125"};
  std::vector<std::string> signals_ggHToWW = {
     // STXS stage 0
     "ggHWW125"};
  std::vector<std::string> signals_qqHToWW = {
     // STXS stage 0
     "qqHWW125"};
  std::vector<std::string> signals_VHToWW = {
      // STXS stage 0
      "WHWW125", "ZHWW125"};
  std::vector<std::string> signals = JoinStr({signals_ggH, signals_qqH, signals_VH});
  std::vector<std::string> signals_HWW = JoinStr({signals_ggHToWW, signals_qqHToWW, signals_VHToWW});

  std::vector<std::string> mssm_ggH_signals = {"ggH_t", "ggH_b", "ggH_i",
                                               "ggh_t", "ggh_b", "ggh_i",
                                               "ggA_t", "ggA_b", "ggA_i",
                                               "ggH2_t", "ggH2_b", "ggH2_i",
                                               "ggH1_t", "ggH1_b", "ggH1_i",
                                               "ggH3_t", "ggH3_b", "ggH3_i"};
  std::vector<std::string> mssm_bbH_signals = {"bbA", "bbH", "bbh", "bbH3", "bbH2", "bbH1", "bbH125"};
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
     ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {""}, 0.996) // for SM bbH125
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
     ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {""}, 0.995) // for SM bbH125
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
     ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {""}, 0.995) // for SM bbH125
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
     ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {"60"}, 1.023)
     ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {"80"}, 1.020)
     ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {"100"}, 1.017)
     ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {"120"}, 1.015)
     ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {"125"}, 1.015)
     ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {""}, 1.015) // for SM bbH125
     ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {"130"}, 1.015)
     ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {"140"}, 1.016)
     ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {"160"}, 1.017)
     ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {"180"}, 1.017)
     ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {"200"}, 1.016)
     ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {"250"}, 1.015)
     ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {"300"}, 1.014)
     ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {"350"}, 1.013)
     ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {"400"}, 1.012)
     ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {"450"}, 1.012)
     ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {"500"}, 1.011)
     ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {"600"}, 1.010)
     ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {"700"}, 1.012)
     ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {"800"}, 1.013)
     ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {"900"}, 1.014)
     ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {"1000"}, 1.016)
     ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {"1200"}, 1.017)
     ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {"1400"}, 1.018)
     ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {"1600"}, 1.018)
     ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {"1800"}, 1.019)
     ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {"2000"}, 1.019)
     ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {"2300"}, 1.020)
     ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {"2600"}, 1.022)
     ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {"2900"}, 1.023)
     ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {"3200"}, 1.025)
     ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {"3500"}, 1.027)
     ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {"60"}, 1.023)
     ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {"80"}, 1.019)
     ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {"100"}, 1.016)
     ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {"120"}, 1.016)
     ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {"125"}, 1.016)
     ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {""}, 1.016) // for SM bbH125
     ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {"130"}, 1.016)
     ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {"140"}, 1.015)
     ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {"160"}, 1.016)
     ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {"180"}, 1.016)
     ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {"200"}, 1.015)
     ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {"250"}, 1.014)
     ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {"300"}, 1.013)
     ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {"350"}, 1.012)
     ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {"400"}, 1.011)
     ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {"450"}, 1.011)
     ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {"500"}, 1.010)
     ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {"600"}, 1.010)
     ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {"700"}, 1.011)
     ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {"800"}, 1.012)
     ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {"900"}, 1.014)
     ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {"1000"}, 1.015)
     ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {"1200"}, 1.015)
     ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {"1400"}, 1.016)
     ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {"1600"}, 1.016)
     ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {"1800"}, 1.016)
     ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {"2000"}, 1.016)
     ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {"2300"}, 1.018)
     ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {"2600"}, 1.019)
     ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {"2900"}, 1.020)
     ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {"3200"}, 1.022)
     ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {"3500"}, 1.024)
     ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {"60"}, 1.022)
     ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {"80"}, 1.018)
     ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {"100"}, 1.015)
     ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {"120"}, 1.016)
     ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {"125"}, 1.016)
     ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {""}, 1.016) // for SM bbH125
     ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {"130"}, 1.015)
     ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {"140"}, 1.014)
     ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {"160"}, 1.016)
     ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {"180"}, 1.016)
     ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {"200"}, 1.015)
     ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {"250"}, 1.014)
     ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {"300"}, 1.012)
     ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {"350"}, 1.011)
     ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {"400"}, 1.010)
     ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {"450"}, 1.010)
     ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {"500"}, 1.010)
     ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {"600"}, 1.009)
     ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {"700"}, 1.010)
     ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {"800"}, 1.012)
     ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {"900"}, 1.013)
     ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {"1000"}, 1.015)
     ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {"1200"}, 1.015)
     ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {"1400"}, 1.016)
     ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {"1600"}, 1.016)
     ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {"1800"}, 1.016)
     ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {"2000"}, 1.016)
     ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {"2300"}, 1.018)
     ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {"2600"}, 1.019)
     ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {"2900"}, 1.020)
     ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {"3200"}, 1.022)
     ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {"3500"}, 1.024));

   cb.cp().process(mssm_bbH_signals).AddSyst(cb, "QCDscaleAndHdamp_bbH_ACCEPT", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
    ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"60"}, 1.010, 0.993)
    ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"80"}, 1.013, 0.992)
    ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"100"}, 1.016, 0.992)
    ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"120"}, 1.009, 0.993)
    ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {"125"}, 1.007, 0.993)
    ({"em","et","mt","tt"}, {"2016"}, nobtag_catagories, {""}, 1.007, 0.993) // for SM bbH125
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
    ({"em","et","mt","tt"}, {"2017"}, nobtag_catagories, {""}, 1.008, 0.992) // for SM bbH125
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
    ({"em","et","mt","tt"}, {"2018"}, nobtag_catagories, {""}, 1.009, 0.992) // for SM bbH125
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
    ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {"60"}, 0.938, 1.042)
    ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {"80"}, 0.932, 1.038)
    ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {"100"}, 0.926, 1.035)
    ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {"120"}, 0.964, 1.028)
    ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {"125"}, 0.974, 1.026)
    ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {""}, 0.974, 1.026) // for SM bbH125
    ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {"130"}, 0.969, 1.024)
    ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {"140"}, 0.958, 1.019)
    ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {"160"}, 0.958, 1.021)
    ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {"180"}, 0.959, 1.020)
    ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {"200"}, 0.961, 1.020)
    ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {"250"}, 0.964, 1.018)
    ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {"300"}, 0.967, 1.017)
    ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {"350"}, 0.971, 1.015)
    ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {"400"}, 0.974, 1.014)
    ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {"450"}, 0.976, 1.017)
    ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {"500"}, 0.978, 1.019)
    ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {"600"}, 0.983, 1.025)
    ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {"700"}, 0.978, 1.023)
    ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {"800"}, 0.974, 1.022)
    ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {"900"}, 0.970, 1.020)
    ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {"1000"}, 0.965, 1.019)
    ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {"1200"}, 0.966, 1.024)
    ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {"1400"}, 0.967, 1.030)
    ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {"1600"}, 0.964, 1.027)
    ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {"1800"}, 0.960, 1.024)
    ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {"2000"}, 0.957, 1.021)
    ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {"2300"}, 0.957, 1.022)
    ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {"2600"}, 0.956, 1.032)
    ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {"2900"}, 0.955, 1.042)
    ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {"3200"}, 0.950, 1.042)
    ({"et","mt","tt","em"}, {"2016"}, btag_catagories, {"3500"}, 0.946, 1.041)
    ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {"60"}, 0.944, 1.043)
    ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {"80"}, 0.932, 1.036)
    ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {"100"}, 0.920, 1.029)
    ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {"120"}, 0.962, 1.026)
    ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {"125"}, 0.972, 1.025)
    ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {""}, 0.972, 1.025) // for SM bbH125
    ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {"130"}, 0.969, 1.023)
    ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {"140"}, 0.962, 1.018)
    ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {"160"}, 0.960, 1.019)
    ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {"180"}, 0.962, 1.018)
    ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {"200"}, 0.963, 1.018)
    ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {"250"}, 0.967, 1.016)
    ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {"300"}, 0.970, 1.015)
    ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {"350"}, 0.974, 1.013)
    ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {"400"}, 0.978, 1.012)
    ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {"450"}, 0.980, 1.015)
    ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {"500"}, 0.981, 1.018)
    ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {"600"}, 0.984, 1.023)
    ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {"700"}, 0.980, 1.021)
    ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {"800"}, 0.976, 1.020)
    ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {"900"}, 0.972, 1.018)
    ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {"1000"}, 0.968, 1.017)
    ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {"1200"}, 0.969, 1.022)
    ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {"1400"}, 0.970, 1.027)
    ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {"1600"}, 0.967, 1.024)
    ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {"1800"}, 0.964, 1.022)
    ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {"2000"}, 0.961, 1.019)
    ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {"2300"}, 0.962, 1.020)
    ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {"2600"}, 0.961, 1.027)
    ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {"2900"}, 0.960, 1.035)
    ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {"3200"}, 0.956, 1.035)
    ({"et","mt","tt","em"}, {"2017"}, btag_catagories, {"3500"}, 0.952, 1.035)
    ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {"60"}, 0.941, 1.042)
    ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {"80"}, 0.930, 1.034)
    ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {"100"}, 0.920, 1.027)
    ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {"120"}, 0.962, 1.026)
    ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {"125"}, 0.972, 1.026)
    ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {""}, 0.972, 1.026) // for SM bbH125
    ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {"130"}, 0.969, 1.023)
    ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {"140"}, 0.962, 1.017)
    ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {"160"}, 0.961, 1.019)
    ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {"180"}, 0.962, 1.018)
    ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {"200"}, 0.964, 1.018)
    ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {"250"}, 0.968, 1.016)
    ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {"300"}, 0.972, 1.015)
    ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {"350"}, 0.975, 1.013)
    ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {"400"}, 0.979, 1.012)
    ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {"450"}, 0.980, 1.015)
    ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {"500"}, 0.982, 1.018)
    ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {"600"}, 0.985, 1.024)
    ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {"700"}, 0.981, 1.022)
    ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {"800"}, 0.978, 1.020)
    ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {"900"}, 0.974, 1.019)
    ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {"1000"}, 0.970, 1.017)
    ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {"1200"}, 0.970, 1.022)
    ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {"1400"}, 0.971, 1.027)
    ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {"1600"}, 0.968, 1.024)
    ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {"1800"}, 0.965, 1.021)
    ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {"2000"}, 0.962, 1.018)
    ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {"2300"}, 0.962, 1.020)
    ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {"2600"}, 0.962, 1.027)
    ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {"2900"}, 0.961, 1.034)
    ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {"3200"}, 0.957, 1.033)
    ({"et","mt","tt","em"}, {"2018"}, btag_catagories, {"3500"}, 0.953, 1.032));

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
  float lumi_unc_corr = 1.0;
  float lumi_unc_1718 = 1.0;
  if (era == 2016) {
      lumi_unc = 1.010;
      lumi_unc_corr = 1.006;
  } else if (era == 2017) {
      lumi_unc = 1.020;
      lumi_unc_corr = 1.009;
      lumi_unc_1718 = 1.006;
  } else if (era == 2018) {
      lumi_unc = 1.015;
      lumi_unc_corr = 1.020;
      lumi_unc_1718 = 1.002;
  }
  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process(mc_processes)
      .AddSyst(cb, "lumi_13TeV_$ERA", "lnN", SystMap<>::init(lumi_unc));
  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process(mc_processes)
      .AddSyst(cb, "lumi_13TeV_correlated", "lnN", SystMap<>::init(lumi_unc_corr));
  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process(mc_processes)
      .AddSyst(cb, "lumi_13TeV_1718", "lnN", SystMap<>::init(lumi_unc_1718));

  // ##########################################################################
  // Uncertainty: ggH Reweighting Hdamp uncertainty
  // References:
  // - 
  // Notes: for Hdamp scales t, b, and i components are decorrelated
  // ##########################################################################

   cb.cp()
       .channel({"et", "mt", "tt", "em"})
       .process({"ggh_i","ggH_i","ggA_i", "ggH1_i", "ggH2_i", "ggH3_i"})
       .AddSyst(cb, "Hdamp_ggH_i_REWEIGHT", "shape", SystMap<>::init(1.00));

   cb.cp()
       .channel({"et", "mt", "tt", "em"})
       .process({"ggh_t","ggH_t","ggA_t", "ggH1_t", "ggH2_t", "ggH3_t"})
       .AddSyst(cb, "Hdamp_ggH_t_REWEIGHT", "shape", SystMap<>::init(1.00));

   cb.cp()
       .channel({"et", "mt", "tt", "em"})
       .process({"ggh_b","ggH_b","ggA_b", "ggH1_b", "ggH2_b", "ggH3_b"})
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
      .process(JoinStr({signals, signals_HWW, mssm_signals, {"EMB", "ZTT", "TTT", "TTL", "VVT", "VVL"}}))
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
        .process(JoinStr({signals, mssm_signals, {"ZTT", "TTT", "TTL", "VVT", "VVL"}}))
        .AddSyst(cb, "CMS_eff_t_"+tauIDbin+"_$ERA", "shape", SystMap<>::init(1.0));
  }
  cb.cp() //second part uncorrelated between channels for IDvsLep
      .channel({"et", "mt"})
      .process(JoinStr({signals, mssm_signals, {"ZTT", "TTT", "TTL", "VVT", "VVL"}}))
      .AddSyst(cb, "CMS_eff_t_$CHANNEL_$ERA", "lnN", SystMap<>::init(1.01));

  // Tau ID: tt with 2 real taus
  for (auto tauIDbin : tauIDdmbins){
    cb.cp()
        .channel({"tt"})
        .process(JoinStr({signals, mssm_signals, {"ZTT", "TTT", "TTL", "VVT", "VVL"}}))
        .AddSyst(cb, "CMS_eff_t_dm"+tauIDbin+"_$ERA", "shape", SystMap<>::init(1.0));
  }
  if (!sm){
      cb.cp()
      .channel({"tt"})
      .process(JoinStr({signals, mssm_signals, {"ZTT", "TTT", "TTL", "VVT", "VVL"}}))
      .AddSyst(cb, "CMS_eff_t_highpT_100-500_$ERA", "shape", SystMap<>::init(1.0));
      cb.cp()
        .channel({"tt"})
        .process(JoinStr({signals, mssm_signals, {"ZTT", "TTT", "TTL", "VVT", "VVL"}}))
        .AddSyst(cb, "CMS_eff_t_highpT_500-inf_$ERA", "shape", SystMap<>::init(1.0));
  }
  cb.cp()
      .channel({"tt"})
      .process(JoinStr({signals, mssm_signals, {"ZTT", "TTT", "TTL", "VVT", "VVL"}}))
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
 
 cb.cp().process(signals_ggH).process({"ggX"}, false).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id>::init
  ({"em"}, {"2016"}, {35}, 0.995, 1.01)
  ({"em"}, {"2017"}, {35}, 0.993, 1.019)
  ({"em"}, {"2018"}, {35}, 0.993, 1.005)
  ({"em"}, {"2016"}, {36}, 0.993, 1.012)
  ({"em"}, {"2017"}, {36}, 0.986, 1.009)
  ({"em"}, {"2018"}, {36}, 0.996, 1.011)
  ({"em"}, {"2016"}, {37}, 0.991, 1.009)
  ({"em"}, {"2017"}, {37}, 0.986, 1.012)
  ({"em"}, {"2018"}, {37}, 0.996, 1.007)
  ({"em"}, {"2016"}, {32}, 0.999, 0.998)
  ({"em"}, {"2017"}, {32}, 1.0, 0.998)
  ({"em"}, {"2018"}, {32}, 1.0, 0.999)
  ({"em"}, {"2016"}, {33}, 1.0, 1.0)
  ({"em"}, {"2017"}, {33}, 1.0, 0.999)
  ({"em"}, {"2018"}, {33}, 1.0, 0.999)
  ({"em"}, {"2016"}, {34}, 1.0, 1.0)
  ({"em"}, {"2017"}, {34}, 1.001, 1.0)
  ({"em"}, {"2018"}, {34}, 1.001, 1.0)
  ({"em"}, {"2016"}, {2}, 1.0, 1.0)
  ({"em"}, {"2017"}, {2}, 1.0, 1.0)
  ({"em"}, {"2018"}, {2}, 1.0, 1.0)
  ({"et"}, {"2016"}, {35}, 0.977, 1.024)
  ({"et"}, {"2017"}, {35}, 0.967, 1.04)
  ({"et"}, {"2018"}, {35}, 0.982, 1.026)
  ({"et"}, {"2016"}, {36}, 0.977, 1.03)
  ({"et"}, {"2017"}, {36}, 0.97, 1.029)
  ({"et"}, {"2018"}, {36}, 0.979, 1.025)
  ({"et"}, {"2016"}, {32}, 1.001, 0.999)
  ({"et"}, {"2017"}, {32}, 1.001, 0.998)
  ({"et"}, {"2018"}, {32}, 1.001, 0.999)
  ({"et"}, {"2016"}, {33}, 1.001, 0.999)
  ({"et"}, {"2017"}, {33}, 1.002, 0.998)
  ({"et"}, {"2018"}, {33}, 1.001, 0.998)
  ({"mt"}, {"2016"}, {35}, 0.977, 1.023)
  ({"mt"}, {"2017"}, {35}, 0.967, 1.036)
  ({"mt"}, {"2018"}, {35}, 0.981, 1.022)
  ({"mt"}, {"2016"}, {36}, 0.977, 1.025)
  ({"mt"}, {"2017"}, {36}, 0.97, 1.035)
  ({"mt"}, {"2018"}, {36}, 0.982, 1.021)
  ({"mt"}, {"2016"}, {32}, 1.001, 0.999)
  ({"mt"}, {"2017"}, {32}, 1.001, 0.998)
  ({"mt"}, {"2018"}, {32}, 1.001, 0.999)
  ({"mt"}, {"2016"}, {33}, 1.001, 0.999)
  ({"mt"}, {"2017"}, {33}, 1.002, 0.998)
  ({"mt"}, {"2018"}, {33}, 1.001, 0.999)
  ({"tt"}, {"2016"}, {32}, 0.978, 1.025)
  ({"tt"}, {"2017"}, {32}, 0.965, 1.039)
  ({"tt"}, {"2018"}, {32}, 0.981, 1.024)
  ({"tt"}, {"2016"}, {35}, 1.001, 0.999)
  ({"tt"}, {"2017"}, {35}, 1.002, 0.998)
  ({"tt"}, {"2018"}, {35}, 1.001, 0.998)
  );

  cb.cp().process(signals_qqH).process({"qqX"}, false).AddSyst(cb, "CMS_htt_eff_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id>::init
  ({"em"}, {"2016"}, {35}, 0.996, 1.0)
  ({"em"}, {"2017"}, {35}, 0.988, 1.01)
  ({"em"}, {"2018"}, {35}, 0.993, 1.002)
  ({"em"}, {"2016"}, {36}, 0.997, 1.007)
  ({"em"}, {"2017"}, {36}, 0.994, 1.015)
  ({"em"}, {"2018"}, {36}, 0.994, 1.002)
  ({"em"}, {"2016"}, {37}, 0.979, 1.014)
  ({"em"}, {"2017"}, {37}, 0.986, 1.017)
  ({"em"}, {"2018"}, {37}, 1.002, 1.016)
  ({"em"}, {"2016"}, {32}, 1.0, 1.0)
  ({"em"}, {"2017"}, {32}, 1.0, 0.998)
  ({"em"}, {"2018"}, {32}, 1.0, 0.999)
  ({"em"}, {"2016"}, {33}, 0.999, 0.999)
  ({"em"}, {"2017"}, {33}, 1.0, 0.999)
  ({"em"}, {"2018"}, {33}, 1.001, 1.0)
  ({"em"}, {"2016"}, {34}, 1.003, 1.001)
  ({"em"}, {"2017"}, {34}, 1.002, 1.0)
  ({"em"}, {"2018"}, {34}, 1.0, 0.999)
  ({"em"}, {"2016"}, {2}, 1.0, 1.0)
  ({"em"}, {"2017"}, {2}, 1.0, 1.0)
  ({"em"}, {"2018"}, {2}, 1.0, 1.0)
  ({"et"}, {"2016"}, {35}, 0.974, 1.026)
  ({"et"}, {"2017"}, {35}, 0.953, 1.045)
  ({"et"}, {"2018"}, {35}, 0.977, 1.025)
  ({"et"}, {"2016"}, {36}, 0.985, 1.021)
  ({"et"}, {"2017"}, {36}, 0.958, 1.038)
  ({"et"}, {"2018"}, {36}, 0.975, 1.021)
  ({"et"}, {"2016"}, {32}, 1.002, 0.998)
  ({"et"}, {"2017"}, {32}, 1.005, 0.995)
  ({"et"}, {"2018"}, {32}, 1.003, 0.997)
  ({"et"}, {"2016"}, {33}, 1.001, 0.998)
  ({"et"}, {"2017"}, {33}, 1.005, 0.995)
  ({"et"}, {"2018"}, {33}, 1.003, 0.997)
  ({"mt"}, {"2016"}, {35}, 0.974, 1.026)
  ({"mt"}, {"2017"}, {35}, 0.954, 1.044)
  ({"mt"}, {"2018"}, {35}, 0.973, 1.028)
  ({"mt"}, {"2016"}, {36}, 0.982, 1.022)
  ({"mt"}, {"2017"}, {36}, 0.957, 1.033)
  ({"mt"}, {"2018"}, {36}, 0.979, 1.026)
  ({"mt"}, {"2016"}, {32}, 1.002, 0.998)
  ({"mt"}, {"2017"}, {32}, 1.005, 0.995)
  ({"mt"}, {"2018"}, {32}, 1.003, 0.997)
  ({"mt"}, {"2016"}, {33}, 1.002, 0.998)
  ({"mt"}, {"2017"}, {33}, 1.005, 0.996)
  ({"mt"}, {"2018"}, {33}, 1.003, 0.997)
  ({"tt"}, {"2016"}, {32}, 0.977, 1.025)
  ({"tt"}, {"2017"}, {32}, 0.956, 1.042)
  ({"tt"}, {"2018"}, {32}, 0.976, 1.026)
  ({"tt"}, {"2016"}, {35}, 1.002, 0.998)
  ({"tt"}, {"2017"}, {35}, 1.005, 0.995)
  ({"tt"}, {"2018"}, {35}, 1.003, 0.997)
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

  cb.cp().process(signals_ggH).process({"ggX"}, false).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id>::init
  ({"em"}, {"2016"}, {35}, 0.905, 1.092)
  ({"em"}, {"2017"}, {35}, 0.892, 1.106)
  ({"em"}, {"2018"}, {35}, 0.867, 1.113)
  ({"em"}, {"2016"}, {36}, 0.914, 1.112)
  ({"em"}, {"2017"}, {36}, 0.874, 1.114)
  ({"em"}, {"2018"}, {36}, 0.858, 1.152)
  ({"em"}, {"2016"}, {37}, 0.91, 1.113)
  ({"em"}, {"2017"}, {37}, 0.882, 1.116)
  ({"em"}, {"2018"}, {37}, 0.849, 1.128)
  ({"em"}, {"2016"}, {32}, 1.006, 0.992)
  ({"em"}, {"2017"}, {32}, 1.008, 0.991)
  ({"em"}, {"2018"}, {32}, 1.014, 0.987)
  ({"em"}, {"2016"}, {33}, 1.002, 0.997)
  ({"em"}, {"2017"}, {33}, 1.004, 0.996)
  ({"em"}, {"2018"}, {33}, 1.006, 0.993)
  ({"em"}, {"2016"}, {34}, 1.003, 0.997)
  ({"em"}, {"2017"}, {34}, 1.005, 0.996)
  ({"em"}, {"2018"}, {34}, 1.008, 0.994)
  ({"em"}, {"2016"}, {2}, 1.0, 1.0)
  ({"em"}, {"2017"}, {2}, 1.0, 1.0)
  ({"em"}, {"2018"}, {2}, 1.0, 1.0)
  ({"et"}, {"2016"}, {35}, 0.939, 1.08)
  ({"et"}, {"2017"}, {35}, 0.921, 1.074)
  ({"et"}, {"2018"}, {35}, 0.882, 1.118)
  ({"et"}, {"2016"}, {36}, 0.937, 1.062)
  ({"et"}, {"2017"}, {36}, 0.926, 1.082)
  ({"et"}, {"2018"}, {36}, 0.894, 1.123)
  ({"et"}, {"2016"}, {32}, 1.002, 0.997)
  ({"et"}, {"2017"}, {32}, 1.004, 0.997)
  ({"et"}, {"2018"}, {32}, 1.007, 0.993)
  ({"et"}, {"2016"}, {33}, 1.003, 0.997)
  ({"et"}, {"2017"}, {33}, 1.004, 0.995)
  ({"et"}, {"2018"}, {33}, 1.007, 0.992)
  ({"mt"}, {"2016"}, {35}, 0.932, 1.08)
  ({"mt"}, {"2017"}, {35}, 0.917, 1.08)
  ({"mt"}, {"2018"}, {35}, 0.878, 1.12)
  ({"mt"}, {"2016"}, {36}, 0.94, 1.074)
  ({"mt"}, {"2017"}, {36}, 0.926, 1.078)
  ({"mt"}, {"2018"}, {36}, 0.886, 1.112)
  ({"mt"}, {"2016"}, {32}, 1.002, 0.997)
  ({"mt"}, {"2017"}, {32}, 1.004, 0.997)
  ({"mt"}, {"2018"}, {32}, 1.006, 0.994)
  ({"mt"}, {"2016"}, {33}, 1.003, 0.997)
  ({"mt"}, {"2017"}, {33}, 1.004, 0.996)
  ({"mt"}, {"2018"}, {33}, 1.008, 0.992)
  ({"tt"}, {"2016"}, {32}, 0.945, 1.074)
  ({"tt"}, {"2017"}, {32}, 0.931, 1.073)
  ({"tt"}, {"2018"}, {32}, 0.892, 1.1)
  ({"tt"}, {"2016"}, {35}, 1.003, 0.996)
  ({"tt"}, {"2017"}, {35}, 1.004, 0.996)
  ({"tt"}, {"2018"}, {35}, 1.007, 0.993)
  );

  cb.cp().process(signals_qqH).process({"qqX"}, false).AddSyst(cb, "CMS_htt_mistag_b_$ERA", "lnN", SystMapAsymm<channel,ch::syst::era,bin_id>::init
  ({"em"}, {"2016"}, {35}, 0.939, 1.09)
  ({"em"}, {"2017"}, {35}, 0.858, 1.103)
  ({"em"}, {"2018"}, {35}, 0.898, 1.125)
  ({"em"}, {"2016"}, {36}, 0.927, 1.083)
  ({"em"}, {"2017"}, {36}, 0.905, 1.117)
  ({"em"}, {"2018"}, {36}, 0.89, 1.103)
  ({"em"}, {"2016"}, {37}, 0.929, 1.068)
  ({"em"}, {"2017"}, {37}, 0.892, 1.101)
  ({"em"}, {"2018"}, {37}, 0.886, 1.119)
  ({"em"}, {"2016"}, {32}, 1.004, 0.993)
  ({"em"}, {"2017"}, {32}, 1.01, 0.991)
  ({"em"}, {"2018"}, {32}, 1.008, 0.988)
  ({"em"}, {"2016"}, {33}, 1.004, 0.994)
  ({"em"}, {"2017"}, {33}, 1.006, 0.992)
  ({"em"}, {"2018"}, {33}, 1.009, 0.992)
  ({"em"}, {"2016"}, {34}, 1.007, 0.997)
  ({"em"}, {"2017"}, {34}, 1.008, 0.994)
  ({"em"}, {"2018"}, {34}, 1.01, 0.99)
  ({"em"}, {"2016"}, {2}, 1.0, 1.0)
  ({"em"}, {"2017"}, {2}, 1.0, 1.0)
  ({"em"}, {"2018"}, {2}, 1.0, 1.0)
  ({"et"}, {"2016"}, {35}, 0.958, 1.047)
  ({"et"}, {"2017"}, {35}, 0.954, 1.046)
  ({"et"}, {"2018"}, {35}, 0.922, 1.073)
  ({"et"}, {"2016"}, {36}, 0.934, 1.038)
  ({"et"}, {"2017"}, {36}, 0.949, 1.04)
  ({"et"}, {"2018"}, {36}, 0.93, 1.08)
  ({"et"}, {"2016"}, {32}, 1.004, 0.996)
  ({"et"}, {"2017"}, {32}, 1.005, 0.995)
  ({"et"}, {"2018"}, {32}, 1.009, 0.991)
  ({"et"}, {"2016"}, {33}, 1.006, 0.996)
  ({"et"}, {"2017"}, {33}, 1.006, 0.995)
  ({"et"}, {"2018"}, {33}, 1.009, 0.99)
  ({"mt"}, {"2016"}, {35}, 0.96, 1.047)
  ({"mt"}, {"2017"}, {35}, 0.952, 1.047)
  ({"mt"}, {"2018"}, {35}, 0.926, 1.072)
  ({"mt"}, {"2016"}, {36}, 0.944, 1.06)
  ({"mt"}, {"2017"}, {36}, 0.953, 1.046)
  ({"mt"}, {"2018"}, {36}, 0.929, 1.074)
  ({"mt"}, {"2016"}, {32}, 1.004, 0.996)
  ({"mt"}, {"2017"}, {32}, 1.005, 0.995)
  ({"mt"}, {"2018"}, {32}, 1.009, 0.991)
  ({"mt"}, {"2016"}, {33}, 1.005, 0.995)
  ({"mt"}, {"2017"}, {33}, 1.005, 0.995)
  ({"mt"}, {"2018"}, {33}, 1.009, 0.991)
  ({"tt"}, {"2016"}, {32}, 0.958, 1.049)
  ({"tt"}, {"2017"}, {32}, 0.956, 1.046)
  ({"tt"}, {"2018"}, {32}, 0.929, 1.072)
  ({"tt"}, {"2016"}, {35}, 1.004, 0.995)
  ({"tt"}, {"2017"}, {35}, 1.005, 0.995)
  ({"tt"}, {"2018"}, {35}, 1.009, 0.991)
  );

  // Do signal b efficiency and mistag uncertainties. Load from json to solve memory issues
  if (not sm){
    Json::Value js;
    string json_file = std::getenv("CMSSW_BASE") + std::string("/src/CombineHarvester/MSSMvsSMRun2Legacy/input/signal_b_efficiencies.json");
    js = ch::ExtractJsonFromFile(json_file);
    std::string systematic_names[2] = {"CMS_htt_eff_b_$ERA","CMS_htt_mistag_b_$ERA"};
    std::string channels[4] = {"em","et","mt","tt"};
    std::string years[3] = {"2016","2017","2018"};
    std::string masses[31] = {"60","80","100","120","125","130","140","160","180","200","250","300","350","400","450","500","600","700","800","900","1000","1200","1400","1600","1800","2000","2300","2600","2900","3200","3500"};
    for (auto s: systematic_names) {
        for (auto c: channels) {
        for (auto yr: years) {
            for (auto m: masses) {
            double bbH_btag_up = js[s]["mssm_bbH_signals"][m]["btag_catagories"][c][yr]["Up"].asDouble();
            double bbH_btag_down = js[s]["mssm_bbH_signals"][m]["btag_catagories"][c][yr]["Down"].asDouble();
            cb.cp().process(mssm_bbH_signals).AddSyst(cb, s, "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init({c}, {yr}, btag_catagories,{m}, bbH_btag_down, bbH_btag_up));
            // Cover SM bbH125 case
            if(m == "125"){
              cb.cp().process({"bbH125"}).AddSyst(cb, s, "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init({c}, {yr}, btag_catagories,{""}, bbH_btag_down, bbH_btag_up));
            }

            double ggH_btag_up = js[s]["mssm_ggH_signals"][m]["btag_catagories"][c][yr]["Up"].asDouble();
            double ggH_btag_down = js[s]["mssm_ggH_signals"][m]["btag_catagories"][c][yr]["Down"].asDouble();
            cb.cp().process(mssm_ggH_signals).AddSyst(cb, s, "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
            ({c}, {yr}, btag_catagories, {m}, ggH_btag_down, ggH_btag_up));

            double bbH_nobtag_up = js[s]["mssm_bbH_signals"][m]["nobtag_catagories"][c][yr]["Up"].asDouble();
            double bbH_nobtag_down = js[s]["mssm_bbH_signals"][m]["nobtag_catagories"][c][yr]["Down"].asDouble();
            cb.cp().process(mssm_bbH_signals).AddSyst(cb, s, "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
            ({c}, {yr}, nobtag_catagories, {m}, bbH_nobtag_down, bbH_nobtag_up));
            // Cover SM bbH125 case
            if(m == "125"){
              cb.cp().process({"bbH125"}).AddSyst(cb, s, "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
              ({c}, {yr}, nobtag_catagories, {""}, bbH_nobtag_down, bbH_nobtag_up));
            }

            double ggH_nobtag_up = js[s]["mssm_ggH_signals"][m]["nobtag_catagories"][c][yr]["Up"].asDouble();
            double ggH_nobtag_down = js[s]["mssm_ggH_signals"][m]["nobtag_catagories"][c][yr]["Down"].asDouble();
            cb.cp().process(mssm_ggH_signals).AddSyst(cb, s, "lnN", SystMapAsymm<channel,ch::syst::era,bin_id,mass>::init
            ({c}, {yr}, nobtag_catagories, {m}, ggH_nobtag_down, ggH_nobtag_up));
            }
        }
        }
    }
  }

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

  // small uncertainty decorrelated by channel to account for statistical uncertainties on corrections, enlarged to cover differences observed between corrections for et and mt channels 
  cb.cp()
      .process({"EMB"})
      .channel({"et", "mt", "tt"})
      .bin_id(mssm_categories)
      .AddSyst(cb, "scale_embed_met_$CHANNEL_$ERA", "shape", SystMap<>::init(0.25)); 

  // the other component of the uncertainty is systematic and correlated between channels (but decorrelated by era) 

  cb.cp()
      .process({"EMB"})
      .bin_id(mssm_categories)
      .channel({"tt","mt","et"})
      .era({"2016"})
      .AddSyst(cb, "scale_embed_met_$ERA", "shape", SystMap<>::init(0.36));
  cb.cp()
      .process({"EMB"})
      .bin_id(mssm_categories)
      .channel({"tt","mt","et"})
      .era({"2017"})
      .AddSyst(cb, "scale_embed_met_$ERA", "shape", SystMap<>::init(0.64));
  cb.cp()
      .process({"EMB"})
      .bin_id(mssm_categories)
      .channel({"tt","mt","et"})
      .era({"2018"})
      .AddSyst(cb, "scale_embed_met_$ERA", "shape", SystMap<>::init(0.14));


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
  if (sm){
        // TT
    cb.cp()
        .channel({"et", "mt", "tt", "em"})
        .process({"TTT", "TTL", "TTJ", "TT"})
        .AddSyst(cb, "CMS_htt_tjXsec", "lnN", SystMap<>::init(1.06));
  }
  else{
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
  }
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
    cb.cp()
       .channel({"em"})
       .process({"QCD"})
       .AddSyst(cb, "subtrMC", "shape", SystMap<>::init(1.00));
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
  //

  // Uncertainty on branching ratio for HTT at 125 GeV
  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process(JoinStr({signals, {"bbH125"}})).process({"qqX", "ggX"}, false)
      .AddSyst(cb, "BR_htt_THU", "lnN", SystMap<>::init(1.017));
  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process(JoinStr({signals, {"bbH125"}})).process({"qqX", "ggX"}, false)
      .AddSyst(cb, "BR_htt_PU_mq", "lnN", SystMap<>::init(1.0099));
  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process(JoinStr({signals, {"bbH125"}})).process({"qqX", "ggX"}, false)
      .AddSyst(cb, "BR_htt_PU_alphas", "lnN", SystMap<>::init(1.0062));

  // 95 GeV samples BR uncertainties 
  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process({"qqX", "ggX"})
      .AddSyst(cb, "BR_htt", "lnN", SystMap<>::init(1.0172));

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

  // QCD scale (no ggH & qqH signals for tautau decay channel) for 125 GeV Higgs
  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process({"ZH125", "Zh", "ZH", "ZH1", "ZHWW125"})
      .AddSyst(cb, "QCDScale_VH", "lnN", SystMap<>::init(1.009));
  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process({"WH125", "Wh", "WH", "WH1", "WHWW125"})
      .AddSyst(cb, "QCDScale_VH", "lnN", SystMap<>::init(1.008));
  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process({"ttH125"})
      .AddSyst(cb, "QCDScale_ttH", "lnN", SystMap<>::init(1.08));
  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process({"bbH125"})
      .AddSyst(cb, "QCDScale_bbH", "lnN", SystMap<>::init(1.22));
  cb.cp()
    .channel({"et", "mt", "tt", "em"})
    .process({signals_ggHToWW})
    .AddSyst(cb, "QCDScale_ggH", "lnN", SystMap<>::init(1.039));
  cb.cp()
    .channel({"et", "mt", "tt", "em"})
    .process({signals_qqHToWW})
    .AddSyst(cb, "QCDScale_qqH", "lnN", SystMap<>::init(1.005));

  // this is the 95 GeV sample
  cb.cp()
    .channel({"et", "mt", "tt", "em"})
    .process({"qqX"})
    .AddSyst(cb, "QCDScale_qqH", "lnN", SystMap<>::init(1.004));

  // PDF for 125 GeV Higgs
  cb.cp()
    .channel({"et", "mt", "tt", "em"})
    .process(JoinStr({signals_ggH,signals_ggHToWW})).process({"ggX"}, false)
    .AddSyst(cb, "pdf_Higgs_gg", "lnN", SystMap<>::init(1.032));
  cb.cp()
    .channel({"et", "mt", "tt", "em"})
    .process(JoinStr({signals_qqH,signals_qqHToWW})).process({"qqX"}, false)
    .AddSyst(cb, "pdf_Higgs_qqbar", "lnN", SystMap<>::init(1.021));
  cb.cp()
    .channel({"et", "mt", "tt", "em"})
    .process({"ZH125", "Zh", "ZH", "ZH1", "ZHWW125"})
    .AddSyst(cb, "pdf_Higgs_VH", "lnN", SystMap<>::init(1.013));
  cb.cp()
    .channel({"et", "mt", "tt", "em"})
    .process({"WH125", "Wh", "WH", "WH1", "WHWW125"})
    .AddSyst(cb, "pdf_Higgs_VH", "lnN", SystMap<>::init(1.018));
  cb.cp()
    .channel({"et", "mt", "tt", "em"})
    .process({"ttH125"})
    .AddSyst(cb, "pdf_Higgs_ttH", "lnN", SystMap<>::init(1.036));

  // No pdf uncertainty for bbH125 available from LHCHWG file

  // 95 GeV sample (use same as 125 GeV for now but update later)
  cb.cp()
    .channel({"et", "mt", "tt", "em"})
    .process({"qqX"})
    .AddSyst(cb, "pdf_Higgs_qqbar", "lnN", SystMap<>::init(1.021));

  if (sm)
  {
    // Gluon-fusion WG1 uncertainty scheme for 125 GeV Higgs
    cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process(signals_ggH).process({"ggX"}, false)
      .AddSyst(cb, "THU_ggH_Mig01", "shape", SystMap<>::init(1.00));
    cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process(signals_ggH).process({"ggX"}, false)
      .AddSyst(cb, "THU_ggH_Mig12", "shape", SystMap<>::init(1.00));
    cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process(signals_ggH).process({"ggX"}, false)
      .AddSyst(cb, "THU_ggH_Mu", "shape", SystMap<>::init(1.00));
    cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process(signals_ggH).process({"ggX"}, false)
      .AddSyst(cb, "THU_ggH_PT120", "shape", SystMap<>::init(1.00));
    cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process(signals_ggH).process({"ggX"}, false)
      .AddSyst(cb, "THU_ggH_PT60", "shape", SystMap<>::init(1.00));
    cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process(signals_ggH).process({"ggX"}, false)
      .AddSyst(cb, "THU_ggH_Res", "shape", SystMap<>::init(1.00));
    cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process(signals_ggH).process({"ggX"}, false)
      .AddSyst(cb, "THU_ggH_VBF2j", "shape", SystMap<>::init(1.00));
    cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process(signals_ggH).process({"ggX"}, false)
      .AddSyst(cb, "THU_ggH_VBF3j", "shape", SystMap<>::init(1.00));
    cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process(signals_ggH).process({"ggX"}, false)
      .AddSyst(cb, "THU_ggH_qmtop", "shape", SystMap<>::init(1.00));
    // VBF WG1 uncertainty scheme for 125 GeV Higgs
   cb.cp()
     .channel({"et", "mt", "tt", "em"})
     .process(signals_qqH).process({"qqX"}, false)
     .AddSyst(cb, "THU_qqH_TOT", "shape", SystMap<>::init(1.00));
   cb.cp()
     .channel({"et", "mt", "tt", "em"})
     .process(signals_qqH).process({"qqX"}, false)
     .AddSyst(cb, "THU_qqH_PTH200", "shape", SystMap<>::init(1.00));
   cb.cp()
     .channel({"et", "mt", "tt", "em"})
     .process(signals_qqH).process({"qqX"}, false)
     .AddSyst(cb, "THU_qqH_Mjj60", "shape", SystMap<>::init(1.00));
   cb.cp()
     .channel({"et", "mt", "tt", "em"})
     .process(signals_qqH).process({"qqX"}, false)
     .AddSyst(cb, "THU_qqH_Mjj120", "shape", SystMap<>::init(1.00));
   cb.cp()
     .channel({"et", "mt", "tt", "em"})
     .process(signals_qqH).process({"qqX"}, false)
     .AddSyst(cb, "THU_qqH_Mjj350", "shape", SystMap<>::init(1.00));
   cb.cp()
     .channel({"et", "mt", "tt", "em"})
     .process(signals_qqH).process({"qqX"}, false)
     .AddSyst(cb, "THU_qqH_Mjj700", "shape", SystMap<>::init(1.00));
   cb.cp()
     .channel({"et", "mt", "tt", "em"})
     .process(signals_qqH).process({"qqX"}, false)
     .AddSyst(cb, "THU_qqH_Mjj1000", "shape", SystMap<>::init(1.00));
   cb.cp()
     .channel({"et", "mt", "tt", "em"})
     .process(signals_qqH).process({"qqX"}, false)
     .AddSyst(cb, "THU_qqH_Mjj1500", "shape", SystMap<>::init(1.00));
   cb.cp()
     .channel({"et", "mt", "tt", "em"})
     .process(signals_qqH).process({"qqX"}, false)
     .AddSyst(cb, "THU_qqH_25", "shape", SystMap<>::init(1.00));
   cb.cp()
     .channel({"et", "mt", "tt", "em"})
     .process(signals_qqH).process({"qqX"}, false)
     .AddSyst(cb, "THU_qqH_JET01", "shape", SystMap<>::init(1.00));
  }
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
