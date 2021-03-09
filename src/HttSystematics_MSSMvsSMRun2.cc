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

  std::vector<std::string> mssm_ggH_t_signals = {"ggh_t", "ggH_t", "ggA_t"};
  std::vector<std::string> mssm_ggH_b_signals = {"ggh_b", "ggH_b", "ggA_b"};
  std::vector<std::string> mssm_ggH_i_signals = {"ggh_i", "ggH_i", "ggA_i"};
  std::vector<std::string> mssm_ggH_signals = JoinStr({mssm_ggH_t_signals, mssm_ggH_b_signals, mssm_ggH_i_signals});
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
  cb.cp().process(mssm_bbH_signals).AddSyst(cb, "QCDScale_QshScale_bbH","lnN", SystMap<channel,bin_id,mass>::init
     ({"em","et","mt","tt"}, nobtag_catagories, {""},   1.04)
     ({"em","et","mt","tt"},   btag_catagories, {""},   0.96)
     ({"em","et","mt","tt"}, nobtag_catagories, {"80"},   1.034)
     ({"em","et","mt","tt"},   btag_catagories, {"80"},   0.827)
     ({"em","et","mt","tt"}, nobtag_catagories, {"90"},   1.036)
     ({"em","et","mt","tt"},   btag_catagories, {"90"},   0.835)
     ({"em","et","mt","tt"}, nobtag_catagories, {"100"},  1.036)
     ({"em","et","mt","tt"},   btag_catagories, {"100"},  0.847)
     ({"em","et","mt","tt"}, nobtag_catagories, {"110"},  1.038)
     ({"em","et","mt","tt"},   btag_catagories, {"110"},  0.853)
     ({"em","et","mt","tt"}, nobtag_catagories, {"120"},  1.038)
     ({"em","et","mt","tt"},   btag_catagories, {"120"},  0.862)
     ({"em","et","mt","tt"}, nobtag_catagories, {"125"},  1.038)
     ({"em","et","mt","tt"},   btag_catagories, {"125"},  0.862)
     ({"em","et","mt","tt"}, nobtag_catagories, {"130"},  1.04 )
     ({"em","et","mt","tt"},   btag_catagories, {"130"},  0.867)
     ({"em","et","mt","tt"}, nobtag_catagories, {"140"},  1.04 )
     ({"em","et","mt","tt"},   btag_catagories, {"140"},  0.872)
     ({"em","et","mt","tt"}, nobtag_catagories, {"160"},  1.038)
     ({"em","et","mt","tt"},   btag_catagories, {"160"},  0.887)
     ({"em","et","mt","tt"}, nobtag_catagories, {"180"},  1.038)
     ({"em","et","mt","tt"},   btag_catagories, {"180"},  0.897)
     ({"em","et","mt","tt"}, nobtag_catagories, {"200"},  1.038)
     ({"em","et","mt","tt"},   btag_catagories, {"200"},  0.902)
     ({"em","et","mt","tt"}, nobtag_catagories, {"250"},  1.035)
     ({"em","et","mt","tt"},   btag_catagories, {"250"},  0.922)
     ({"em","et","mt","tt"}, nobtag_catagories, {"300"},  1.035)
     ({"em","et","mt","tt"},   btag_catagories, {"300"},  0.922)
     ({"em","et","mt","tt"}, nobtag_catagories, {"350"},  1.033)
     ({"em","et","mt","tt"},   btag_catagories, {"350"},  0.939)
     ({"em","et","mt","tt"}, nobtag_catagories, {"400"},  1.036)
     ({"em","et","mt","tt"},   btag_catagories, {"400"},  0.934)
     ({"em","et","mt","tt"}, nobtag_catagories, {"450"},  1.035)
     ({"em","et","mt","tt"},   btag_catagories, {"450"},  0.94 )
     ({"em","et","mt","tt"}, nobtag_catagories, {"500"},  1.032)
     ({"em","et","mt","tt"},   btag_catagories, {"500"},  0.95 )
     ({"em","et","mt","tt"}, nobtag_catagories, {"600"},  1.034)
     ({"em","et","mt","tt"},   btag_catagories, {"600"},  0.948)
     ({"em","et","mt","tt"}, nobtag_catagories, {"700"},  1.034)
     ({"em","et","mt","tt"},   btag_catagories, {"700"},  0.952)
     ({"em","et","mt","tt"}, nobtag_catagories, {"800"},  1.036)
     ({"em","et","mt","tt"},   btag_catagories, {"800"},  0.948)
     ({"em","et","mt","tt"}, nobtag_catagories, {"900"},  1.036)
     ({"em","et","mt","tt"},   btag_catagories, {"900"},  0.95 )
     ({"em","et","mt","tt"}, nobtag_catagories, {"1000"}, 1.037)
     ({"em","et","mt","tt"},   btag_catagories, {"1000"}, 0.949)
     ({"em","et","mt","tt"}, nobtag_catagories, {"1200"}, 1.037)
     ({"em","et","mt","tt"},   btag_catagories, {"1200"}, 0.951)
     ({"em","et","mt","tt"}, nobtag_catagories, {"1400"}, 1.034)
     ({"em","et","mt","tt"},   btag_catagories, {"1400"}, 0.957)
     ({"em","et","mt","tt"}, nobtag_catagories, {"1600"}, 1.041)
     ({"em","et","mt","tt"},   btag_catagories, {"1600"}, 0.943)
     ({"em","et","mt","tt"}, nobtag_catagories, {"1800"}, 1.037)
     ({"em","et","mt","tt"},   btag_catagories, {"1800"}, 0.952)
     ({"em","et","mt","tt"}, nobtag_catagories, {"2000"}, 1.035)
     ({"em","et","mt","tt"},   btag_catagories, {"2000"}, 0.956)
     ({"em","et","mt","tt"}, nobtag_catagories, {"2300"}, 1.035)
     ({"em","et","mt","tt"},   btag_catagories, {"2300"}, 0.956)
     ({"em","et","mt","tt"}, nobtag_catagories, {"2600"}, 1.039)
     ({"em","et","mt","tt"},   btag_catagories, {"2600"}, 0.95 )
     ({"em","et","mt","tt"}, nobtag_catagories, {"2900"}, 1.036)
     ({"em","et","mt","tt"},   btag_catagories, {"2900"}, 0.954)
     ({"em","et","mt","tt"}, nobtag_catagories, {"3200"}, 1.034)
     ({"em","et","mt","tt"},   btag_catagories, {"3200"}, 0.957));

  // // ##########################################################################
  // // Uncertainty: b tagging acceptance uncertainties for pdf and scale and hdamp variations.
  // // References:
  // // - "Talk in MSSM HTT Meeting by Danny"
  // //   (https://indico.cern.ch/event/990440/contributions/4167707/attachments/2167087/3659678/MSSM_signal.pdf)
  // //   (file source: https://cernbox.cern.ch/index.php/s/9cgjdpaTeqYEaFZ
  // // Notes:
  // // - FIXME: Uncomment for new powheg bbH samples
  // // ##########################################################################
  // cb.cp().process(mssm_bbH_signals).AddSyst(cb, "pdf_bbH_ACCEPT", "lnN", SystMap<era,bin_id,mass>::init
  //    ({"2016"}, nobtag_catagories, {"60"}, 0.996)
  //    ({"2016"},   btag_catagories, {"60"}, 1.023)
  //    ({"2017"}, nobtag_catagories, {"60"}, 0.996)
  //    ({"2017"},   btag_catagories, {"60"}, 1.023)
  //    ({"2018"}, nobtag_catagories, {"60"}, 0.996)
  //    ({"2018"},   btag_catagories, {"60"}, 1.022)
  //    ({"2016"}, nobtag_catagories, {"100"}, 0.996)
  //    ({"2016"},   btag_catagories, {"100"}, 1.017)
  //    ({"2017"}, nobtag_catagories, {"100"}, 0.996)
  //    ({"2017"},   btag_catagories, {"100"}, 1.016)
  //    ({"2018"}, nobtag_catagories, {"100"}, 0.996)
  //    ({"2018"},   btag_catagories, {"100"}, 1.015)
  //    ({"2016"}, nobtag_catagories, {"125"}, 0.996)
  //    ({"2016"},   btag_catagories, {"125"}, 1.015)
  //    ({"2017"}, nobtag_catagories, {"125"}, 0.995)
  //    ({"2017"},   btag_catagories, {"125"}, 1.016)
  //    ({"2018"}, nobtag_catagories, {"125"}, 0.995)
  //    ({"2018"},   btag_catagories, {"125"}, 1.016)
  //    ({"2016"}, nobtag_catagories, {"140"}, 0.995)
  //    ({"2016"},   btag_catagories, {"140"}, 1.016)
  //    ({"2017"}, nobtag_catagories, {"140"}, 0.995)
  //    ({"2017"},   btag_catagories, {"140"}, 1.015)
  //    ({"2018"}, nobtag_catagories, {"140"}, 0.995)
  //    ({"2018"},   btag_catagories, {"140"}, 1.014)
  //    ({"2016"}, nobtag_catagories, {"160"}, 0.995)
  //    ({"2016"},   btag_catagories, {"160"}, 1.017)
  //    ({"2017"}, nobtag_catagories, {"160"}, 0.994)
  //    ({"2017"},   btag_catagories, {"160"}, 1.016)
  //    ({"2018"}, nobtag_catagories, {"160"}, 0.994)
  //    ({"2018"},   btag_catagories, {"160"}, 1.016)
  //    ({"2016"}, nobtag_catagories, {"400"}, 0.993)
  //    ({"2016"},   btag_catagories, {"400"}, 1.012)
  //    ({"2017"}, nobtag_catagories, {"400"}, 0.993)
  //    ({"2017"},   btag_catagories, {"400"}, 1.011)
  //    ({"2018"}, nobtag_catagories, {"400"}, 0.993)
  //    ({"2018"},   btag_catagories, {"400"}, 1.010)
  //    ({"2016"}, nobtag_catagories, {"600"}, 0.993)
  //    ({"2016"},   btag_catagories, {"600"}, 1.010)
  //    ({"2017"}, nobtag_catagories, {"600"}, 0.993)
  //    ({"2017"},   btag_catagories, {"600"}, 1.010)
  //    ({"2018"}, nobtag_catagories, {"600"}, 0.993)
  //    ({"2018"},   btag_catagories, {"600"}, 1.009)
  //    ({"2016"}, nobtag_catagories, {"1000"}, 0.987)
  //    ({"2016"},   btag_catagories, {"1000"}, 1.016)
  //    ({"2017"}, nobtag_catagories, {"1000"}, 0.987)
  //    ({"2017"},   btag_catagories, {"1000"}, 1.015)
  //    ({"2018"}, nobtag_catagories, {"1000"}, 0.986)
  //    ({"2018"},   btag_catagories, {"1000"}, 1.015)
  //    ({"2016"}, nobtag_catagories, {"1400"}, 0.985)
  //    ({"2016"},   btag_catagories, {"1400"}, 1.018)
  //    ({"2017"}, nobtag_catagories, {"1400"}, 0.985)
  //    ({"2017"},   btag_catagories, {"1400"}, 1.016)
  //    ({"2018"}, nobtag_catagories, {"1400"}, 0.985)
  //    ({"2018"},   btag_catagories, {"1400"}, 1.016)
  //    ({"2016"}, nobtag_catagories, {"2000"}, 0.983)
  //    ({"2016"},   btag_catagories, {"2000"}, 1.019)
  //    ({"2017"}, nobtag_catagories, {"2000"}, 0.983)
  //    ({"2017"},   btag_catagories, {"2000"}, 1.016)
  //    ({"2018"}, nobtag_catagories, {"2000"}, 0.983)
  //    ({"2018"},   btag_catagories, {"2000"}, 1.016)
  //    ({"2016"}, nobtag_catagories, {"2300"}, 0.982)
  //    ({"2016"},   btag_catagories, {"2300"}, 1.020)
  //    ({"2017"}, nobtag_catagories, {"2300"}, 0.981)
  //    ({"2017"},   btag_catagories, {"2300"}, 1.018)
  //    ({"2018"}, nobtag_catagories, {"2300"}, 0.981)
  //    ({"2018"},   btag_catagories, {"2300"}, 1.018)
  //    ({"2016"}, nobtag_catagories, {"2900"}, 0.979)
  //    ({"2016"},   btag_catagories, {"2900"}, 1.023)
  //    ({"2017"}, nobtag_catagories, {"2900"}, 0.979)
  //    ({"2017"},   btag_catagories, {"2900"}, 1.020)
  //    ({"2018"}, nobtag_catagories, {"2900"}, 0.978)
  //    ({"2018"},   btag_catagories, {"2900"}, 1.020)
  //    ({"2016"}, nobtag_catagories, {"3500"}, 0.974)
  //    ({"2016"},   btag_catagories, {"3500"}, 1.027)
  //    ({"2017"}, nobtag_catagories, {"3500"}, 0.974)
  //    ({"2017"},   btag_catagories, {"3500"}, 1.024)
  //    ({"2018"}, nobtag_catagories, {"3500"}, 0.973)
  //    ({"2018"},   btag_catagories, {"3500"}, 1.024)
  // );
  // cb.cp().process(mssm_bbH_signals).AddSyst(cb, "QCDscaleAndHdamp_bbH_ACCEPT", "lnN", SystMapAsymm<era,bin_id,mass>::init
  //    ({"2016"}, nobtag_catagories, {"60"}, 1.010, 0.993)
  //    ({"2016"},   btag_catagories, {"60"}, 0.938, 1.042)
  //    ({"2017"}, nobtag_catagories, {"60"}, 1.010, 0.992)
  //    ({"2017"},   btag_catagories, {"60"}, 0.944, 1.043)
  //    ({"2018"}, nobtag_catagories, {"60"}, 1.011, 0.992)
  //    ({"2018"},   btag_catagories, {"60"}, 0.941, 1.042)
  //    ({"2016"}, nobtag_catagories, {"100"}, 1.016, 0.992)
  //    ({"2016"},   btag_catagories, {"100"}, 0.926, 1.035)
  //    ({"2017"}, nobtag_catagories, {"100"}, 1.020, 0.993)
  //    ({"2017"},   btag_catagories, {"100"}, 0.920, 1.029)
  //    ({"2018"}, nobtag_catagories, {"100"}, 1.021, 0.993)
  //    ({"2018"},   btag_catagories, {"100"}, 0.920, 1.027)
  //    ({"2016"}, nobtag_catagories, {"125"}, 1.007, 0.993)
  //    ({"2016"},   btag_catagories, {"125"}, 0.974, 1.026)
  //    ({"2017"}, nobtag_catagories, {"125"}, 1.008, 0.992)
  //    ({"2017"},   btag_catagories, {"125"}, 0.972, 1.025)
  //    ({"2018"}, nobtag_catagories, {"125"}, 1.009, 0.992)
  //    ({"2018"},   btag_catagories, {"125"}, 0.972, 1.026)
  //    ({"2016"}, nobtag_catagories, {"140"}, 1.012, 0.995)
  //    ({"2016"},   btag_catagories, {"140"}, 0.958, 1.019)
  //    ({"2017"}, nobtag_catagories, {"140"}, 1.012, 0.994)
  //    ({"2017"},   btag_catagories, {"140"}, 0.962, 1.018)
  //    ({"2018"}, nobtag_catagories, {"140"}, 1.013, 0.994)
  //    ({"2018"},   btag_catagories, {"140"}, 0.962, 1.017)
  //    ({"2016"}, nobtag_catagories, {"160"}, 1.013, 0.994)
  //    ({"2016"},   btag_catagories, {"160"}, 0.958, 1.021)
  //    ({"2017"}, nobtag_catagories, {"160"}, 1.014, 0.993)
  //    ({"2017"},   btag_catagories, {"160"}, 0.960, 1.019)
  //    ({"2018"}, nobtag_catagories, {"160"}, 1.015, 0.993)
  //    ({"2018"},   btag_catagories, {"160"}, 0.961, 1.019)
  //    ({"2016"}, nobtag_catagories, {"400"}, 1.014, 0.993)
  //    ({"2016"},   btag_catagories, {"400"}, 0.974, 1.014)
  //    ({"2017"}, nobtag_catagories, {"400"}, 1.014, 0.993)
  //    ({"2017"},   btag_catagories, {"400"}, 0.978, 1.012)
  //    ({"2018"}, nobtag_catagories, {"400"}, 1.014, 0.992)
  //    ({"2018"},   btag_catagories, {"400"}, 0.979, 1.012)
  //    ({"2016"}, nobtag_catagories, {"600"}, 1.012, 0.983)
  //    ({"2016"},   btag_catagories, {"600"}, 0.983, 1.025)
  //    ({"2017"}, nobtag_catagories, {"600"}, 1.012, 0.982)
  //    ({"2017"},   btag_catagories, {"600"}, 0.984, 1.023)
  //    ({"2018"}, nobtag_catagories, {"600"}, 1.012, 0.981)
  //    ({"2018"},   btag_catagories, {"600"}, 0.985, 1.024)
  //    ({"2016"}, nobtag_catagories, {"1000"}, 1.028, 0.985)
  //    ({"2016"},   btag_catagories, {"1000"}, 0.965, 1.019)
  //    ({"2017"}, nobtag_catagories, {"1000"}, 1.028, 0.985)
  //    ({"2017"},   btag_catagories, {"1000"}, 0.968, 1.017)
  //    ({"2018"}, nobtag_catagories, {"1000"}, 1.028, 0.984)
  //    ({"2018"},   btag_catagories, {"1000"}, 0.970, 1.017)
  //    ({"2016"}, nobtag_catagories, {"1400"}, 1.027, 0.975)
  //    ({"2016"},   btag_catagories, {"1400"}, 0.967, 1.030)
  //    ({"2017"}, nobtag_catagories, {"1400"}, 1.028, 0.975)
  //    ({"2017"},   btag_catagories, {"1400"}, 0.970, 1.027)
  //    ({"2018"}, nobtag_catagories, {"1400"}, 1.029, 0.973)
  //    ({"2018"},   btag_catagories, {"1400"}, 0.971, 1.027)
  //    ({"2016"}, nobtag_catagories, {"2000"}, 1.038, 0.981)
  //    ({"2016"},   btag_catagories, {"2000"}, 0.957, 1.021)
  //    ({"2017"}, nobtag_catagories, {"2000"}, 1.040, 0.981)
  //    ({"2017"},   btag_catagories, {"2000"}, 0.961, 1.019)
  //    ({"2018"}, nobtag_catagories, {"2000"}, 1.040, 0.980)
  //    ({"2018"},   btag_catagories, {"2000"}, 0.962, 1.018)
  //    ({"2016"}, nobtag_catagories, {"2300"}, 1.040, 0.979)
  //    ({"2016"},   btag_catagories, {"2300"}, 0.957, 1.022)
  //    ({"2017"}, nobtag_catagories, {"2300"}, 1.040, 0.979)
  //    ({"2017"},   btag_catagories, {"2300"}, 0.962, 1.020)
  //    ({"2018"}, nobtag_catagories, {"2300"}, 1.042, 0.978)
  //    ({"2018"},   btag_catagories, {"2300"}, 0.962, 1.020)
  //    ({"2016"}, nobtag_catagories, {"2900"}, 1.041, 0.961)
  //    ({"2016"},   btag_catagories, {"2900"}, 0.955, 1.042)
  //    ({"2017"}, nobtag_catagories, {"2900"}, 1.041, 0.964)
  //    ({"2017"},   btag_catagories, {"2900"}, 0.960, 1.035)
  //    ({"2018"}, nobtag_catagories, {"2900"}, 1.043, 0.962)
  //    ({"2018"},   btag_catagories, {"2900"}, 0.961, 1.034)
  //    ({"2016"}, nobtag_catagories, {"3500"}, 1.051, 0.962)
  //    ({"2016"},   btag_catagories, {"3500"}, 0.946, 1.041)
  //    ({"2017"}, nobtag_catagories, {"3500"}, 1.050, 0.963)
  //    ({"2017"},   btag_catagories, {"3500"}, 0.952, 1.035)
  //    ({"2018"}, nobtag_catagories, {"3500"}, 1.053, 0.964)
  //    ({"2018"},   btag_catagories, {"3500"}, 0.953, 1.032)
  // );

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
  // Notes:
  // ##########################################################################

  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process(mssm_ggH_t_signals)
      .AddSyst(cb, "Hdamp_ggH_t_REWEIGHT", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process(mssm_ggH_b_signals)
      .AddSyst(cb, "Hdamp_ggH_b_REWEIGHT", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process(mssm_ggH_i_signals)
      .AddSyst(cb, "Hdamp_ggH_i_REWEIGHT", "shape", SystMap<>::init(1.00));

  // ##########################################################################
  // Uncertainty: ggH Reweighting QCDscale uncertainty
  // References:
  // - 
  // Notes:
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

  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_htt_eff_b_$ERA", "shape", SystMap<>::init(1.00));

  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_htt_mistag_b_$ERA", "shape", SystMap<>::init(1.00));

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

  cb.cp()
      .channel({"em", "et"})
      .process(mc_processes)
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

  // Embedded MET correction uncertainties.
  // Source: https://indico.cern.ch/event/1005631/contributions/4221822/attachments/2185217/3693362/Embed_MET_correction.pdf
  if (not sm) {
      if (era == 2016) {
          // Statistical component
          cb.cp()
              .channel({"tt", "et", "mt", "em"})
              .process({"EMB"})
              .AddSyst(cb, "CMS_scale_embed_met_stat_$CHANNEL_$ERA", "shape", SystMap<channel>::init
                          ({"tt"}, 0.10)
                          ({"et"}, 0.07)
                          ({"mt"}, 0.05)
                          ({"em"}, 0.25));
          // Systematic component from MC.
          cb.cp()
              .channel({"tt", "et", "mt", "em"})
              .process({"EMB"})
              .AddSyst(cb, "CMS_scale_embed_met_syst_$CHANNEL_$ERA", "shape", SystMap<channel>::init
                          ({"tt"}, 0.48)
                          ({"et"}, 0.95)
                          ({"mt"}, 1.25)
                          ({"em"}, 1.48));
          // Additional fractional correction for determination of correction.
          cb.cp()
              .channel({"tt", "et", "mt", "em"})
              .process({"EMB"})
              .AddSyst(cb, "CMS_scale_embed_met_frac_$CHANNEL_$ERA", "shape", SystMap<>::init(0.20));
      }
      else if (era == 2017) {
          // Statistical component
          cb.cp()
              .channel({"tt", "et", "mt", "em"})
              .process({"EMB"})
              .AddSyst(cb, "CMS_scale_embed_met_stat_$CHANNEL_$ERA", "shape", SystMap<channel>::init
                          ({"tt"}, 0.06)
                          ({"et"}, 0.08)
                          ({"mt"}, 0.04)
                          ({"em"}, 0.05));
          // Systematic component from MC.
          cb.cp()
              .channel({"tt", "et", "mt", "em"})
              .process({"EMB"})
              .AddSyst(cb, "CMS_scale_embed_met_syst_$CHANNEL_$ERA", "shape", SystMap<channel>::init
                          ({"tt"}, 0.15)
                          ({"et"}, 0.47)
                          ({"mt"}, 0.29)
                          ({"em"}, 0.36));
          // Additional fractional correction for determination of correction.
          cb.cp()
              .channel({"tt", "et", "mt", "em"})
              .process({"EMB"})
              .AddSyst(cb, "CMS_scale_embed_met_frac_$CHANNEL_$ERA", "shape", SystMap<>::init(0.20));
      }
      else if (era == 2018) {
          // Statistical component
          cb.cp()
              .channel({"tt", "et", "mt", "em"})
              .process({"EMB"})
              .AddSyst(cb, "CMS_scale_embed_met_stat_$CHANNEL_$ERA", "shape", SystMap<channel>::init
                          ({"tt"}, 0.05)
                          ({"et"}, 0.05)
                          ({"mt"}, 0.03)
                          ({"em"}, 0.05));
          // Systematic component from MC.
          cb.cp()
              .channel({"tt", "et", "mt", "em"})
              .process({"EMB"})
              .AddSyst(cb, "CMS_scale_embed_met_syst_$CHANNEL_$ERA", "shape", SystMap<channel>::init
                          ({"tt"}, 0.26)
                          ({"et"}, 0.40)
                          ({"mt"}, 0.38)
                          ({"em"}, 0.68));
          // Additional fractional correction for determination of correction.
          cb.cp()
              .channel({"tt", "et", "mt", "em"})
              .process({"EMB"})
              .AddSyst(cb, "CMS_scale_embed_met_frac_$CHANNEL_$ERA", "shape", SystMap<>::init(0.20));
      }
  }

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
  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process({"TTT", "TTL", "TTJ", "TT"})
      .AddSyst(cb, "CMS_htt_tjXsec", "lnN", SystMap<>::init(1.06));

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
      .AddSyst(cb, "CMS_htt_ttbarShape", "shape", SystMap<>::init(1.00));

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
  cb.cp()
      .channel({"et"})
      .process({"ZL"})
      .AddSyst(cb, "CMS_ZLShape_$CHANNEL_1prong_barrel_$ERA", "shape",
               SystMap<>::init(1.00));

  cb.cp()
      .channel({"et"})
      .process({"ZL"})
      .AddSyst(cb, "CMS_ZLShape_$CHANNEL_1prong1pizero_barrel_$ERA", "shape",
               SystMap<>::init(1.00));
  cb.cp()
      .channel({"et"})
      .process({"ZL"})
      .AddSyst(cb, "CMS_ZLShape_$CHANNEL_1prong_endcap_$ERA", "shape",
               SystMap<>::init(1.00));

  cb.cp()
      .channel({"et"})
      .process({"ZL"})
      .AddSyst(cb, "CMS_ZLShape_$CHANNEL_1prong1pizero_endcap_$ERA", "shape",
               SystMap<>::init(1.00));

  // Electron fakes
  cb.cp()
      .channel({"et"})
      .process({"ZL"})
      .AddSyst(cb, "CMS_fake_e_BA_$ERA", "shape",
               SystMap<>::init(1.00));
  cb.cp()
      .channel({"et"})
      .process({"ZL"})
      .AddSyst(cb, "CMS_fake_e_EC_$ERA", "shape",
               SystMap<>::init(1.00));

  // Muon fakes
  cb.cp()
      .channel({"mt"})
      .process({"ZL"})
      .AddSyst(cb, "CMS_fake_m_WH1_$ERA", "shape",
               SystMap<>::init(1.00));
  cb.cp()
      .channel({"mt"})
      .process({"ZL"})
      .AddSyst(cb, "CMS_fake_m_WH2_$ERA", "shape",
               SystMap<>::init(1.00));
  cb.cp()
      .channel({"mt"})
      .process({"ZL"})
      .AddSyst(cb, "CMS_fake_m_WH3_$ERA", "shape",
               SystMap<>::init(1.00));
  cb.cp()
      .channel({"mt"})
      .process({"ZL"})
      .AddSyst(cb, "CMS_fake_m_WH4_$ERA", "shape",
               SystMap<>::init(1.00));
  cb.cp()
      .channel({"mt"})
      .process({"ZL"})
      .AddSyst(cb, "CMS_fake_m_WH5_$ERA", "shape",
               SystMap<>::init(1.00));

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
        for (auto reg: unc_regions) {
            for (auto unc: qcd_tt_uncs) {
                cb.cp()
                    .channel({"et", "mt", "tt"})
                    .process({"jetFakes"})
                    .AddSyst(cb, "CMS_ff_total_qcd_stat_"+njet+"_jet_pt_"+reg+"_"+unc+"_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));
            }
        }
    }

    // W shape stat.
    for (auto njet: jet_bins) {
        for (auto reg: unc_regions) {
            for (auto unc: wjets_uncs) {
                cb.cp()
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
            cb.cp()
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
            .AddSyst(cb, "CMS_ff_total_wjets_stat_extrap_"+unc+"_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));

        cb.cp()
            .channel({"et", "mt"})
            .process({"jetFakes"})
            .AddSyst(cb, "CMS_ff_total_ttbar_stat_met_"+unc+"_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));

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
        .channel({"et", "mt"})
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
