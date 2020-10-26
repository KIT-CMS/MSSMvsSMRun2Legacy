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

void AddMSSMvsSMRun2Systematics(CombineHarvester &cb, bool jetfakes, bool embedding, bool regional_jec, bool ggh_wg1, bool qqh_wg1, int era, std::string heavy_mass, std::string light_mass) {

  // ##########################################################################
  // Define groups of signal processes
  // ##########################################################################

  std::vector<std::string> signals_ggH = {
      // STXS stage 0
      "ggH125",
      // STXS stage 1.1
      "ggH_GG2H_FWDH125",
      "ggH_GG2H_PTH_GT200125",
      "ggH_GG2H_0J_PTH_0_10125",
      "ggH_GG2H_0J_PTH_GT10125",
      "ggH_GG2H_1J_PTH_0_60125",
      "ggH_GG2H_1J_PTH_60_120125",
      "ggH_GG2H_1J_PTH_120_200125",
      "ggH_GG2H_GE2J_MJJ_0_350_PTH_0_60125",
      "ggH_GG2H_GE2J_MJJ_0_350_PTH_60_120125",
      "ggH_GG2H_GE2J_MJJ_0_350_PTH_120_200125",
      "ggH_GG2H_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_0_25125",
      "ggH_GG2H_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_GT25125",
      "ggH_GG2H_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_0_25125",
      "ggH_GG2H_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_GT25125",
      };
  std::vector<std::string> signals_qqH = {
      // STXS stage 0
      "qqH125",
      // STXS stage 1
      "qqH_QQ2HQQ_FWDH125",
      "qqH_QQ2HQQ_0J125",
      "qqH_QQ2HQQ_1J125",
      "qqH_QQ2HQQ_GE2J_MJJ_0_60125",
      "qqH_QQ2HQQ_GE2J_MJJ_60_120125",
      "qqH_QQ2HQQ_GE2J_MJJ_120_350125",
      "qqH_QQ2HQQ_GE2J_MJJ_GT350_PTH_GT200125",
      "qqH_QQ2HQQ_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_0_25125",
      "qqH_QQ2HQQ_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_GT25125",
      "qqH_QQ2HQQ_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_0_25125",
      "qqH_QQ2HQQ_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_GT25"
      };
  std::vector<std::string> signals_VH = {
      // STXS stage 0
      "WH125", "ZH125", "ttH125", "VH125"};
  std::vector<std::string> signals_ggHToWW = {
     // STXS stage 0
     "ggHWW125"};
  std::vector<std::string> signals_qqHToWW = {
     // STXS stage 0
     "qqHWW125"};
  std::vector<std::string> signals = JoinStr({signals_ggH, signals_qqH, signals_VH, {"qqh"}});

  std::vector<std::string> nmssm_signals = {"NMSSM_"+heavy_mass+"_125_"+light_mass};

  // All processes being taken from simulation
  std::vector<std::string> mc_processes =
      JoinStr({
              signals,
              signals_ggHToWW,
              signals_qqHToWW,
              nmssm_signals,
              {"ZTT", "TT", "TTT", "TTL", "TTJ", "ZJ", "ZL", "VV", "VVT", "VVL", "VVJ", "ST", "W"}
              });

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
      // mt cross trigger
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
          .AddSyst(cb, "CMS_eff_trigger_tt_dm"+tauTriggerbin+"_$ERA", "shape", SystMap<>::init(1.00));
          
      cb.cp()
          .channel({"tt"})
          .process({"EMB"})
          .AddSyst(cb, "CMS_eff_trigger_emb_tt_dm"+tauTriggerbin+"_$ERA", "shape", SystMap<>::init(0.866));

      // Correlated component acting on Embedded
      cb.cp()
          .channel({"tt"})
          .process({"EMB"})
          .AddSyst(cb, "CMS_eff_trigger_tt_dm"+tauTriggerbin+"_$ERA", "shape", SystMap<>::init(0.5));
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

  std::string tauIDptbins[5] = {"30-35", "35-40", "40-500", "500-1000", "1000-inf"};
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
        .process(JoinStr({signals, nmssm_signals, {"ZTT", "TTT", "TTL", "VVT", "VVL"}}))
        .AddSyst(cb, "CMS_eff_t_"+tauIDbin+"_$ERA", "shape", SystMap<>::init(1.0));
  }
  cb.cp() //second part uncorrelated between channels for IDvsLep
      .channel({"et", "mt"})
      .process(JoinStr({signals, nmssm_signals, {"ZTT", "TTT", "TTL", "VVT", "VVL"}}))
      .AddSyst(cb, "CMS_eff_t_$CHANNEL_$ERA", "lnN", SystMap<>::init(1.01));

  // Tau ID: tt with 2 real taus
  for (auto tauIDbin : tauIDdmbins){
    cb.cp()
        .channel({"tt"})
        .process(JoinStr({signals, nmssm_signals, {"ZTT", "TTT", "TTL", "VVT", "VVL"}}))
        .AddSyst(cb, "CMS_eff_t_dm"+tauIDbin+"_$ERA", "shape", SystMap<>::init(1.0));
  }
  cb.cp()
      .channel({"tt"})
      .process(JoinStr({signals, nmssm_signals, {"ZTT", "TTT", "TTL", "VVT", "VVL"}}))
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
      .AddSyst(cb, "CMS_eff_t_emb_$CHANNEL_$ERA", "lnN", SystMap<>::init(1.005));

  // Tau ID: tt with 2 real taus
  for (auto tauIDbin : tauIDdmbins){
    cb.cp()
        .channel({"tt"})
        .process({"EMB"})
        .AddSyst(cb, "CMS_eff_t_emb_dm"+tauIDbin+"_$ERA", "shape", SystMap<>::init(0.866));
  }
  cb.cp()
      .channel({"tt"})
      .process({"EMB"})
      .AddSyst(cb, "CMS_eff_t_emb_$CHANNEL_$ERA", "lnN", SystMap<>::init(1.007));


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
      .AddSyst(cb, "CMS_eff_t_$CHANNEL_$ERA", "lnN", SystMap<>::init(1.0087));

  // Tau ID: tt with 2 real taus
  for (auto tauIDbin : tauIDdmbins){
    cb.cp()
        .channel({"tt"})
        .process({"EMB"})
        .AddSyst(cb, "CMS_eff_t_dm"+tauIDbin+"_$ERA", "shape", SystMap<>::init(0.5));
  }
  cb.cp()
      .channel({"tt"})
      .process({"EMB"})
      .AddSyst(cb, "CMS_eff_t_$CHANNEL_$ERA", "lnN", SystMap<>::init(1.012));

  // Tau ID: tt with 1 real taus and 1 jet fake
  cb.cp()
      .channel({"tt"})
      .process({"W", "ZJ", "TTJ", "VVJ"})
      .AddSyst(cb, "CMS_eff_t_$ERA", "lnN", SystMap<>::init(1.06));

  cb.cp()
      .channel({"tt"})
      .process({"W", "ZJ", "TTJ", "VVJ"})
      .AddSyst(cb, "CMS_eff_t_$CHANNEL_$ERA", "lnN", SystMap<>::init(1.02));



  // ##########################################################################
  // Uncertainty: b-tag and mistag efficiency
  // References:
  // Notes:
  // - FIXME: References?
  // ##########################################################################

  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_htt_eff_b_$ERA", "shape", SystMap<>::init(0.71));
  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_htt_eff_b_$ERA_$CHANNEL", "shape", SystMap<>::init(0.71));

  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_htt_mistag_b_$ERA", "shape", SystMap<>::init(0.71));
  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_htt_mistag_b_$ERA_$CHANNEL", "shape", SystMap<>::init(0.71));

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
      //.AddSyst(cb, "CMS_scale_mc_e", "shape", SystMap<>::init(0.71));

  cb.cp()
      .channel({"em", "et"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_res_e", "shape", SystMap<>::init(1.00));
      //.AddSyst(cb, "CMS_scale_mc_e", "shape", SystMap<>::init(0.71));

  // Embedded uncorrelated uncertainty

  cb.cp()
      .channel({"em", "et"})
      .process({"EMB"})
      .AddSyst(cb, "CMS_scale_e_emb", "shape", SystMap<>::init(1.00));
      //.AddSyst(cb, "CMS_scale_emb_e", "shape", SystMap<>::init(0.71));


  // ##########################################################################
  // Uncertainty: Tau energy scale
  // References:
  // Notes:
  // - Tau energy scale is split by decay mode.
  // - FIXME: References?
  // - FIXME: Need it for H->WW in mt, et, (and tt)?
  // ##########################################################################


  // Common component acting on MC

  cb.cp()
      .channel({"et", "mt", "tt"})
      .process(JoinStr({signals, nmssm_signals, {"ZTT", "TTT", "TTL", "VVT", "VVL", "jetFakes"}}))
      .AddSyst(cb, "CMS_scale_t_1prong_$ERA", "shape", SystMap<>::init(1.0));

  cb.cp()
      .channel({"et", "mt", "tt"})
      .process(JoinStr({signals, nmssm_signals, {"ZTT", "TTT", "TTL", "VVT", "VVL", "jetFakes"}}))
      .AddSyst(cb, "CMS_scale_t_1prong1pizero_$ERA", "shape",
               SystMap<>::init(1.0));

  cb.cp()
      .channel({"et", "mt", "tt"})
      .process(JoinStr({signals, nmssm_signals, {"ZTT", "TTT", "TTL", "VVT", "VVL", "jetFakes"}}))
      .AddSyst(cb, "CMS_scale_t_3prong_$ERA", "shape", SystMap<>::init(1.0));

  cb.cp()
      .channel({"et", "mt", "tt"})
      .process(JoinStr({signals, nmssm_signals, {"ZTT", "TTT", "TTL", "VVT", "VVL", "jetFakes"}}))
      .AddSyst(cb, "CMS_scale_t_3prong1pizero_$ERA", "shape",
               SystMap<>::init(1.0));

  // Component for EMB only
  if(embedding){
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"EMB", "jetFakes"})
      .AddSyst(cb, "CMS_scale_t_emb_1prong_$ERA", "shape", SystMap<>::init(0.866));

  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"EMB", "jetFakes"})
      .AddSyst(cb, "CMS_scale_t_emb_1prong1pizero_$ERA", "shape", SystMap<>::init(0.866));

  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"EMB", "jetFakes"})
      .AddSyst(cb, "CMS_scale_t_emb_3prong_$ERA", "shape", SystMap<>::init(0.866));

  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"EMB", "jetFakes"})
      .AddSyst(cb, "CMS_scale_t_emb_3prong1pizero_$ERA", "shape", SystMap<>::init(0.866));
  }
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

  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process({"ZTT", "TT", "TTT", "TTL", "TTJ", "ZJ", "ZL", "VV", "VVT", "VVL", "VVJ", "ST", "W"})  //Z and W processes are only included due to the EWK fraction. Make sure that there is no contribution to the shift from the DY or Wjets samples.
      .AddSyst(cb, "CMS_scale_met_unclustered", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process(JoinStr({signals, nmssm_signals, signals_ggHToWW, signals_qqHToWW, {"WHWW125", "ZHWW125"}, {"ZTT", "ZL", "ZJ", "W"}}))
      .AddSyst(cb, "CMS_htt_boson_scale_met_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process(JoinStr({signals, nmssm_signals, signals_ggHToWW, signals_qqHToWW, {"WHWW125", "ZHWW125"}, {"ZTT", "ZL", "ZJ", "W"}}))
      .AddSyst(cb, "CMS_htt_boson_res_met_$ERA", "shape", SystMap<>::init(1.00));


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

  // QCD
  cb.cp()
      .channel({"et"})
      .process({"QCD"})
      .AddSyst(cb, "CMS_ExtrapSSOS_$CHANNEL_$ERA", "lnN", SystMap<>::init(1.05));
  cb.cp()
      .channel({"mt"})
      .process({"QCD"})
      .AddSyst(cb, "CMS_ExtrapSSOS_$CHANNEL_$ERA", "lnN", SystMap<>::init(1.03));
  cb.cp()
      .channel({"tt"})
      .process({"QCD"})
      .AddSyst(cb, "CMS_ExtrapABCD_$CHANNEL_$ERA", "lnN", SystMap<>::init(1.03));

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
      .AddSyst(cb, "CMS_htt_ttbarShape", "shape", SystMap<>::init(0.71));
  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process({"TTT", "TTL", "TTJ", "TT"})
      .AddSyst(cb, "CMS_htt_ttbarShape_$ERA", "shape", SystMap<>::init(0.71));

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
  //cb.cp()
  //    .channel({"et"})
  //    .process({"ZL"})
  //    .AddSyst(cb, "CMS_fake_e_$ERA", "lnN", SystMap<>::init(1.15));

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
  //cb.cp()
  //    .channel({"mt"})
  //    .process({"ZL"})
  //    .AddSyst(cb, "CMS_fake_m_$ERA", "lnN", SystMap<>::init(1.25));

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
  // Uncertainty: Jet to tau fakes
  // References:
  // Notes:
  // - FIXME: Adapt for fake factor and embedding
  // - FIXME: References?
  // ##########################################################################

  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"W", "TTJ", "ZJ", "VVJ"})
      .AddSyst(cb, "CMS_htt_fake_j_$ERA", "shape", SystMap<>::init(1.00));

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


//  cb.cp()
//      .channel({"et", "mt", "tt", "em"})
//      .process(nmssm_signals)
//      .AddSyst(cb, "MG_scale_choice", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process(nmssm_signals)
      .AddSyst(cb, "MG_scale_norm", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process(nmssm_signals)
      .AddSyst(cb, "PDF_scale", "shape", SystMap<>::init(1.00));


  // Uncertainty on branching ratio for HTT at 125 GeV
  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process(signals)
      .AddSyst(cb, "BR_htt", "lnN", SystMap<>::init(1.021));

  // Uncertainty on branching ratio for HWW at 125 GeV
  cb.cp()
      .channel({"et", "mt", "tt", "em"})
     .process(JoinStr({signals_ggHToWW,signals_qqHToWW}))
     .AddSyst(cb, "BR_hww_THU", "lnN", SystMap<>::init(1.0099));
  cb.cp()
      .channel({"et", "mt", "tt", "em"})
     .process(JoinStr({signals_ggHToWW,signals_qqHToWW}))
     .AddSyst(cb, "BR_hww_PU_mq", "lnN", SystMap<>::init(1.0099));
  cb.cp()
      .channel({"et", "mt", "tt", "em"})
     .process(JoinStr({signals_ggHToWW,signals_qqHToWW}))
     .AddSyst(cb, "BR_hww_PU_alphas", "lnN", SystMap<>::init(1.0066));
  // QCD scale
  if (!ggh_wg1) {
  cb.cp()
      .channel({"et", "mt", "tt", "em"})
     .process(JoinStr({signals_ggH,signals_ggHToWW}))
      .AddSyst(cb, "QCDScale_ggH", "lnN", SystMap<>::init(1.039));
  }
  if (!qqh_wg1) {
      cb.cp()
          .channel({"et", "mt", "tt", "em"})
         .process(JoinStr({signals_qqH,signals_qqHToWW, {"qqh"}}))
          .AddSyst(cb, "QCDScale_qqH", "lnN", SystMap<>::init(1.005));
  }
  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process({"ZH125","VH125"})
      .AddSyst(cb, "QCDScale_VH", "lnN", SystMap<>::init(1.01));
  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process({"ttH125"})
      .AddSyst(cb, "QCDScale_ttH", "lnN", SystMap<>::init(1.08));

  // PDF
  cb.cp()
      .channel({"et", "mt", "tt", "em"})
     .process(JoinStr({signals_ggH,signals_ggHToWW}))
      .AddSyst(cb, "pdf_Higgs_gg", "lnN", SystMap<>::init(1.032));
  cb.cp()
     .channel({"et", "mt", "tt", "em"})
     .process(JoinStr({signals_qqH,signals_qqHToWW, {"qqh"}}))
     .AddSyst(cb, "pdf_Higgs_qqbar", "lnN", SystMap<>::init(1.021));
  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process({"ZH125"})
      .AddSyst(cb, "pdf_Higgs_VH", "lnN", SystMap<>::init(1.013));
  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process({"WH125","VH125"})
      .AddSyst(cb, "pdf_Higgs_VH", "lnN", SystMap<>::init(1.018));
  cb.cp()
      .channel({"et", "mt", "tt", "em"})
      .process({"ttH125"})
      .AddSyst(cb, "pdf_Higgs_ttH", "lnN", SystMap<>::init(1.036));

  // Gluon-fusion WG1 uncertainty scheme
  if (ggh_wg1) {
    cb.cp()
      .channel({"et", "mt", "tt", "em"})
       .process({signals_ggH})
      .AddSyst(cb, "THU_ggH_Mig01", "shape", SystMap<>::init(1.00));
    cb.cp()
      .channel({"et", "mt", "tt", "em"})
       .process({signals_ggH})
      .AddSyst(cb, "THU_ggH_Mig12", "shape", SystMap<>::init(1.00));
    cb.cp()
      .channel({"et", "mt", "tt", "em"})
       .process({signals_ggH})
      .AddSyst(cb, "THU_ggH_Mu", "shape", SystMap<>::init(1.00));
    cb.cp()
      .channel({"et", "mt", "tt", "em"})
       .process({signals_ggH})
      .AddSyst(cb, "THU_ggH_PT120", "shape", SystMap<>::init(1.00));
    cb.cp()
      .channel({"et", "mt", "tt", "em"})
       .process({signals_ggH})
      .AddSyst(cb, "THU_ggH_PT60", "shape", SystMap<>::init(1.00));
    cb.cp()
      .channel({"et", "mt", "tt", "em"})
       .process({signals_ggH})
      .AddSyst(cb, "THU_ggH_Res", "shape", SystMap<>::init(1.00));
    cb.cp()
      .channel({"et", "mt", "tt", "em"})
       .process({signals_ggH})
      .AddSyst(cb, "THU_ggH_VBF2j", "shape", SystMap<>::init(1.00));
    cb.cp()
      .channel({"et", "mt", "tt", "em"})
       .process({signals_ggH})
      .AddSyst(cb, "THU_ggH_VBF3j", "shape", SystMap<>::init(1.00));
    cb.cp()
      .channel({"et", "mt", "tt", "em"})
       .process({signals_ggH})
      .AddSyst(cb, "THU_ggH_qmtop", "shape", SystMap<>::init(1.00));
    cb.cp()
      .channel({"et", "mt", "tt", "em"})
     .process({signals_ggHToWW})
      .AddSyst(cb, "QCDScale_ggHWW", "lnN", SystMap<>::init(1.039));
  }
 // VBF WG1 uncertainty scheme
 if (qqh_wg1) {
   cb.cp()
     .channel({"et", "mt", "tt", "em"})
     .process({signals_qqH})
     .AddSyst(cb, "THU_qqH_TOT", "shape", SystMap<>::init(1.00));
   cb.cp()
     .channel({"et", "mt", "tt", "em"})
     .process({signals_qqH})
     .AddSyst(cb, "THU_qqH_PTH200", "shape", SystMap<>::init(1.00));
   cb.cp()
     .channel({"et", "mt", "tt", "em"})
     .process({signals_qqH})
     .AddSyst(cb, "THU_qqH_Mjj60", "shape", SystMap<>::init(1.00));
   cb.cp()
     .channel({"et", "mt", "tt", "em"})
     .process({signals_qqH})
     .AddSyst(cb, "THU_qqH_Mjj120", "shape", SystMap<>::init(1.00));
   cb.cp()
     .channel({"et", "mt", "tt", "em"})
     .process({signals_qqH})
     .AddSyst(cb, "THU_qqH_Mjj350", "shape", SystMap<>::init(1.00));
   cb.cp()
     .channel({"et", "mt", "tt", "em"})
     .process({signals_qqH})
     .AddSyst(cb, "THU_qqH_Mjj700", "shape", SystMap<>::init(1.00));
   cb.cp()
     .channel({"et", "mt", "tt", "em"})
     .process({signals_qqH})
     .AddSyst(cb, "THU_qqH_Mjj1000", "shape", SystMap<>::init(1.00));
   cb.cp()
     .channel({"et", "mt", "tt", "em"})
     .process({signals_qqH})
     .AddSyst(cb, "THU_qqH_Mjj1500", "shape", SystMap<>::init(1.00));
   cb.cp()
     .channel({"et", "mt", "tt", "em"})
     .process({signals_qqH})
     .AddSyst(cb, "THU_qqH_25", "shape", SystMap<>::init(1.00));
   cb.cp()
     .channel({"et", "mt", "tt", "em"})
     .process({signals_qqH})
     .AddSyst(cb, "THU_qqH_JET01", "shape", SystMap<>::init(1.00));
   cb.cp()
     .channel({"et", "mt", "tt", "em"})
      .process({signals_qqHToWW})
     .AddSyst(cb, "QCDScale_qqH", "lnN", SystMap<>::init(1.005));
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
    .AddSyst(cb, "CMS_htt_emb_ttbar_$ERA", "shape", SystMap<>::init(1.0));

  // Uncertainty of hadronic tau track efficiency correction
  // uncorrelated between eras
  cb.cp()
    .channel({"et", "mt", "tt"})
    .process({"EMB"})
    .AddSyst(cb, "CMS_3ProngEff_$ERA", "shape", SystMap<>::init(0.71));

  cb.cp()
    .channel({"et", "mt", "tt"})
    .process({"EMB"})
    .AddSyst(cb, "CMS_1ProngPi0Eff_$ERA", "shape", SystMap<>::init(0.71));
  // correlated between eras
  cb.cp()
    .channel({"et", "mt", "tt"})
    .process({"EMB"})
    .AddSyst(cb, "CMS_3ProngEff", "shape", SystMap<>::init(0.71));

  cb.cp()
    .channel({"et", "mt", "tt"})
    .process({"EMB"})
    .AddSyst(cb, "CMS_1ProngPi0Eff", "shape", SystMap<>::init(0.71));

  // ##########################################################################
  // Uncertainty: Jet fakes
  // References:
  // - https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauJet2TauFakes
  // Notes:
  // - FIXME: add 2017 norm uncertainties, and properly correlate across years
  // ##########################################################################
  
  // QCD shape stat.
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_qcd_njet0_morphed_stat_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_qcd_njet1_morphed_stat_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_qcd_njet2_morphed_stat_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));


  // W shape stat.
  cb.cp()
      .channel({"et", "mt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_w_njet0_morphed_stat_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"et", "mt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_w_njet1_morphed_stat_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"et", "mt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_w_njet2_morphed_stat_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));



  // TT shape stat.
  cb.cp()
      .channel({"et", "mt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_tt_njet0_morphed_stat_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));

  cb.cp()
      .channel({"et", "mt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_tt_njet1_morphed_stat_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));
  
  // MC subtraction uncertainty
  // uncorrelated between eras
  cb.cp()
      .channel({"et", "mt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_w_mc_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_qcd_mc_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"et", "mt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_frac_w_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));

      
  // Shape syst. of different contributions (QCD/W/tt)
  // uncorrelated between eras
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_qcd_mvis_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_qcd_mvis_osss_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_corr_qcd_mvis_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_corr_qcd_mvis_osss_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"et", "mt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_qcd_muiso_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"et", "mt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_corr_qcd_muiso_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_qcd_tau2_pt_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_corr_qcd_tau2_pt_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_w_syst_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_tt_syst_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"et", "mt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_tt_morphed_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"et", "mt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_corr_tt_syst_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
  //cb.cp()
  //    .channel({"et", "mt"})
  //    .process({"jetFakes"})
  //    .AddSyst(cb, "CMS_ff_tt_sf_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"et", "mt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_w_lepPt_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"et", "mt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_corr_w_lepPt_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"et", "mt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_w_mt_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"et", "mt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_corr_w_mt_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
            

  //below: jetFakes norm uncertainties. Current values are for 2016, which are probably a good approx. for 2017. To be updated.
          


//   // Stat. norm (uncorrelated across years)
//   cb.cp()
//       .channel({"et", "mt", "tt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_norm_stat_$CHANNEL_$ERA", "lnN", SystMap<channel>::init
// 	       ({"mt"}, 1.026) //incl
// 	       ({"et"}, 1.046) //incl
// 	       ({"tt"}, 1.029) //incl
// 	       );
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_norm_stat_$CHANNEL_$BIN_$ERA", "lnN", SystMap<channel, bin_id>::init
	       ({"mt"}, {1},  1.026)
	       ({"mt"}, {2},  1.026) 
	       ({"mt"}, {3},  1.026)
	       ({"mt"}, {4},  1.026)
	       ({"mt"}, {5},  1.026) 
	       ({"mt"}, {6},  1.026)
	       ({"mt"}, {7},  1.026) 
	       ({"mt"}, {8},  1.026)
	       ({"mt"}, {9},  1.026)
	       ({"mt"}, {10},  1.026)
	       ({"mt"}, {11},  1.026)
	       ({"mt"}, {12},  1.026)
	       ({"mt"}, {300}, 1.026)
	       ({"et"}, {1},  1.046)
	       ({"et"}, {2},  1.046) 
	       ({"et"}, {3},  1.046)
	       ({"et"}, {4},  1.046)
	       ({"et"}, {5},  1.046) 
	       ({"et"}, {6},  1.046)
	       ({"et"}, {7},  1.046) 
	       ({"et"}, {8},  1.046)
	       ({"et"}, {9},  1.046)
	       ({"et"}, {10},  1.046)
	       ({"et"}, {11},  1.046)
	       ({"et"}, {12},  1.046)
	       ({"et"}, {300}, 1.046)
	       ({"tt"}, {1},  1.029)
	       ({"tt"}, {2},  1.029) 
	       ({"tt"}, {3},  1.029)
	       ({"tt"}, {4},  1.029)
	       ({"tt"}, {5},  1.029) 
	       ({"tt"}, {6},  1.029)
	       ({"tt"}, {7},  1.029) 
	       ({"tt"}, {8},  1.029)
	       ({"tt"}, {9},  1.029)
	       ({"tt"}, {10},  1.029)
	       ({"tt"}, {11},  1.029)
	       ({"tt"}, {12},  1.029)
	       ({"tt"}, {300}, 1.029)
	       );
  // Syst. norm: Bin-correlated
  // uncorrelated between eras
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_norm_syst_$CHANNEL_$ERA", "lnN", SystMap<channel>::init
	       ({"mt"}, 1.042) //incl
	       ({"et"}, 1.042) //incl
	       ({"tt"}, 1.067) //incl
	       );
  // correlated between eras
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_norm_syst_$CHANNEL", "lnN", SystMap<channel>::init
	       ({"mt"}, 1.042) //incl
	       ({"et"}, 1.042) //incl
	       ({"tt"}, 1.067) //incl
	       );
}
} // namespace ch
