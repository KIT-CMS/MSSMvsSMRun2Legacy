#include <iterator>
#include <algorithm>
#include <string>
#include <map>
#include <set>
#include <iostream>
#include <utility>
#include <vector>
#include <cstdlib>
#include "boost/algorithm/string/predicate.hpp"
#include "boost/program_options.hpp"
#include "boost/lexical_cast.hpp"
#include "CombineHarvester/CombineTools/interface/CombineHarvester.h"
#include "CombineHarvester/CombineTools/interface/Observation.h"
#include "CombineHarvester/CombineTools/interface/Process.h"
#include "CombineHarvester/CombineTools/interface/Utilities.h"
#include "CombineHarvester/CombineTools/interface/Systematics.h"
#include "CombineHarvester/CombineTools/interface/BinByBin.h"
#include "CombineHarvester/CombineTools/interface/AutoRebin.h"
#include "CombineHarvester/CombinePdfs/interface/MorphFunctions.h"
#include "CombineHarvester/CombineTools/interface/HttSystematics.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "TH2.h"

using namespace std;
using boost::starts_with;
using ch::syst::SystMap;
using ch::syst::SystMapAsymm;
using ch::syst::era;
using ch::syst::channel;
using ch::syst::bin_id;
using ch::syst::process;
using ch::JoinStr;
using namespace ch;
namespace po = boost::program_options;

typedef vector<string> VString;
typedef vector<pair<int, string>> Categories;

bool debug = false;


void dout() {
    if (debug) std::cout << std::endl;
}
template <typename Head, typename... Tail>
void dout(Head H, Tail... T) {
    if (debug) std::cout << H << ' ';
    dout(T...);
}
template<typename T>
void printVector(const T& t) {
    std::copy(t.cbegin(), t.cend(), std::ostream_iterator<typename T::value_type>(std::cout, ", "));
}
template<typename T>
void printVectorInVector(const T& t) {
    std::for_each(t.cbegin(), t.cend(), printVector<typename T::value_type>);
}
template <class T>
bool contains(std::vector<T> const &v, T const &x) {
    return ! (v.empty() && std::find(v.begin(), v.end(), x) == v.end());
}


int main(int argc, char** argv)
{
    // Parsing command line arguments
    string SM125= "";
    string mass = "mA";
    string output_folder = "output/mssm_run2";
    string input_folder = "/nfs/dust/cms/user/glusheno/Combine/MSSM-Full-2016/KIT/";
    string input_file = "/afs/desy.de/user/g/glusheno/RWTH/KIT/Shapes/ES-subanalysis/converted_shapes/htt_et.inputs-etFes.root";
    string postfix = "-mttot";
    bool newera = false;
    bool decorelate_emb = false;
    bool fakefactors = false;
    bool embedding = false;
    int binsize = 2;
    vector<float> vbins = {};
    po::variables_map vm;
    po::options_description config("configuration");
    VString channels = {"et"};
    VString catteg = {"et_0jet_dm0", "et_0jet_dm1", "et_0jet_dm10"};
    config.add_options()  // Declare the supported options.
        ("mass,m",         po::value<string>(&mass)->default_value(mass))
        ("binsize",         po::value<int>(&binsize)->default_value(binsize))
        ("vbins",         po::value<vector<float>>(&vbins)->multitoken(), "bins")
        ("input-folder",   po::value<string>(&input_folder)->default_value(input_folder))
        ("input-file,i",   po::value<string>(&input_file)->default_value(input_file))
        ("output-folder",  po::value<string>(&output_folder)->default_value(output_folder))
        ("postfix",        po::value<string>(&postfix)->default_value(postfix))
        ("newera",         po::bool_switch(&newera), "2017 syst unc model")
        ("decorelate-emb",         po::bool_switch(&decorelate_emb), "use the values for decorelated with embedded parts")
        ("embedding",         po::bool_switch(&embedding), "used embedding")
        ("fakefactors",         po::bool_switch(&fakefactors), "used fakefactors")
        ("SM125,h",        po::value<string>(&SM125)->default_value(SM125))
        ("help",           "produce help message")
        ("catteg",         po::value<vector<string>>(&catteg)->multitoken(), "cattegories: et_nojets_alldm et_0jet_dm0 et_0jet_dm1 et_0jet_dm10 ")
        ("channels,c",     po::value<vector<string>>(&channels)->multitoken(), "channels.")
        ("debug,d",        po::bool_switch(&debug), "debug printout");
    po::store(po::parse_command_line(argc, argv, config), vm);
    po::notify(vm);
    if (vm.count("help"))
    {
        cout << config << "\n";
        return 0;
    }
    if (debug)
    {
        std::cout << "checnnels: ";
        for (auto i = channels.begin(); i != channels.end(); ++i)
        std::cout << *i << ' ';
        std::cout << '\n';
    }
    string input_dir = input_folder;
    string output_dir = output_folder + "/";
    if (!boost::filesystem::exists(output_dir))
        boost::filesystem::create_directories(output_dir);
    // Prepering the combine run
    RooRealVar faketaues("faketaues", "faketaues", 0.0, -4.0, 18.0);
    vector<string> energy_scales = {
        "-4","-3","-2","-1.75","-1.5","-1.25","-1","-0.75","-0.5","-0.25",
        "0",
        "0.25","0.5","0.75","1","1.25","1.5","1.75","2","2.25","2.5",
        "2.75","3","3.25","3.5","3.75","4","5","6","7","8",
        "9","10","11","12","13","14","15","16","17","18"
    };
    map<string, VString> bkg_procs;
    map<string, Categories> cattegories;
    map<string, vector<double> > binning;
    // m_vis original binning: [0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,155,160]
    // vector<double> _nojets = {70,75,80,85,90,95,100,105,110,115,120,125,130,135,140}; // fit of m_vis > 70
    // vector<double> _nojets = {64,72,80,88,96,104,112,120,128,136,144};
    // vector<double> _nojets = {66,68,70,72,74,76,78,80,82,84,86,88,90,92,94,96,98,100,102,104,106,108,110,112,114,116,118,120,122,124,126,128,130,132,134,136,138,140};
    vector<double> _nojets;
    if (vbins.size() > 0)
    {
        dout("the binning set manually:", vbins.size());
        // for(auto v: vbins)
        // {
        //     dout(v, "->", vbins[v]);
        // }
        for (const auto v: vbins) dout(v);
        _nojets.insert(_nojets.end(), std::begin(vbins), std::end(vbins));
    }
    else if (binsize == 2)
    {
        // double a [] = {70,72,74,76,78,80,82,84,86,88,90,92,94,96,98,100,102,104,106,108,110,112,114,116,118,120,122,124,126,128,130,132,134,136,138,140};
        double a [] = {70,72,74,76,78,80,82,84,86,88,90,92,94,96,98,100,102,104,106,108,110,112,114,116,118,120};
        _nojets.insert(_nojets.end(), std::begin(a), std::end(a));
    }
    else if (binsize == 5)
    {
        double a [] = {70,75,80,85,90,95,100,105,110,115,120,125,130,135,140}; // fit of m_vis > 70
        _nojets.insert(_nojets.end(), std::begin(a), std::end(a));
    }
    else
    {
        double a [] = {72,80,88,96,104,112,120,128,136,144};
        _nojets.insert(_nojets.end(), std::begin(a), std::end(a));
    }

    // TODO std::set https://wandbox.org/permlink/8TFTqyVoize91zcW
    bkg_procs["et"] = {"TTL", "VVL"};
    if (fakefactors)
        bkg_procs["et"].push_back("jetFakes");
    else
    {
        std::vector<string> newElements = {"ZJ", "TTJ", "VVJ", "W", "QCD"};
        std::copy(begin(newElements), end(newElements), std::back_inserter(bkg_procs["et"]));
    }

    if (embedding)
        bkg_procs["et"].push_back("EMB");
    else
    {
        std::vector<string> newElements = {"ZTT", "TTT", "VVT"};
        std::copy(begin(newElements), end(newElements), std::back_inserter(bkg_procs["et"]));
    }

    cattegories["et_13TeV"] = {};


    std::map<string, int> cattegories_bins;
    // cattegories_bins["et_nojets_alldm"]=1;
    // cattegories_bins["et_njet0_alldm"]=1;
    // cattegories_bins["et_0jet_dm0"]=2;
    // cattegories_bins["et_njetN_dm0"]=2;
    // cattegories_bins["et_nojets_dm1"]=3;
    // cattegories_bins["et_njet0_dm1"]=3;
    // cattegories_bins["et_0jet_dm10"]=4;
    // cattegories_bins["et_njetN_dm10"]=4;
    // new categories
        // inclusive eta
        cattegories_bins["et_inc_eta_2_njet0_alldm"] = 1;
        cattegories_bins["et_inc_eta_2_njet0_dm0"] = 2;
        cattegories_bins["et_inc_eta_2_njet0_dm1"] = 3;
        cattegories_bins["et_inc_eta_2_njetN_alldm"] = 4;
        cattegories_bins["et_inc_eta_2_njetN_dm0"] = 5;
        cattegories_bins["et_inc_eta_2_njetN_dm1"] = 6;
        // barel
        cattegories_bins["et_eta_2_barel_real_njet0_alldm"] = 7;
        cattegories_bins["et_eta_2_barel_real_njet0_dm0"] = 8;
        cattegories_bins["et_eta_2_barel_real_njet0_dm1"] = 9;
        cattegories_bins["et_eta_2_barel_real_njetN_alldm"] = 10;
        cattegories_bins["et_eta_2_barel_real_njetN_dm0"] = 11;
        cattegories_bins["et_eta_2_barel_real_njetN_dm1"] = 12;
        //gap
        cattegories_bins["et_eta_2_gap_njet0_alldm"] = 13;
        cattegories_bins["et_eta_2_gap_njet0_dm0"] = 14;
        cattegories_bins["et_eta_2_gap_njet0_dm1"] = 15;
        cattegories_bins["et_eta_2_gap_njetN_alldm"] = 16;
        cattegories_bins["et_eta_2_gap_njetN_dm0"] = 17;
        cattegories_bins["et_eta_2_gap_njetN_dm1"] = 18;
        //endcap
        cattegories_bins["et_eta_2_endcap_real_njet0_alldm"] = 19;
        cattegories_bins["et_eta_2_endcap_real_njet0_dm0"] = 20;
        cattegories_bins["et_eta_2_endcap_real_njet0_dm1"] = 21;
        cattegories_bins["et_eta_2_endcap_real_njetN_alldm"] = 22;
        cattegories_bins["et_eta_2_endcap_real_njetN_dm0"] = 23;
        cattegories_bins["et_eta_2_endcap_real_njetN_dm1"] = 24;


    for(auto c: catteg)
    {
        dout(c, "->", cattegories_bins[c]);
        cattegories["et_13TeV"].push_back({cattegories_bins[c], c});
        binning[c] = _nojets;
    }
    dout("categzise", cattegories["et_13TeV"].size());


    // Create an empty CombineHarvester instance that will hold all of the
    // datacard configuration and histograms etc.
    dout("Create an empty CombineHarvester instance");
    ch::CombineHarvester cb;
    for(auto channel : channels)
    {
        dout(channel);
        cb.AddObservations(  // data_obs
            {"*"},  // masspoint
            {"htt"},  // analysis
            {"13TeV"},  // era
            {channel},
            cattegories[channel + "_13TeV"]  // bin
        );
        cb.AddProcesses(
            {"*"},
            {"htt"},
            {"13TeV"},
            {channel},
            bkg_procs[channel],
            cattegories[channel + "_13TeV"],
            false  // background
        );
        cb.AddProcesses(  // Example : TH1 <channel>/<proc>_<masspoint>
            {energy_scales},
            {"htt"},
            {"13TeV"},
            {channel},
            {"ZL"},
            cattegories[channel + "_13TeV"],
            true
        );
    }

    auto signal = Set2Vec(cb.cp().signals().SetFromProcs(std::mem_fn(&Process::process)));
    if (debug)
    {
        std::cout << "signal(" << signal.size() << ") : ";
        for (auto i = signal.begin(); i != signal.end(); ++i)
            std::cout << *i << ' ';
        std::cout << '\n';

        cb.PrintAll();
    }

    // Example path : TH1 et_nobtag_tight/ZTT_CMS_scale_t_et_13TeVUp

    // TES
    // cb.cp().process({"ZTT", "TTT", "VVT"}).AddSyst(cb,"CMS_scale_t_3prong_13TeV", "shape", SystMap<>::init(1.00));
    // cb.cp().process({"ZTT", "TTT", "VVT"}).AddSyst(cb,"CMS_scale_t_1prong_13TeV", "shape", SystMap<>::init(1.00));
    // cb.cp().process({"ZTT", "TTT", "VVT"}).AddSyst(cb,"CMS_scale_t_1prong1pizero_13TeV", "shape", SystMap<>::init(1.00));

    /* missing shapes unc.
    db.add_shape_systematic("CMS_htt_ttbarShape_13TeV", 1.0, channels,["TTT", "TTJ"])
    db.add_shape_systematic("CMS_htt_jetToTauFake_13TeV", 1.0, channels,["ZJ", "W", "TTJ", "VVJ"])
    db.add_shape_systematic("CMS_eFakeTau_1prong_13TeV", 1.0, "et", ["ZL"])
    db.add_shape_systematic("CMS_eFakeTau_1prong1pizero_13TeV", 1.0, "et", ["ZL"])
    db.add_shape_systematic("CMS_ZLShape_et_1prong_13TeV", 1.0, "et", ["ZL"])
    db.add_shape_systematic("CMS_ZLShape_et_1prong1pizero_13TeV", 1.0, "et", ["ZL"])
    */
    // ?  CMS_eff_trigger_et: ZJ, TTJ, VVL, VVK, W, QCD
    if (!newera) // 2016
    {
        cb.cp().process({"ZL"}).AddSyst(cb,"CMS_eFakeTau_13TeV", "lnN",SystMap<>::init(1.30)); // tight wp of anti-e discrim
        cb.cp().process({"ZTT", "ZL", "ZJ", "TTT", "TTL", "TTJ", "VVT", "VVL", "VVJ", "W"}).AddSyst(cb, "CMS_eff_e_13TeV", "lnN", SystMap<>::init(1.02));
        cb.cp().process({"ZTT", "ZL", "ZJ", "TTT", "TTL", "TTJ", "VVT", "VVL", "VVJ", "W"}).AddSyst(cb, "lumi_13TeV", "lnN", SystMap<>::init(1.062));
        cb.cp().process({"ZTT", "ZL", "ZJ", "TTT", "TTL", "TTJ", "VVT", "VVL", "VVJ", "W"}).AddSyst(cb, "CMS_eff_trigger_et_13TeV", "lnN", SystMap<>::init(1.02));
        cb.cp().process({"ZTT", "ZL",       "TTT",               "VVT"}).AddSyst(cb, "CMS_eff_t_et_13TeV", "lnN", SystMap<>::init(1.08));

        cb.cp().process({"ZL",  "ZJ", "ZTT"}).AddSyst(cb, "CMS_htt_zttNorm_13TeV", "lnN", SystMap<>::init(1.20));
        cb.cp().process({"VVT", "VVJ", "VVL"}).AddSyst(cb, "CMS_htt_VVNorm_13TeV", "lnN", SystMap<>::init(1.10));
        cb.cp().process({"TTT", "TTJ", "TTL"}).AddSyst(cb, "CMS_htt_TTNorm_13TeV", "lnN", SystMap<>::init(1.10));
        cb.cp().process({"W"}).AddSyst(cb, "CMS_htt_WNorm_13TeV", "lnN", SystMap<>::init(1.10));

        cb.cp().process({"QCD", "W"}).AddSyst(cb, "CMS_ExtrapSSOS_13TeV", "lnN", SystMap<>::init(1.10));
    }
    else // 2017 : https://github.com/KIT-CMS/CombineHarvester/blob/SMHTT2017-dev/HTTSM2017/src/HttSystematics_SMRun2.cc#L144-L153
    {
        cb.cp().process({"ZL"}).AddSyst(cb,"CMS_eFakeTau_13TeV", "lnN",SystMap<>::init(1.16));  // https://github.com/KIT-CMS/CombineHarvester/blame/ce6a51be625aca1fa270b9f67c5dd0758ebf9262/HTTSM2017/src/HttSystematics_SMRun2.cc#L553
        cb.cp().process({"ZTT", "ZJ", "ZL", "TTT", "TTL", "TTJ", "VVT", "VVL", "VVJ", "W", "QCD"}).AddSyst(cb, "CMS_eff_e_13TeV", "lnN", SystMap<>::init(1.02));  // https://github.com/KIT-CMS/CombineHarvester/blame/62ed47b40ca2c0a69c22abb9eb6f99e3b0e62ebd/HTTSM2017/src/HttSystematics_SMRun2.cc#L142
        cb.cp().process({"ZTT", "ZJ", "ZL", "TTT", "TTL", "TTJ", "VVT", "VVL", "VVJ", "W", "QCD"}).AddSyst(cb, "lumi_13TeV", "lnN", SystMap<>::init(1.023));  // https://github.com/KIT-CMS/CombineHarvester/blame/29942200e5df5aef905335f8e1b52b4d22b31e68/HTTSM2017/src/HttSystematics_SMRun2.cc#L74
        cb.cp().process({"ZTT", "ZL", "ZJ", "TTT", "TTL", "TTJ", "VVT", "VVL", "VVJ", "W", "QCD"}).AddSyst(cb, "CMS_eff_trigger_et_13TeV", "lnN", SystMap<>::init(1.02));  // https://github.com/KIT-CMS/CombineHarvester/blame/ad54b4d0d626cb14d988088c5efaf90922e2e079/HTTSM2017/src/HttSystematics_SMRun2.cc#L91
        cb.cp().process({"ZTT", "ZL",       "TTT",               "VVT",                         }).AddSyst(cb, "CMS_eff_t_13TeV", "lnN", SystMap<>::init(1.045));  // https://github.com/KIT-CMS/CombineHarvester/blame/62ed47b40ca2c0a69c22abb9eb6f99e3b0e62ebd/HTTSM2017/src/HttSystematics_SMRun2.cc#L154
        cb.cp().process({"ZTT", "ZL",       "TTT",               "VVT",                         }).AddSyst(cb, "CMS_eff_t_et_13TeV", "lnN", SystMap<>::init(1.02));  // https://github.com/KIT-CMS/CombineHarvester/blame/62ed47b40ca2c0a69c22abb9eb6f99e3b0e62ebd/HTTSM2017/src/HttSystematics_SMRun2.cc#L159

    if (fakefactors)
    {
        // https://github.com/KIT-CMS/CombineHarvester/blame/ce6a51be625aca1fa270b9f67c5dd0758ebf9262/HTTSM2017/src/HttSystematics_SMRun2.cc
        cb.cp().channel({"et"}).process({"jetFakes"}).AddSyst(cb, "CMS_ff_qcd_dm0_njet0_$CHANNEL_stat_$ERA", "shape", SystMap<>::init(1.00));
        cb.cp().channel({"et"}).process({"jetFakes"}).AddSyst(cb, "CMS_ff_qcd_dm0_njet1_$CHANNEL_stat_$ERA", "shape", SystMap<>::init(1.00));
        // W shape stat.
        cb.cp().channel({"et"}).process({"jetFakes"}).AddSyst(cb, "CMS_ff_w_dm0_njet0_$CHANNEL_stat_$ERA", "shape", SystMap<>::init(1.00));
        cb.cp().channel({"et"}).process({"jetFakes"}).AddSyst(cb, "CMS_ff_w_dm0_njet1_$CHANNEL_stat_$ERA", "shape", SystMap<>::init(1.00));
        // TT shape stat.
        cb.cp().channel({"et"}).process({"jetFakes"}).AddSyst(cb, "CMS_ff_tt_dm0_njet0_stat_$ERA", "shape", SystMap<>::init(1.00));
        cb.cp().channel({"et"}).process({"jetFakes"}).AddSyst(cb, "CMS_ff_tt_dm0_njet1_stat_$ERA", "shape", SystMap<>::init(1.00));
        // Shape syst. of different contributions (QCD/W/tt)
        cb.cp().channel({"et"}).process({"jetFakes"}).AddSyst(cb, "CMS_ff_qcd_$CHANNEL_syst_$ERA", "shape", SystMap<>::init(1.00));
        cb.cp().channel({"et"}).process({"jetFakes"}).AddSyst(cb, "CMS_ff_w_syst_$ERA", "shape", SystMap<>::init(1.00));
        cb.cp().channel({"et"}).process({"jetFakes"}).AddSyst(cb, "CMS_ff_tt_syst_$ERA", "shape", SystMap<>::init(1.00));
        // Shape syst. in tautau due to using QCD FF also for W / tt
        // PASS

        // https://github.com/KIT-CMS/CombineHarvester/blame/79a6f18c55eefe19c9021029fa504831071a326b/HTTSM2017/src/HttSystematics_SMRun2.cc#L781
        // Stat. norm (uncorrelated across years)
        cb.cp().channel({"et"}).process({"jetFakes"}).AddSyst(cb, "CMS_ff_norm_stat_$CHANNEL_$BIN_$ERA", "lnN", SystMap<>::init(1.05));
        // Syst. norm: Bin-correlated
        cb.cp().channel({"et"}).process({"jetFakes"}).AddSyst(cb, "CMS_ff_norm_syst_$CHANNEL_$ERA", "lnN", SystMap<>::init(1.05));
        // Syst. norm: Bin-dependent ?: what are the bins?
        cb.cp().channel({"et"}).process({"jetFakes"}).AddSyst(cb, "CMS_ff_sub_syst_$CHANNEL_$BIN_$ERA", "lnN", SystMap<>::init(1.05));
    }

    // Cross-sections
    // https://github.com/KIT-CMS/CombineHarvester/blame/29942200e5df5aef905335f8e1b52b4d22b31e68/HTTSM2017/src/HttSystematics_SMRun2.cc#L687
    cb.cp().process({"ZTT", "ZL", "ZJ"}).AddSyst(cb, "CMS_htt_zjXsec_13TeV", "lnN", SystMap<>::init(1.04));
    cb.cp().process({"VVT", "VVJ", "VVL", "VV", "ST"}).AddSyst(cb, "CMS_htt_vvXsec_13TeV", "lnN", SystMap<>::init(1.05));
    cb.cp().process({"TTT", "TTL", "TTJ", "TT"}).AddSyst(cb, "CMS_htt_tjXsec_13TeV", "lnN", SystMap<>::init(1.06));
    cb.cp().process({"W"}).AddSyst(cb, "CMS_htt_wjXsec_13TeV", "lnN", SystMap<>::init(1.04));

    // Method-dependant
    cb.cp().process({"QCD"}).AddSyst(cb, "CMS_ExtrapSSOS_13TeV", "lnN", SystMap<>::init(1.05));
    //
    // recoil correction : https://github.com/KIT-CMS/CombineHarvester/commit/78bd81120ffde206febbe46eeeae74d835cb409d#diff-28d11c464d5b138c4ab8c83502ca62ef

    if (embedding)  // https://github.com/cms-analysis/CombineHarvester/blob/SMHTT2017-dev/HTTSM2017/src/HttSystematics_SMRun2.cc#L505
    {
        // uncertainty of hadronic tau track efficiency correction
        // cb.cp().process({"ZL", "TTT", "VVT", "EWKZ",  "ZJ", "W", "TTJ", "VVJ", "QCD"}).AddSyst(cb,"CMS_3ProngEff_13TeV", "shape", SystMap<>::init(1.00));
        // cb.cp().process({"ZL", "TTT", "VVT", "EWKZ",  "ZJ", "W", "TTJ", "VVJ", "QCD"}).AddSyst(cb,"CMS_1ProngPi0Eff_13TeV", "shape", SystMap<>::init(1.00));
        // cb.cp().process({}).AddSyst(cb,"", "shape", SystMap<>::init(1.00));
        cb.cp().process({"EMB"}).AddSyst(cb, "CMS_eff_trigger_et_2017", "lnN", SystMap<>::init(1.02));
        cb.cp().process({"EMB"}).AddSyst(cb, "CMS_eff_trigger_emb_et", "lnN", SystMap<>::init(1.02));

        // if args.embedding:
        // db.add_shape_systematic("CMS_htt_dyShape_13TeV", 1.0, channels, ["ZL", "ZJ"])
        // db.add_shape_systematic("CMS_htt_mistag_b_13TeV", 1.0, channels, ["ZL", "TTT", "VVT", "EWKZ" + "ZJ", "W", "TTJ", "VVJ", "QCD")
        // db.add_shape_systematic( "CMS_scale_j_13TeV", 1.0, channels, ["ZL", "TTT", "VVT", "EWKZ"] +  "ZJ", "W", "TTJ", "VVJ", "QCD")
        // db.add_shape_systematic( "CMS_scale_met_unclustered_13TeV", 1.0, channels, ["ZL", "TTT", "VVT", "EWKZ"] +  "ZJ", "W", "TTJ", "VVJ", "QCD")
    }

    if (decorelate_emb)
    {
        // https://github.com/KIT-CMS/CombineHarvester/blob/29942200e5df5aef905335f8e1b52b4d22b31e68/HTTSM2017/bin/MorphingSM2017.cpp
    }
    }

    //! [part7]
    for (auto channel: channels)
    {
        dout("Extracting shapes for channel ", channel);
        dout("Backgrounds file...", input_file);
        dout(input_file);

        // Problem: TH1 et_nojets_alldm/data_obs not found in /nfs/dust/cms/user/glusheno/Combine/MSSM-Full-2016/KIT/htt_et.inputs-mssm-13TeV-mttot.root
        cb.cp().channel({channel}).backgrounds().ExtractShapes(
        input_file,
        // input_dir + "htt_" + channel + ".inputs-mssm-13TeV" + postfix + ".root",
        "$BIN/$PROCESS",
        "$BIN/$PROCESS_$SYSTEMATIC");

        dout("process...");
        cb.cp().channel({channel}).process({"ZL"}).ExtractShapes(
        input_file,
        // input_dir + "htt_" + channel + ".inputs-mssm-13TeV" + postfix + ".root",
        "$BIN/ZL_$MASS",
        "$BIN/ZL_$MASS_$SYSTEMATIC");
    }
    if (debug) cb.PrintAll();
    // exit(1);
    // {
    // auto bins = cb.bin_set();
    //     for (auto b : bins)
    //     {
    //       auto procs = cb.cp().bin({b}).signals().process_set();
    //       for (auto p : procs)
    //       {
    //         dout(p);
    //         dout(cb.cp().bin({b}).process({"ZL"}).mass({"0"}).GetRate());
    //       }};exit(1);



    //   }
    dout("bin_set...");
    auto bins = cb.cp().bin_set();
    bool manual_rebin = true;
    bool auto_rebin = false;
    std::map<std::string, TH1F> before_rebin;
    std::map<std::string, TH1F> after_rebin;
    std::map<std::string, TH1F> after_rebin_neg;
    auto rebin = ch::AutoRebin()
    .SetBinThreshold(0.)
    .SetBinUncertFraction(0.9)
    .SetRebinMode(1)
    .SetPerformRebin(true)
    .SetVerbosity(1);
    if (auto_rebin) rebin.Rebin(cb, cb);

    if (manual_rebin)
    {
        for(auto b : bins)
        {
            std::cout << "Rebinning by hand for bin: " << b <<  std::endl;
            cb.cp().bin({b}).VariableRebin(binning[b]);
        }
    }

    cout << "Generating bbb uncertainties...";
    auto bbb = ch::BinByBinFactory()
    .SetAddThreshold(0.)
    .SetMergeThreshold(0.4)
    .SetFixNorm(false);
    bbb.MergeAndAdd(cb.cp().process({"ZTT", "ZJ", "TTT", "TTL", "TTJ", "VVT", "VVL", "VVJ", "W", "QCD", "EMB", "jetFakes"}), cb);
    // auto bbbsig = ch::BinByBinFactory()
    //   .SetAddThreshold(0.)
    //   .SetMergeThreshold(0.4)
    //   .SetFixNorm(false);
    // bbbsig.MergeAndAdd(cb.cp().process({"QCD"}), cb);
    cout << " done\n";

    // This function modifies every entry to have a standardised bin name of
    // the form: {analysis}_{channel}_{bin_id}_{era}
    // which is commonly used in the htt analyses
    ch::SetStandardBinNames(cb);

    //! [part9]
    dout("First we generate a set of bin names:");
    RooWorkspace ws("htt", "htt");
    string demo_file = output_dir + "htt_mssm_demo.root";
    TFile demo(demo_file.c_str(), "RECREATE");
    bool do_morphing = true;
    if (do_morphing)
    {
        auto bins = cb.bin_set();
        for (auto b : bins)
        {
            auto procs = cb.cp().bin({b}).signals().process_set();
            for (auto p : procs)
            {
                // ws, cb, std::string const& bin, std::string const& process, RooAbsReal& mass_var, std::string norm_postfix, bool allow_morph, bool verbose, bool force_template_limit, TFile * file) {
                ch::BuildRooMorphing(ws, cb, b, p, faketaues, "norm", true, true, false, NULL);
            }
        }
    }

    ws.var("faketaues")->setVal(0.0);
    // ws.var("CMS_th1x_htt_et_2_13TeV")->setVal(0.0);

    demo.Close();
    dout("AddWorkspace:");
    cb.AddWorkspace(ws);
    dout("ExtractPdfs:");
    cb.cp().process({"ZL"}).ExtractPdfs(cb, "htt", "$BIN_$PROCESS_morph");
    dout("cb.PrintAll():");
    cb.PrintAll();

    string foldercmb = output_dir + "cmb";
    boost::filesystem::create_directories(foldercmb);

    //Write out datacards. Naming convention important for rest of workflow. We
    //make one directory per channel-cat, one per channel and cmb. In this code we only
    //store the individual datacards for each directory to be combined later, but
    //note that it's also possible to write out the full combined card with CH

    cout << "Writing datacards ...";

    //Individual channel-cattegories
    for (string channel : channels)
    {
        string folderchn = output_dir + channel;
        auto bins = cb.cp().channel({channel}).bin_set();
        for (auto b : bins)
        {
            dout("±±±±±±±±±±±±±±±±±±±±±±±±± \n\t channel", channel, "bin", b);
            string folderchncat = output_dir + b;
            boost::filesystem::create_directories(folderchn);
            boost::filesystem::create_directories(folderchncat);

            TFile output((foldercmb + "/" + b + "_input.root").c_str(), "RECREATE");
            TFile outputchn((folderchn + "/" + b + "_input.root").c_str(), "RECREATE");
            TFile outputchncat((folderchncat + "/" + b + "_input.root").c_str(), "RECREATE");
            cb.cp().channel({channel}).bin({b}).mass({"*"}).PrintAll();
            cb.cp().channel({channel}).bin({b}).mass({"*"}).WriteDatacard(folderchn + "/" + b + ".txt", outputchn);
            cb.cp().channel({channel}).bin({b}).mass({"*"}).WriteDatacard(folderchncat + "/" + b + ".txt", outputchncat);
            cb.cp().channel({channel}).bin({b}).mass({"*"}).WriteDatacard(foldercmb + "/" + b + ".txt", output);

            output.Close();
            outputchn.Close();
            outputchncat.Close();
        }
    }
    dout("done");
}
