#include "CombineHarvester/MSSMvsSMRun2Legacy/interface/BinomialBinByBin.h"
#include <iostream>
#include <string>
#include <vector>
#include "boost/format.hpp"
#include "boost/lexical_cast.hpp"
#include "Math/QuantFunc.h"

namespace ch {

BinomialBinByBinFactory::BinomialBinByBinFactory() : BinByBinFactory() {}

void BinomialBinByBinFactory::AddBinomialBinByBin(CombineHarvester &src, CombineHarvester &dest) {
  std::vector<Process *> procs;
  src.ForEachProc([&](Process *p) { 
    procs.push_back(p);
  });
  for (unsigned i = 0; i < procs.size(); ++i) {
    if (!procs[i]->shape()) continue;
    TH1 const* h = procs[i]->shape();
    if (h->GetSumw2N() == 0) {
      std::cout << "Process " << procs[i]->process()
                << " does not continue the weights information needed for "
                   "valid errors, skipping\n";
      continue;
    }
    unsigned n_pop_bins = 0;
    for (int j = 1; j <= h->GetNbinsX(); ++j) {
      if (h->GetBinContent(j) > 0.0) ++n_pop_bins;
    }
    if (n_pop_bins <= 1 && BinomialBinByBinFactory::GetFixNorm()) {
      if (BinomialBinByBinFactory::GetVerbosity() >= 1) {
        std::cout << "Requested fixed_norm but template has <= 1 populated "
                     "bins, skipping\n";
        std::cout << Process::PrintHeader << *(procs[i]) << "\n";
      }
      continue;
    }
    for (int j = 1; j <= h->GetNbinsX(); ++j) {
      double val = h->GetBinContent(j);
      double err = std::sqrt(val*binomial_p_/binomial_n_*(1.-binomial_p_));

      ch::Systematic sys;
      ch::SetProperties(&sys, procs[i]);
      sys.set_type("shape");
      std::string name = BinomialBinByBinFactory::GetPattern();
      boost::replace_all(name, "$ANALYSIS", sys.analysis());
      boost::replace_all(name, "$CHANNEL", sys.channel());
      boost::replace_all(name, "$BIN", sys.bin());
      boost::replace_all(name, "$BINID", boost::lexical_cast<std::string>(sys.bin_id()));
      boost::replace_all(name, "$ERA", sys.era());
      boost::replace_all(name, "$PROCESS", sys.process());
      boost::replace_all(name, "$MASS", sys.mass());
      boost::replace_all(name, "$#", boost::lexical_cast<std::string>(j));
      sys.set_name(name);
      sys.set_asymm(true);
      std::unique_ptr<TH1> h_d(static_cast<TH1 *>(h->Clone()));
      std::unique_ptr<TH1> h_u(static_cast<TH1 *>(h->Clone()));
      h_d->SetBinContent(j, val - err);
      if (h_d->GetBinContent(j) < 0.) h_d->SetBinContent(j, 0.);
      if (!(h_d->Integral() > 0.)) h_d->SetBinContent(j,0.00001*h->Integral());
      h_u->SetBinContent(j, val + err);
      if (BinomialBinByBinFactory::GetFixNorm()) {
        sys.set_value_d(1.0);
        sys.set_value_u(1.0);
      } else {
        sys.set_value_d(h_d->Integral()/h->Integral());
        sys.set_value_u(h_u->Integral()/h->Integral());
      }
      sys.set_shapes(std::move(h_u), std::move(h_d), nullptr);
      dest.CreateParameterIfEmpty(sys.name());
      dest.InsertSystematic(sys);
    }
  }
}

}
