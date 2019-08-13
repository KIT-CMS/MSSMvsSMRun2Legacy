#ifndef MSSMvsSMRun2Legacy_BinomialBinByBin_h
#define MSSMvsSMRun2Legacy_BinomialBinByBin_h
#include <string>
#include "CombineHarvester/CombineTools/interface/CombineHarvester.h"
#include "CombineHarvester/CombineTools/interface/BinByBin.h"

namespace ch {
/**
 * Merges bin uncertainties and creates binomial bin-by-bin statistical uncertainties
 *
 * Typical usage:
 * 
 *     auto binomialbbb = ch::BinomialBinByBinFactory()
 *         .SetAddThreshold(0.1)
 *         .SetMergeThreshold(0.5)
 *         .SetFixNorm(true);
 *     bbb.MergeBinErrors(cb.cp().backgrounds());
 *     bbb.AddBinomialBinByBin(cb.cp().backgrounds(), cb);
 * 
 * See below for details on each class method.
 */
class BinomialBinByBinFactory: public BinByBinFactory {
 public:
  BinomialBinByBinFactory();

  /**
   * Create binomial bin-by-bin shape uncertainties for every process in **src**, and
   * add these to **dest**
   *
   * The behaviour of this function is controlled by three parameters:
   * 
   *   * The first is the uncertainty threshold which determines whether a 
   *     systematic will be created for a given histogram bin
   *     (\ref SetAddThreshold). For a bin with content \f$x\f$ and bin error
   *     \f$e\f$, this threshold is thefractional error \f$e/x\f$.
   *   * The second parameters (\ref SetFixNorm) is a flag that when set to
   *     `true` will re-normalize the up and down shape templates to have the
   *     same integral as the nominal. This means the nuisance parameter for
   *     this bin is only able to change the shape of the distribution: if a
   *     bin is shifted up, all other bins are shifted down slightly to
   *     preserve the normalisation. This feature is useful when the
   *     statistical uncertainty in the process yield is considered as an
   *     independent nuisance parameter.
   *   * The third parameters is a pattern-string which prescribes how the
   *     bin-by-bin uncertainties should be named. By default this is
   *     `CMS_$ANALYSIS_$CHANNEL_$BIN_$ERA_$PROCESS_bin_$#`, but can be
   *     changed with \ref SetPattern. Note that the special term `$#`
   *     should always be included as this is replaced with the bin index.
   */
  void AddBinomialBinByBin(CombineHarvester &src, CombineHarvester &dest);

  /**
   * By default this class only produces output on the screen when an error
   * occurs, set to a value greater than zero for more verbose output
   */
  inline BinomialBinByBinFactory& SetVerbosity(unsigned verbosity) {
    return static_cast<BinomialBinByBinFactory&> (BinByBinFactory::SetVerbosity(verbosity));
  }

  /**
   * Set the fractional bin error threshold for bin-by-bin creation and
   * for participation in the merging algorithm
   */
  inline BinomialBinByBinFactory& SetAddThreshold(double val) {
    return static_cast<BinomialBinByBinFactory&> (BinByBinFactory::SetAddThreshold(val));
  }

  /**
   * The threshold for the merging algorithm
   */
  inline BinomialBinByBinFactory& SetMergeThreshold(double val) {
    return static_cast<BinomialBinByBinFactory&> (BinByBinFactory::SetMergeThreshold(val));
  }

  /**
   * The pattern-string for the systematic naming convention
   */
  inline BinomialBinByBinFactory& SetPattern(std::string const& pattern) {
    return static_cast<BinomialBinByBinFactory&> (BinByBinFactory::SetPattern(pattern));
  }

  /**
   * Whether or not the bin-by-bin systematics are allowed to vary the process
   * normalisation
   */
  inline BinomialBinByBinFactory& SetFixNorm(bool fix) {
    return static_cast<BinomialBinByBinFactory&> (BinByBinFactory::SetFixNorm(fix));
  }

  /**
   * Construct approximate Poisson uncertainties instead of default Gaussian
   */
  inline BinomialBinByBinFactory& SetPoissonErrors(bool poisson_errors) {
    return static_cast<BinomialBinByBinFactory&> (BinByBinFactory::SetPoissonErrors(poisson_errors));
  }

  /**
   * Set whether bins with zero content can participate in the merging procedure
   */
  inline BinomialBinByBinFactory& SetMergeZeroBins(bool merge) {
    return static_cast<BinomialBinByBinFactory&> (BinByBinFactory::SetMergeZeroBins(merge));
  }

  /**
   * Set whether bins with error >= content participate in the merging procedure
   */
  inline BinomialBinByBinFactory& SetMergeSaturatedBins(bool merge) {
    return static_cast<BinomialBinByBinFactory&> (BinByBinFactory::SetMergeSaturatedBins(merge));
  }

  /**
   * Set value of probability p for event passing in embedded kinematic filtering
   */
  inline BinomialBinByBinFactory& SetBinomialP(double val) {
    binomial_p_ = val;
    return *this;
  }

  /**
   * Set value of number of tries n for event passing in embedded kinematic filtering
   */
  inline BinomialBinByBinFactory& SetBinomialN(double val) {
    binomial_n_ = val;
    return *this;
  }


 private:
  double binomial_p_;
  double binomial_n_;
};
}

#endif
