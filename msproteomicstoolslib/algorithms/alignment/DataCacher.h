
#include <string>
#include <vector>
#include <stdexcept>

#include <iostream>

struct c_data_cache {

  std::vector< std::vector<double> > fdr_;
  std::vector< std::vector<double> > rt_;

public:

  c_data_cache() {};
  void addValues(std::vector<double> & fdr, std::vector<double> & rt)
  {
    fdr_.push_back(fdr);
    rt_.push_back(rt);
  }
  void retrieveValues(std::vector<double> & rt1, std::vector<double> & rt2, int run1, int run2)
  {
    for (size_t k = 0; k < fdr_.size(); k++)
    {
      if (fdr_[k][run1] >= 0 && fdr_[k][run2] >= 0)
      {
        rt1.push_back( rt_[k][run1] );
        rt2.push_back( rt_[k][run2] );
      }
    }
  }
};
