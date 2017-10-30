#ifndef MSPROTEOMICSTOOLSLIB_PEAKGROUP_H
#define MSPROTEOMICSTOOLSLIB_PEAKGROUP_H

#include <string>
#include <vector>
#include <stdexcept>

#include <iostream>

// forward decl
struct c_precursor;

struct c_peakgroup {

public:
  double fdr_score;
  double normalized_retentiontime;
  std::string internal_id_;
  double intensity_;
  double dscore_;
  int cluster_id_;

  c_precursor* precursor;

  c_peakgroup() {};
  c_precursor* getPeptide() {return precursor;};

};

#endif
