
#include <string>

struct c_peakgroup {

public:
  double fdr_score;
  double normalized_retentiontime;
  std::string internal_id_;
  double intensity_;
  double dscore_;
  int cluster_id_;

  c_peakgroup() {};
};

