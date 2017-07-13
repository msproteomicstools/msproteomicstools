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

  // copy constructor
  c_peakgroup( const c_peakgroup &rhs)
  {
    // std::cout << " copy costr: c_peakgroup" << std::endl;
    fdr_score   = rhs.fdr_score;
    normalized_retentiontime = rhs.normalized_retentiontime;
    internal_id_= rhs.internal_id_;  
    intensity_  =  rhs.intensity_;    
    dscore_     =  rhs.dscore_;       
    cluster_id_ =  rhs.cluster_id_;   
    precursor   =  rhs.precursor;     
  }

  c_peakgroup& operator = (const c_peakgroup &rhs)
  {
    // std::cout << " assignment operator : c_peakgroup" << std::endl;
    fdr_score   = rhs.fdr_score;
    normalized_retentiontime = rhs.normalized_retentiontime;
    internal_id_= rhs.internal_id_;  
    intensity_  =  rhs.intensity_;    
    dscore_     =  rhs.dscore_;       
    cluster_id_ =  rhs.cluster_id_;   
    precursor   =  rhs.precursor;     
  }


};

#endif
