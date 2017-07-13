
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

struct c_precursor {

public:
  bool decoy;
  std::vector<c_peakgroup> peakgroups;
  std::string curr_id_;
  std::string protein_name_;
  std::string sequence_;
  std::string run_id_;
  std::string precursor_group_id; // need that?

  c_precursor() {};
  c_precursor(std::string id, std::string run_id) : curr_id_(id), run_id_(run_id) {}

  std::string getRunId() {return run_id_;}

  std::string get_id() {return curr_id_;}
  void add_peakgroup_tpl(c_peakgroup & pg, std::string tpl_id, int cluster_id=-1)
  {
    peakgroups.push_back(pg);
  }

  void setClusterID(std::string this_id, int cl_id)
  {
    _setClusterID(this_id, cl_id);
  }

  void _setClusterID(std::string this_id, int cl_id)
  {
    int nr_hit = 0;
    for (std::vector<c_peakgroup>::iterator it = peakgroups.begin(); it != peakgroups.end(); it++)
    {
      if (it->internal_id_ == this_id)
      {
        it->cluster_id_ = cl_id;
        nr_hit++;
      }
    }
    if (nr_hit != 1) {throw std::invalid_argument("Did not find pg with specified id."); }
  }

};

