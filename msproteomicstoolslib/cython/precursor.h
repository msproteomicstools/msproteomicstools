#ifndef MSPROTEOMICSTOOLSLIB_PRECURSOR_H
#define MSPROTEOMICSTOOLSLIB_PRECURSOR_H

#include <string>
#include <vector>
#include <stdexcept>

#include <iostream>

#include "peakgroup.h"

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

#endif
