
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

class c_linear_interpolate {
    typedef std::vector<double>::iterator itT;
public:
  c_linear_interpolate() {}
  c_linear_interpolate(std::vector<double>& x, std::vector<double>& y, double abs_err)
  {
    if (x.size() < 2) {throw std::invalid_argument("Needs at least 2 values."); }
    if (x.size() != y.size()) {throw std::invalid_argument("Needs equal size for x any y."); }

    // remove duplicate entries
    x_.push_back(x[0]);
    y_.push_back(y[0]);
    for (size_t k = 1; k < x.size(); k++)
    {
      if (x[k] > x_.back() ) 
      {
        x_.push_back(x[k]);
        y_.push_back(y[k]);
      }
    }

    if (!isSorted()) {throw std::invalid_argument("Needs sorted arrays.");}
  }

  bool isSorted()
  {
    for (size_t k = 1; k < x_.size(); k++)
    {
      if( x_[k] < x_[k-1]) {return false;}
    }
    return true;
  }

  double predict(double xnew)
  {
    // binary search
    itT x_lo = std::lower_bound(x_.begin(), x_.end(), xnew);
    itT x_hi = x_lo;

    itT y_lo = y_.begin();
    std::iterator_traits< itT >::difference_type iterator_pos = std::distance((itT)x_.begin(), x_lo);
    std::advance(y_lo, iterator_pos);
    itT y_hi = y_lo;

    // At this point, all iterators should point to the same position, namely
    // the std::lower_bound which is the first value *higher* than the xnew
    // value. We now need to decrease the x_lo and y_low iterator by one in
    // order to make the point to the last value *lower* the xnew value.
    if (x_lo == x_.begin())
    {
      // Special case, x_lo and y_lo are already at the beginning (the xnew
      // value is lower than the first value we have, just linearly interpolate
      // from the first two values)
      x_hi++;
      y_hi++;
    }
    else if (x_lo == x_.end())
    {
      // Special case, x_lo and y_lo are at the end (the xnew value is larger
      // than the last value we have, just linearly interpolate from the last
      // two values)
      x_lo--; x_lo--;
      y_lo--; y_lo--;
      x_hi--;
      y_hi--;
    }
    else
    {
      // Normal case, we lower the x_lo and y_lo iterator to point to the last
      // value *lower* the xnew value.
      x_lo--;
      y_lo--;
    }

    // deal with duplicates
    // while (x_hi != x_.end() && (*x_hi- *x_lo) == 0)
    // {
    //   x_hi++;
    //   y_hi++;
    // }

    double slope = (*y_hi - *y_lo) / (*x_hi- *x_lo);
    double y_new = slope * (xnew - *x_lo) + *y_lo;
    return y_new;
  }

private:

  std::vector<double> x_;
  std::vector<double> y_;

};
