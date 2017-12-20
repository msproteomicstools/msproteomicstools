#include <string>
#include <vector>
#include <stdexcept>

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
