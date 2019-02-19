from libcpp.vector cimport vector as libcpp_vector

cdef extern from "_linear_interpol.h":
    cdef cppclass c_linear_interpolate:
        c_linear_interpolate(libcpp_vector[double] & x, libcpp_vector[double] & y, double abs_err)
        double predict(double xnew)

