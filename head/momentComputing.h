#ifndef FIELDCOMPUTING_H_INCLUDED
#define FIELDCOMPUTING_H_INCLUDED

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <typeinfo>
#include <vector>
#include <algorithm>

// Does dot product between 2 vectors of equal length
// param a first vector
// param b second vector
// return dot product of a and b
template <typename T>
T innerProduct
(
    const std::vector<T> &vec1,
    const std::vector<T> &vec2
)
{
    if(vec1.size() != vec2.size()) throw std::runtime_error("vector size mismatch");
    T result = 0.0;
    auto it_vec2 = begin(vec2);
    for(auto elem_vec1 : vec1) result += elem_vec1 * (*it_vec2++);
    return result;
}

// Does multiplication between matrix and vector of equal length
// param a first matrix
// param b second vector
// return vector product of a and b
template <typename T, typename U>
U MVmultiplyProduct
(
    const T &matrix1,
    const U &vector2
)
{
    if(matrix1.at(0).size() != vector2.size()) throw std::runtime_error("matrix and vector size mismatch");
    //auto e = matrix1.at(0).size();
    auto e = vector2.size();
    U result(e, 0);
    for (auto i = 0; i < e; ++i) result[i] = innerProduct(matrix1.at(i), vector2);
    return result;
}

// Calculates zeroth moment of a node based on formula in LBIntro
// param node distribution function node containing nine discrete velocity vectors
// return zeroth moment of node. It actually is the pressure.
template <typename T>
T getZerothMoment
(
    const std::vector<T> &df
)
{
  T result = 0.0;
  for(auto i : df) result += i;
  return result;
}

// Calculates first moment of a node based on formula in LBIntro
// param node distribution function node containing nine discrete velocity vectors
// return first moment of node. 2D should have dim_e = 2
template <typename T>
T getFirstMoment
(
    const T &df,
    const std::vector<T> &e
)
{
    auto dim_e = e.at(0).size();
    T result(dim_e, 0);
    for(auto d = 0u; d < dim_e; ++d)
    {
        auto it_df = begin(df);
        for(auto dir : e) result[d] += (*it_df++) * dir[d];
    }  // d
    return result;
}

// Checks if velocity (u, v) is in steady state based on the following formula
// |u_curr - u_prev|/|u_curr| <= tol && |v_curr - v_prev|/|v_curr| <= tol
// This function knows that the model is 2D. "Dry" boundary nodes such as full-
// way bounceback nodes need to be excluded from the velocity vectors being passed in.
// param u_prev velocity lattice from the previous time step
// param u_curr velocity lattice from the current time step
// param tolerance absolute tolerance value for steady state checking
// return TRUE steady state reached
//        FALSE steady state not reached
template <typename T, typename U>
bool checkSteadyState
(
    const T &u_prev,
    const T &u_curr,
    U tolerance
)
{
    const auto nn = u_prev.size();
    const auto ii = u_prev.at(0).size();
    std::vector<U> diff_sum(ii, 0);
    std::vector<U> sum(ii, 0);
    std::vector<U> error(ii, 0);
    for (auto n = 0u; n < nn; ++n)
    {
        for (auto i = 0u; i < ii; ++i)
        {
            diff_sum[i] += fabs(u_curr[n][i] - u_prev[n][i]);
            sum[i] += fabs(u_curr[n][i]);
            error[i] = diff_sum[i] / sum[i];
        }
    }
    return (*std::max_element(std::begin(error), std::end(error)) < tolerance);
}

// Checks if velocity (u, v) is in steady state based on the following formula
// |u_curr - u_prev|/|u_curr| <= tol && |v_curr - v_prev|/|v_curr| <= tol
// This function knows that the model is 2D. "Dry" boundary nodes such as full-
// way bounceback nodes need to be excluded from the velocity vectors being passed in.
// param u_prev velocity lattice from the previous time step
// param u_curr velocity lattice from the current time step
// param tolerance absolute tolerance value for steady state checking
// return TRUE steady state reached
//        FALSE steady state not reached
template <typename T>
double checkError
(
    const T &u_prev,
    const T &u_curr
)
{
    const auto nn = u_prev.size();
    const auto ii = u_prev.at(0).size();
    std::vector<double> diff_sum(ii, 0.0);
    std::vector<double> sum(ii, 0.0);
    std::vector<double> error(ii, 0.0);
    for (auto n = 0u; n < nn; ++n)
    {
        for (auto i = 0u; i < ii; ++i)
        {
            diff_sum[i] += fabs(u_curr[n][i] - u_prev[n][i]);
            sum[i] += fabs(u_curr[n][i]);
            error[i] = diff_sum[i] / sum[i];
        }
    }
    return *std::max_element(std::begin(error), std::end(error));
}

#endif // FIELDCOMPUTING_H_INCLUDED
