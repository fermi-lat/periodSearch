/** \file ChiSquaredProb.cxx
    \brief Implmementation for ChiSquaredProb class.
    \author Masaharu Hirayama, GSSC
*/
#include <stdexcept>

#include "ChiSquaredProb.h"

namespace periodSearch {

  ChiSquaredProb::ChiSquaredProb(int dof, double min_pdf): m_dof(dof), m_dof_minus_2(dof-2), m_lognorm(0.0) {
    // check if dof and min_pdf are positive or not
    if (dof <= 0) throw std::logic_error("ChiSquaredProb::ChiSquaredProb: non-positive number of degrees of freedom");
    if (min_pdf <= 0.) throw std::logic_error("ChiSquaredProb::ChiSquaredProb: non-positive minimum number of probability density function");

    // pre-compute normalization
    double nterm;
    if (dof % 2) {
      m_lognorm = - 0.5 * log(M_PI);
      nterm = dof / 2;
    } else {
      m_lognorm = 0.0;
      nterm = dof / 2 - 2;
    }
    for (int ii=1; ii<=nterm; ii++) {
      m_lognorm -= log(m_dof / 2.0 - double(ii));
    }
    m_lognorm -= m_dof / 2.0 * M_LN2;

    // find the maximum chi-squared value to compute
    // NOTE: in the following search for m_max_chisq, the initial value
    // of x_cur must be > m_dof_minus_2 and > 0.
    double x_cur = m_dof + 1.0;
    if (pdf(x_cur) > min_pdf) {
      double log_min_pdf = (log(min_pdf) - m_lognorm) * 2.0;
      const double MY_EPSILON = 1.0e-10;
      const double min_diff = 2.0 * log(1.0 + MY_EPSILON);
      double log_diff;
      while (1) {
        log_diff = m_dof_minus_2 * log(x_cur) - x_cur - log_min_pdf;
        if (fabs(log_diff) > min_diff) x_cur += log_diff;
        else break;
      }
    }
    m_max_chisq = x_cur;
  }

  std::pair<double,double> ChiSquaredProb::operator() (double chisq, double precision, int max_iteration) const {
    // return 1.0 in trivial cases
    if (chisq <= 0.0) return std::make_pair(1.0, 1.0);

    // return the maximum residual if chisq >= m_max_chisq
    if (chisq >= m_max_chisq) {
      if (m_max_chisq > m_dof) {
        double max_residual = 2.0 * m_max_chisq * pdf(m_max_chisq)
  	/ (m_max_chisq - m_dof);
        return std::make_pair(0.0, (1.0 < max_residual ? 1.0 : max_residual));
      } else {
        return std::make_pair(0.0, 1.0);
      }
    }

    // prepare for computing step size
    double p_value = precision / 2.0;
    double x_offset = 0.0;
    if (m_dof > 2.1) {
      x_offset = 2.0 * p_value + sqrt(2.0 * p_value * m_dof_minus_2);
    }

    // set initial values for numerical integration
    double x_cur = chisq;
    double f_cur = pdf(x_cur);
    double Qmin = 0.0;
    double Qmax = 0.0;
    double Rn = 1.0;

    // helper variables for numerical integration
    double x_step, x_next, f_next;
    int num_iter = 0;

    // integration where pdf(x) increases with x
    if (m_dof > 2.1) {
      for (; (num_iter < max_iteration) && (x_cur < m_dof_minus_2)
  	   && (x_cur < m_max_chisq);
  	 num_iter++, x_cur = x_next, f_cur = f_next) {

        // next sampling point
        x_step = p_value * fabs(eta(x_cur));
        x_next = x_cur + x_step;
        if (x_next > m_dof_minus_2) {
  	x_next = m_dof_minus_2;
  	x_step = x_next - x_cur;
        }
        f_next = pdf(x_next);

        // accumulate estimators
        Qmax += f_next * x_step;
        Qmin += f_cur  * x_step;
      }
    }

    // integration where pdf(x) decreases with x
    for (; (num_iter < max_iteration) && (x_cur < m_max_chisq);
         num_iter++, x_cur = x_next, f_cur = f_next) {
      // next sampling point
      x_step = p_value * fabs(eta(x_cur + x_offset));
      x_next = x_cur + x_step;
      if (x_next > m_max_chisq) {
        x_next = m_max_chisq;
        x_step = x_next - x_cur;
      }
      f_next = pdf(x_next);

      // accumulate estimators
      Qmax += f_cur  * x_step;
      Qmin += f_next * x_step;

      if (x_next > m_dof) {
        // compute residual
        Rn = 2.0 * x_next * f_next / (x_next - m_dof);

        // terminate numerical integration
        if (Rn <= Qmin* p_value) break;
      }
    }

    // add residual (upper limit)
    Qmax += Rn;
    Qmax = (Qmax < 1.0 ? Qmax : 1.0);

    // return integration and upper limit for residual
    return std::make_pair(Qmin, Qmax);
  }

}
