/** \file ChiSquaredProb.h
    \brief Declaration for ChiSquaredProb class.
    \author Masaharu Hirayama, GSSC
*/
#ifndef periodSearch_ChiSquaredProb_h
#define periodSearch_ChiSquaredProb_h

#include <cmath>
#include <utility>

/** \class ChiSquaredProb
    \brief Functor for estimating upper/lower limits of the cumulative chi-squared probability by numerical intgration.
*/
class ChiSquaredProb {
  public:
    /** \brief Create ChiSquaredProb object.
        \param dof Number of degrees of freedom.
        \param min_pdf Smallest probability density function (pdf) to be computed.
    */
    ChiSquaredProb(int dof, double min_pdf=1.0e-99);

    /** \brief Compute lower/upper limit of chi squared probability.
        \param chisq The chi-squared value.
        \param precision Fractional precision of chi squared probability desired.
        \param max_iteration Maximum number of iterations in numerical integration.
    */
    std::pair<double,double> operator() (double chisq, double precision=1.0e-4, int max_iteration=1000000) const;

  private:
   // Degrees of freedom
    double m_dof;
    double m_dof_minus_2; // = m_dof - 2

    // Log of normalization of chisqare distribution.
    double m_lognorm;

    // Maximum chi-squared value to compute.
    double m_max_chisq;

    // Chi-squared probability distribution function.
    inline double pdf(double x) const {
      return exp(m_lognorm + (m_dof_minus_2 * log(x) - x) / 2.0);
    }

    // Function for computing step size.
    inline double eta(double x) const {
      return 2.0 / (1.0 - m_dof_minus_2 / x);
    }
};

#endif
