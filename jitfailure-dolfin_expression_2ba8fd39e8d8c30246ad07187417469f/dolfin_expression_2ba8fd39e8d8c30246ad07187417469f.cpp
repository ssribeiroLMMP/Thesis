
// Based on https://gcc.gnu.org/wiki/Visibility
#if defined _WIN32 || defined __CYGWIN__
    #ifdef __GNUC__
        #define DLL_EXPORT __attribute__ ((dllexport))
    #else
        #define DLL_EXPORT __declspec(dllexport)
    #endif
#else
    #define DLL_EXPORT __attribute__ ((visibility ("default")))
#endif

#include <dolfin/function/Expression.h>
#include <dolfin/math/basic.h>
#include <Eigen/Dense>


// cmath functions
using std::cos;
using std::sin;
using std::tan;
using std::acos;
using std::asin;
using std::atan;
using std::atan2;
using std::cosh;
using std::sinh;
using std::tanh;
using std::exp;
using std::frexp;
using std::ldexp;
using std::log;
using std::log10;
using std::modf;
using std::pow;
using std::sqrt;
using std::ceil;
using std::fabs;
using std::floor;
using std::fmod;
using std::max;
using std::min;

const double pi = DOLFIN_PI;


namespace dolfin
{
  class dolfin_expression_2ba8fd39e8d8c30246ad07187417469f : public Expression
  {
     public:
       double etaInf;
double eta0;
double K;
double nPow;
double ts;
double eps;


       dolfin_expression_2ba8fd39e8d8c30246ad07187417469f()
       {
            
       }

       void eval(Eigen::Ref<Eigen::VectorXd> values, Eigen::Ref<const Eigen::VectorXd> x) const override
       {
          values[0] = (1 - exp(-eta0*(gammaDot)/tauY_t))*                                         (tauY_t/(gammaDot+eps) +                                         K*(pow((abs(gammaDot)+eps),nPow)/(gammaDot+eps)))                                         + etaInf;

       }

       void set_property(std::string name, double _value) override
       {
          if (name == "etaInf") { etaInf = _value; return; }          if (name == "eta0") { eta0 = _value; return; }          if (name == "K") { K = _value; return; }          if (name == "nPow") { nPow = _value; return; }          if (name == "ts") { ts = _value; return; }          if (name == "eps") { eps = _value; return; }
       throw std::runtime_error("No such property");
       }

       double get_property(std::string name) const override
       {
          if (name == "etaInf") return etaInf;          if (name == "eta0") return eta0;          if (name == "K") return K;          if (name == "nPow") return nPow;          if (name == "ts") return ts;          if (name == "eps") return eps;
       throw std::runtime_error("No such property");
       return 0.0;
       }

       void set_generic_function(std::string name, std::shared_ptr<dolfin::GenericFunction> _value) override
       {

       throw std::runtime_error("No such property");
       }

       std::shared_ptr<dolfin::GenericFunction> get_generic_function(std::string name) const override
       {

       throw std::runtime_error("No such property");
       }

  };
}

extern "C" DLL_EXPORT dolfin::Expression * create_dolfin_expression_2ba8fd39e8d8c30246ad07187417469f()
{
  return new dolfin::dolfin_expression_2ba8fd39e8d8c30246ad07187417469f;
}

