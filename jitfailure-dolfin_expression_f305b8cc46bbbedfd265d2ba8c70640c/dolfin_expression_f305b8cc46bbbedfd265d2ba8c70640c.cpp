
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
  class dolfin_expression_f305b8cc46bbbedfd265d2ba8c70640c : public Expression
  {
     public:
       double CMax;
double CMin;
double k;
std::shared_ptr<dolfin::GenericFunction> generic_function_cAdv;
double x0;


       dolfin_expression_f305b8cc46bbbedfd265d2ba8c70640c()
       {
            
       }

       void eval(Eigen::Ref<Eigen::VectorXd> values, Eigen::Ref<const Eigen::VectorXd> x) const override
       {
          double cAdv;
            generic_function_cAdv->eval(Eigen::Map<Eigen::Matrix<double, 1, 1>>(&cAdv), x);
          values[0] = (CMax-CMin)/(1+exp(IntIncl*(-cAdv+x0)))+CMin;

       }

       void set_property(std::string name, double _value) override
       {
          if (name == "CMax") { CMax = _value; return; }          if (name == "CMin") { CMin = _value; return; }          if (name == "k") { k = _value; return; }          if (name == "x0") { x0 = _value; return; }
       throw std::runtime_error("No such property");
       }

       double get_property(std::string name) const override
       {
          if (name == "CMax") return CMax;          if (name == "CMin") return CMin;          if (name == "k") return k;          if (name == "x0") return x0;
       throw std::runtime_error("No such property");
       return 0.0;
       }

       void set_generic_function(std::string name, std::shared_ptr<dolfin::GenericFunction> _value) override
       {
          if (name == "cAdv") { generic_function_cAdv = _value; return; }
       throw std::runtime_error("No such property");
       }

       std::shared_ptr<dolfin::GenericFunction> get_generic_function(std::string name) const override
       {
          if (name == "cAdv") return generic_function_cAdv;
       throw std::runtime_error("No such property");
       }

  };
}

extern "C" DLL_EXPORT dolfin::Expression * create_dolfin_expression_f305b8cc46bbbedfd265d2ba8c70640c()
{
  return new dolfin::dolfin_expression_f305b8cc46bbbedfd265d2ba8c70640c;
}

