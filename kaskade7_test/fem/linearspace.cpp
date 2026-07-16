#include <boost/fusion/container/generation/make_vector.hpp>
#include <boost/fusion/include/make_vector.hpp>


#include "fem/linearspace.hh"

using namespace Kaskade;

#ifdef UNITTEST

int main()
{
  using namespace boost::fusion;
  using Vector = Dune::FieldVector<double,2>;
  using Seq = vector<Vector>;
  LinearProductSpace<double,Seq> xs(make_vector(Vector(2)));
  LinearProductSpace<double,Seq> ys(make_vector(Vector(3)));

  auto sp = xs*ys;
  std::cout << "scalar product is " << sp << "\n";

  return 0;
};

#endif
