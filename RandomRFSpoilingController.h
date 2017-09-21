#ifndef RandomRFSpoilingController_h
#define RandomRFSpoilingController_h

#include "mkl_vsl.h"

namespace RandomGenerators {
  static const int niederreiter = VSL_BRNG_NIEDERR;
  static const int sobol = VSL_BRNG_SOBOL;
  static const int pseudorandom = VSL_BRNG_SFMT19937;
}

template <typename value_type, int brng, size_t generatorSize = 1024>
class RandomRFSpoilingController {
  public:

    RandomRFSpoilingController():
      curPhase(0)
    {
      vslNewStream(&stream, brng, 1);
      vsRngUniform(VSL_RNG_METHOD_UNIFORM_STD,
        stream, generatorSize, angles.data(),
        value_type(0), value_type(2.0 * M_PI)
      );
    }

    value_type nextPhase() {
      if(generatorSize == curPhase) {
        vsRngUniform(VSL_RNG_METHOD_UNIFORM_STD,
          stream, generatorSize, angles.data(),
          value_type(0), value_type(2.0 * M_PI)
        );

        curPhase = 0;
      }

      value_type ret = angles[curPhase];

      curPhase++;

      return ret;
    }

  protected:
    size_t curPhase;
    VSLStreamStatePtr stream;
    std::array<value_type, generatorSize> angles;
};

#endif
