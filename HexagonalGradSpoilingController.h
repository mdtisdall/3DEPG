#ifndef HexagonalGradSpoilingController_h
#define HexagonalGradSpoilingController_h

#include <array>

template <typename GradType>
class HexagonalGradSpoilingController {
  public:
    HexagonalGradSpoilingController() :
      gradSchedule(),
      curIndex(0)
    {
      gradSchedule[0] = GradType({1, 0, 1});
      gradSchedule[1] = GradType({0, 1, 1});
      gradSchedule[2] = GradType({-1, 1, 1});
      gradSchedule[3] = GradType({-1, 0, 1});
      gradSchedule[4] = GradType({0, -1, 1});
      gradSchedule[5] = GradType({1, -1, 1});
    }

    GradType nextGrad() {
      GradType ret = gradSchedule[curIndex];
      curIndex++;
      curIndex = curIndex % 6;

      return ret; 
    }

  protected:
    std::array<GradType, 6> gradSchedule;
    size_t curIndex;
};

#endif
