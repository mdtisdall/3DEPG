#ifndef ResettingExponentialRFSpoilingController_h
#define ResettingExponentialRFSpoilingController_h

template <typename value_type>
class ResettingExponentialRFSpoilingController {
  public:
    ResettingExponentialRFSpoilingController(
      const value_type increment,
      const size_t runLength):
      increment(increment),
      runLength(runLength),
      curIndex(0)
    {}

    value_type nextPhase() {
      value_type ret = increment * value_type(curIndex + 1) *
        value_type(curIndex) * value_type(0.5);

      curIndex++;

      curIndex = curIndex % runLength;

      return ret;
    }

  protected:
    value_type increment;
    size_t runLength;
    size_t curIndex;
};

#endif
