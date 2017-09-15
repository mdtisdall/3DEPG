#ifndef ExponentialRFSpoilingController_h
#define ExponentialRFSpoilingController_h

template <typename value_type>
class ExponentialRFSpoilingController {
  public:
    ExponentialRFSpoilingController(const value_type increment):
      increment(increment),
      curIndex(0)
    {}

    value_type nextPhase() {
      value_type ret = increment * value_type(curIndex + 1) *
        value_type(curIndex) * value_type(0.5);

      curIndex++;

      return ret;
    }

  protected:
    value_type increment;
    size_t curIndex;
};

#endif
