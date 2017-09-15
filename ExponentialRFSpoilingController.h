#ifndef ExponentialRFSpoilingController_h
#define ExponentialRFSpoilingController_h

template <typename value_type>
class ExponentialRFSpoilingController {
  public:
    ExponentialRFSpoilingController(const value_type increment):
      baseIncrement(increment),
      curIncrement(0),
      curPhase(0)
    {}

    value_type nextPhase() {
      value_type ret = curPhase;
      curIncrement += baseIncrement;
      curIncrement = fmod(curIncrement, value_type(2) * value_type(M_PI));
      curPhase += curIncrement;
      curPhase = fmod(curPhase, value_type(2) * value_type(M_PI));

      return ret;
    }

  protected:
    value_type baseIncrement;
    value_type curIncrement;
    value_type curPhase;
};

#endif
