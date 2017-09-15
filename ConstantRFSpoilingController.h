#ifndef ConstantRFSpoilingController_h
#define ConstantRFSpoilingController_h

template <typename value_type>
class ConstantRFSpoilingController {
  public:
    ConstantRFSpoilingController() {}

    value_type nextPhase() {
      return value_type(0);
    }
};

#endif
