#ifndef ConstantGradSpoilingController_h
#define ConstantGradSpoilingController_h

template <typename GradType>
class ConstantGradSpoilingController {
  public:
    GradType nextGrad() {
      return GradType({0,0,1}); 
    }
};

#endif
