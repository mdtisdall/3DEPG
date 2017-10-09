#ifndef FLASHTR_h
#define FLASHTR_h

template <typename SpinStates>
class FLASHTR {
  public:
    typedef typename SpinStates::value_type value_type;

    FLASHTR(
      const value_type flipAngle,
      const value_type flipPhase,
      const typename SpinStates::StateIndex spoilingGrad,
      const value_type tr,
      const value_type t1,
      const value_type t2,
      const value_type threshold = value_type(0)
      ) :
      excitation(flipAngle, flipPhase),
      spoiling(spoilingGrad, threshold),
      relaxation(tr, t1, t2),
      invPhase(std::polar(value_type(1.0),
        -flipPhase + value_type(M_PI) / value_type(2)))
    {}

    typename SpinStates::cvalue_type operator()(SpinStates *states) {
      excitation(states);
      typename SpinStates::cvalue_type ret =
        states->getState(typename SpinStates::StateIndex({{0,0,0}})).fPlus * 
        invPhase;
      relaxation(states); 
      spoiling(states); 

      return ret;
    }


  protected:
    typename SpinStates::Excitation excitation;
    typename SpinStates::Spoiling spoiling;
    typename SpinStates::Relaxation relaxation;
    typename SpinStates::cvalue_type invPhase;
};

#endif
