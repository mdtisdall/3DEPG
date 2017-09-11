#ifndef SpinStates_h
#define SpinStates_h

#include <map>
#include <tuple>
#include <vector>

#include "blas_local.h"

template <typename index_type, typename value_type>
class SpinStates {
  public:
    typedef std::tuple<index_type, index_type, index_type> StateIndex;

    //struct StateIndex {
    //  index_type z, y, x;
    //};
   
    struct SpinState {
      value_type fPlus, fMinus, z; 
    };

    SpinStates() :
      curStates(0),
      maxStates(1024),
      stateValues(maxStates * 3)
    {
      SpinState initState = {0, 0, 1};
      StateIndex initIndex = {0, 0, 0}; 

      insertState(initIndex, initState);
    }

  protected:
    void insertState(const StateIndex &index, const SpinState &state) {
      if(curStates == maxStates) {
        expandStates();
      }

      stateValues[curStates] = state.fPlus;
      stateValues[curStates + maxStates] = state.fMinus;
      stateValues[curStates + 2 * maxStates] = state.z;

      stateMap[index] = curStates;
    }

    void expandStates() {
      size_t newMaxStates = maxStates * 2;
      std::vector<value_type> newStateValues(newMaxStates * 3);
     
      BLAS::copy(maxStates,
          stateValues.data(), 1,
          newStateValues.data(), 1);

      BLAS::copy(maxStates,
          stateValues.data() + maxStates, 1,
          newStateValues.data() + newMaxStates, 1);

      BLAS::copy(maxStates,
          stateValues.data() + 2 * maxStates, 1,
          newStateValues.data() + 2 * newMaxStates, 1);

      maxStates = newMaxStates;
      stateValues = newStateValues;
    }

    typedef std::map<StateIndex, size_t> StateMap;

    StateMap stateMap;
    size_t curStates;
    size_t maxStates;
    
    std::vector<value_type> stateValues;

};

#endif

