#ifndef SpinStates_h
#define SpinStates_h

#include <map>
#include <array>
#include <vector>

#include "blas_local.h"
template <typename index_type, typename _value_type>
class SpinStates {
  public:
    typedef _value_type value_type;
    typedef std::complex<value_type> cvalue_type;
    typedef std::array<index_type, 3> StateIndex;

    //struct StateIndex {
    //  index_type z, y, x;
    //};
   
    struct SpinState {
      cvalue_type fPlus, fMinus, z; 
    };

    SpinStates() :
      curStates(0),
      maxStates(4),
      stateValues(maxStates * 3),
      stateValuesBuffer(maxStates * 3)
    {
      SpinState initState = {0, 0, 1};
      StateIndex initIndex = {0, 0, 0}; 

      insertState(initIndex, initState);
    }

    SpinState getState(StateIndex index) {
      SpinState ret;

      typename StateMapT::const_iterator it = stateMap.find(index);
      
      if(stateMap.end() == it) {
        ret = {0,0,0}; 
      }
      else {
        size_t offset = it->second;

        ret.fPlus = stateValues[offset]; 
        ret.fMinus = stateValues[offset + maxStates]; 
        ret.z = stateValues[offset + 2 * maxStates]; 
      }

      return ret;
    }

    size_t size() {
      return curStates; 
    }
    
    class Excitation {
      public:
        Excitation(value_type flipAngle, value_type phase) :
          opMatrixTrans(9) {
          value_type cosAlpha = cos(flipAngle);
          value_type sinAlpha = sin(flipAngle);
          value_type alphaDiv2 = flipAngle/value_type(2);
          value_type cosAlphaDiv2 = cos(alphaDiv2);
          value_type cosAlphaDiv2Sq = cosAlphaDiv2 * cosAlphaDiv2;
          value_type sinAlphaDiv2 = sin(alphaDiv2);
          value_type sinAlphaDiv2Sq = sinAlphaDiv2 * sinAlphaDiv2;
          cvalue_type expPhase = std::polar(value_type(1.0), phase);
          cvalue_type expPhaseSq =
            std::polar(value_type(1.0), value_type(2.0) * phase);
          cvalue_type expPhaseInv =
            std::polar(value_type(1.0), - phase);
          cvalue_type expPhaseSqInv =
            std::polar(value_type(1.0), value_type(-2.0) * phase);

          opMatrixTrans[0] = cosAlphaDiv2Sq;
          opMatrixTrans[1] = expPhaseSq * sinAlphaDiv2Sq;
          opMatrixTrans[2] = cvalue_type(value_type(0), value_type(-1.0)) *
            expPhase * sinAlpha;
          opMatrixTrans[3] = expPhaseSqInv * sinAlphaDiv2Sq;
          opMatrixTrans[4] = cosAlphaDiv2Sq;
          opMatrixTrans[5] = cvalue_type(value_type(0), value_type(1.0)) *
            expPhaseInv * sinAlpha;
          opMatrixTrans[6] = cvalue_type(value_type(0), value_type(-0.5)) *
            expPhaseInv * sinAlpha;
          opMatrixTrans[7] = cvalue_type(value_type(0), value_type(0.5)) *
            expPhase * sinAlpha;
          opMatrixTrans[8] = cosAlpha;
        }

        void operator()(SpinStates<index_type, value_type> *states) {
          const unsigned int curStatesLocal = states->curStates;
          const unsigned int maxStatesLocal = states->maxStates;
          BLAS::gemm(
            BLAS::ColMajor, BLAS::NoTrans, BLAS::NoTrans,
            curStatesLocal, 3, 3,
            cvalue_type(1.0), states->stateValues.data(), maxStatesLocal,
            opMatrixTrans.data(), 3,
            cvalue_type(0), states->stateValuesBuffer.data(), maxStatesLocal);

          // move data back from the temporary buffer
          // this could be eliminated with proper double-buffering...
          BLAS::copy(curStatesLocal, states->stateValuesBuffer.data(),
            1, states->stateValues.data(), 1);
          BLAS::copy(
            curStatesLocal, states->stateValuesBuffer.data() + maxStatesLocal,
            1, states->stateValues.data() + maxStatesLocal, 1);
          BLAS::copy(
            curStatesLocal,
            states->stateValuesBuffer.data() + 2 * maxStatesLocal,
            1, states->stateValues.data() + 2 * maxStatesLocal, 1);
        }

      protected:
        std::vector< std::complex<value_type> > opMatrixTrans;
    };

    class Relaxation {
      public:
        Relaxation(value_type duration, value_type t1, value_type t2Star) {

          relaxationTerm[0] = relaxationTerm[1] = exp(-duration/t2Star);
          relaxationTerm[2] = exp(-duration/t1);
          relaxationTerm[3] = value_type(1) - relaxationTerm[2];
        }

        void operator()(SpinStates<index_type, value_type> *states) {
          const unsigned int curStatesLocal = states->curStates;
          const unsigned int maxStatesLocal = states->maxStates;

          BLAS::scal(
            curStatesLocal, relaxationTerm, states->stateValues.data(), 1);
          
          BLAS::scal(
            curStatesLocal,
            relaxationTerm + 1,
            states->stateValues.data() + maxStatesLocal, 1);
          
          BLAS::scal(
            curStatesLocal,
            relaxationTerm + 2,
            states->stateValues.data() + 2 * maxStatesLocal, 1);

          StateIndex originIndex({0,0,0});
          size_t originOffset = states->stateMap[originIndex];

          states->stateValues[originOffset + 2*maxStatesLocal] +=
            relaxationTerm[3];
        }

      protected:
        cvalue_type relaxationTerm[4];
    };
    
    class Spoiling {
      public:
        Spoiling(StateIndex spoilGrad, value_type trimThreshold = 0.0) :
          spoilGrad(spoilGrad),
          trimThreshold(trimThreshold){}

        void operator()(SpinStates<index_type, value_type> *states) {
          const unsigned int curStatesLocal = states->curStates;
          const unsigned int maxStatesLocal = states->maxStates;
          const unsigned int maxStatesBuffer = maxStatesLocal;

          BLAS::copy(curStatesLocal,
            states->stateValues.data(), 1,
            states->stateValuesBuffer.data(), 1);
          
          BLAS::copy(curStatesLocal,
            states->stateValues.data() + maxStatesLocal, 1,
            states->stateValuesBuffer.data() + maxStatesLocal, 1);

          memset(states->stateValues.data(), 0,
            curStatesLocal * sizeof(cvalue_type));
          
          memset(states->stateValues.data() + maxStatesLocal, 0,
            curStatesLocal * sizeof(cvalue_type));
         
          typename StateMapT::const_iterator mapIt = states->stateMap.begin();

          // TODO: consider whether emplace() removes need for copy
          StateMapT newStateMap = states->stateMap;

          for(; mapIt != states->stateMap.end(); mapIt++) {
            StateIndex curIndex = mapIt->first;
            size_t curOffset = mapIt->second;
 
            StateIndex curfPlusIndex;
            StateIndex curfMinusIndex;

            for(unsigned int i = 0; i < 3; i++) {
              curfPlusIndex[i] = curIndex[i] + spoilGrad[i];
              curfMinusIndex[i] = curIndex[i] - spoilGrad[i];
            }

            {
              typename StateMapT::iterator fPlusIt =
                newStateMap.find(curfPlusIndex);

              // if the element where fPlus is landing doesn't yet exist,
              // make it
              if(newStateMap.end() == fPlusIt) {
                // if the state list is full, expand it
                if(states->curStates == states->maxStates) {
                  states->expandStates(); 
                }

                newStateMap[curfPlusIndex] = states->curStates;

                states->stateValues[states->curStates] =
                  states->stateValuesBuffer[curOffset];
                
                states->stateValues[states->curStates + states->maxStates]
                  = cvalue_type(0);
                
                states->stateValues[states->curStates + 2 * states->maxStates]
                  = cvalue_type(0);

                states->curStates++;
              }
              // if the element does exist
              else {
                states->stateValues[fPlusIt->second] = 
                  states->stateValuesBuffer[curOffset];
              }
            }
            
            {
              typename StateMapT::iterator fMinusIt =
                newStateMap.find(curfMinusIndex);

              // if the element where fMinus is landing doesn't yet exist,
              // make it
              if(newStateMap.end() == fMinusIt) {
                // if the state list is full, expand it
                if(states->curStates == states->maxStates) {
                  states->expandStates(); 
                }

                newStateMap[curfMinusIndex] = states->curStates;
                
                states->stateValues[states->curStates]
                  = cvalue_type(0);

                states->stateValues[states->curStates + states->maxStates] = 
                  states->stateValuesBuffer[curOffset + maxStatesBuffer];
                
                states->stateValues[states->curStates + 2 * states->maxStates]
                  = cvalue_type(0);

                states->curStates++;
              }
              // if the element does exist
              else {
                states->stateValues[fMinusIt->second + states->maxStates] = 
                  states->stateValuesBuffer[curOffset + maxStatesBuffer];
              }
            }
          }

          states->stateMap = newStateMap;

          states->trimStates(trimThreshold);
        }

      protected:
        StateIndex spoilGrad;
        const value_type trimThreshold;
    };

  protected:
    void insertState(const StateIndex &index, const SpinState &state) {
      if(curStates == maxStates) {
        expandStates();
      }

      stateValues[curStates] = state.fPlus;
      stateValues[curStates + maxStates] = state.fMinus;
      stateValues[curStates + 2 * maxStates] = state.z;

      stateMap[index] = curStates;
      curStates++;
    }

    void expandStates() {
      size_t newMaxStates = maxStates * 2;
      std::vector<cvalue_type> newStateValues(newMaxStates * 3);
     
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
      stateValuesBuffer.resize(newMaxStates * 3);
    }

    void trimStates(const value_type threshold) {
      StateMapT newMap;
      size_t newCurStates = 0;

      typename StateMapT::const_iterator mapIt = stateMap.begin();

      for(; stateMap.end() != mapIt; mapIt++) {
        value_type nrm =
          BLAS::nrm2(3, stateValues.data() + mapIt->second, maxStates);

        //only copy states larger than threshold
        if(nrm > threshold) {
          newMap[mapIt->first] = newCurStates;
          stateValuesBuffer[newCurStates] =
            stateValues[mapIt->second];
          stateValuesBuffer[newCurStates + maxStates] =
            stateValues[mapIt->second + maxStates];
          stateValuesBuffer[newCurStates + 2 * maxStates] =
            stateValues[mapIt->second + 2 * maxStates];

          newCurStates++;
        }
      }

      stateMap = newMap;
      curStates = newCurStates;
      // move data back from the temporary buffer
      // this could be eliminated with proper double-buffering...
      BLAS::copy(curStates, stateValuesBuffer.data(),
        1, stateValues.data(), 1);
      BLAS::copy(
        curStates, stateValuesBuffer.data() + maxStates,
        1, stateValues.data() + maxStates, 1);
      BLAS::copy(
        curStates,
        stateValuesBuffer.data() + 2 * maxStates,
        1, stateValues.data() + 2 * maxStates, 1);
    }

    typedef std::map<StateIndex, size_t> StateMapT;

    StateMapT stateMap;
    size_t curStates;
    size_t maxStates;
    
    std::vector< cvalue_type > stateValues;
    std::vector< cvalue_type > stateValuesBuffer;

};

#endif

