#ifndef SpinStates_h
#define SpinStates_h

#include <map>
#include <array>
#include <vector>

#include "blas_local.h"

#include "mkl_ops.h"

#include "MKLAllocator/MKLAllocator.h"


template <typename index_type, typename _value_type>
class SpinStates {
  public:
    typedef _value_type value_type;
    typedef std::complex<value_type> cvalue_type;
  
    typedef std::array<index_type, 3> StateIndex;

  protected:

    typedef std::vector<cvalue_type, mkl_allocator<cvalue_type> > BufferVec;
    typedef std::unique_ptr<BufferVec> BufferPtr;

    static value_type stateIndexNormSq(const StateIndex &index) {
      return value_type(index[0] * index[0]) +
        value_type(index[1] * index[1]) +
        value_type(index[2] * index[2]);
    }
    
    static value_type stateIndexNorm(const StateIndex &index) {
      return std::sqrt(stateIndexNormSq(index));
    }

  public:

    struct SpinState {
      cvalue_type fPlus, fMinus, z; 
    }; 
    
    SpinStates(
      const value_type state_increment_inv_mm, 
      const value_type diffusionCoeff_mm_sq_per_s
      ) :
      diffusionCoeff_mm_sq_per_s(diffusionCoeff_mm_sq_per_s),
      state_increment_inv_mm(state_increment_inv_mm),
      curRadii({{1,1,1}}),
      curSideLengths({{1,1,1}}),
      curStates(curSideLengths[0] * curSideLengths[1] * curSideLengths[2]),
      frontBuffer(std::make_unique<BufferVec>(curStates * 3)),
      backBuffer(std::make_unique<BufferVec>(curStates * 3)),
      indexNormSqBuffer(
        std::make_unique<BufferVec>(curStates)),
      indexNormBuffer(
        std::make_unique<BufferVec>(curStates))
    {
      //initialize signal at origin
      (*frontBuffer)[0] = 0;
      (*frontBuffer)[1] = 0;
      (*frontBuffer)[2] = 1;

      fillIndexNormBuffer();
    }

    SpinStates(
      const value_type state_increment_gradient_area_mT_s_per_m,
      const value_type larmour_freq_Hz_per_mT, 
      const value_type diffusionCoeff_mm_sq_per_s
      ) :
      SpinStates(
        state_increment_gradient_area_mT_s_per_m * larmour_freq_Hz_per_mT *
        0.001,
        diffusionCoeff_mm_sq_per_s)
    {}

    SpinState getState(StateIndex index) {
      SpinState ret;

      size_t offset = stateIndexToOffset(index);

      ret.fPlus = (*frontBuffer)[offset]; 
      ret.fMinus = (*frontBuffer)[offset + curStates]; 
      ret.z = (*frontBuffer)[offset + 2 * curStates]; 

      return ret;
    }
   

    size_t size() const {
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
          BLAS::gemm(
            BLAS::ColMajor, BLAS::NoTrans, BLAS::NoTrans,
            curStatesLocal, 3, 3,
            cvalue_type(1.0), states->frontBuffer->data(), curStatesLocal,
            opMatrixTrans.data(), 3,
            cvalue_type(0), states->backBuffer->data(), curStatesLocal);
      
          // swap the front and back buffers to our updates are visible
          std::swap(states->frontBuffer, states->backBuffer); 
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

          BLAS::scal(
            curStatesLocal, relaxationTerm, states->frontBuffer->data(), 1);
          
          BLAS::scal(
            curStatesLocal,
            relaxationTerm + 1,
            states->frontBuffer->data() + curStatesLocal, 1);
          
          BLAS::scal(
            curStatesLocal,
            relaxationTerm + 2,
            states->frontBuffer->data() + 2 * curStatesLocal, 1);

          StateIndex originIndex({{0,0,0}});
          size_t originOffset = states->stateIndexToOffset(originIndex);

          states->frontBuffer->data()[originOffset + 2 * curStatesLocal] +=
            relaxationTerm[3];
        }

      protected:
        cvalue_type relaxationTerm[4];
    };
    
    
    class Spoiling {
      public:
        Spoiling(
          StateIndex spoilGrad,
          value_type TR_ms,
          value_type trimThreshold = 0.0) :
          spoilGrad(spoilGrad),
          spoilGradNorm(stateIndexNorm(spoilGrad)),
          TR_ms(TR_ms),
          trimThreshold(trimThreshold)
        {}

        void operator()(SpinStates<index_type, value_type> *states)
        {
          const unsigned int curStatesLocal = states->curStates;
          
          // compute the common scaling factor used in all diffusion
          // spoiling calculations
          const cvalue_type negCommonDiffScale =
            states->diffusionCoeff_mm_sq_per_s *
            states->state_increment_inv_mm * states->state_increment_inv_mm *
            TR_ms * -0.001;
    
          // zero out the first third of the back buffer
          std::fill(
            states->backBuffer->begin(),
            states->backBuffer->begin() + curStatesLocal, 0);

          // calculate the scale factor for longitudinal states
          // and store it in the first third of the back buffer
          BLAS::axpy(
            curStatesLocal,
            &negCommonDiffScale, states->indexNormSqBuffer->data(), 1,
            states->backBuffer->data(), 1
          );

          // exponentiate the factor to get the actual weighting factor
          // and store it in the second third of the back buffer
          MKL::vExp(
            curStatesLocal,
            states->backBuffer->data(),
            states->backBuffer->data() + curStatesLocal
          );

          // scale all the z states by this weighting factor
          MKL::vMul(
            curStatesLocal,
            states->backBuffer->data() + curStatesLocal,
            states->frontBuffer->data() + 2 * curStatesLocal,
            states->frontBuffer->data() + 2 * curStatesLocal
          );
          
          // update the scale factor with additional terms for
          // the transverse states and store it in the first third of
          // the back buffer
          BLAS::axpy(
            curStatesLocal,
            &negCommonDiffScale, states->indexNormBuffer->data(), 1,
            states->backBuffer->data(), 1
          );

          typename BufferVec::iterator it = states->backBuffer->begin();
          typename BufferVec::iterator itEnd = it + curStatesLocal;

          value_type additionalTerm =
            std::abs(negCommonDiffScale) / value_type(3.0);

          while(it != itEnd) {
            *it *= additionalTerm;
            it++;
          }

          // exponentiate the transfer factor to get the actual
          // transverse weighting factor and store it in the second third of
          // the back buffer
          MKL::vExp(
            curStatesLocal,
            states->backBuffer->data(),
            states->backBuffer->data() + curStatesLocal
          );
         
          // scale the fPlus states by this factor 
          MKL::vMul(
            curStatesLocal,
            states->backBuffer->data() + curStatesLocal,
            states->frontBuffer->data(),
            states->frontBuffer->data()
          );
          
          // scale the fMinus states by this factor 
          MKL::vMul(
            curStatesLocal,
            states->backBuffer->data() + curStatesLocal,
            states->frontBuffer->data() + curStatesLocal,
            states->frontBuffer->data() + curStatesLocal
          );

          // zero out the first 2/3 of the back buffer
          std::fill(
            states->backBuffer->begin(),
            states->backBuffer->begin() + 2 * curStatesLocal, 0
          );

          // copy the z values since they don't shift
          BLAS::copy(curStatesLocal,
              states->frontBuffer->data() + 2 * curStatesLocal, 1,
              states->backBuffer->data() + 2 * curStatesLocal, 1
          );

          // iterate over each position in the front buffer
          size_t curOffset = 0;
          StateIndex curRadii = states->curRadii;
          StateIndex curIndex =
            StateIndex({{
              1 - curRadii[0],
              1 - curRadii[1],
              1 - curRadii[2]
            }});

          // we'll keep extra indicies for the backBuffer in case
          // we need to expand it
          StateIndex backRadii = states->curRadii;
          StateIndex backSideLengths = states->curSideLengths;
          size_t backStates = states->curStates;


          //const value_type oneThird = value_type(1.0)/value_type(3.0);

          for(;curIndex[0] < curRadii[0]; curIndex[0]++)
          {
            for(
              curIndex[1] = 1 - curRadii[1];
              curIndex[1] < curRadii[1];
              curIndex[1]++)
            {

              for(
                curIndex[2] = 1 - curRadii[2];
                curIndex[2] < curRadii[2];
                curIndex[2]++, curOffset++)
              {
                
                // move fPlus component to destination state
                {
                  StateIndex curfPlusIndex;

                  // compute the destination states for the components 
                  for(unsigned int i = 0; i < 3; i++)
                  {
                    curfPlusIndex[i] = curIndex[i] + spoilGrad[i];
                  }

                  bool needsExpansion = false;
                  StateIndex neededRadii = backRadii;

                  for(unsigned int i = 0; i < 3; i++)
                  {
                    while(std::abs(curfPlusIndex[i]) >= neededRadii[i])
                    {
                      neededRadii[i] *= 2;
                      needsExpansion = true;
                    }
                  }

                  // if we don't need to expand the back buffer,
                  // or we do need to expand but the current fPlus
                  // component is big enough to warrant buffer expansion
                  if((!needsExpansion) ||
                      (std::abs(states->frontBuffer->data()[curOffset]) >=
                        trimThreshold))
                  {

                    if(needsExpansion) {
                      states->expandBackBuffer(
                        neededRadii, &backRadii, &backSideLengths, &backStates);
                    }

                    size_t destOffset =
                      stateIndexToOffset(
                        curfPlusIndex, backRadii, backSideLengths);
                    
                    states->backBuffer->data()[destOffset] = 
                      states->frontBuffer->data()[curOffset];
                  }
                }
                
                // move fMinus component to destination state
                {
                  StateIndex curfMinusIndex;

                  // compute the destination states for the components 
                  for(unsigned int i = 0; i < 3; i++)
                  {
                    curfMinusIndex[i] = curIndex[i] - spoilGrad[i];
                  }

                  bool needsExpansion = false;
                  StateIndex neededRadii = backRadii;

                  for(unsigned int i = 0; i < 3; i++)
                  {
                    while(std::abs(curfMinusIndex[i]) >= neededRadii[i])
                    {
                      neededRadii[i] *= 2;
                      needsExpansion = true;
                    }
                  }

                  // if we don't need to expand the back buffer,
                  // or we do need to expand but the current fMinus
                  // component is big enough to warrant buffer expansion
                  if((!needsExpansion) ||
                      std::abs(
                        states->frontBuffer->data()[
                          curOffset + curStatesLocal]) >=
                      trimThreshold)
                  {
                    if(needsExpansion) {
                      states->expandBackBuffer(
                        neededRadii, &backRadii, &backSideLengths, &backStates);
                    }

                    size_t destOffset =
                      stateIndexToOffset(
                        curfMinusIndex, backRadii, backSideLengths);
                    
                    states->backBuffer->data()[destOffset + backStates] = 
                      states->frontBuffer->data()[curOffset + curStatesLocal];
                  }
                }
              }
            }
          }

          // if we expanded the back buffer, we need to make a new front
          // buffer before we do the swap
          if(curStatesLocal != backStates){
            states->curStates = backStates;
            states->curSideLengths = backSideLengths;
            states->curRadii = backRadii;

            BufferPtr newFrontBuffer(
              std::make_unique<BufferVec>(states->curStates * 3));

            std::swap(states->frontBuffer, newFrontBuffer);

            // we also need to allocate and refill the indexNorm buffers
            BufferPtr newIndexNormSqBuffer(
              std::make_unique<BufferVec>(states->curStates));
            BufferPtr newIndexNormBuffer(
              std::make_unique<BufferVec>(states->curStates));

            std::swap(states->indexNormBuffer, newIndexNormBuffer);
            std::swap(states->indexNormSqBuffer, newIndexNormSqBuffer);

            states->fillIndexNormBuffer();
          }

          // swap the front and back buffers so our updates are visible
          std::swap(states->frontBuffer, states->backBuffer);

          // erase the back buffer
          //std::fill(states->backBuffer->begin(), states->backBuffer->end(), 0);
        }
      
      protected:
        StateIndex spoilGrad;
        const value_type spoilGradNorm;
        const value_type TR_ms;
        const value_type trimThreshold;
    };

  protected:
    
    static size_t stateIndexToOffset(
      const StateIndex index,
      const StateIndex radii,
      const StateIndex sideLengths)
    {
      return (
          (index[0] + radii[0] - 1) * sideLengths[1] +
          (index[1] + radii[1] - 1)
        ) * sideLengths[2] + index[2] + radii[2] - 1;
    }

    size_t stateIndexToOffset(const StateIndex index) const
    {
      return stateIndexToOffset(index, curRadii, curSideLengths); 
    }

    void fillIndexNormBuffer() {
      // iterate over each position in the front buffer
      size_t curOffset = 0;
      StateIndex curIndex =
        StateIndex({{
          1 - curRadii[0],
          1 - curRadii[1],
          1 - curRadii[2]
        }});
          
      for(;curIndex[0] < curRadii[0]; curIndex[0]++)
      {
        for(
          curIndex[1] = 1 - curRadii[1];
          curIndex[1] < curRadii[1];
          curIndex[1]++)
        {

          for(
            curIndex[2] = 1 - curRadii[2];
            curIndex[2] < curRadii[2];
            curIndex[2]++, curOffset++)
          {

            indexNormSqBuffer->data()[curOffset] = stateIndexNormSq(curIndex);
          }
        }
      }

      MKL::vSqrt(curStates, indexNormSqBuffer->data(), indexNormBuffer->data());
    }

    void expandBackBuffer(
      const StateIndex neededRadii,
      StateIndex *backRadii,
      StateIndex *backSideLengths,
      size_t *backStates
    )
    {
      StateIndex newSideLengths;

      size_t newStates = 1;

      for(unsigned int i = 0; i < 3; i++)
      {
        newSideLengths[i] = neededRadii[i] * 2 - 1;
        newStates *= newSideLengths[i];
      }

      // We'll allocate the new buffer here, but then swap the ptr by the end
      // of this method, so the new buffer will be retained
      BufferPtr newBuffer(std::make_unique<BufferVec>(newStates * 3));
      std::fill(newBuffer->begin(), newBuffer->end(), 0);

      // Perform the copy into the new buffer.
      // We do one copy per inner-most line of the buffer (e.g., if we think
      // of coords as z-y-x, then we do one copy of complete x-line
      // per z/y pair)
      cvalue_type *readHead = backBuffer->data();
      
      std::array<index_type, 2> curIndicies({{0,0}});

      std::array<size_t, 3> writeOffsets;

      writeOffsets[0] = (newSideLengths[0] - (*backSideLengths)[0]) / 2;

      writeOffsets[1] = (newSideLengths[1] - (*backSideLengths)[1]) / 2;

      writeOffsets[2] = (newSideLengths[2] - (*backSideLengths)[2]) / 2;

      for(; curIndicies[0] < (*backSideLengths)[0]; curIndicies[0]++)
      {
        size_t outerOffset = (writeOffsets[0] + curIndicies[0]) *
          newSideLengths[1] * newSideLengths[2];

        for(curIndicies[1] = 0; curIndicies[1] < (*backSideLengths)[1]; curIndicies[1]++)
        {
          size_t middleOffset =  (writeOffsets[1] + curIndicies[1]) *
            newSideLengths[2];

          cvalue_type *writeHead =
            newBuffer->data() + outerOffset +
              middleOffset + writeOffsets[2];
          
          BLAS::copy((*backSideLengths)[2], readHead, 1, writeHead, 1);
          BLAS::copy((*backSideLengths)[2],
            readHead + (*backStates), 1,
            writeHead + newStates, 1);
          BLAS::copy((*backSideLengths)[2],
            readHead + 2 * (*backStates), 1,
            writeHead + 2 * newStates, 1);

          readHead += (*backSideLengths)[2]; 
        }
      }

      std::swap(newBuffer, backBuffer); 

      *backStates = newStates;
      *backSideLengths = newSideLengths;
      *backRadii = neededRadii;
    }

    value_type diffusionCoeff_mm_sq_per_s;
    value_type state_increment_inv_mm;

    std::array<index_type, 3> curRadii;
    std::array<index_type, 3> curSideLengths; 
    index_type curStates;

    BufferPtr frontBuffer;
    BufferPtr backBuffer;
    BufferPtr indexNormSqBuffer;
    BufferPtr indexNormBuffer;

};

#endif

