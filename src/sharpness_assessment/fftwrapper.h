#ifndef FFTWRAPPER_H
#define FFTWRAPPER_H
#include "s3tools.h"

#ifdef USE_FFTW
#include "../fftw/fftw3.h"
#endif
/**
 * @brief The FFTWrapper class - wrapper for DFT implementation
 */
class FFTWrapper
{
public:
    FFTWrapper();
    /**
     * @brief absFFT2D_LL - get abs value of 2D FFT inplce
     * @param data
     * @param rows
     * @param columns
     */
    void absFFT2D_LL(s3real* data, int rows, int columns);

#ifdef USE_FFTW
public:
   void absFFT2D(float* data, int rows, int columns);
   void absFFT2D(double* data, int rows, int columns);
protected:
   fftwf_complex* scratchf;
   fftw_complex* scratch;
   unsigned int scratchSize;
#endif


};

#endif // FFTWRAPPER_H
