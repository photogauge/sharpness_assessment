#include "fftwrapper.h"
#include <cmath>
//#include <QDebug>
#include "low_level_fft.h"
#include <iostream>
#include <string.h>
FFTWrapper::FFTWrapper()
{
#ifdef USE_FFTW
    scratchf=0;
    scratch=0;
#endif
}
void FFTWrapper::absFFT2D_LL(s3real *data, int rows, int columns)
{
    static std::vector<double> tmp;
    static std::vector<double> line;
    unsigned int cnt=rows*columns;
    if(tmp.size()!=cnt)
    {
        tmp.resize(cnt*2);
    }
    if(line.size()!=static_cast<unsigned int>(rows*2))
        line.resize(rows*2);

    memset(tmp.data(),0,sizeof(double)*cnt*2);
    for(unsigned int i=0,k=0;i<cnt;i++,k+=2)
    {
        tmp[k]=data[i];
    }

    for(int row=0;row<rows;row++)
    {
        LowLevelFFT::cdft(columns*2,-1,(double*)(tmp.data()+columns*row*2));
    }
    for(int col=0;col<columns;col++)
    {
        for(int i=0,i2=0;i<rows;i++,i2+=2)
        {
            line[i2]=tmp[(col+i*columns)*2];
            line[i2+1]=tmp[(col+i*columns)*2+1];
        }
        LowLevelFFT::cdft(rows*2,-1,(double*)(line.data()));
        for(int i=0,i2=0;i<rows;i++,i2+=2)
         {
            tmp[(col+i*columns)*2]=line[i2];
            tmp[(col+i*columns)*2+1]=line[i2+1];
        }
    }
    double* d=tmp.data();
    for(unsigned int i=0;i<cnt;i++)
    {
        double v=(*d)*(*d);
        d++;
        v+=(*d)*(*d);
        d++;
        data[i]=sqrt(v);
    }
}
#ifdef USE_FFTW

void FFTWrapper::absFFT2D(float *data, int rows, int columns)
{
    if(scratchSize!=rows*columns)
    {
        scratchSize=rows*columns;
        if(scratchf!=0)
            fftwf_free(scratchf);
        scratchf = fftwf_alloc_complex(scratchSize);
    }
    fftwf_complex a;
    int count=rows*columns;
    float* inputData=data;
    for(int i=0;i<count;i++)
    {
        scratchf[i][0]=*inputData++;
        scratchf[i][1]=0;
    }
    fftwf_execute(fftwf_plan_dft_2d(rows,columns,scratchf,scratchf,FFTW_FORWARD,FFTW_ESTIMATE));

    inputData=data;
    for(int i=0;i<count;i++)
    {
        inputData[i]=sqrt(scratchf[i][0]*scratchf[i][0] + scratchf[i][1]*scratchf[i][1]);
    }
}

void FFTWrapper::absFFT2D(double *data, int rows, int columns)
{
    if(scratchSize!=rows*columns)
    {
        scratchSize=rows*columns;
        if(scratch!=0)
            fftw_free(scratch);
        scratch = fftw_alloc_complex(scratchSize);
    }
    fftwf_complex a;
    int count=rows*columns;
    double* inputData=data;
    for(int i=0;i<count;i++)
    {
        scratch[i][0]=*inputData++;
        scratch[i][1]=0;
    }
    fftw_execute(fftw_plan_dft_2d(rows,columns,scratch,scratch,FFTW_FORWARD,FFTW_ESTIMATE));
    for(int i=0;i<count;i++)
    {
        data[i]=sqrt(scratch[i][0]*scratch[i][0] + scratch[i][1]*scratch[i][1]);
    }
}
#endif
