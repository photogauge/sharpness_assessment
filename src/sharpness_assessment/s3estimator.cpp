#include "s3estimator.h"
#include <iostream>
#include <cmath>
#include "fftwrapper.h"
#include "windowing.h"
//#include "qttools.h"
#include <cstring>

S3Estimator::S3Estimator()
{

}
void S3Estimator::s3_map(Matrix& img)
{
    m_s0=img;
    m_s1=spectral_map(img,16);
    Matrix s_map21=spatial_map(img,8);
    Matrix s_map22=spatial_map(img,4);
    m_s2=s_map21;

    m_s3 = Matrix(img.rows(),img.columns());
    s3real* s1=m_s1.data();
    s3real* s2=m_s2.data();
    s3real* s3=m_s3.data();
    s3real* s22=s_map22.data();
    for(int i=0;i<m_s2.count();i++)
    {
        s2[i]=std::max(s2[i],s22[i]);
        s3[i]=sqrt(s1[i]*s2[i]);
    }
}
Matrix S3Estimator::spatial_map(Matrix& img, int padding)
{
    Matrix input;
    if(padding)
    {
        input = pad(img,padding);
    }
    else
    {
        input=img;
    }
    Matrix res(input.rows(),input.columns());
    res.zero();

    const int blk_size = 8;
    int num_rows = input.rows();
    int num_cols = input.columns();
    for(int r=blk_size/2;r<=(num_rows-blk_size/2); r+=blk_size)
    {
        for(int c=blk_size/2;c<=(num_cols-blk_size/2); c+=blk_size)
        {
            Matrix gry_blk = input.block(r-blk_size/2, c-blk_size/2, blk_size, blk_size);
            s3real maxVal=0;
            for(int i=0;i<blk_size-1;i++)
            {
                for(int j=0;j<blk_size-1;j++)
                {
                    s3real tmp=0;
                    tmp += std::abs(gry_blk.at(i,j) - gry_blk.at(i,j+1));
                    tmp += std::abs(gry_blk.at(i,j) - gry_blk.at(i+1,j));
                    tmp += std::abs(gry_blk.at(i,j) - gry_blk.at(i+1,j+1));
                    tmp += std::abs(gry_blk.at(i+1,j+1) - gry_blk.at(i+1,j));
                    tmp += std::abs(gry_blk.at(i+1,j) - gry_blk.at(i,j+1));
                    tmp += std::abs(gry_blk.at(i+1,j+1) - gry_blk.at(i,j+1));
                    maxVal = std::max(maxVal,tmp);
                }
            }
            s3real tv_max=maxVal/4./255.;

            for(int ir=0;ir<blk_size;ir++)
            {
                for(int ic=0;ic<blk_size;ic++)
                {
                    res.at(ir+r-blk_size/2,ic+c-blk_size/2)=tv_max;
                }
            }
        }
    }

    if(padding>0)
    {
        return res.block(padding,padding,res.rows()-2*padding,res.columns()-2*padding);
    }
    else
    {
        return res;
    }
}
Matrix S3Estimator::pad(Matrix& src, int padding)
{
    Matrix out(src.rows()+padding*2,src.columns()+padding*2);
    out.zero();
    for(int r=0;r<src.rows();r++)
    {
        memcpy(out.dataAt(r+padding,padding),src.dataAt(r,0),sizeof(s3real)*src.columns());
        for(int i=0;i<padding;i++)
        {
            out.at(r+padding,i)=src.at(r,padding-i-1);
            out.at(r+padding,out.columns()-i-1)=src.at(r,src.columns()-(padding-i));
        }
    }
    for(int r=0;r<padding;r++)
    {
         memcpy(out.dataAt(padding-r-1,0),out.dataAt(r+padding,0),sizeof(s3real)*out.columns());
         memcpy(out.dataAt(out.rows()-1-r,0),out.dataAt(out.rows()-2*padding+r,0),sizeof(s3real)*out.columns());
    }
    return out;
}

s3real S3Estimator::blk_amp_spec_slope_eo_toy(Matrix& blk)
{
    if(m_window.rows()!=blk.rows())
    {
        int cnt=blk.rows();
        std::vector<s3real> window(cnt);
        Windowing<s3real>::hanning(window.data(),cnt);
        m_window=Matrix(cnt,cnt);
        //m_window.name="Window";
        for(int i=0;i<cnt;i++)
        {
            for(int j=0;j<cnt;j++)
            {
                m_window.at(i,j)=window[i]*window[j];
            }
        }
        //m_window.print();
    }

    s3real* blkData=blk.data();
    s3real* windowData=m_window.data();
    for(int i=0;i<blk.count();i++)
    {
        (*blkData++)*=*windowData++;
    }

    fft.absFFT2D_LL(blk.data(),blk.rows(),blk.columns());
    std::vector<s3real> s = eo_polaraverage(blk);
    static std::vector<s3real> f;
    f=eo_frequencies(s.size());
    f.erase(f.begin());
    s.erase(s.begin());

    for(unsigned int i=0;i<s.size();i++)
    {
        s[i]=log(s[i]);
        f[i]=log(f[i]);
    }
    s3real slope=fitLine(f,s);
    return -slope;
}
Matrix S3Estimator::spectral_map(Matrix& img, int pad_len)
{
    Matrix input;
    if(pad_len>0)
    {
        input=pad(img,pad_len);
    }
    else
    {
        input=img;
    }

    int blk_size=32;
    int d_blk=blk_size/4;
    Matrix res(input.rows(), input.columns());
    Matrix luminance(input);
    applyLuminance(luminance);
    //std::fill(res.begin(),res.end(),-100);
    res.zero();
    s3real contrast_threshold=0;

    for(int r = blk_size/2; r<=(input.rows()-blk_size/2);r+=d_blk)
    {
        for(int c = blk_size/2; c<=(input.columns()-blk_size/2);c+=d_blk)
        {
            Matrix gry_blk = input.block(r-blk_size/2,c-blk_size/2,blk_size,blk_size);
            Matrix tmp=luminance.block(r-blk_size/2,c-blk_size/2,blk_size,blk_size);
            Matrix contrastMap = contrast_map_overlap(tmp);
            s3real val=0;
            if(*std::max_element(contrastMap.begin(), contrastMap.end())>contrast_threshold)
            {
                val=blk_amp_spec_slope_eo_toy(gry_blk);

                val=1.-1./(1.+exp(-3*(val-2)));
            }

            for(int ir=0;ir<d_blk;ir++)
            {
                for(int ic=0;ic<d_blk;ic++)
                {
                    res.at(ir+r-d_blk/2,ic+c-d_blk/2)=val;
                }
            }

        }
    }
    if(pad_len>0)
    {
        return res.block(pad_len,pad_len,res.rows()-2*pad_len,res.columns()-2*pad_len);
    }
    else
    {
        return res;
    }
}

s3real S3Estimator::fitLine(std::vector<s3real>& x,std::vector<s3real>& y)
{
    s3real sumX=0, sumY=0, sumXY=0, sumX2=0,slope;
    //, intercept;
    int nPoints=x.size();
    for(int i=0; i<nPoints; i++)
    {
        sumX += x[i];
        sumY += y[i];
        sumXY += x[i] * y[i];
        sumX2 += x[i] * x[i];
    }

    s3real xMean = sumX / nPoints;
    s3real yMean = sumY / nPoints;
    s3real denominator = sumX2 - sumX * xMean;
    if( std::abs(denominator) == 0. ) {
        slope=0;
//        intercept=0;
    }
    else
    {
        slope = (sumXY - sumX * yMean) / denominator;
//        intercept  = yMean - slope * xMean;
    }
    return slope;
}
void S3Estimator::printVector(std::vector<s3real>& v)
{
    for(unsigned int i=0;i<v.size();i++)
    {
        std::cout<<v[i]<<" ";
    }
     std::cout<<std::endl;
}

std::vector<s3real> S3Estimator::eo_frequencies(int count)
{
     std::vector<s3real> f(count);
    for(int i=0;i<count;i++)
    {
        f[i]=(i*0.5/(count-1));
    }
    return f;
}
void S3Estimator::eo_generate_map(int block_size)
{
    m_eo_map.clear();
    s3real rr=M_PI*2;
    int dr=360;
    int n=block_size;
     int cc=block_size;
    for(int r=1;r<=n/2;r++)
    {
        for(int ith=0; ith<dr; ith++)
        {
            s3real th=(s3real)ith/dr;
            s3real x = r*sin(th*rr);
            s3real y = r*cos(th*rr);

            int x1 = copysign(std::floor(std::abs(x)),x);
            int x2 = copysign(std::ceil(std::abs(x)),x);
            int y1 = copysign(std::floor(std::abs(y)),y);
            int y2 = copysign(std::ceil(std::abs(y)),y);

            s3real ex=std::abs(x-x1);
            s3real ey=std::abs(y-y1);

            if(x2<0.)
            {
                ex=std::abs(x-x2);
                if(x1<0.)
                    x1=n+x1;
                x2=n+x2;
            }
            if(y2<0.)
            {
                ey = std::abs(y - y2);
                if(y1<0.)
                    y1 = n + y1;
                y2 = n + y2;
            }
            eo_map_struct zzz;
            zzz.ex=ex;
            zzz.ey=ey;
            zzz.f11 =x1*cc+ y1;
            zzz.f12 = x1*cc+ y2;
            zzz.f21 = x2*cc+y1;
            zzz.f22 = x2*cc+y2;
            m_eo_map.push_back(zzz);
        }
    }
//      qDebug()<<this<<"generated map size"<<m_eo_map.size();
}
std::vector<s3real> S3Estimator::eo_polaraverage(Matrix& data)
{
    unsigned int dr=360;
    unsigned int n=std::max(data.rows(),data.columns());
    if(n/2*dr!=m_eo_map.size())
    {
        eo_generate_map(n);
    }
    data.at(0,0)=(data.at(1,0) + data.at(0,1))*0.5;
    std::vector<s3real> s(n/2+1);
    s[0]=0;

    eo_map_struct* mapData=m_eo_map.data();
    s3real* dataPointer=data.data();
  //  int cc=data.columns();
    for(unsigned int r=1;r<=n/2;r++)
    {
        s3real zs=0;
        for(unsigned int ith=0; ith<dr; ith++)
        {
            const s3real f11 = dataPointer[mapData->f11];
            s3real f12 = *(dataPointer+ mapData->f12);
            s3real f21 = *(dataPointer+ mapData->f21);
            s3real f22 = *(dataPointer+ mapData->f22);

            s3real z=(f21-f11)*mapData->ex*(1.-mapData->ey) + (f12-f11)*(1.-mapData->ex)*mapData->ey + (f22-f11)*mapData->ex*mapData->ey + f11;
            zs = zs + z;
            mapData++;
        }
        s[r]=zs/dr;
    }
    return s;
}
void S3Estimator::applyLuminance(Matrix& m)
{
    s3real* d=m.data();
    for(int i=0;i<m.count();i++)
    {
        s3real val=0.7656 + 0.0364*(*d);
        //*d=val*val;
        *d = std::pow(val,2.2);
        d++;
    }
}

#define OPTIMIZE
Matrix S3Estimator::contrast_map_overlap(Matrix& img)
{
    Matrix out(img.rows(),img.columns());
    out.zero();
    const int blk_size=8;
    const int d_blk=blk_size/2;
    for(int r=0;r<img.rows()-d_blk;r+=d_blk)
    {
        for(int c=0;c<img.columns()-d_blk;c+=d_blk)
        {
            Matrix block=img.block(r,c,blk_size,blk_size);
#ifdef OPTIMIZE
            s3real* d=block.data();
            s3real min_lum=*d;
            s3real max_lum=*d;
            s3real m_lum=0.;
            for(int i=0;i<block.count();i++)
            {
                min_lum=std::min(min_lum,*d);
                max_lum=std::max(max_lum,*d);
                m_lum+=*d;
                d++;
            }
            m_lum/=block.count();

#else
            s3real min_lum=(*std::min_element(block.begin(),block.end()));
            s3real max_lum=(*std::max_element(block.begin(),block.end()));
            s3real m_lum=mean2(block);
#endif

            if(m_lum>127.5)
            {
                s3real* d=block.data();
                for(int i=0;i<block.count();i++)
                {
                    *d=255.-(*d);
                    d++;
                }
                m_lum=255.-m_lum;
            }

            //std::max_element

            s3real contrast=0.;
            if(m_lum>2. && (max_lum-min_lum)>5.)
            {
                contrast=block.std2(m_lum)/m_lum;
            }
            contrast=std::min(contrast,static_cast<s3real>(5.))/5.;

            for(int rs1=0;rs1<d_blk;rs1++)
            {
                s3real* d=out.dataAt(rs1+r,c);
                for(int cs1=0;cs1<d_blk;cs1++)
                {
                    *d++=contrast;
//                    out.at(rs1+r,cs1+c)=contrast;
                }
            }
        }

    }
    return out;
}
