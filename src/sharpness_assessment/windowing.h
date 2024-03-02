/*
a = [0.35875 0.48829 0.14128 0.01168];
window(i+1) = a(1) - a(2)*cos(2*pi*i/N) + a(3)*cos(4*pi*i/N) - a(4)*cos(6*pi*i/N); %blackman-harris
window(i+1) = 0.42 - 0.5*cos(2*pi*i/N) + 0.08*cos(4*pi*i/N);%blackman
window(i+1) = 0.54 - 0.46*cos(2*pi*i/N);% hamming
window(i+1) = 0.5*(1-cos(2*pi*i/N));% hann
*/

#ifndef WINDOWING_H
#define WINDOWING_H
#include <cmath>
#ifndef M_PI
#define M_PI       3.14159265358979323846
#endif
template <typename T>
/**
 * @brief The Windowing class - some windoming functions
 */
class Windowing
{
public:

	static void rectangular(T* out, int length)
	{
		for(int i=0;i<length;i++)
			out[i]=1.;
	}
    static void blackman(T* out,int length)
	{
		int n=length-1;
		for(int i=0;i<length;i++)
		{
			out[i] = 0.42 - 0.5*cos(2*M_PI*(double)i/(double)n) + 0.08*cos(4*M_PI*(double)i/(double)n);
		}
	}
	static void hamming(T* out,int length)
	{
		int n=length-1;
		for(int i=0;i<length;i++)
		{
			out[i] = 0.54 - 0.46*cos(2*M_PI*(double)i/(double)n);
		}

	}
	static void hann(T* out,int length)
	{
		int n=length-1;
		for(int i=0;i<length;i++)
		{
			out[i] = 0.5*(1.-cos(2*M_PI*(double)i/(double)n));
		}

	}
    static void hanning(T* out,int length)
    {
        int n=length-1;
        for(int i=0;i<length;i++)
        {
            out[i] = 0.5*(1.-cos(2*M_PI*(double)(i+1)/(double)(n+2)));
        }

    }
	static void blackmanharris(T* out,int length)
	{
		int n=length-1;
		for(int i=0;i<length;i++)
		{
			out[i]	= 0.35875 - 0.48829*cos(2*M_PI*(double)i/(double(n)))
					+ 0.14128*cos(4*M_PI*(double)i/(double(n)))
					- 0.01168*cos(6*M_PI*(double)i/(double(n)));
		}

	}

	static void nuttall(T* out,int length)
	{
		int n=length-1;
		for(int i=0;i<length;i++)
		{
			out[i]	= 0.3635819 - 0.4891775*cos(2*M_PI*(double)i/(double(n)))
				+ 0.1365995*cos(4*M_PI*(double)i/(double(n)))
				- 0.0106411*cos(6*M_PI*(double)i/(double(n)));
		}

	}
};

#endif
