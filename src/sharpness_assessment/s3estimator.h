#ifndef S3ESTIMATOR_H
#define S3ESTIMATOR_H
#include "s3tools.h"
#include <vector>
#include <algorithm>
#include "fftwrapper.h"
using namespace S3Tools;
class S3Estimator
{
public:
    S3Estimator();
    /**
     * @brief s3_map - calculte sharpness maps
     * @param img - input matrix
     */
    void s3_map(Matrix& img);
    /**
     * @brief m_s0 - input matrix
     */
    Matrix m_s0;
    /**
     * @brief m_s1 - spectral map
     */
    Matrix m_s1;
    /**
     * @brief m_s2 - spatial map
     */
    Matrix m_s2;
    /**
     * @brief m_s3 - combined sharpness map
     */
    Matrix m_s3;

protected:
    /**
     * @brief contrast_map_overlap
     * @param img - input matrix
     * @return
     */
    Matrix contrast_map_overlap(Matrix& img);

    /**
     * @brief eo_polaraverage
     * @param blk - input matrix
     * @return
     */
    std::vector<s3real> eo_polaraverage(Matrix& blk);

    /**
     * @brief eo_frequencies - generate vector of freuqency values for "blk_amp_spec_slope_eo_toy"
     * @param count
     * @return
     */
    std::vector<s3real> eo_frequencies(int count);

    /**
     * @brief eo_generate_map - precompute some data to optimize "eo_polaraverage"
     * @param block_size
     */
    void eo_generate_map(int block_size);

    /**
     * @brief applyLuminance - convert grayscale value to luminance
     * @param m - matrix to process
     */
    void applyLuminance(Matrix& m);

    /**
     * @brief fitLine - fit line in set of points and return slope
     * @param x
     * @param y
     * @return
     */
    s3real fitLine(std::vector<s3real>& x,std::vector<s3real>& y);

    s3real blk_amp_spec_slope_eo_toy(Matrix& blk);
    /**
     * @brief spectral_map - calculate spectral map
     * @param img - source matrix
     * @param padding - patting size
     * @return
     */
    Matrix spectral_map(Matrix& img, int padding);

    /**
     * @brief spatial_map - calculate spatial map
     * @param img - source matrix
     * @param padding - padding size
     * @return
     */
    Matrix spatial_map(Matrix& img, int padding);

    /**
     * @brief pad - create matrix padded with mirrored rows/columns
     * @param src - source matrix
     * @param padding - padding size
     * @return
     */
    Matrix pad(Matrix& src,int padding);

    /**
     * @brief printVector - print vector to stdout
     * @param v
     */
    void printVector(std::vector<s3real>& v);

    /**
     * @brief The eo_map_struct - structure for holding eo_polaraverage grid
     */
    struct eo_map_struct
    {
        int f11,f12,f21,f22;
        s3real ex,ey;
    };

    /**
     * @brief m_eo_map - array of precomputed points for optimization of "eo_polaraverage"
     */
    std::vector<eo_map_struct> m_eo_map;

    /**
     * @brief fft - fft implementation
     */
    FFTWrapper fft;

    /**
     * @brief m_window - window to apply before fft (hanning)
     */
    Matrix m_window;
};

#endif // S3ESTIMATOR_H
