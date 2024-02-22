#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
#include "s3tools.h"
#include "s3estimator.h"
#include "iostream"
#include <string>


#define LOG_TAG "ImageQuality-Native"

using namespace std;
using namespace cv;

#ifdef __cplusplus
extern "C" {
#endif

// JNIEXPORT void JNICALL Java_com_photogauge_opencv_OpenCV_imageQualityNew
//         (JNIEnv *, jclass, jstring, jstring, jdouble, jdouble);

#ifdef __cplusplus
}
#endif

S3Estimator m_estimator;

Matrix matToMatrix(Mat image)
{
    Matrix img(image.rows, image.cols);
    for(int r=0;r<img.rows();r++)
    {
        for(int c=0;c<img.columns();c++)
        {
            img.at(r,c) = image.at<uchar>(r,c);
        }
    }
    return img;
}

Mat matrixToMat(Matrix image)
{
    Mat img(image.rows(), image.columns(), CV_8UC1);
    for(int r=0; r<img.rows; r++)
    {
        for(int c=0; c<img.cols;c++)
        {
            float v = std::max((s3real)0.,image.at(r,c))*255.;
            img.at<uchar>(r,c) = (unsigned char)v;
        }
    }

    //img = 255 - img;//reverting the image here
    return img;
}

Mat merge_image(Mat img, Mat mask)
{
    //Mat img(5, 5, CV_64FC3); // declare three channels image
    Mat ch1, ch2, ch3; // declare three matrices
    // "channels" is a vector of 3 Mat arrays:
    vector<Mat> channels(3);
    // split img:
    split(img, channels);
    // get the channels (follow BGR order in OpenCV)
    add(channels[1], 255, channels[1], mask = mask);
    //ch2 = channels[1];
    //ch3 = channels[2];
    // modify channel// then merge

    merge(channels, img);

    return img;
}

Mat matrixToMatS1(Matrix image)
{
    Mat img(image.rows(), image.columns(), CV_8UC1);
    unsigned char redColor = (unsigned char) 16711680;
    for(int r=0; r<img.rows; r++)
    {
        for(int c=0; c<img.cols;c++)
        {
            float v = std::max((s3real)0.,image.at(r,c))*255.;
            if((int)v < 128){
                img.at<uchar>(r,c) = (unsigned char)(redColor);
            } else {
                img.at<uchar>(r,c) = (unsigned char)v;
            }
        }
    }

    return img;
}

void getSharpnessMap(Matrix src_mat, Mat* output_mat){
    Mat s_mat = matrixToMat(src_mat);
    // Mat norm_image;
    // double min;
    // double max;
    // minMaxLoc(s_mat, &min, &max);
    // cout<<"min "<< min<<endl;
    // cout<<"max "<< max<<endl;
    // normalize(s_mat, norm_image, 0, 255, NORM_MINMAX);

    // Refer this link for various colormaps: https://learnopencv.com/applycolormap-for-pseudocoloring-in-opencv-c-python/

    applyColorMap(s_mat, *output_mat, COLORMAP_JET);


}
void matrixStats(const Matrix& m)
{
    s3real mean=m.mean2();
    //LOGE("mean: %lf std: %lf diag sum: %lf", mean, m.std2(mean), m.diagSum());
    printf("Photogauge Sharpness map mean=%lf,std=%lf,sum=%lf", mean, m.std2(mean), m.diagSum());
}

// Resize input image from max height
Mat resize_with_max_height(Mat InputImage, int MaxHeight){
    double resizeVal = (double)MaxHeight/(double)InputImage.rows;
    Mat OutputImage;
    if (resizeVal > 1.0){
        resize(InputImage, OutputImage, Size(), resizeVal, resizeVal, INTER_LINEAR);
    }
    else{
        resize(InputImage, OutputImage, Size(), resizeVal, resizeVal, INTER_AREA);
    }
    return OutputImage;

}

string getFileName(const string& s) {

	char sep = '/';

#ifdef _WIN32
	sep = '\\';
#endif

	size_t i = s.rfind(sep, s.length());
	if (i != string::npos) {
		return(s.substr(i + 1, s.length() - i));
	}

	return("");
}

void assess_sharpness(const char* imagePath, const char* outPath, float scale_factor) {
	cout << "..........." << imagePath << "\n" << "........" << outPath;

    // read input image
    Mat src_color;
	vector<String> filenames;
	glob(imagePath, filenames);
	//cout << "............" << filenames.size();
	for (int i = 0; i < filenames.size(); i++)
	{
		src_color = imread(filenames[i]);
		size_t lastdot = getFileName(filenames[i]).find_last_of(".");
		string write_name = getFileName(filenames[i]).substr(0, lastdot);

		cout << write_name;
		// Resize input image
		Mat dst_color;
		//dst_color = resize_with_max_height(src_color, scale_factor);
		resize(src_color, dst_color, Size(0, 0), scale_factor, scale_factor, cv::INTER_LANCZOS4);

		// COnvert to grayscale
		Mat  dst_gray;
		cv::cvtColor(dst_color, dst_gray, cv::COLOR_BGR2GRAY);

		// Get the sharpness scores
		Matrix img = matToMatrix(dst_gray);
		m_estimator.s3_map(img);

		// Get S1 score
		Mat s1;
		cout << "Getting S1 map" << endl;
		getSharpnessMap(m_estimator.m_s1, &s1);
		string s1_path = outPath + write_name + "_s1.jpg";
		imwrite(s1_path, s1);

		// Get S2 score
		Mat s2;
		cout << "Getting S2 map" << endl;
		getSharpnessMap(m_estimator.m_s2, &s2);
		string s2_path = outPath + write_name + "_s2.jpg";
		imwrite(s2_path, s2);

		// Get S3 score
		Mat s3;
		string s3_path = outPath + write_name + "_s3.jpg";
		cout << "Getting S3 map" << endl;
		getSharpnessMap(m_estimator.m_s3, &s3);
		imwrite(s3_path, s3);
	}
}