#ifndef S3TOOLS_H
#define S3TOOLS_H


#include <vector>
#include <string>

/**
 * @brief s3real data type used for calulations float/double
 */
typedef double s3real;

namespace S3Tools
{
/**
     * @brief The Matrix class
     */
    class Matrix: public std::vector<s3real>
    {
    public:
        Matrix(int rowCount=0,int columnCount=0);
        Matrix(const Matrix& src);
        Matrix& operator=(const Matrix&);
        ~Matrix();
        /**
         * @brief zero - init matrix with zeros
         */
        void zero();
        /**
         * @brief at - get reference to value at (row,col)
         * @param row
         * @param col
         * @return
         */
        s3real& at(int row, int col);
        /**
         * @brief at - get const reference to value at (row,col)
         * @param row
         * @param col
         * @return
         */
        const s3real& at(int row, int col) const;
        /**
         * @brief print - print matrix contents to stdout
         */
        void print();
        /**
         * @brief dataAt - get pointer to data at (row,col)
         * @param row
         * @param col
         * @return
         */
        s3real* dataAt(int row, int col);
        /**
         * @brief block - create matrix as part of existing one
         * @param row
         * @param col
         * @param rowCount
         * @param colCount
         * @return
         */
        Matrix block(int row, int col, int rowCount, int colCount) const;
        int rows() const;
        int columns() const;
        /**
         * @brief count - get elements count
         * @return
         */
        int count() const;
        /**
         * @brief mean2 - get matix mean value
         * @return
         */
        s3real mean2() const;
        /**
         * @brief mean2 - get matrix std, using supplied mean value
         * @param mean
         * @return
         */
        s3real std2(s3real mean) const;
        /**
         * @brief std2 - get matrix std
         * @return
         */
        s3real std2() const;
        /**
         * @brief max - matrix max value
         * @return
         */
        s3real max() const;
        /**
         * @brief eye - init eye matrix
         */
        void eye();
        /**
         * @brief diagSum - calculate sum of main diagonal
         * @return
         */
        s3real diagSum() const;
        /**
         * @brief name - matrix name used in debug output
         */
        std::string name;
    protected:
        /**
         * @brief m_rows  - rows count
         */
        int m_rows;
        /**
         * @brief m_columns - column count
         */
        int m_columns;
    };
}

#endif // S3TOOLS_H
