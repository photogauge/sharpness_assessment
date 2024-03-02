#include "s3tools.h"

#include <cstring>
#include <cassert>
#include <iostream>
#include <cmath>
#include <sstream>
#include <fstream>
#include <iostream>
#include <math.h>
#include <string>
#include <filesystem>
using namespace S3Tools;

Matrix::Matrix(int rowCount,int columnCount):
    vector<s3real>(static_cast<unsigned long>(rowCount*columnCount)),
    name("matrix"),
    m_rows(rowCount),
    m_columns(columnCount)
{
}
Matrix::Matrix(const Matrix& src):vector<s3real>(src),name(src.name),m_rows(src.m_rows),m_columns(src.m_columns)
{

}
Matrix& Matrix::operator=(const Matrix& m)
{
    if (this == &m) {
        return *this;
    }
    name=m.name;
    resize(static_cast<unsigned long>(m.count()));
    memcpy(this->data(),m.data(),static_cast<unsigned long>(m.count())*sizeof(s3real));
    this->m_rows=m.m_rows;
    this->m_columns=m.m_columns;
    return *this;
}
Matrix::~Matrix()
{
  //  std::cout<<name<<" deleted"<<std::endl;
}
s3real Matrix::mean2() const
{
    s3real sum=0;
    int cnt=count();
    s3real* d=const_cast<s3real*>(data());
    for(int i=0;i<cnt;i++)
    {
        sum+=*d++;
    }
    return sum/=cnt;
}
s3real Matrix::std2(s3real mean) const
{
    s3real sum=0;
    s3real* d=const_cast<s3real*>(data());
    int cnt=count();
    for(int i=0;i<cnt;i++)
    {
        s3real tmp=*d++-mean;
        sum+=tmp*tmp;
    }
    return static_cast<s3real>(sqrt(sum/(cnt-1)));
}
s3real Matrix::max() const
{
    if(!count())
        return -1;
    s3real* d=const_cast<s3real*>(data());
    s3real m=d[0];

    int cnt=count();
    for(int i=1;i<cnt;i++)
    {
        m=std::max(m,d[i]);
    }
    return m;
}
s3real Matrix::std2() const
{
    return std2(mean2());
}

s3real Matrix::diagSum() const
{
    int len=std::min(m_rows,m_columns);
    s3real sum=0;
    for(int i=0;i<len;i++)
    {
        sum+=at(i,i);
    }
    return sum;
}
void  Matrix::zero()
{
    memset(data(),0,sizeof(s3real)*static_cast<unsigned long>(m_rows*m_columns));
}
void Matrix::eye()
{
    zero();
    int s=std::min(m_rows,m_columns);
    for(int i=0;i<s;i++)
    {
        at(i,i)=1;
    }
}
s3real&  Matrix::at(int row, int col)
{
    assert(row<m_rows && col<m_columns);
    return std::vector<s3real>::at(m_columns*row+col);
}
const s3real&  Matrix::at(int row, int col) const
{
    assert(row<m_rows && col<m_columns);
    return std::vector<s3real>::at(m_columns*row+col);
}
void  Matrix::print()
{
    std::cout<<name<<" "<<m_rows<<"x"<<m_columns<<std::endl;
    s3real* d=data();
    for(int r=0;r<m_rows;r++)
    {
        for(int c=0;c<m_columns;c++)
        {
            std::cout<<*d++<<" ";
        }
        std::cout<<std::endl;
    }
}
s3real*  Matrix::dataAt(int row, int col)
{
    assert(row<m_rows && col<m_columns);
    return data() + m_columns*row+col;
}
Matrix   Matrix::Matrix::block(int row, int col, int rowCount, int colCount) const
{
    assert((row+rowCount)<=m_rows);
    assert((col+colCount)<=m_columns);

    Matrix out(rowCount,colCount);
    for(int i=0;i<rowCount;i++)
    {
        s3real* d=const_cast<Matrix*>(this)->dataAt(i+row,col);
        memcpy(out.dataAt(i,0),d,sizeof(s3real)*static_cast<unsigned long>(colCount));
    }
    return out;
}
int Matrix::count() const
{
    return m_rows*m_columns;
}

int Matrix::rows() const
{
    return m_rows;
}
int Matrix::columns() const
{
    return m_columns;
}
