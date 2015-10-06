#ifndef EXTRZEROCLASS_H
#define EXTRZEROCLASS_H

#include <QDebug>

#define EXTR_TOTAL_MAX_NUM 3

class ExtrZeroClass
{
public:
    ExtrZeroClass();
//    ExtrZeroClass(const double * buf_x, const double * buf_y, const int buf_in_len);
    ~ExtrZeroClass();
    int getExtrZeroError();

    int getExtrMaxNum();
    int getExtrMinNum();
    int getExtrTotalNum();
    int getZeroNum();
    bool isBufMonotonic(const double * buf_x, const int buf_in_len);
    bool isExtrNumOK();
    bool isExtrNumOtherZero();
    bool isExtrZeroIMF();

    void showAllExtrData();
    void showAllZeroData();
    void computeExtrResult(const double * buf_y, const int buf_in_len);
    void computeZeroResult(const double * buf_y, const int buf_in_len);
    int *extr_min_x;
    int *extr_max_x;
    int *zero_x;
private:
    int extrzero_error;
    int extr_min_x_num;
    int extr_max_x_num;
    int zero_x_num;
};

#endif // EXTRZEROCLASS_H
