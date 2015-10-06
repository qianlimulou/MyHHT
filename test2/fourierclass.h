#ifndef FOURIERCLASS_H
#define FOURIERCLASS_H

#include <QDebug>
#include <math.h>
#define MY_PI 3.14159265358979323846

class FourierClass
{
public:
    FourierClass();
    int getFourierError();
    bool isLenSuitFft(const int buf_len);
    void computeInverse(double * buf, const int buf_len);
    void computeDIFFft(double *buf_r, double *buf_i, const int buf_len, const bool inv_flag);
    void computeDFT(double *buf_r, double *buf_i, const int buf_len, const bool inv_flag);
private:
    int fourier_error;

};

#endif // FOURIERCLASS_H
