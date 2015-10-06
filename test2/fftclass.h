#ifndef FFTCLASS_H
#define FFTCLASS_H

#include <QDebug>
#include <math.h>
//#include <complex>

#define pi 3.14159265358979323846

class fftClass
{
public:
    fftClass();
    fftClass(const int buf_len);
    int getFftError();
    int getFftLen();
    void computeInverse(double * buf);
    void computeDIFFft(double *buf_r, double *buf_i, const int inv_flag);
private:
    int fft_error;
    int fft_len;

};

#endif // FFTCLASS_H
