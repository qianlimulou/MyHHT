#ifndef HILBERTCLASS_H
#define HILBERTCLASS_H

#include <QDebug>

#include "fourierclass.h"

class HilbertClass: public FourierClass
{
public:
    HilbertClass();
    HilbertClass(const int buf_len);
    ~HilbertClass();
    int getHilbertError();
    int getHilbertLen();
    void showHilbertData();
    void computeHilbert(const double *buf_in, double * buf_out);
    void computeInstantFrequency(const double *buf_in, double * buf_out);

private:
    void buildHilbertFactor();
    int hilbert_error;
    int hilbert_len;
    double * hilbert_buf_in;
    double * hilbert_factor;
};

#endif // HILBERTCLASS_H
