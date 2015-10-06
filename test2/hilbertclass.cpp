#include "hilbertclass.h"

HilbertClass::HilbertClass()
{
    hilbert_error = -2;
}

HilbertClass::HilbertClass(const int buf_len):FourierClass()
{
//    if (buf_in == nullptr)
//    {
//        hilbert_error = -1;
//        return;
//    }
    hilbert_len = buf_len;
    if (hilbert_len <= 0)
    {
        hilbert_error = -2;
    }
    else
    {
        hilbert_error = 0;
    }

    hilbert_buf_in = new double [hilbert_len];
    hilbert_factor = new double [hilbert_len];
    //Ï£ÊÏÒò×ÓµÄÆµÆ×1<=i<=N/2-1,H[i]=2,f>0;H[i]=0,f<0;


//    for (int i = 0; i < hilbert_len; i++)
//    {
//        hilbert_buf_in[i] = buf_in[i];
//        hilbert_factor[i] = 0;
////        if ((i >= 1) && (i <= hilbert_len / 2 - 1))
////        {
////            hilbert_factor[i] = 2;
////        }
////        else if (i == hilbert_len / 2)
////        {
////            hilbert_factor[i] = 1;
////        }
////        else if (i > hilbert_len / 2)
////        {
////            hilbert_factor[0] = 0;
////        }
////        qDebug() << "o" <<hilbert_factor[i];
//    }
//    hilbert_factor[0] = 1;
    buildHilbertFactor();
//    showHilbertData();
}


HilbertClass::~HilbertClass()
{
    delete [] hilbert_buf_in;
    hilbert_buf_in = nullptr;
    delete [] hilbert_factor;
    hilbert_factor = nullptr;
}


int HilbertClass::getHilbertError()
{
    return hilbert_error;
}

int HilbertClass::getHilbertLen()
{
    return hilbert_len;
}

void HilbertClass::showHilbertData()
{
    for (int i = 0; i < hilbert_len; i++)
    {
        qDebug() << "hilbert_buf_in[" << i << "]=" << hilbert_buf_in[i];
    }
    for (int i = 0; i < hilbert_len; i++)
    {
        qDebug() << "hilbert_factor[" << i << "]=" << hilbert_factor[i];
    }
}


void HilbertClass::computeHilbert(const double *buf_in, double *buf_out)
{
    if ((buf_in == nullptr)||(buf_out == nullptr))
    {
        hilbert_error = -1;
        return;
    }
    for (int i = 0; i < hilbert_len; i++)
    {
        hilbert_buf_in[i] = buf_in[i];
    }
//    if (buf_out_len != hilbert_len)
//    {
//        hilbert_error = -4;
//        return;
//    }

    for (int i = 0; i < hilbert_len; i++)
    {
        buf_out[i] = 0;
    }
    if (isLenSuitFft(hilbert_len))
    {
//        qDebug() << "why";
        computeDIFFft(hilbert_buf_in, buf_out, hilbert_len, false);//FFT
        for (int i = 0; i < hilbert_len; i++)
        {
            hilbert_buf_in[i] = hilbert_buf_in[i] * hilbert_factor[i];
            buf_out[i] = buf_out[i] * hilbert_factor[i];
        }
        computeDIFFft(hilbert_buf_in, buf_out, hilbert_len, true);//IFFT
    }
    else
    {
//        qDebug() << "what";
        computeDFT(hilbert_buf_in, buf_out, hilbert_len, false);//DFT

//        for(int i = 0; i < 8; i++)
//        {
//            qDebug() << "hilbert_buf_in[" << i << "]=:" <<hilbert_buf_in[i];
//        }
//        for(int i = 0; i < 8; i++)
//        {
//            qDebug() << "buf_out[" << i << "]=:" <<buf_out[i];
//        }
        for (int i = 0; i < hilbert_len; i++)
        {
            hilbert_buf_in[i] = hilbert_buf_in[i] * hilbert_factor[i];
            buf_out[i] = buf_out[i] * hilbert_factor[i];
        }
        computeDFT(hilbert_buf_in, buf_out, hilbert_len, true);//IDFT
    }

}


void HilbertClass::buildHilbertFactor()
{
    if (hilbert_len == hilbert_len / 2 * 2)
    {
        for (int i = 0; i < hilbert_len; i++)
        {
            if (i == 0)
            {
                hilbert_factor[i] = 1;
            }
            if ((i >= 1) && (i <= hilbert_len / 2 - 1))
            {
                hilbert_factor[i] = 2;
            }
            else if (i == hilbert_len / 2)
            {
                hilbert_factor[i] = 1;
            }
            else if (i > hilbert_len / 2)
            {
                hilbert_factor[0] = 0;
            }
//            qDebug() << "o" <<hilbert_factor[i];
        }
        hilbert_factor[0] = 1;

    }
    else
    {
        for (int i = 0; i < hilbert_len; i++)
        {
            if (i == 0)
            {
                hilbert_factor[i] = 1;
            }
            if ((i >= 1) && (i <= (hilbert_len - 1) / 2))
            {
                hilbert_factor[i] = 2;
            }
//            else if (i == hilbert_len / 2)
//            {
//                hilbert_factor[i] = 1;
//            }
            else if (i > (hilbert_len - 1) / 2)
            {
                hilbert_factor[0] = 0;
            }
//            qDebug() << "o" <<hilbert_factor[i];
        }
        hilbert_factor[0] = 1;
    }
}


void HilbertClass::computeInstantFrequency(const double *buf_in, double * buf_out)
{
    if ((buf_in == nullptr)||(buf_out == nullptr))
    {
        hilbert_error = -1;
        return;
    }
    double * buf_i = new double [hilbert_len];
    computeHilbert(buf_in, buf_i);
    double re = 0;
    double im = 0;
    for (int i = 1; i < hilbert_len - 1; i++)
    {
        re = hilbert_buf_in[i + 1] * hilbert_buf_in[i - 1] + buf_i[i + 1] * buf_i[i - 1];
        im = hilbert_buf_in[i + 1] * buf_i[i - 1] - hilbert_buf_in[i - 1] * buf_i[i + 1];
//        double a = im / re;
        buf_out[i] = 250 * std::atan2(im ,re) / MY_PI;
    }

    delete [] buf_i;
    buf_i = nullptr;
}
