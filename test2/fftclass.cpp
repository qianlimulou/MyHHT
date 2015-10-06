#include "fftclass.h"

fftClass::fftClass()
{
    fft_error = 0;
    fft_len = 0;
}

fftClass::fftClass(const int buf_len)
{
    fft_len = buf_len;
    if (fft_len <= 0)
    {
        fft_error = -2;
    }
    else
    {
        fft_error = 0;
    }
    int my_len;
    if ((fft_len & (fft_len - 1)) != 0)
    {
        while (my_len <= 32)
        {
            if ((1 << my_len) > fft_len)
            {
                break;
            }
            my_len++;
        }
    }
    fft_len = my_len;

}

int fftClass::getFftError()
{
    return fft_error;
}

int fftClass::getFftLen()
{
    return fft_len;
}

/****************************************************************************************

    为DIF FFT服务的倒序子程序
    example:
        input:
            buf_in[8] = {0,1,2,3,4,5,6,7}
        output:
            buf_in[8] = (0,4,2,6,1,5,3,7)

*****************************************************************************************/
void fftClass::computeInverse(double * buf)//TEST PASS
{
    if (buf == nullptr)
    {
        fft_error = -1;
        return;
    }
    if (fft_len == 0)
    {
        fft_error = -2;
        return;
    }


    int j = fft_len / 2;
    int k = 0;
    double temp = 0;

    for(int i = 1; i <= fft_len -2; i++)
    {
        if (i < j)
        {
            temp = buf[i];
            buf[i] = buf[j];
            buf[j] = temp;
        }
        k = fft_len / 2;

        if (j < k)
        {
            j = j + k;
        }
        else
        {
            while (j >= k)
            {
                j = j - k;
                k = k / 2;
            }
            j = j + k;
        }
    }
}

//buf_r为实部，buf_i为虚部，inv_flag =0为fft; =1为ifft;
void fftClass::computeDIFFft(double * buf_r,double * buf_i, const int inv_flag)
{
    if ((buf_r == nullptr) || (buf_i == nullptr))
    {
        fft_error = -1;
        return;
    }

    if (fft_len == 0)
    {
        fft_error = -2;
        return;
    }


//    if (fft_len & (fft_len - 1) != 0)
//    {

//    }
    //fft_len为fft的点数,M为2的次方数
    int M = 1, nn = fft_len;
    int B = 0;
    double P = 0;
    double Wr = 0, Wi = 0;
    double Tr = 0, Ti = 0;

    while(nn / 2 != 1)
    {
        M = M + 1;
        nn = nn / 2;

    }         //求m的值
//    qDebug() << M << "MMMM";
    if (inv_flag)
    {
        for(int i = 0; i < fft_len; i++)
        {
            buf_i[i] = -buf_i[i];           //如果inv=-1,取共轭
        }
    }
    computeInverse(buf_r);       //调用倒序子函数Inverse,倒序

    computeInverse(buf_i);

    for (int i = 1; i <= M; i++)                     //DIF-FFT Decimation-In-Time 基2-FFT
    {
//        B = 2 << (i - 1);
//        B = (int)std::pow((double)2, (i - 1));
        B = 1 << (i - 1);
        for (int j = 0; j <= B - 1; j++)
        {
            int k = 0;
//            P = j * std::pow((double)2, (M - i));

            P = j * (1 << (M - i));     //计算旋转因子
//            qDebug() << M-i << "aboutMMMM";
//            qDebug() << j << ":" << P;
            Wr = std::cos(2 * pi * (double)(P / fft_len));
            Wi = std::sin(2 * pi * (double)(P / fft_len));

//            qDebug() << "Wr[" << j << "]=" << Wr;
//            qDebug() << "Wi[" << j << "]=" << Wi;
            //sin2PI and cos 1.5PI 很奇怪 做滤波
//            if ((Wr > -2.5e-16) && (Wr < 2.5e-16))
//            {
//                Wr = 0;
//            }
//            if ((Wi > -2.5e-16) && (Wi < 2.5e-16))
//            {
//                Wi = 0;
//            }
//            qDebug() << i << ":" << P;
//            qDebug() << "oooooWr[" << j << "]=" << Wr;
//            qDebug() << "oooooWi[" << j << "]=" << Wi;
            for (k = j; k < fft_len - 1; k = k + (1 << i))//(int)std::pow((double)2, i))//k = k + (2 << i))
            {
                Tr = buf_r[k + B] * Wr + buf_i[k + B] * Wi;

                Ti = buf_i[k + B] * Wr - buf_r[k + B] * Wi;
//                qDebug() << "i" << i;
//                qDebug() << "B" << B;
//                qDebug() << "j" << j;
//                qDebug() << "k" << k;
//                qDebug() << "Tr[" << k + B << "]=" << Tr;
//                qDebug() << "Ti[" << k + B << "]=" << Ti;
                //Xr[k]为实部,Xi[k]为虚部
                buf_r[k + B] = buf_r[k] - Tr;
                buf_i[k + B] = buf_i[k] - Ti;

                buf_r[k] = buf_r[k] + Tr;
                buf_i[k] = buf_i[k] + Ti;
            }
        }
    }
//    qDebug() << "fft_why1";
//    for (int i = 0; i < fft_len; i++)
//    {
//        qDebug() << "buf_r[" << i << "] =" <<  buf_r[i];
//    }
//    for (int i = 0; i < fft_len; i++)
//    {
//        qDebug() << "buf_i[" << i << "] =" <<  buf_i[i];
//    }

    if(inv_flag)
    {
        for(int i = 0; i < fft_len; i++)
        {
//            buf_i[i] = -buf_i[i];           //如果inv=-1,取共轭
            //如果inv=-1,要在X[N]前乘以1/N
            buf_r[i] = (1.0 / fft_len) * buf_r[i];
            buf_i[i] = (1.0 / fft_len) * buf_i[i];
        }
    }
}


