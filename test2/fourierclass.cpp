#include "fourierclass.h"

FourierClass::FourierClass()
{
    fourier_error = 0;
}

int FourierClass::getFourierError()
{
    return fourier_error;
}

bool FourierClass::isLenSuitFft(const int buf_len)
{
    int len = buf_len;
    if (len <= 0)
    {
        fourier_error = -2;
        return false;
    }
    if ((len & (len - 1)) != 0)
    {
        return false;
    }
    else
    {
        return true;
    }
}

/****************************************************************************************

    ΪDIF FFT����ĵ����ӳ���
    example:
        input:
            buf_in[8] = {0,1,2,3,4,5,6,7}
        output:
            buf_in[8] = (0,4,2,6,1,5,3,7)

*****************************************************************************************/
void FourierClass::computeInverse(double *buf, const int buf_len)
{
    if (buf == nullptr)
    {
        fourier_error = -1;
        return;
    }
    if (buf_len <= 0)
    {
        fourier_error = -2;
        return;
    }


    int j = buf_len / 2;
    int k = 0;
    double temp = 0;

    for(int i = 1; i <= buf_len -2; i++)
    {
        if (i < j)
        {
            temp = buf[i];
            buf[i] = buf[j];
            buf[j] = temp;
        }
        k = buf_len / 2;

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


//inv_flag = false ����FFT���任 = true����IFFT���任
void FourierClass::computeDIFFft(double *buf_r, double *buf_i, const int buf_len, const bool inv_flag)
{
    if ((buf_r == nullptr) || (buf_i == nullptr))
    {
        fourier_error = -1;
        return;
    }

    if (buf_len <= 0)
    {
        fourier_error = -2;
        return;
    }

    if (!isLenSuitFft(buf_len))
    {
        fourier_error = -3;
        return;
    }
//    if (fft_len & (fft_len - 1) != 0)
//    {

//    }
    //fft_lenΪfft�ĵ���,MΪ2�Ĵη���
    int M = 1, nn = buf_len;
    int B = 0;
    double P = 0;
    double Wr = 0, Wi = 0;
    double Tr = 0, Ti = 0;

    while(nn / 2 != 1)
    {
        M = M + 1;
        nn = nn / 2;

    }         //��m��ֵ
//    qDebug() << M << "MMMM";
    if (inv_flag)
    {
        for(int i = 0; i < buf_len; i++)
        {
            buf_i[i] = -buf_i[i];           //���inv=-1,ȡ����
        }
    }
    computeInverse(buf_r, buf_len);       //���õ����Ӻ���Inverse,����

    computeInverse(buf_i, buf_len);

    for (int i = 1; i <= M; i++)                     //DIF-FFT Decimation-In-Time ��2-FFT
    {
//        B = 2 << (i - 1);
//        B = (int)std::pow((double)2, (i - 1));
        B = 1 << (i - 1);
        for (int j = 0; j <= B - 1; j++)
        {
            int k = 0;
//            P = j * std::pow((double)2, (M - i));

            P = j * (1 << (M - i));     //������ת����
//            qDebug() << M-i << "aboutMMMM";
//            qDebug() << j << ":" << P;
            Wr = std::cos(2 * MY_PI * (double)(P / buf_len));
            Wi = std::sin(2 * MY_PI * (double)(P / buf_len));

//            qDebug() << "Wr[" << j << "]=" << Wr;
//            qDebug() << "Wi[" << j << "]=" << Wi;
            //sin2PI and cos 1.5PI ����� ���˲�
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
            for (k = j; k < buf_len - 1; k = k + (1 << i))//(int)std::pow((double)2, i))//k = k + (2 << i))
            {
                Tr = buf_r[k + B] * Wr + buf_i[k + B] * Wi;

                Ti = buf_i[k + B] * Wr - buf_r[k + B] * Wi;
//                qDebug() << "i" << i;
//                qDebug() << "B" << B;
//                qDebug() << "j" << j;
//                qDebug() << "k" << k;
//                qDebug() << "Tr[" << k + B << "]=" << Tr;
//                qDebug() << "Ti[" << k + B << "]=" << Ti;
                //Xr[k]Ϊʵ��,Xi[k]Ϊ�鲿
                buf_r[k + B] = buf_r[k] - Tr;
                buf_i[k + B] = buf_i[k] - Ti;

                buf_r[k] = buf_r[k] + Tr;
                buf_i[k] = buf_i[k] + Ti;
            }
        }
    }
    if(inv_flag)
    {
        for(int i = 0; i < buf_len; i++)
        {
//            buf_i[i] = -buf_i[i];           //���inv=-1,ȡ����
            //���inv=-1,Ҫ��X[N]ǰ����1/N
            buf_r[i] = (1.0 / buf_len) * buf_r[i];
            buf_i[i] = (1.0 / buf_len) * buf_i[i];
        }
    }
}

//inv_flag = false ����DFT���任 = true����IDFT���任
void FourierClass::computeDFT(double *buf_r, double *buf_i, const int buf_len, const bool inv_flag)
{
    if ((buf_r == nullptr) || (buf_i == nullptr))
    {
        fourier_error = -1;
        return;
    }

    if (buf_len <= 0)
    {
        fourier_error = -2;
        return;
    }

    double Num1, Num2;
    double *buf_temp_r = new double[buf_len];
    double *buf_temp_i = new double[buf_len];

    for(int i = 0; i < buf_len; i++)
    {
        buf_temp_r[i] = 0;
        buf_temp_i[i] = 0;
    }
    if (!inv_flag)
    {
        for(int i = 0; i < buf_len; i++)
        {
            Num1 = 0.0;
            Num2 = 0.0;
            for(int j = 0; j < buf_len; j++)
            {
//                if (i == 1)
//                {
//                    qDebug() << cos(2 * (MY_PI / buf_len) * i * j) << endl
//                             << "Num1" << Num1;
//                }
                //����Ҷ�任
                Num1 += cos(2 * (MY_PI / buf_len) * i * j) * buf_r[j];
                Num2 += -sin(2 * (MY_PI / buf_len) * i * j) * buf_r[j];
//                if (i == 1)
//                {
//                    qDebug() << cos(2 * (MY_PI / buf_len) * i * j) << endl
//                             << "Num1" << Num1;
//                }
            }
            buf_temp_r[i] = Num1;
            buf_temp_i[i] = Num2;
        }
        for(int i = 0; i < buf_len; i++)
        {
            buf_r[i] = buf_temp_r[i];
            buf_i[i] = buf_temp_i[i];
//            buf_i[i] = -buf_i[i];           //���inv=-1,ȡ����
        }
    }
    else
    {
        for(int i = 0; i < buf_len; i++)
        {
            Num1 = 0.0;
            Num2 = 0.0;
            for(int j = 0; j < buf_len; j++)
            {
                //����Ҷ���任
                Num1 += buf_r[j] * cos(2 * MY_PI / buf_len * i * j)
                        - buf_i[j] * sin(2 * MY_PI / buf_len * i * j);
                Num2 += buf_r[j] * sin(2 * MY_PI / buf_len * i * j)
                        + buf_i[j] * cos(2 * MY_PI / buf_len * i * j);
            }
            buf_temp_r[i] = Num1 / buf_len;
            buf_temp_i[i] = Num2 / buf_len;
        }
        for(int i = 0; i < buf_len; i++)
        {
            buf_r[i] = buf_temp_r[i];
            buf_i[i] = buf_temp_i[i];
            buf_i[i] = -buf_i[i];           //���inv=-1,ȡ����
        }
    }



    delete [] buf_temp_r;
    buf_temp_r = nullptr;
    delete [] buf_temp_i;
    buf_temp_i = nullptr;
}
