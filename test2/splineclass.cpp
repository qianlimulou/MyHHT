#include "splineclass.h"

SplineClass::SplineClass()
{
    spline_out_len = 0;
//    spline_buf_out_x = new double [1];
//    spline_buf_out_y = new double [1];
    spline_con = NOTaKNOT;
    spline_boudary_A = 0;
    spline_boudary_B = 0;
    spline_error = 0;
}
/*
SplineClass::SplineClass(double *buf_in_x, double *buf_in_y, int in_len,
                         double *buf_out_x, double *buf_out_y, int out_len)
{
    if ((in_len <= 1)||(out_len <=1 ))
    {
        qDebug() << "in_len or out_len must > 1";
        return ;
    }
    spline_buf_in_x = new double [in_len];
    spline_buf_in_x = buf_in_x;
    spline_buf_in_y = new double [in_len];
    spline_buf_in_y = buf_in_y;
    spline_in_len = in_len;
    spline_buf_out_x = new double [out_len];
    spline_buf_out_x = buf_out_x;
    spline_out_len = out_len;
    spline_buf_out_y = new double [out_len];
    spline_buf_out_y = buf_out_y;
    spline_con = NOTaKNOT;
}

SplineClass::SplineClass(double *buf_in_x, double *buf_in_y, int in_len,
                         double *buf_out_x, double *buf_out_y, int out_len,
                         SplineClass::spline_boudary_condition con, int a, int b)
{
    if ((in_len <= 1)||(out_len <=1 ))
    {
        qDebug() << "in_len or out_len must > 1";
        return ;
    }
    spline_buf_in_x = new double [in_len];
    spline_buf_in_x = buf_in_x;
    spline_buf_in_y = new double [in_len];
    spline_buf_in_y = buf_in_y;
    spline_in_len = in_len;
    spline_buf_out_x = new double [out_len];
    spline_buf_out_x = buf_out_x;
//    spline_buf_out_y = buf_out_y;
    spline_out_len = out_len;
    spline_buf_out_y = new double [out_len];
    buf_out_y = spline_buf_out_y;
    spline_con = con;
    spline_boudary_A = a;
    spline_boudary_B = b;
}
*/
SplineClass::~SplineClass()
{
//    delete [] spline_buf_out_x; spline_buf_out_x = nullptr;
//    delete [] spline_buf_out_y; spline_buf_out_y = nullptr;
}

int SplineClass::getSplineError()
{
    return spline_error;
}

int SplineClass::getSplineOutLen()
{
    return spline_out_len;
}

//void SplineClass::showAllSplineData()
//{
//    for (int i = 0; i < spline_out_len; i++)
//    {
////        qDebug() << "spline_out_x[" << i << "]=:" << spline_buf_out_x[i];
//    }

//    for (int i = 0; i < spline_out_len; i++)
//    {
////        qDebug() << "spline_out_y[" << i << "]=:" << spline_buf_out_y[i];
//    }
//}

void SplineClass::setThreeSplineCondition(SplineClass::spline_boudary_condition con, int a, int b)
{
    spline_con = con;
    spline_boudary_A = a;
    spline_boudary_B = b;
}

void SplineClass::computeThreeSpline(const int *buf_in_x, const double *buf_in_y, const int in_len, const int start_pos,
                                     /*double *buf_out_x, */double *buf_out_y, const int out_len)
{
    if ((in_len <= 1) || (out_len <=1 ) || (in_len >= out_len))
    {
        spline_error = -2;
        return;
    }
    if ((start_pos < 0) || (start_pos > out_len))
    {
        spline_error = -3;
        return;
    }
    if ((buf_in_x == nullptr) || (buf_in_y == nullptr) || (buf_out_y == nullptr))
    {
        spline_error = -1;
        return;
    }
    if (in_len <= 3)
    {
        for (int i = 0; i < out_len; i ++)
        {
            buf_out_y[i] = buf_in_y[0];
        }
        return;
    }
    spline_out_len = out_len;
//// need thinking
//    if ((buf_in_x[in_len - 1] - buf_in_x[0]) > out_len)
//    {
//        spline_error = -2;
//        return;
//    }
//    qDebug() << "what!!!";
//    delete [] spline_buf_out_x; //spline_buf_out_x = nullptr;
//    delete [] spline_buf_out_y; //spline_buf_out_y = nullptr;
//    qDebug() << "why";
//    spline_buf_out_x = new double [out_len];
//    spline_buf_out_x = buf_out_x;

//    spline_buf_out_y = new double [out_len];
//    spline_buf_out_y = buf_out_y;
//    qDebug() << "why???";

    int n = in_len - 1;//插值曲线个数
    int i, j;//循环计数
    double a0, b1, c2, d3;//三次曲线系数
    double temp; //交换临时存数
    double *a = new double [n + 1];
    double *b = new double [n + 1];
    double *c = new double [n + 1];
    double *d = new double [n + 1];
    double *f = new double [n + 1];
    double *m = new double [n + 1];
    double *h = new double [n + 1];

//    double interval_x ;
//    interval_x = (buf_in_x[in_len - 1] - buf_in_x[0]) / out_len;

//    for (i = 0; i < out_len; i ++)
//    {
//        spline_buf_out_x[i] = buf_in_x[0] + i * interval_x;
//    }

    /* 计算 h[] d[] */
    for (i = 0; i < n; i++)
    {
        h[i] = buf_in_x[i + 1] - buf_in_x[i];
        d[i] = (buf_in_y[i + 1] - buf_in_y[i]) / h[i];
        /* printf("%f\t%f\n", h[i], d[i]); */
    }

    /* 计算 a[] b[] d[] f[] */
    switch(spline_con)
    {
    case NATURAL:
    {
        f[0] = 0;
        f[n] = 0;
        a[0] = 0;
        c[n] = 0;
        c[0] = 0;
        a[n] = 0;
        b[0] = 1;
        b[n] = 1;
        break;
    }
    case CLAMPED:
    {
        f[0] = 6 * (d[0] - spline_boudary_A);
        f[n] = 6 * (spline_boudary_B - d[n - 1]);
        a[0] = 0;
        c[n] = 0;
        c[0] = h[0];
        a[n] = h[n - 1];
        b[0] = 2 * h[0];
        b[n] = 2 * h[n - 1];
        break;
    }
    case NOTaKNOT:
    {
        f[0] = 0;
        f[n] = 0;
        a[0] = -h[0];
        c[n] = -h[n - 1];
        c[0] = h[0] + h[1];
        a[n] = h[n - 2] + h[n - 1];
        b[0] = -h[1];
        b[n] = -h[n - 2];
        break;
    }
    }

    for (i = 1; i < n; i++)
    {
        a[i] = h[i - 1];
        b[i] = 2 * (h[i - 1] + h[i]);
        c[i] = h[i];
        f[i] = 6 * (d[i] - d[i - 1]);
    }

    /*************
    ** 高斯消元
    *************/
    /* 第0行到第(n-3)行的消元 */
    for(i = 0; i <= n - 3; i++)
    {
        /* 选主元 */
        if( fabs(a[i + 1]) > fabs(b[i]) )
        {
            temp = a[i + 1]; a[i + 1] = b[i]; b[i] = temp;
            temp = b[i + 1]; b[i + 1] = c[i]; c[i] = temp;
            temp = c[i + 1]; c[i + 1] = a[i]; a[i] = temp;
            temp = f[i + 1]; f[i + 1] = f[i]; f[i] = temp;
        }
        temp = a[i + 1] / b[i];
        a[i + 1] = 0;
        b[i + 1] = b[i + 1] - temp * c[i];
        c[i + 1] = c[i + 1] - temp * a[i];
        f[i + 1] = f[i + 1] - temp * f[i];
    }

    /* 最后3行的消元 */
    /* 选主元 */
    if( fabs(a[n - 1]) > fabs(b[n - 2]) )
    {
        temp = a[n - 1]; a[n - 1] = b[n - 2]; b[n - 2] = temp;
        temp = b[n - 1]; b[n - 1] = c[n - 2]; c[n - 2] = temp;
        temp = c[n - 1]; c[n - 1] = a[n - 2]; a[n - 2] = temp;
        temp = f[n - 1]; f[n - 1] = f[n - 2]; f[n - 2] = temp;
    }
    /* 选主元 */
    if( fabs(c[n]) > fabs(b[n - 2]) )
    {
        temp = c[n]; c[n] = b[n - 2]; b[n - 2] = temp;
        temp = a[n]; a[n] = c[n - 2]; c[n - 2] = temp;
        temp = b[n]; b[n] = a[n - 2]; a[n - 2] = temp;
        temp = f[n]; f[n] = f[n - 2]; f[n - 2] = temp;
    }
    /* 第(n-1)行消元 */
    temp = a[n - 1] / b[n - 2];
    a[n - 1] = 0;
    b[n - 1] = b[n - 1] - temp * c[n - 2];
    c[n - 1] = c[n - 1] - temp * a[n - 2];
    f[n - 1] = f[n - 1] - temp * f[n - 2];
    /* 第n行消元 */
    temp = c[n] / b[n - 2];
    c[n] = 0;
    a[n] = a[n] - temp * c[n - 2];
    b[n] = b[n] - temp * a[n - 2];
    f[n] = f[n] - temp * f[n - 2];
    /* 选主元 */
    if( fabs(a[n]) > fabs(b[n - 1]) )
    {
        temp = a[n]; a[n] = b[n - 1]; b[n - 1] = temp;
        temp = b[n]; b[n] = c[n - 1]; c[n - 1] = temp;
        temp = f[n]; f[n] = f[n - 1]; f[n - 1] = temp;
    }
    /* 最后一次消元 */
    temp = a[n] / b[n-1];
    a[n] = 0;
    b[n] = b[n] - temp * c[n - 1];
    f[n] = f[n] - temp * f[n - 1];

    /* 回代求解 m[] */
    m[n] = f[n] / b[n];
    m[n - 1] = (f[n - 1] - c[n - 1] * m[n]) / b[n-1];
    for(i = n - 2; i >= 0; i--)
    {
        m[i] = ( f[i] - (m[i + 2] * a[i] + m[i + 1] * c[i]) ) / b[i];
    }

    int start_pos_for = start_pos;
    for(i = start_pos_for; i < spline_out_len + start_pos_for; i++)
    {
        for(j = 0; j < n; j++)
        {
            if(i >= buf_in_x[j] && i <= buf_in_x[j + 1])
            {
                d3 = (m[j + 1] - m[j]) / (6 * h[j]);
                c2 = m[j] / 2;
                b1 = d[j] - (h[j] / 6) * (2 * m[j] + m[j + 1]);
                a0 = buf_in_y[j];

                buf_out_y[i - start_pos_for] =
                            d3 * (i - buf_in_x[j]) * (i - buf_in_x[j]) * (i - buf_in_x[j])
                          + c2 * (i - buf_in_x[j]) * (i - buf_in_x[j])
                          + b1 * (i - buf_in_x[j])
                          + a0;
//                qDebug() << "buf_out_y[" << i - start_pos_for << "]=" << buf_out_y[i - start_pos_for];
                //printf("mySpline%d=%f\n",i,(mySpline[i]));
            }
        }
    }
//    qDebug() <<"what???";
    delete [] a; a = nullptr;
//    qDebug() <<"why";
    delete [] b; b = nullptr;
    delete [] c; c = nullptr;
    delete [] d; d = nullptr;
    delete [] f; f = nullptr;
    delete [] m; m = nullptr;
    delete [] h; h = nullptr;

}
