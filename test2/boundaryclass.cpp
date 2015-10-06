#include "boundaryclass.h"

BoundaryClass::BoundaryClass()
{
    boundary_x_min = new int [1];
    boundary_x_max = new int [1];
    boundary_y_min = new double [1];
    boundary_y_max = new double [1];
    boundary_min_num = 0;
    boundary_max_num = 0;
    boundary_error = 0;
    boundary_min_start = 0;
    boundary_max_start = 0;
}

BoundaryClass::~BoundaryClass()
{
    delete [] boundary_x_min;
    boundary_x_min = nullptr;
    delete [] boundary_x_max;
    boundary_x_max = nullptr;
    delete [] boundary_y_min;
    boundary_y_min = nullptr;
    delete [] boundary_y_max;
    boundary_y_max = nullptr;
}

int BoundaryClass::getBoundaryError()
{
    return boundary_error;
}

int BoundaryClass::getBoundaryMaxNum()
{
    return boundary_max_num;
}

int BoundaryClass::getBoundaryMinNum()
{
    return boundary_min_num;
}

int BoundaryClass::getBoundaryMaxStart()
{
    return boundary_max_start;
}

int BoundaryClass::getBoundaryMinStart()
{
    return boundary_min_start;
}

int BoundaryClass::getBoundaryTotalNum()
{
    return (boundary_max_num + boundary_min_num);
}

void BoundaryClass::mixBoundaryExtr(int *mix_extr_x, double *mix_extr_y)
{
    if ((mix_extr_x == nullptr) || (mix_extr_y == nullptr))
    {
        boundary_error = -1;
        return;
    }
    for (int i = 0; i < boundary_min_num; i++)
    {
        mix_extr_x[i] = boundary_x_min[i];
        mix_extr_y[i] = boundary_y_min[i];
    }
    for (int i = 0; i < boundary_max_num; i++)
    {
        mix_extr_x[i + boundary_min_num] = boundary_x_max[i];
        mix_extr_y[i + boundary_min_num] = boundary_y_max[i];
    }

    for(int i = 0; i < boundary_min_num + boundary_max_num - 1; i++)
    {
        for(int j = 0; j < boundary_min_num + boundary_max_num - 1; j++)
        {
            int temp;
            double temp_double;
            if(mix_extr_x[j] > mix_extr_x[j + 1])
            {
                temp = mix_extr_x[j];
                temp_double = mix_extr_y[j];
                mix_extr_x[j] = mix_extr_x[j + 1];
                mix_extr_y[j] = mix_extr_y[j + 1];
                mix_extr_x[j+1] = temp;
                mix_extr_y[j+1] = temp_double;
            }
        }
    }

}

void BoundaryClass::showAllBoundaryData()
{
    for (int i = 0; i < boundary_min_num; i++)
    {
        qDebug() << "boundary_x_min[" << i << "]=" << boundary_x_min[i];
    }

    for (int i = 0; i < boundary_max_num; i++)
    {
        qDebug() << "boundary_x_max[" << i << "]=" << boundary_x_max[i];
    }

    for (int i = 0; i < boundary_min_num; i++)
    {
        qDebug() << "boundary_y_min[" << i << "]=" << boundary_y_min[i];
    }

    for (int i = 0; i < boundary_max_num; i++)
    {
        qDebug() << "boundary_y_max[" << i << "]=" << boundary_y_max[i];
    }
}

void BoundaryClass::computeBoundaryMainlyMirrorSymmetryCondition(double *buf_in, const int buf_in_len,
                                         const int *max_pos, const int max_num,
                                         const int *min_pos, const int min_num)
{

    int * t = new int [buf_in_len];
    int i = 0;
    int lmax[2];
    int lmin[2];
    int lsym = 0;
    int rmax[2];
    int rmin[2];
    int rsym = 0;
    int tlmin[2];
    int tlmax[2];
    int trmin[2];
    int trmax[2];
    double zlmax[2];
    double zrmin[2];
    double zrmax[2];
    double zlmin[2];
    int lminCot = 0;
    int lmaxCot = 0;
    int rminCot = 0;
    int rmaxCot = 0;
    int maxTemp = 0;
    int mintemp = 0;

    for(i = 0; i < buf_in_len; i++)
    {
        t[i] = i;
    }

    //求信号起点边界
    if(max_pos[0] < min_pos[0])
    {
        if(buf_in[0] > buf_in[min_pos[0]])
        {
            if(max_num > 2)
            {
                lmax[0] = max_pos[2];
                lmax[1] = max_pos[1];
                lmaxCot = 2;
            }
            else
            {
                lmax[0] = max_pos[1];
                lmaxCot = 1;
            }
            if(min_num >= 2)
            {
                lmin[0] = min_pos[1];
                lmin[1] = min_pos[0];
                lminCot = 2;
            }
            else
            {
                lmin[0] = min_pos[0];
                lminCot = 1;
            }
            lsym = max_pos[0];
        }
        else
        {
            if(max_num >= 2)
            {
                lmax[0] = max_pos[1];
                lmax[1] = max_pos[0];
                lmaxCot = 2;
            }
            else
            {
                lmax[0] = max_pos[0];
                lmaxCot = 1;
            }
            lmin[0] = min_pos[0];
            lmin[1] = 0;//1 before modify
            lminCot = 2;
            lsym = 0;
        }
    }
    else
    {
        if(buf_in[0] < buf_in[max_pos[0]])
        {
            if(max_num >= 2)
            {
                lmax[0] = max_pos[1];
                lmax[1] = max_pos[0];
                lmaxCot = 2;
            }
            else
            {
                lmax[0] = max_pos[0];
                lmaxCot = 1;
            }
            if(min_num > 2)
            {
                lmin[0] = min_pos[2];
                lmin[1] = min_pos[1];
                lminCot = 2;
            }
            else
            {
                lmin[0] = min_pos[1];
                lminCot = 1;
            }
            lsym = min_pos[0];
        }
        else
        {
            lmax[0] = max_pos[0];
            lmax[1] = 0;//1 before modify
            lmaxCot = 2;
            if(min_num >= 2)
            {
                lmin[0] = min_pos[1];
                lmin[1] = min_pos[0];
                lminCot = 2;
            }
            else
            {
                lmin[0] = min_pos[0];
                lminCot = 1;
            }
            lsym = 0;
        }
    }
    //求信号终点边界
    if(max_pos[max_num - 1] < min_pos[min_num - 1])
    {
        if(buf_in[buf_in_len - 1] < buf_in[max_pos[max_num - 1]])
        {
            if(max_num >= 2)
            {
                rmax[0] = max_pos[max_num - 1];
                rmax[1] = max_pos[max_num - 2];
                rmaxCot = 2;
            }
            else
            {
                rmax[0] = max_pos[max_num - 1];
                rmaxCot = 1;
            }
            if(min_num > 2)
            {
                rmin[0] = min_pos[min_num - 2];
                rmin[1] = min_pos[min_num - 3];
                rminCot = 2;
            }
            else
            {
                rmin[0] = min_pos[min_num - 2];
                rminCot = 1;
            }
            rsym = min_pos[min_num - 1];
        }
        else
        {
            rmax[0] = buf_in_len - 1;
            rmax[1] = max_pos[max_num - 1];
            rmaxCot = 2;
            if(min_num >= 2)
            {
                rmin[0] = min_pos[min_num - 1];
                rmin[1] = min_pos[min_num - 2];
                rminCot = 2;
            }
            else
            {
                rmin[0] = min_pos[min_num - 1];
                rminCot = 1;
            }
            rsym = buf_in_len - 1;
        }
    }
    else
    {
        if(buf_in[buf_in_len - 1] > buf_in[min_pos[min_num - 1]])
        {
            if(max_num > 2)
            {
                rmax[0] = max_pos[max_num - 2];
                rmax[1] = max_pos[max_num - 3];
                rmaxCot = 2;
            }
            else
            {
                rmax[0] = max_pos[max_num - 2];
                rmaxCot = 1;
            }
            if(min_num >= 2)
            {
                rmin[0] = min_pos[min_num - 1];
                rmin[1] = min_pos[min_num - 2];
                rminCot = 2;
            }
            else
            {
                rmin[0] = min_pos[min_num - 1];
                rminCot = 1;
            }
            rsym = max_pos[max_num - 1];
        }
        else
        {
            if(max_num >= 2)
            {
                rmax[0] = max_pos[max_num - 1];
                rmax[1] = max_pos[max_num - 2];
                rmaxCot = 2;
            }
            else
            {
                rmax[0] = max_pos[max_num - 1];
                rmaxCot = 1;
            }
            rmin[0] = buf_in_len - 1;
            rmin[1] = min_pos[min_num - 1];
            rminCot = 2;
            rsym = buf_in_len - 1;
        }
    }

    tlmin[0] = 2 * t[lsym] - t[lmin[0]];
    if(lminCot > 1)
    {
        tlmin[1] = 2 * t[lsym] - t[lmin[1]];
    }
    tlmax[0] = 2 * t[lsym] - t[lmax[0]];
    if(lmaxCot > 1)
    {
        tlmax[1] = 2 * t[lsym] - t[lmax[1]];
    }
    trmin[0] = 2 * t[rsym] - t[rmin[0]];
    if(rminCot > 1)
    {
        trmin[1] = 2 * t[rsym] - t[rmin[1]];
    }
    trmax[0] = 2 * t[rsym] - t[rmax[0]];
    if(rmaxCot > 1)
    {
        trmax[1] = 2 * t[rsym] - t[rmax[1]];
    }

    if((tlmin[0] > t[0]) || (tlmax[0] > t[0]))
    {
        if(lsym == max_pos[0])
        {
            if(max_num >= 2)
            {
                lmax[0] = max_pos[1];
                lmax[1] = max_pos[0];
                lmaxCot = 2;
            }
            else
            {
                lmax[0] = max_pos[0];
                lmaxCot = 1;
            }
        }
        else
        {
            if(min_num >= 2)
            {
                lmin[0] = min_pos[1];
                lmin[1] = min_pos[0];
                lminCot = 2;
            }
            else
            {
                lmin[0] = min_pos[0];
                lminCot = 1;
            }
        }
        if(lsym == 1)
        {
            boundary_error = -5;
            return;
        }
        lsym = 1;
        tlmin[0] = 2 * t[lsym] - t[lmin[0]];
        if(lminCot > 1)
        {
            tlmin[1] = 2 * t[lsym] - t[lmin[1]];
        }
        tlmax[0]= 2 * t[lsym] - t[lmax[0]];
        if(lmaxCot > 1)
        {
            tlmax[1] = 2 * t[lsym] - t[lmax[1]];
        }
    }
    maxTemp = rmaxCot;
    mintemp = rminCot;
    if((trmin[mintemp - 1] < t[buf_in_len - 1]) || (trmax[maxTemp - 1] < t[buf_in_len - 1]))
    {
        if(rsym == max_pos[max_num - 1])
        {
            if(max_num >= 2)
            {
                rmax[0] = max_pos[max_num - 1];
                rmax[1] = max_pos[max_num - 2];
                rmaxCot = 2;
            }
            else
            {
                rmax[0] = max_pos[max_num - 1];
                rmaxCot = 1;
            }
        }
        else
        {
            if(min_num >= 2)
            {
                rmin[0] = min_pos[min_num - 1];
                rmin[1] = min_pos[min_num - 2];
                rminCot = 2;
            }
            else
            {
                rmin[0] = min_pos[min_num - 1];
                rminCot = 1;
            }
        }
        if(rsym == buf_in_len - 1)
        {
            boundary_error = -6;
            return;
        }
        rsym = buf_in_len - 1;
        trmin[0] = 2 * t[rsym] - t[rmin[0]];
        if(rminCot > 1)
        {
            trmin[1] = 2 * t[rsym] - t[rmin[1]];
        }
        trmax[0] = 2 * t[rsym] - t[rmax[0]];
        if(rmaxCot > 1)
        {
            trmax[1] = 2 * t[rsym] - t[rmax[1]];
        }
    }

    //获取相应的y值
    zlmax[0] = buf_in[lmax[0]];
    if(lmaxCot > 1)
    {
        zlmax[1] = buf_in[lmax[1]];
    }
    zlmin[0] = buf_in[lmin[0]];
    if(lminCot > 1)
    {
        zlmin[1] = buf_in[lmin[1]];
    }
    zrmax[0] = buf_in[rmax[0]];
    if(rmaxCot > 1)
    {
        zrmax[1] = buf_in[rmax[1]];
    }
    zrmin[0] = buf_in[rmin[0]];
    if(rminCot > 1)
    {
        zrmin[1] = buf_in[rmin[1]];
    }

    //将结果写入
    delete [] boundary_x_min;
    delete [] boundary_x_max;
    delete [] boundary_y_min;
    delete [] boundary_y_max;

    boundary_min_num = lminCot + min_num + rminCot;

    boundary_x_min = new int [boundary_min_num];
    boundary_y_min = new double [boundary_min_num];
    boundary_x_min[0] = tlmin[0];
    boundary_y_min[0] = zlmin[0];
    if (lminCot > 1)
    {
        boundary_x_min[1] = tlmin[1];
        boundary_y_min[1] = zlmin[1];
    }
    for (i = 0; i < min_num; i++)
    {
        boundary_x_min[i + lminCot] = min_pos[i];
        boundary_y_min[i + lminCot] = buf_in[min_pos[i]];
    }
    boundary_x_min[lminCot + min_num] = trmin[0];
    boundary_y_min[lminCot + min_num] = zrmin[0];
    if (rminCot > 1)
    {
        boundary_x_min[lminCot + min_num + 1] = trmin[1];
        boundary_y_min[lminCot + min_num + 1] = zrmin[1];
    }

    boundary_max_num = lmaxCot + max_num + rmaxCot;

    boundary_x_max = new int [boundary_max_num];
    boundary_y_max = new double [boundary_max_num];
    boundary_x_max[0] = tlmax[0];
    boundary_y_max[0] = zlmax[0];
    if (lmaxCot > 1)
    {
        boundary_x_max[1] = tlmax[1];
        boundary_y_max[1] = zlmax[1];
    }
    for (i = 0; i < max_num; i++)
    {
        boundary_x_max[i + lmaxCot] = max_pos[i];
        boundary_y_max[i + lmaxCot] = buf_in[max_pos[i]];
    }
    boundary_x_max[lmaxCot + max_num] = trmax[0];
    boundary_y_max[lmaxCot + max_num] = zrmax[0];
    if (rmaxCot > 1)
    {
        boundary_x_max[lmaxCot + max_num + 1] = trmax[1];
        boundary_y_max[lmaxCot + max_num + 1] = zrmax[1];
    }
//    boundary_max_start = boundary_x_max[lmaxCot];
//    boundary_min_start = boundary_x_min[lminCot];
    delete [] t;
}

void BoundaryClass::computeBoundaryCharacteristicWaveCondition(double *buf_in, const int buf_in_len,
                                                               const int *max_pos, const int max_num,
                                                               const int *min_pos, const int min_num)
{
    if ((buf_in == nullptr) || (max_pos == nullptr) || (min_pos == nullptr))
    {
         boundary_error = -1;
         return;
    }
    if ((buf_in_len <= 0) || (max_num < 0) || (min_num < 0))
    {
        boundary_error = -2;
        return;
    }
    int i;
    int lminCot = 1;
    int rminCot = 1;
    int lmaxCot = 1;
    int rmaxCot = 1;
    int tlmin[2];
    int tlmax[2];
    int trmin[2];
    int trmax[2];

    double zlmin[2];
    double zlmax[2];
    double zrmin[2];
    double zrmax[2];

    int Tlmax = 0;
    int Trmax = 0;
    int Tlmin = 0;
    int Trmin = 0;

    if (max_num >= 4)
    {
        Tlmax = (max_pos[3] - max_pos[0]) / 3;
        Trmax = (max_pos[max_num - 1] - max_pos[max_num - 4]) / 3;
        zlmax[0] = (buf_in[max_pos[2]] + buf_in[max_pos[1]] + buf_in[max_pos[0]]) / 3;
        zrmax[0] = (buf_in[max_pos[max_num - 3]]
                + buf_in[max_pos[max_num - 2]]
                + buf_in[max_pos[max_num - 1]]) / 3;
    }
    else if(max_num == 3)
    {
        Tlmax = (max_pos[2] - max_pos[0]) / 2;
        Trmax = (max_pos[max_num - 1] - max_pos[max_num - 3]) / 2;
        zlmax[0] = (buf_in[max_pos[2]] + buf_in[max_pos[1]] + buf_in[max_pos[0]]) / 3;
        zrmax[0] = (buf_in[max_pos[max_num - 3]]
                + buf_in[max_pos[max_num - 2]]
                + buf_in[max_pos[max_num - 1]]) / 3;
    }
    else if(max_num == 2)
    {
        Tlmax = (max_pos[1] - max_pos[0]);
        Trmax = (max_pos[max_num - 1] - max_pos[max_num - 2]);
        zlmax[0] = (buf_in[max_pos[1]] + buf_in[max_pos[0]]) / 2;
        zrmax[0] = ( buf_in[max_pos[max_num - 2]]
                + buf_in[max_pos[max_num - 1]]) / 2;
    }
    else
    {
        boundary_error = -3;
        return;
    }

    if (min_num >= 4)
    {
        Tlmin = (min_pos[3] - min_pos[0]) / 3;
        Trmin = (min_pos[min_num - 1] - min_pos[min_num - 4]) / 3;
        zlmin[0] = (buf_in[min_pos[2]] + buf_in[min_pos[1]] + buf_in[min_pos[0]]) / 3;
        zrmin[0] = (buf_in[min_pos[min_num - 3]]
                + buf_in[min_pos[min_num - 2]]
                + buf_in[min_pos[min_num - 1]]) / 3;
    }
    else if(min_num == 3)
    {
        Tlmin = (min_pos[2] - min_pos[0]) / 2;
        Trmin = (min_pos[min_num - 1] - min_pos[min_num - 3]) / 2;
        zlmin[0] = (buf_in[min_pos[2]] + buf_in[min_pos[1]] + buf_in[min_pos[0]]) / 3;
        zrmin[0] = (buf_in[min_pos[min_num - 3]]
                + buf_in[min_pos[min_num - 2]]
                + buf_in[min_pos[min_num - 1]]) / 3;
    }
    else if(min_num == 2)
    {
        Tlmin = (min_pos[1] - min_pos[0]) ;
        Trmin = (min_pos[min_num - 1] - min_pos[min_num - 2]) ;
        zlmin[0] = (buf_in[min_pos[1]] + buf_in[min_pos[0]]) / 2;
        zrmin[0] = (buf_in[min_pos[min_num - 2]]
                + buf_in[min_pos[min_num - 1]]) / 2;
    }
    else
    {
        boundary_error = -4;
        return;
    }

    tlmax[0] = max_pos[0] - Tlmax;
    tlmax[1] = 1;
    trmax[0] = max_pos[max_num - 1] + Trmax;
    trmax[1] = max_num - 2;

    tlmin[0] = min_pos[0] - Tlmin;
    tlmin[1] = 1;
    trmin[0] = min_pos[min_num - 1] + Trmin;
    trmin[1] = min_num - 2;

//    zlmin[0] = (buf_in[min_pos[2]] + buf_in[min_pos[1]] + buf_in[min_pos[0]]) / 3;
    zlmin[1] = 1;
//    zlmax[0] = (buf_in[max_pos[2]] + buf_in[max_pos[1]] + buf_in[max_pos[0]]) / 3;
    zlmax[1] = 1;

//    zrmin[0] = (buf_in[min_pos[min_num - 3]]
//            + buf_in[min_pos[min_num - 2]]
//            + buf_in[min_pos[min_num - 1]]) / 3;
    zrmin[1] = 1;

//    zrmax[0] = (buf_in[max_pos[max_num - 3]]
//            + buf_in[max_pos[max_num - 2]]
//            + buf_in[max_pos[max_num - 1]]) / 3;
    zrmax[1] = 1;
    if (DEBUG_FOR_BOUNDARY2)
    {
        qDebug() << "after" << endl
                 << tlmin[0] << endl
                 << trmin[0] << endl
                 << tlmax[0] << endl
                 << trmax[0] << endl;
    }
    if (tlmin[0] < 0)
    {
        boundary_min_start = -tlmin[0];
        lminCot = 1;
    }
    else
    {
        boundary_min_start = 0;
        tlmin[0] = 0;
        zlmin[0] = buf_in[0];
        lminCot = 0;
    }

    if ((trmin[0] + boundary_min_start) > (buf_in_len - 1))
    {
//        rminCot = trmin[0] - buf_in_len + 1;
        rminCot = 1;
    }
    else
    {
        trmin[0] = buf_in_len - 1;
        zrmin[0] = buf_in[buf_in_len - 1];
        rminCot = 0;
    }

    if (tlmax[0] < 0)
    {
        boundary_max_start = -tlmax[0];
        lmaxCot = 1;
    }
    else
    {
        boundary_max_start = 0;
        tlmax[0] = 0;
        zlmax[0] = buf_in[0];
        lmaxCot = 0;
    }

    if ((trmax[0] + boundary_max_start) > (buf_in_len - 1))
    {
//        rmaxCot = trmax[0] - buf_in_len + 1;
        rmaxCot = 1;
    }
    else
    {
        trmax[0] = buf_in_len - 1;
        zrmax[0] = buf_in[buf_in_len - 1];
        rmaxCot = 0;
    }

    if (DEBUG_FOR_BOUNDARY2)
    {
        qDebug() <<"about tl tr max/min" << endl
                << tlmin[0] << endl
                << trmin[0] << endl
                << tlmax[0] << endl
                << trmax[0] << endl;
        qDebug() <<"about l r max/min Cot" << endl
                << lminCot << endl
                << rminCot << endl
                << lmaxCot << endl
                << rmaxCot << endl;
    }

    //将结果写入
    delete [] boundary_x_min;
    delete [] boundary_x_max;
    delete [] boundary_y_min;
    delete [] boundary_y_max;

    boundary_min_num = lminCot + min_num + rminCot;

    boundary_x_min = new int [boundary_min_num];
    boundary_y_min = new double [boundary_min_num];
    if (lminCot == 1)
    {
        boundary_x_min[0] = tlmin[0] + boundary_min_start;
        boundary_y_min[0] = zlmin[0];
    }
    if (0)//lminCot > 1)
    {
        boundary_x_min[1] = tlmin[1];
        boundary_y_min[1] = zlmin[1];
    }
    for (i = 0; i < min_num; i++)
    {
        boundary_x_min[i + lminCot] = min_pos[i] + boundary_min_start;
        boundary_y_min[i + lminCot] = buf_in[min_pos[i]];
    }
    if (rminCot == 1)
    {
        boundary_x_min[lminCot + min_num] = trmin[0] + boundary_min_start;
        boundary_y_min[lminCot + min_num] = zrmin[0];
    }

    if (0)//rminCot > 1)
    {
        boundary_x_min[lminCot + min_num + 1] = trmin[1];
        boundary_y_min[lminCot + min_num + 1] = zrmin[1];
    }

    boundary_max_num = lmaxCot + max_num + rmaxCot;

    boundary_x_max = new int [boundary_max_num];
    boundary_y_max = new double [boundary_max_num];
    if (lmaxCot == 1)
    {
        boundary_x_max[0] = tlmax[0] + boundary_max_start;
        boundary_y_max[0] = zlmax[0];
    }
    if (0)//lmaxCot > 1)
    {
        boundary_x_max[1] = tlmax[1];
        boundary_y_max[1] = zlmax[1];
    }
    for (i = 0; i < max_num; i++)
    {
        boundary_x_max[i + lmaxCot] = max_pos[i] + boundary_max_start;
        boundary_y_max[i + lmaxCot] = buf_in[max_pos[i]];
    }
    if (rmaxCot == 1)
    {
        boundary_x_max[lmaxCot + max_num] = trmax[0] + boundary_max_start;
        boundary_y_max[lmaxCot + max_num] = zrmax[0];
    }
    if (0)//rmaxCot > 1)
    {
        boundary_x_max[lmaxCot + max_num + 1] = trmax[1];
        boundary_y_max[lmaxCot + max_num + 1] = zrmax[1];
    }
    if (DEBUG_FOR_BOUNDARY2)
    {
        qDebug() <<"about max/min" << endl
                << min_num << endl
                << max_num << endl
                << boundary_min_num << endl
                << boundary_max_num << endl;
    }
}
