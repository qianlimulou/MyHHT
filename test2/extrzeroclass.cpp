#include "extrzeroclass.h"

ExtrZeroClass::ExtrZeroClass()
{
    extr_min_x = new int [1];
    extr_max_x = new int [1];
    zero_x = new int [1];
    extr_min_x_num = 0;
    extr_max_x_num = 0;
    zero_x_num = 0;
    extrzero_error = 0;
}

ExtrZeroClass::~ExtrZeroClass()
{
    delete [] extr_min_x;
    extr_min_x = nullptr;
    delete [] extr_max_x;
    extr_max_x = nullptr;
    delete [] zero_x;
    zero_x = nullptr;
}

int ExtrZeroClass::getExtrZeroError()
{
    return extrzero_error;
}


int ExtrZeroClass::getExtrMaxNum()
{
    return extr_max_x_num;
}

int ExtrZeroClass::getExtrMinNum()
{
    return extr_min_x_num;
}

int ExtrZeroClass::getExtrTotalNum()
{
    return (getExtrMaxNum() + getExtrMinNum());
}

int ExtrZeroClass::getZeroNum()
{
    return zero_x_num;
}

bool ExtrZeroClass::isBufMonotonic(const double *buf_x, const int buf_in_len)
{
    if (buf_x == nullptr)
    {
        extrzero_error = -1;
        return false;
    }
    if (buf_in_len <= 0)
    {
        extrzero_error = -2;
        return false;
    }

    for (int i = 0; i < buf_in_len -1; i++)
    {
        if (buf_x[i] > buf_x[i + 1])
        {
            return false;
        }
    }
    return true;
}

void ExtrZeroClass::showAllExtrData()
{
    if (extr_min_x_num <= 0)
    {
        qDebug() << "extr_min_x_num=:" << extr_min_x_num << " <= 0 ";
    }
    for (int i = 0; i < extr_min_x_num; i++)
    {
        qDebug() << "extr_min_x[" << i << "]=" << extr_min_x[i];
    }
    if (extr_max_x_num <= 0)
    {
        qDebug() << "extr_max_x_num=:" << extr_max_x_num << " <= 0 ";
    }
    for (int i = 0; i < extr_max_x_num; i++)
    {
        qDebug() << "extr_max_x[" << i << "]=" << extr_max_x[i];
    }
}

void ExtrZeroClass::showAllZeroData()
{
    if (zero_x_num <= 0)
    {
        qDebug() << "zero_x_num=:" << zero_x_num << " <= 0 ";
    }
    for (int i = 0; i < zero_x_num; i++)
    {
        qDebug() << "zero_x[" << i << "]=" << zero_x[i];
    }
}

bool ExtrZeroClass::isExtrNumOK()
{
    if (extr_min_x_num + extr_max_x_num < EXTR_TOTAL_MAX_NUM)
    {
        return true;
    }
    else
    {
        return false;
    }
//    if ((extr_min_x_num == 0) || (extr_max_x_num == 0))
//    {
//        return false;
//    }

}

bool ExtrZeroClass::isExtrNumOtherZero()
{
    if ((extr_min_x_num == 0) || (extr_max_x_num == 0))
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool ExtrZeroClass::isExtrZeroIMF()
{
    int ret = extr_min_x_num + extr_max_x_num - zero_x_num;
    if ((ret < -1) || (ret > 1))
    {
        return false;
    }
    else
    {
        return true;
    }
}



void ExtrZeroClass::computeExtrResult(const double *buf_y, const int buf_in_len)
{
    if (buf_y == nullptr)
    {
        extrzero_error = -1;
        return;
    }
    if (buf_in_len <= 0)
    {
        extrzero_error = -2;
        return;
    }
    int i, j;//for循环使用
    int indminTempCot = 0;//临时极小值个数
    int indmaxTempCot = 0;//临时极大值个数
    int mmTempCot = 0;//临时幅值相同值个数
    int maxTempCot=0;//临时相同幅值极大值个数
    int minTempCot=0;//临时相同幅值极小值个数
    int *indminTempPos = new int [buf_in_len];//临时极小值位置
    int *indmaxTempPos = new int [buf_in_len];//临时极大值位置
    int *mmTempPos = new int [buf_in_len];//临时相同幅值的位置
    int *minTempPos = new int [buf_in_len];//临时相同幅值极大值位置
    int *maxTempPos = new int [buf_in_len];//临时相同幅值极小值位置
    double *buf_in_k = new double [buf_in_len - 1];//buf_in一次导数值
    //求极值
    //求一次导数
    for(i = 0; i < buf_in_len - 1; i++)
    {
        buf_in_k[i] = buf_y[i + 1] - buf_y[i];
    }
    //求一般情况极值位置
    for(i = 0; i < buf_in_len - 2; i++)
    {
        if((buf_in_k[i] * buf_in_k[i + 1] < 0) && (buf_in_k[i] < 0))
        {
            indminTempPos[indminTempCot++] = i + 1;
        }
        if((buf_in_k[i] * buf_in_k[i + 1] < 0) && (buf_in_k[i] > 0))
        {
            indmaxTempPos[indmaxTempCot++] = i + 1;
        }
    }

    ////////此处同样要考虑幅值连续相等的情况
    for(i = 0; i < buf_in_len - 1; i++)
    {
        //if(x[i]==0)
       //     continue;
        int start_for_same_Extremum = 0;
        int end_for_same_Extremum = 0;
        if(buf_y[i] == buf_y[i + 1])
        {
            start_for_same_Extremum = i;
        }
        else
        {
            continue;
        }
        //modify by wang
        //for(j = 0; j < buf_in_len; j++)
        for(j = 0; j < buf_in_len - i - 1; j++)
        //end modify
        {
            if(buf_y[i + j] != buf_y[i + j + 1])
            {
                end_for_same_Extremum = i + j;
                mmTempPos[2 * mmTempCot] = start_for_same_Extremum;
                mmTempPos[2 * mmTempCot + 1] = end_for_same_Extremum;
                mmTempCot++;
                i = end_for_same_Extremum;
                break;
            }
        }
    }
    //求相同幅值极值位置
    for(i = 0; i < mmTempCot; i++)
    {
        int temp;
        double dTemp;
        if(mmTempPos[2 * i]==0)
        {
            continue;////////////////////////////////////////端点算不算极值点
            if(buf_y[mmTempPos[2 * i + 1] + 1] > buf_y[mmTempPos[2 * i + 1]])
            {
                minTempPos[minTempCot++] = (int)(mmTempPos[2 * i]
                        + mmTempPos[2 * i + 1]) / 2 + 1;
            }
            else
            {
                maxTempPos[maxTempCot++] = (int)(mmTempPos[2 * i]
                        + mmTempPos[2 * i + 1]) / 2 + 1;
            }
        }
        else if(mmTempPos[2 * i + 1] == buf_in_len - 1)
        {
            continue;////////////////////////////////////////端点算不算极值点
            if(buf_y[mmTempPos[2 * i]] < buf_y[mmTempPos[2 * i] - 1])
            {
                minTempPos[minTempCot++] = (int)(mmTempPos[2 * i]
                        + mmTempPos[2 * i + 1]) / 2 + 1;
            }
            else
            {
                maxTempPos[maxTempCot++] = (int)(mmTempPos[2 * i]
                        + mmTempPos[2 * i + 1]) / 2 + 1;
            }
        }
        else
        {
            if((buf_y[mmTempPos[2 * i]] > buf_y[mmTempPos[2 * i] - 1])
                    && (buf_y[mmTempPos[2 * i + 1]] > buf_y[mmTempPos[2 * i + 1] + 1]))
            {
                dTemp = (mmTempPos[2 * i] + mmTempPos[2 * i + 1]) / 2.0;
                temp = (int)dTemp;
                if(dTemp - temp > 0)
                {
                    temp = temp + 1;
                }
                else
                {
                    temp = (int)dTemp;
                }
                maxTempPos[maxTempCot++] = temp;
            }
            if((buf_y[mmTempPos[2 * i]] < buf_y[mmTempPos[2 * i] - 1])
                    && (buf_y[mmTempPos[2 * i + 1]] < buf_y[mmTempPos[2 * i + 1] + 1]))
            {
                dTemp = (mmTempPos[2 * i] + mmTempPos[2 * i + 1]) / 2.0;
                temp = (int)dTemp;
                if(dTemp - temp > 0)
                {
                    temp = temp + 1;
                }
                else
                {
                    temp = (int)dTemp;
                }
                minTempPos[minTempCot++] = temp;
            }
        }
    }

    //结合一般情况和相同幅值情况极值位置和个数
    for(i = 0; i < minTempCot; i++)
    {
        indminTempPos[indminTempCot++] = minTempPos[i];
    }
    for(i = 0; i < maxTempCot; i++)
    {
        indmaxTempPos[indmaxTempCot++] = maxTempPos[i];
    }
    //根据位置排序
    for(i = 0; i < indminTempCot - 1; i++)
    {
        for(j = 0; j < indminTempCot - 1; j++)
        {
            int temp;
            if(indminTempPos[j] > indminTempPos[j + 1])
            {
                temp = indminTempPos[j];
                indminTempPos[j] = indminTempPos[j + 1];
                indminTempPos[j+1] = temp;
            }
        }
    }

    for(i = 0; i < indmaxTempCot - 1; i++)
    {
        for(j = 0; j < indmaxTempCot - 1; j++)
        {
            int temp;
            if(indmaxTempPos[j] > indmaxTempPos[j + 1])
            {
                temp = indmaxTempPos[j];
                indmaxTempPos[j] = indmaxTempPos[j + 1];
                indmaxTempPos[j + 1] = temp;
            }
        }
    }

    //写入结果
    extr_min_x_num = indminTempCot;
    extr_max_x_num = indmaxTempCot;
    delete [] extr_min_x;
    extr_min_x = nullptr;
    delete [] extr_max_x;
    extr_max_x = nullptr;
    extr_min_x = new int [indminTempCot];
    extr_max_x = new int [indmaxTempCot];
    for(i = 0; i < indminTempCot; i++)
    {
        extr_min_x[i] = indminTempPos[i];
    }

    for(i = 0; i < indmaxTempCot; i++)
    {
        extr_max_x[i] = indmaxTempPos[i];
    }

    delete [] minTempPos;
    minTempPos = nullptr;
    delete [] maxTempPos;
    maxTempPos = nullptr;
    delete [] mmTempPos;
    mmTempPos = nullptr;
    delete [] indminTempPos;
    indminTempPos = nullptr;
    delete [] indmaxTempPos;
    indmaxTempPos = nullptr;
    delete [] buf_in_k;
    buf_in_k = nullptr;

}

void ExtrZeroClass::computeZeroResult(const double * buf_y, const int buf_in_len)
{
    if (buf_y == nullptr)
    {
        extrzero_error = -1;
        return;
    }
    if (buf_in_len <= 0)
    {
        extrzero_error = -2;
        return;
    }
    int i, j;
    int conNum = 0;//值为连续零点个数
    int indZeroTempCot = 0;//临时零点个数
    int conZeroCot = 0;//值为0的个数
    int *indZeroTempPos = new int [buf_in_len];//临时零点位置
    int *conZeroPos = new int [buf_in_len];//值为0的位置

    //求一般情况零点位置
    for(i = 0; i < buf_in_len - 1; i++)
    {
        if(buf_y[i] * buf_y[i + 1] < 0)
        {
            indZeroTempPos[indZeroTempCot++] = i;
        }
    }
    ////////此处要考虑连零的情况
    ///

    for(i = 0; i < buf_in_len - 1; i++)
    {
        if(buf_y[i] == 0)
        {
//            int conNum = 1;
            while((i < buf_in_len - 1) && (buf_y[i + 1] == 0))
            {
                i++;
                conNum++;
            }
            conZeroPos[conZeroCot++] = i - (int)(conNum/2);
        }
    }
    //modify by wang but is wrong
//    for(i = 0; i < buf_in_len; i++)
//    {
//        if(buf_y[i] == 0)
//        {
//            conZeroPos[conZeroCot++] = i;
//        }
//    }
    //end modify
    //结合一般情况和连续零情况位置和个数
    for(i = 0; i < conZeroCot; i++)
    {
        indZeroTempPos[indZeroTempCot ++] = conZeroPos[i];
    }

    //根据位置排序
    for(i = 0; i < indZeroTempCot - 1; i++)
    {
        for(j = 0; j < indZeroTempCot - 1; j++)
        {
            int temp;
            if(indZeroTempPos[j] > indZeroTempPos[j+1])
            {
                temp = indZeroTempPos[j];
                indZeroTempPos[j] = indZeroTempPos[j + 1];
                indZeroTempPos[j + 1] = temp;
            }
        }
    }

    zero_x_num = indZeroTempCot;
    delete [] zero_x;
    zero_x = nullptr;
    zero_x = new int [indZeroTempCot];
    for(i = 0; i < indZeroTempCot; i++)
    {
        zero_x[i] = indZeroTempPos[i];
    }

    delete [] indZeroTempPos;
    indZeroTempPos = nullptr;
    delete [] conZeroPos;
    conZeroPos = nullptr;

}
