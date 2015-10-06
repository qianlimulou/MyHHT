#include "myemd.h"

void showBoundaryValue(boundary * boundary)
{
    for(int i = 0; i < boundary->lentmax; i++)
    {
        qDebug() << "tmax[" << i << "]=" << boundary->tmax[i];
    }

    for(int i = 0; i < boundary->lentmin; i++)
    {
        qDebug() << "tmin[" << i << "]=" << boundary->tmin[i];
    }

    for(int i = 0; i < boundary->lenzmax; i++)
    {
        qDebug() << "zmax[" << i << "]=" << boundary->zmax[i];
    }

    for(int i = 0; i < boundary->lenzmin; i++)
    {
        qDebug() << "zmin[" << i << "]=" << boundary->zmin[i];
    }
}

void computeExtremumAndZeroPoints(threePoints * threePoints, const double *buf_in, const int buf_in_len)
{
    int i, j;//for循环使用
    int indzerTempCot = 0;//临时零点个数
    int conZerCot = 0;//值为0的个数
    int indminTempCot = 0;//临时极小值个数
    int indmaxTempCot = 0;//临时极大值个数
    int mmTempCot = 0;//临时幅值相同值个数
    int maxTempCot=0;//临时相同幅值极大值个数
    int minTempCot=0;//临时相同幅值极小值个数
    int *indzerTempPos = new int [buf_in_len];//临时零点位置
    int *conZerpos = new int [buf_in_len];//值为0的位置
    int *indminTempPos = new int [buf_in_len];//临时极小值位置
    int *indmaxTempPos = new int [buf_in_len];//临时极大值位置
    int *mmTempPos = new int [buf_in_len];//临时相同幅值的位置
    int *minTempPos = new int [buf_in_len];//临时相同幅值极大值位置
    int *maxTempPos = new int [buf_in_len];//临时相同幅值极小值位置
    double *buf_in_k = new double [buf_in_len - 1];//buf_in一次导数值

    //求一般情况零点位置
    for(i = 0; i < buf_in_len - 1; i++)
    {
        if(buf_in[i] * buf_in[i + 1] < 0)
        {
            indzerTempPos[indzerTempCot++] = i;
        }
    }
    ////////此处要考虑连零的情况
    ///
    //modify by wang
//    for(i=0;i<length-1;i++)
//    {
//        if(x[i]==0)
//        {
//            //int conNum=1;
//            while(i<length-1&&x[i+1]==0)
//            {
//                i++;
//                conNum++;
//            }
//            conZerpos[conZerCot++]=i-(int)(conNum/2);
//        }
//    }

    for(i = 0; i < buf_in_len; i++)
    {
        if(buf_in[i] == 0)
        {
            conZerpos[conZerCot++] = i;
        }
    }
    //end modify
    //结合一般情况和连续零情况位置和个数
    for(i = 0; i < conZerCot; i++)
    {
        indzerTempPos[indzerTempCot ++] = conZerpos[i];
    }

    //根据位置排序
    for(i = 0; i < indzerTempCot - 1; i++)
    {
        for(j = 0; j < indzerTempCot - 1; j++)
        {
            int temp;
            if(indzerTempPos[j] > indzerTempPos[j+1])
            {
                temp = indzerTempPos[j];
                indzerTempPos[j] = indzerTempPos[j + 1];
                indzerTempPos[j + 1] = temp;
            }
        }
    }

    //求极值
    //求一次导数
    for(i = 0; i < buf_in_len - 1; i++)
    {
        buf_in_k[i] = buf_in[i + 1] - buf_in[i];
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
        if(buf_in[i] == buf_in[i + 1])
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
            if(buf_in[i + j] != buf_in[i + j + 1])
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
            if(buf_in[mmTempPos[2 * i + 1] + 1] > buf_in[mmTempPos[2 * i + 1]])
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
            if(buf_in[mmTempPos[2 * i]] < buf_in[mmTempPos[2 * i] - 1])
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
            if((buf_in[mmTempPos[2 * i]] > buf_in[mmTempPos[2 * i] - 1])
                    && (buf_in[mmTempPos[2 * i + 1]] > buf_in[mmTempPos[2 * i + 1] + 1]))
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
            if((buf_in[mmTempPos[2 * i]] < buf_in[mmTempPos[2 * i] - 1])
                    && (buf_in[mmTempPos[2 * i + 1]] < buf_in[mmTempPos[2 * i + 1] + 1]))
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

    //将结果写入结构体
    threePoints->count_time++;
    threePoints->lenPointMax = indmaxTempCot;
    threePoints->lenPointMin = indminTempCot;
    threePoints->lenPointZero = indzerTempCot;
    for(i = 0; i < indmaxTempCot; i++)
    {
        threePoints->posPointMax[i] = indmaxTempPos[i];
    }
    for(i = 0; i < indminTempCot; i++)
    {
        threePoints->posPointMin[i] = indminTempPos[i];
    }
    for(i = 0; i < indzerTempCot; i++)
    {
        threePoints->posPointZero[i] = indzerTempPos[i];
    }
    delete [] minTempPos;
    delete [] maxTempPos;
    delete [] mmTempPos;
    delete [] indzerTempPos;
    delete [] conZerpos;
    delete [] indminTempPos;
    delete [] indmaxTempPos;
    delete [] buf_in_k;
}

bool judgeExtremumAndZeroPointsIsSatisfyIMF(threePoints * threePoints)
{
    if (threePoints->count_time == 0)
    {
        qDebug() << "No exist any Extremum and Zero Points!!!";
        return false;
    }
    int Extremum_num = threePoints->lenPointMax + threePoints->lenPointMin;
    int Zero_num = threePoints->lenPointZero;
    if (Extremum_num > Zero_num)
    {
        if ((Extremum_num - Zero_num) <= 1)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        if ((Zero_num - Extremum_num) <= 1)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
}

bool judgeIsMonotonic()//const double *buf_in, const int buf_in_len)
{
    return true;
}

void setBoundaryCondition(threePoints * threePoints, boundary * boundary, const double *buf_in, const int buf_in_len)
{
    int * t = new int [6000];
    int i=0;
    int lmax[2];
    int lmin[2];
    int lsym=0;
    int rmax[2];
    int rmin[2];
    int rsym=0;
    int tlmin[2];
    int tlmax[2];
    int trmin[2];
    int trmax[2];
    double zlmax[2];
    double zrmin[2];
    double zrmax[2];
    double zlmin[2];
    int lminCot=0;
    int lmaxCot=0;
    int rminCot=0;
    int rmaxCot=0;
    int maxTemp=0;
    int mintemp=0;

    for(i = 0; i < buf_in_len; i++)
    {
        t[i] = i;
    }
    if (threePoints->count_time == 0)
    {
        qDebug() << "No exist any Extremum and Zero Points!!!";
        return;
    }

    if((threePoints->lenPointMax + threePoints->lenPointMin) < 3)
    {
        qDebug() << "All Extremum == 0!!!";
        threePoints->error = 1;
        boundary->error = 1;
        return;
    }

    if(((threePoints->lenPointMax == 0) || (threePoints->lenPointMin)) ==0)
    {
        qDebug() << "Extremum of Max or Min is == 0!!!";
        threePoints->error = 1;
        boundary->error = 1;
        return;
    }
    //求信号起点边界
    if(threePoints->posPointMax[0] < threePoints->posPointMin[0])
    {
        if(buf_in[0] > buf_in[threePoints->posPointMin[0]])
        {
            if(threePoints->lenPointMax > 2)
            {
                lmax[0] = threePoints->posPointMax[2];
                lmax[1] = threePoints->posPointMax[1];
                lmaxCot = 2;
            }
            else
            {
                lmax[0] = threePoints->posPointMax[1];
                lmaxCot = 1;
            }
            if(threePoints->lenPointMin >= 2)
            {
                lmin[0] = threePoints->posPointMin[1];
                lmin[1] = threePoints->posPointMin[0];
                lminCot = 2;
            }
            else
            {
                lmin[0] = threePoints->posPointMin[0];
                lminCot = 1;
            }
            lsym = threePoints->posPointMax[0];
        }
        else
        {
            if(threePoints->lenPointMax >= 2)
            {
                lmax[0] = threePoints->posPointMax[1];
                lmax[1] = threePoints->posPointMax[0];
                lmaxCot = 2;
            }
            else
            {
                lmax[0] = threePoints->posPointMax[0];
                lmaxCot = 1;
            }
            lmin[0] = threePoints->posPointMin[0];
            lmin[1] = 1;
            lminCot = 2;
            lsym = 0;
        }
    }
    else
    {
        if(buf_in[0] < buf_in[threePoints->posPointMax[0]])
        {
            if(threePoints->lenPointMax >= 2)
            {
                lmax[0] = threePoints->posPointMax[1];
                lmax[1] = threePoints->posPointMax[0];
                lmaxCot = 2;
            }
            else
            {
                lmax[0] = threePoints->posPointMax[0];
                lmaxCot = 1;
            }
            if(threePoints->lenPointMin > 2)
            {
                lmin[0] = threePoints->posPointMin[2];
                lmin[1] = threePoints->posPointMin[1];
                lminCot = 2;
            }
            else
            {
                lmin[0] = threePoints->posPointMin[1];
                lminCot = 1;
            }
            lsym = threePoints->posPointMin[0];
        }
        else
        {
            lmax[0] = threePoints->posPointMax[0];
            lmax[1] = 1;
            lmaxCot = 2;
            if(threePoints->lenPointMin >= 2)
            {
                lmin[0] = threePoints->posPointMin[1];
                lmin[1] = threePoints->posPointMin[0];
                lminCot = 2;
            }
            else
            {
                lmin[0] = threePoints->posPointMin[0];
                lminCot = 1;
            }
            lsym = 0;
        }
    }
    //求信号终点边界
    if(threePoints->posPointMax[threePoints->lenPointMax - 1]
            < threePoints->posPointMin[threePoints->lenPointMin - 1])
    {
        if(buf_in[buf_in_len-1] < buf_in[threePoints->posPointMax[threePoints->lenPointMax-1]])
        {
            if(threePoints->lenPointMax >= 2)
            {
                rmax[0] = threePoints->posPointMax[threePoints->lenPointMax - 1];
                rmax[1] = threePoints->posPointMax[threePoints->lenPointMax - 2];
                rmaxCot = 2;
            }
            else
            {
                rmax[0] = threePoints->posPointMax[threePoints->lenPointMax - 1];
                rmaxCot = 1;
            }
            if(threePoints->lenPointMin > 2)
            {
                rmin[0] = threePoints->posPointMin[threePoints->lenPointMin - 2];
                rmin[1] = threePoints->posPointMin[threePoints->lenPointMin - 3];
                rminCot = 2;
            }
            else
            {
                rmin[0] = threePoints->posPointMin[threePoints->lenPointMin - 2];
                rminCot = 1;
            }
            rsym = threePoints->posPointMin[threePoints->lenPointMin - 1];
        }
        else
        {
            rmax[0] = buf_in_len - 1;
            rmax[1] = threePoints->posPointMax[threePoints->lenPointMax - 1];
            rmaxCot = 2;
            if(threePoints->lenPointMin >= 2)
            {
                rmin[0] = threePoints->posPointMin[threePoints->lenPointMin - 1];
                rmin[1] = threePoints->posPointMin[threePoints->lenPointMin - 2];
                rminCot = 2;
            }
            else
            {
                rmin[0] = threePoints->posPointMin[threePoints->lenPointMin - 1];
                rminCot = 1;
            }
            rsym = buf_in_len - 1;
        }
    }
    else
    {
        if(buf_in[buf_in_len - 1] > buf_in[threePoints->posPointMin[threePoints->lenPointMin - 1]])
        {
            if(threePoints->lenPointMax > 2)
            {
                rmax[0] = threePoints->posPointMax[threePoints->lenPointMax - 2];
                rmax[1] = threePoints->posPointMax[threePoints->lenPointMax - 3];
                rmaxCot = 2;
            }
            else
            {
                rmax[0] = threePoints->posPointMax[threePoints->lenPointMax - 2];
                rmaxCot = 1;
            }
            if(threePoints->lenPointMin >= 2)
            {
                rmin[0] = threePoints->posPointMin[threePoints->lenPointMin - 1];
                rmin[1] = threePoints->posPointMin[threePoints->lenPointMin - 2];
                rminCot = 2;
            }
            else
            {
                rmin[0] = threePoints->posPointMin[threePoints->lenPointMin - 1];
                rminCot = 1;
            }
            rsym = threePoints->posPointMax[threePoints->lenPointMax - 1];
        }
        else
        {
            if(threePoints->lenPointMax >= 2)
            {
                rmax[0] = threePoints->posPointMax[threePoints->lenPointMax - 1];
                rmax[1] = threePoints->posPointMax[threePoints->lenPointMax - 2];
                rmaxCot = 2;
            }
            else
            {
                rmax[0] = threePoints->posPointMax[threePoints->lenPointMax - 1];
                rmaxCot = 1;
            }
            rmin[0] = buf_in_len - 1;
            rmin[1] = threePoints->posPointMin[threePoints->lenPointMin - 1];
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
        if(lsym == threePoints->posPointMax[0])
        {
            if(threePoints->lenPointMax >= 2)
            {
                lmax[0] = threePoints->posPointMax[1];
                lmax[1] = threePoints->posPointMax[0];
                lmaxCot = 2;
            }
            else
            {
                lmax[0] = threePoints->posPointMax[0];
                lmaxCot = 1;
            }
        }
        else
        {
            if(threePoints->lenPointMin >= 2)
            {
                lmin[0] = threePoints->posPointMin[1];
                lmin[1] = threePoints->posPointMin[0];
                lminCot = 2;
            }
            else
            {
                lmin[0] = threePoints->posPointMin[0];
                lminCot = 1;
            }
        }
        if(lsym == 1)
        {
            qDebug() << "------------------ERROR!------------------";
            threePoints->error = 1;
            boundary->error = 1;
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
        if(rsym == threePoints->posPointMax[threePoints->lenPointMax - 1])
        {
            if(threePoints->lenPointMax >= 2)
            {
                rmax[0] = threePoints->posPointMax[threePoints->lenPointMax - 1];
                rmax[1] = threePoints->posPointMax[threePoints->lenPointMax - 2];
                rmaxCot = 2;
            }
            else
            {
                rmax[0] = threePoints->posPointMax[threePoints->lenPointMax - 1];
                rmaxCot = 1;
            }
        }
        else
        {
            if(threePoints->lenPointMin >= 2)
            {
                rmin[0] = threePoints->posPointMin[threePoints->lenPointMin - 1];
                rmin[1] = threePoints->posPointMin[threePoints->lenPointMin - 2];
                rminCot = 2;
            }
            else
            {
                rmin[0] = threePoints->posPointMin[threePoints->lenPointMin - 1];
                rminCot = 1;
            }
        }
        if(rsym == buf_in_len - 1)
        {
            threePoints->error = 1;
            boundary->error = 1;
            return;
        }
        rsym = buf_in_len-1;
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


    boundary->tmin[0] = tlmin[0];
    if(lminCot > 1)
    {
        boundary->tmin[1] = tlmin[1];
    }
    for(i = 0; i < threePoints->lenPointMin; i++)
    {
        boundary->tmin[i + lminCot] = threePoints->posPointMin[i];
    }

    boundary->tmin[threePoints->lenPointMin + lminCot] = trmin[0];
    if(rminCot > 1)
    {
        boundary->tmin[threePoints->lenPointMin + lminCot + 1] = trmin[1];
    }
    boundary->lentmin = threePoints->lenPointMin + lminCot + rminCot;

    boundary->tmax[0] = tlmax[0];
    if(lmaxCot > 1)
    {
        boundary->tmax[1] = tlmax[1];
    }
    for(i = 0; i < threePoints->lenPointMax; i++)
    {
        boundary->tmax[i + lmaxCot] = threePoints->posPointMax[i];
    }
    boundary->tmax[threePoints->lenPointMax + lmaxCot] = trmax[0];
    if(rmaxCot > 1)
    {
        boundary->tmax[threePoints->lenPointMax + lmaxCot + 1] = trmax[1];
    }
    boundary->lentmax = threePoints->lenPointMax + lmaxCot + rmaxCot;

    boundary->zmin[0] = zlmin[0];
    if(lminCot > 1)
    {
        boundary->zmin[1] = zlmin[1];
    }
    for(i = 0; i < threePoints->lenPointMin; i++)
    {
        boundary->zmin[i + lminCot] = buf_in[threePoints->posPointMin[i]];
    }
    boundary->zmin[threePoints->lenPointMin + lminCot] = zrmin[0];
    if(rminCot > 1)
    {
        boundary->zmin[threePoints->lenPointMin + lminCot + 1] = zrmin[1];
    }
    boundary->lenzmin = threePoints->lenPointMin + lminCot + rminCot;

    boundary->zmax[0] = zlmax[0];
    if(lmaxCot > 1)
    {
        boundary->zmax[1] = zlmax[1];
    }
    for(i = 0; i < threePoints->lenPointMax; i++)
    {
        boundary->zmax[i + lmaxCot] = buf_in[threePoints->posPointMax[i]];
    }
    boundary->zmax[threePoints->lenPointMax + lmaxCot] = zrmax[0];
    if(rmaxCot > 1)
    {
        boundary->zmax[threePoints->lenPointMax + lmaxCot + 1] = zrmax[1];
    }
    boundary->lenzmax = threePoints->lenPointMax + rmaxCot + lmaxCot;

    delete [] t;

}
void computeSplineResult(const double * x, const double * y, const int xy_len,
                                double * result,  const int buf_in_len, const condition_spline con,
                                const int A, const int B)
{
    int n = xy_len - 1;//插值曲线个数
    int i, j;//循环计数
    double a0, b1, c2, d3;//三次曲线系数
    double temp; //交换临时存数
    double *a = new double [n];
    double *b = new double [n];
    double *c = new double [n];
    double *d = new double [n];
    double *f = new double [n];
    double *m = new double [n];
    double *h = new double [n];

    /* 计算 h[] d[] */
    for (i = 0; i < n; i++)
    {
        h[i] = x[i + 1] - x[i];
        d[i] = (y[i + 1] - y[i]) / h[i];
        /* printf("%f\t%f\n", h[i], d[i]); */
    }

    /* 计算 a[] b[] d[] f[] */
    switch(con)
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
        f[0] = 6 * (d[0] - A);
        f[n] = 6 * (B - d[n - 1]);
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

    for(i = 0; i < buf_in_len; i++)
    {
        for(j = 0; j < n - 1; j++)
        {
            if(i >= x[j] && i < x[j + 1])
            {
                d3 = (m[j + 1] - m[j]) / (6 * h[j]);
                c2 = m[j] / 2;
                b1 = d[j] - (h[j] / 6) * (2 * m[j] + m[j + 1]);
                a0 = y[j];
                result[i] = d3 * (i - x[j]) * (i - x[j]) * (i - x[j])
                          + c2 * (i - x[j]) * (i - x[j])
                          + b1 * (i - x[j])
                          + a0;
                //printf("mySpline%d=%f\n",i,(mySpline[i]));
            }
        }
    }
    delete [] a;
    delete [] b;
    delete [] c;
    delete [] d;
    delete [] f;
    delete [] m;
    delete [] h;


}

void proceedEmdAndJudgeStopRequirement()//double *buf_in, int buf_in_len)
{

}

