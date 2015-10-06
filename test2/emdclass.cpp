#include "emdclass.h"

EmdClass::EmdClass()
{
    MyExtrZero = new ExtrZeroClass;
    MyBoundary = new BoundaryClass;
    MySpline = new SplineClass;
    emd_envmoy = new double[1];
    emd_amp = new double[1];
    emd_imf = new double[1];
    emd_residue = new double[1];
    emd_imf_len = 0;
    emd_imf_num = 0;
    emd_sd = 0.05;
    emd_sd2 = 0.5;
    emd_tol = 0.05;
    emd_imf_max_num = 10;
    emd_iteration = 0;
    emd_max_iteration = 2000;
    emd_error = 0;
    test_a = 0;
}

EmdClass::~EmdClass()
{
    delete MyExtrZero;
    delete MyBoundary;
    delete MySpline;
    delete [] emd_envmoy;
    delete [] emd_amp;
    delete [] emd_imf;
    delete [] emd_residue;
}

void EmdClass::showAllEnvmoyAndAmp()
{
    if (emd_imf_len <= 0)
    {
        qDebug() << "emd_imf_len=:" << emd_imf_len << " <= 0 ";
    }
    for(int i = 0; i < emd_imf_len; i++)
    {
        qDebug() << "Envmoy[" << i << "]=:" << emd_envmoy[i];
    }
    for(int i = 0; i < emd_imf_len; i++)
    {
        qDebug() << "Amp[" << i << "]=:" << emd_amp[i];
    }
}

void EmdClass::showAllEmdImfData()
{
    if ((emd_imf_len <= 0) || (emd_imf_num <= 0))
    {
        qDebug() << "emd_imf_len=:" << emd_imf_len << " <= 0 ";
    }
    for(int i = 0; i < emd_imf_len * emd_imf_num; i++)
    {
        qDebug() << "emd_imf[" << i << "]=:" << emd_imf[i];
    }
}

void EmdClass::showAllEmdResidueData()
{
    if (emd_imf_len <= 0)
    {
        qDebug() << "emd_imf_len=:" << emd_imf_len << " <= 0 ";
    }
    for(int i = 0; i < emd_imf_len; i++)
    {
        qDebug() << "emd_residue[" << i << "]=:" << emd_residue[i];
    }
}

int EmdClass::getEmdImfLen()
{
    return emd_imf_len;
}

int EmdClass::getEmdImfNum()
{
    return emd_imf_num;
}

int EmdClass::getEmdImfSetMaxNum()
{
    return emd_imf_max_num;
}

int EmdClass::getEmdIteration()
{
    return emd_iteration;
}

int EmdClass::getEmdImfSetMaxIteration()
{
    return emd_max_iteration;
}

int EmdClass::getEmdError()
{
    return emd_error;
}

double EmdClass::getEmdSdValue()
{
    return emd_sd;
}

double EmdClass::getEmdSd2Value()
{
    return emd_sd2;
}

double EmdClass::getEmdTolValue()
{
    return emd_tol;
}

void EmdClass::shiftSPMS(const int *buf_extr_x, const double *buf_extr_y, const int buf_extr_len,
                         double *buf_envmoy, double *buf_amp, bool choose)
{
    if ((buf_extr_x == nullptr) || (buf_extr_y == nullptr) || (buf_envmoy == nullptr))
    {
        emd_error = -1;
        return;
    }
    if (buf_extr_len <= 1)
    {
        emd_error = -2;
        return;
    }
    double a, b, c;
    if (choose)//TPDS
    {
        for (int i = 1; i < buf_extr_len - 1; i++)
        {
            a = 0.5 * ((buf_extr_x[i + 1] - buf_extr_x[i]) / (buf_extr_x[i + 1] - buf_extr_x[i  - 1]));
            b = 0.5;
            c = 0.5 * ((buf_extr_x[i] - buf_extr_x[i - 1]) / (buf_extr_x[i + 1] - buf_extr_x[i  - 1]));
            buf_envmoy[i - 1] = a * buf_extr_y[i -1] + b * buf_extr_y[i] + c * buf_extr_y[i + 1];
            buf_amp[i - 1] = a * buf_extr_y[i -1] - b * buf_extr_y[i] + c * buf_extr_y[i + 1];
        }
    }
    else//TPSS
    {
        for (int i = 1; i < buf_extr_len - 1; i++)
        {
            a = 0.5 * (1 + 2 * ((buf_extr_x[i] - buf_extr_x[i - 1]) / (buf_extr_x[i + 1] - buf_extr_x[i  - 1])))
                    * ((buf_extr_x[i + 1] - buf_extr_x[i]) / (buf_extr_x[i + 1] - buf_extr_x[i  - 1]))
                    * ((buf_extr_x[i + 1] - buf_extr_x[i]) / (buf_extr_x[i + 1] - buf_extr_x[i  - 1]));
            b = 0.5;
            c = 0.5 * (1 + 2 * ((buf_extr_x[i] - buf_extr_x[i + 1]) / (buf_extr_x[i + 1] - buf_extr_x[i  - 1])))
                    * ((buf_extr_x[i - 1] - buf_extr_x[i]) / (buf_extr_x[i + 1] - buf_extr_x[i  - 1]))
                    * ((buf_extr_x[i - 1] - buf_extr_x[i]) / (buf_extr_x[i + 1] - buf_extr_x[i  - 1]));
            buf_envmoy[i - 1] = a * buf_extr_y[i -1] + b * buf_extr_y[i] + c * buf_extr_y[i + 1];
            buf_amp[i - 1] = a * buf_extr_y[i -1] - b * buf_extr_y[i] + c * buf_extr_y[i + 1];
        }
    }
}

void EmdClass::setEmdStopShiftPara(const int sd, const int sd2, const int tol,
                                   const int imf_max_num, const int max_iteration)
{
    emd_sd = sd;
    emd_sd2 = sd2;
    emd_tol = tol;
    emd_imf_max_num = imf_max_num;
    emd_max_iteration = max_iteration;
}



bool EmdClass::isStopEmd(const double * buf_in_y, const int buf_in_len)
{
    if (buf_in_y == nullptr)
    {
         emd_error = -1;
         return true;
    }
    if (buf_in_len <= 0)
    {
        emd_error = -2;
        return true;
    }
    MyExtrZero->computeExtrResult(buf_in_y, buf_in_len);

    if (MyExtrZero->isExtrNumOK())
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool EmdClass::isStopEmd2(const double *buf_in_a, const double *buf_in_b, const int buf_in_len)
{
    if ((buf_in_a == nullptr) || (buf_in_b == nullptr))
    {
         emd_error = -1;
         return true;
    }
    if (buf_in_len <= 0)
    {
        emd_error = -2;
        return true;
    }
    double maxa = -5000;
    double maxb = -5000;
    for (int i = 0; i < buf_in_len; i++)
    {
        if(buf_in_a[i] > maxa)
        {
            maxa = buf_in_a[i];
        }
        if(buf_in_b[i] > maxb)
        {
            maxb = buf_in_b[i];
        }
    }
    if(maxa < maxb /10000.0)
    {
        return true;
    }
    else
    {
        return false;
    }
}

int EmdClass::getShiftStopFlag(const double * buf_in_x, const double * buf_in_y, const int buf_in_len)
{
    if ((buf_in_x == nullptr) || (buf_in_y == nullptr))
    {
         emd_error = -1;
         return -1;
    }
    if (buf_in_len <= 0)
    {
        emd_error = -2;
        return -2;
    }
    computeMeanAndAmp(buf_in_x, buf_in_y, buf_in_len);
    int ret = MyExtrZero->getExtrZeroError() + MyBoundary->getBoundaryError() + MySpline->getSplineError();
    if ( 0 != ret)
    {
        emd_error = -3;
        qDebug() << "everything ERROR!!";
        return -3;
    }

    double * test = new double [buf_in_len];
    double compTemp = 0;
    int compTemp2 = 0;
    for(int i = 0; i < buf_in_len; i++)
    {
        if(emd_envmoy[i] < 0)
        {
            test[i] = -emd_envmoy[i] / emd_amp[i];
        }
        else
        {
            test[i] = emd_envmoy[i] / emd_amp[i];
        }
//        if (emd_imf_num == TEST_IMF_NUM)
//        {
//            qDebug() << "++++++++++++" << test[i];
//        }
        if (test[i] > emd_sd)
        {
            compTemp++;
        }
        else
        {
            compTemp += 0;
        }
        if (test[i] > emd_sd2)
        {
            compTemp2 = 1;
        }
    }

    compTemp  = compTemp / buf_in_len;
    delete [] test;
    if (DEBUG_FLAG)
    {
        qDebug() << "compTemp = " << compTemp << endl
                 << "Total Extr:" << MyExtrZero->getExtrTotalNum()
                 << "Zero:" << MyExtrZero->getZeroNum();
    }
    if(!(((compTemp > emd_tol) || (compTemp2 == 1)) && (MyExtrZero->getExtrTotalNum() > 2)) && MyExtrZero->isExtrZeroIMF())
    {
        if (DEBUG_FLAG)
        {
            qDebug() << "getShiftStopFlag should stop shifting!" ;
        }
        return 1;
    }
    else
    {
        return 0;
    }
}


void EmdClass::computeMeanAndAmp(const double * buf_in_x, const double * buf_in_y, const int buf_in_len)
{
    if ((buf_in_x == nullptr) || (buf_in_y == nullptr))
    {
         emd_error = -1;
         return;
    }
    if (buf_in_len <= 0)
    {
        emd_error = -2;
        return;
    }
    emd_imf_len = buf_in_len;
    if (!MyExtrZero->isBufMonotonic(buf_in_x, buf_in_len))
    {
        emd_error = -3;
        return;
    }
    MyExtrZero->computeExtrResult(buf_in_y, buf_in_len);
    MyExtrZero->computeZeroResult(buf_in_y, buf_in_len);

    if (MyExtrZero->isExtrNumOK())
    {
        emd_error = -4;
    }
    else if (MyExtrZero->isExtrNumOtherZero())
    {
        emd_error = -5;
    }
    else
    {
        MyBoundary->computeBoundaryMainlyMirrorSymmetryCondition((double *)buf_in_y, buf_in_len,
                                                                 MyExtrZero->extr_max_x, MyExtrZero->getExtrMaxNum(),
                                                                 MyExtrZero->extr_min_x, MyExtrZero->getExtrMinNum());
    }

    if (MyBoundary->getBoundaryError() !=  0)
    {
        emd_error = -6;
        return;
    }
    /*
    int * mix_extr_x = new int [MyBoundary->getBoundaryTotalNum()];
    double * mix_extr_y = new double [MyBoundary->getBoundaryTotalNum()];
    MyBoundary->mixBoundaryExtr(mix_extr_x, mix_extr_y);

    double * mix_envmoy = new double [MyBoundary->getBoundaryTotalNum()];
    double * mix_amp = new double [MyBoundary->getBoundaryTotalNum()];
    shiftSPMS(mix_extr_x, mix_extr_y, MyBoundary->getBoundaryTotalNum(), mix_envmoy, mix_amp, true);
    if((emd_imf_num == TEST_IMF_NUM) && ((emd_iteration == TEST_ITERATION) || (emd_iteration == TEST_ITERATION + 1)))//&&(test_a == 0))
    {
        for (int i = 0; i < MyBoundary->getBoundaryTotalNum(); i++)
        {
            qDebug() << "x:" << mix_extr_x[i] << "y:" << mix_extr_y[i] << "env:" << mix_envmoy[i] << "amp:" << mix_amp[i];
        }
    }
    delete [] emd_envmoy;
    delete [] emd_amp;
    emd_envmoy = new double [buf_in_len];
    emd_amp = new double [buf_in_len];
    MySpline->computeThreeSpline(&mix_extr_x[1], mix_envmoy,
                                 MyBoundary->getBoundaryTotalNum() - 1, 0,
                                 emd_envmoy, buf_in_len);
    delete [] mix_extr_x;
    mix_extr_x = nullptr;
    delete [] mix_extr_y;
    mix_extr_y = nullptr;
    delete [] mix_envmoy;
    mix_envmoy = nullptr;
    */

    double *emd_envmoy_min = new double [buf_in_len];
    MySpline->computeThreeSpline(MyBoundary->boundary_x_min, MyBoundary->boundary_y_min,
                                 MyBoundary->getBoundaryMinNum(),MyBoundary->getBoundaryMinStart(),
                                 emd_envmoy_min, buf_in_len);

    double *emd_envmoy_max = new double [buf_in_len];
    MySpline->computeThreeSpline(MyBoundary->boundary_x_max, MyBoundary->boundary_y_max,
                                 MyBoundary->getBoundaryMaxNum(), MyBoundary->getBoundaryMaxStart(),
                                 emd_envmoy_max, buf_in_len);




    delete [] emd_envmoy;
    delete [] emd_amp;
    emd_envmoy = new double [buf_in_len];
    emd_amp = new double [buf_in_len];
    for(int i = 0; i < buf_in_len; i++)
    {
        emd_envmoy[i] = (emd_envmoy_min[i] + emd_envmoy_max[i]) / 2;
        emd_amp[i] = (emd_envmoy_min[i] - emd_envmoy_max[i])/2;
        if(emd_amp[i] < 0)
        {
            emd_amp[i] = -emd_amp[i];
        }
    }

    ////for test
//    if((emd_imf_num == TEST_IMF_NUM) && ((emd_iteration == TEST_ITERATION) || (emd_iteration == TEST_ITERATION + 1)))//&&(test_a == 0))
//    {

//    if (emd_imf_num == TEST_IMF_NUM)
//    {
//        qDebug() << MyExtrZero->getExtrMaxNum();
//        qDebug() << MyExtrZero->getExtrMinNum();
//        qDebug() << MyExtrZero->getZeroNum();
//        MyExtrZero->showAllExtrData();
//        MyExtrZero->showAllZeroData();
//        qDebug() << MyBoundary->getBoundaryMaxNum();
//        qDebug() << MyBoundary->getBoundaryMinNum();
//        MyBoundary->showAllBoundaryData();
//        qDebug() << "****************************************************************************************************************";
//        qDebug() << "-------------------------------max spline";
//        for (int i = 0; i < buf_in_len; i ++ )
//        {
//            qDebug() << "max[" << i << "]=" << emd_envmoy_max[i];
//        }
//        qDebug() << "-------------------------------min spline";
//        for (int i = 0; i < buf_in_len; i ++ )
//        {
//            qDebug() << "min[" << i << "]=" << emd_envmoy_min[i];
//        }
//        qDebug() << "-------------------------------env spline";
//        for (int i = 0; i < buf_in_len; i ++ )
//        {
//            qDebug() << "envmoy[" << i << "]=" << emd_envmoy[i];
//        }
//        qDebug() << "****************************************************************************************************************";
//    }
    delete [] emd_envmoy_max;
    emd_envmoy_max = nullptr;
    delete [] emd_envmoy_min;
    emd_envmoy_min = nullptr;

}



void EmdClass::computeEmd(const double * buf_in_x, const double * buf_in_y, const int buf_in_len)
{
    if ((buf_in_x == nullptr) || (buf_in_y == nullptr))
    {
         emd_error = -1;
         return;
    }
    if (buf_in_len <= 0)
    {
        emd_error = -2;
        return;
    }
    else
    {
        emd_imf_len = buf_in_len;
    }
    double * buf_y_temp = new double [buf_in_len];
    delete [] emd_residue;
    emd_residue = new double [buf_in_len];
    for (int i = 0; i < buf_in_len; i++)
    {
        buf_y_temp[i] = buf_in_y[i];
        emd_residue[i] = buf_in_y[i];
    }
    delete [] emd_imf;
    emd_imf = new double [emd_imf_max_num * buf_in_len];
    while (!isStopEmd(emd_residue, buf_in_len) && (emd_imf_num < emd_imf_max_num || emd_imf_max_num == 0))
    {
        if (isStopEmd2(buf_in_y, buf_y_temp, buf_in_len))
        {
            qDebug() << "forced stop of EMD : too small amplitude";
            break;
        }
        emd_iteration = 0;
//        qDebug() << "out+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++";
        while ((getShiftStopFlag(buf_in_x, buf_y_temp, buf_in_len) == 0) && (emd_iteration < emd_max_iteration))
        {
            if (DEBUG_FLAG)
            {
                qDebug() << "in----------------------------------------------------------------------------------------";
            }
            for (int i = 0; i < buf_in_len; i++)
            {
                buf_y_temp[i] = buf_y_temp[i] - emd_envmoy[i];
//                qDebug() << buf_y_temp;
//                if ((emd_imf_num == TEST_IMF_NUM) && ((emd_iteration == TEST_ITERATION) || (emd_iteration == TEST_ITERATION + 1)))
//                {
//                    qDebug() << "m++[" << i << "]=" << QString::number(buf_y_temp[i], 'g', 8) ;
////                    qDebug() << buf_y_temp[i];
//                }
            }
//            if (emd_imf_num == TEST_IMF_NUM)
//            {
//                qDebug() << "m[" << 4 << "]=" << QString::number(buf_y_temp[4], 'g', 8) ;
//            }
//            ret_for_stop_shifting = getShiftStopFlag(buf_in_x, buf_y_temp, buf_in_len);
            if (DEBUG_FLAG)
            {
                qDebug() << "emd_iteration£º" << emd_iteration << endl
                         << "in_error :::::" << emd_error << endl;
            }
            emd_iteration++;
        }
        for (int i = 0; i < buf_in_len; i++)
        {

            emd_imf[emd_imf_num * buf_in_len + i] = buf_y_temp[i];
//            if (emd_imf_num == TEST_IMF_NUM)
//            {
//                qDebug() << "imf" << TEST_IMF_NUM + 1 << "[" << i << "]=" << QString::number(buf_y_temp[i], 'g', 8) ;
//            }
            emd_residue[i] = emd_residue[i] - buf_y_temp[i];
            buf_y_temp[i] = emd_residue[i];
//            if (emd_imf_num == TEST_IMF_NUM)
//            {
//                qDebug() << "R" << TEST_IMF_NUM + 1 << "[" << i << "]=" << QString::number(emd_residue[i], 'g', 8) ;
//            }
        }
        emd_imf_num++;
        if (DEBUG_FLAG)
        {
            qDebug() << "Emd_imf_num=:" << emd_imf_num;/// << test_a;
            qDebug() << "emd_iteration_num=:" << emd_iteration;
        }
    }
    qDebug() << "Emd_end!";
    delete [] buf_y_temp;
    buf_y_temp = nullptr;
}
