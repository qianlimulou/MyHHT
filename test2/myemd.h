#ifndef MYEMD_H
#define MYEMD_H
#define LENGTH 6000
#include <QDebug>
#include <cmath>
#include <QWidget>
//#include "spline.h"

typedef struct _threePoints   //����ֵ/��Сֵ/��� λ��/����
{
    int posPointMin[LENGTH];
    int posPointMax[LENGTH];
    int posPointZero[LENGTH];
    int lenPointMin;
    int lenPointMax;
    int lenPointZero;
    int count_time;
    int error;
}threePoints;

typedef struct _boundary    //�߽�����
{
    double tmin[LENGTH];
    int     lentmin;
    double tmax[LENGTH];
    int     lentmax;
    double zmin[LENGTH];
    int     lenzmin;
    double zmax[LENGTH];
    int     lenzmax;
    int error;
}boundary;

typedef struct _envelopeLine //��ֵ֮����
{
    int length;
    double  envmin[LENGTH];
    double  envmax[LENGTH];
    double  envmoy[LENGTH];
    double  amp[LENGTH];
    int nem;
    int nzm;
    int error;
}envelopeLine;

typedef struct _stopSiftingState  //ֹͣɸѡ״̬
{
    int stop;
    double s;
    double  envmoy[LENGTH];
    int error;
}stopSiftingState;

typedef struct _IMFS    //����ģ̬��������
{
    double  *imf;
    int length;
    int piles;
    int error;
}IMFS;

typedef enum _condition //����������ֵ�߽�ѡ��
{
    NATURAL,
    CLAMPED,
    NOTaKNOT
}condition_spline;

//    int getThreePointsMinNum(threePoints * threePoints){return threePoints->lenPointMin;}
//    int getThreePointsMaxNum(threePoints * threePoints){return threePoints->lenPointMax;}
//    int getThreePointsZeroNum(threePoints * threePoints){return threePoints->lenPointZero;}
//    int getThreePointsError(threePoints * threePoints){return threePoints->error;}
//    int getThreePointsCount(threePoints * threePoints){return threePoints->count_time;}
//    int getBoundaryTminNum(boundary * boundary){return boundary->lentmin;}
//    int getBoundaryTmanNum(boundary * boundary){return boundary->lentmax;}
//    int getBoundaryZminNum(boundary * boundary){return boundary->lenzmin;}
//    int getBoundaryZmaxNum(boundary * boundary){return boundary->lenzmax;}
//    int getBoundaryError(boundary * boundary){return boundary->error;}
    void showBoundaryValue(boundary * boundary);
//    int getEnvelopeLineNem(envelopeLine * envelopeLine){return envelopeLine->nem;}
//    int getEnvelopeLineNzm(envelopeLine * envelopeLine){return envelopeLine->nzm;}
//    int getEnvelopeLineError(envelopeLine * envelopeLine){return envelopeLine->error;}
    void computeExtremumAndZeroPoints(threePoints * threePoints, const double * buf_in, const int buf_in_len);
    bool judgeExtremumAndZeroPointsIsSatisfyIMF(threePoints * threePoints);
    bool judgeIsMonotonic();//const double * buf_in, const int buf_in_len);
    void setBoundaryCondition(threePoints * threePoints, boundary * boundary, const double * buf_in, const int buf_in_len);
    void computeSplineResult(const double * x, const double * y, const int xy_len,
                             double * result,  const int buf_in_len, const condition_spline con,
                             const int A, const int B);
    void proceedEmdAndJudgeStopRequirement();//double * buf_in, int buf_in_len);


#endif // MYEMD_H
