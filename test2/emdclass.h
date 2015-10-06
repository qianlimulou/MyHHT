#ifndef EMDCLASS_H
#define EMDCLASS_H

#include <QDebug>

#include "extrzeroclass.h"
#include "boundaryclass.h"
#include "splineclass.h"

#define TEST_IMF_NUM 0
#define TEST_ITERATION 0
#define DEBUG_FLAG 0
class EmdClass
{
public:
    EmdClass();
    ~EmdClass();
    void showAllEnvmoyAndAmp();
    void showAllEmdResidueData();
    void showAllEmdImfData();
    int getEmdImfLen();
    int getEmdImfNum();
    int getEmdImfSetMaxNum();
    int getEmdIteration();
    int getEmdImfSetMaxIteration();
    int getEmdError();
    double getEmdSdValue();
    double getEmdSd2Value();
    double getEmdTolValue();
    //CHOOSE: TRUE = TPDS; FALSE = TPSS
    //buf_envmoy 至少长度为buf_extr_len - 1
    void shiftSPMS(const int * buf_extr_x, const double * buf_extr_y, const int buf_extr_len, double * buf_envmoy, double *buf_amp, bool choose);
    void setEmdStopShiftPara(const int sd, const int sd2, const int tol, const int imf_max_num, const int max_iteration);
    bool isStopEmd(const double * buf_in_y, const int buf_in_len);
    bool isStopEmd2(const double * buf_in_a, const double * buf_in_b, const int buf_in_len);
    int getShiftStopFlag(const double * buf_in_x, const double * buf_in_y, const int buf_in_len);
    void computeMeanAndAmp(const double * buf_in_x, const double * buf_in_y, const int buf_in_len);
    void computeEmd(const double * buf_in_x, const double * buf_in_y, const int buf_in_len);
    double * emd_imf;
    double * emd_residue;
private:
    ExtrZeroClass *MyExtrZero;
    BoundaryClass *MyBoundary;
    SplineClass *MySpline;
    double * emd_envmoy;
    double * emd_amp;
    int emd_imf_len;
    int emd_imf_num;
    int emd_imf_max_num;
    int emd_iteration;
    int emd_max_iteration;
    int emd_error;
    double emd_sd;
    double emd_sd2;
    double emd_tol;
    int test_a;

};

#endif // EMDCLASS_H
