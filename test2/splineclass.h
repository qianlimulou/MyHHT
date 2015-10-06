#ifndef SPLINECLASS_H
#define SPLINECLASS_H
//TEST OVER
#include <QDebug>
class SplineClass
{
public:
    typedef enum _condition
    {
        NATURAL,
        CLAMPED,
        NOTaKNOT
    }spline_boudary_condition;
    SplineClass();
//    SplineClass(double *buf_in_x, double *buf_in_y, int in_len,
//                double *buf_out_x, double *buf_out_y, int out_len);
//    SplineClass(double * buf_in_x, double * buf_in_y, int in_len,
//                double * buf_out_x, double * buf_out_y, int out_len,
//                SplineClass::spline_boudary_condition con, int a, int b);
    ~SplineClass();
    int getSplineError();
    int getSplineOutLen();
//    void showAllSplineData();
    void setThreeSplineCondition(SplineClass::spline_boudary_condition con,
                                 int a, int b);
    void computeThreeSpline(const int *buf_in_x, const double *buf_in_y, const int in_len, const int start_pos,
                            /*double *buf_out_x, */double *buf_out_y, const int out_len);
private:
    spline_boudary_condition spline_con;
    int spline_out_len;
//    double * spline_buf_out_x;
//    double * spline_buf_out_y;
    int spline_boudary_A;
    int spline_boudary_B;
    int spline_error;

};

#endif // SPLINECLASS_H
