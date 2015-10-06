#ifndef BOUNDARYCLASS_H
#define BOUNDARYCLASS_H

#include <QDebug>

#define DEBUG_FOR_BOUNDARY2 0

class BoundaryClass
{
public:
    BoundaryClass();
//    BoundaryClass();
    ~BoundaryClass();
    int getBoundaryError();
    int getBoundaryMinNum();
    int getBoundaryMaxNum();
    int getBoundaryMinStart();
    int getBoundaryMaxStart();
    int getBoundaryTotalNum();
    void showAllBoundaryData();
    void mixBoundaryExtr(int * mix_extr_x, double * mix_extr_y);//保证mix_extr_x的大小至少为极值点之和
    void computeBoundaryMainlyMirrorSymmetryCondition(double * buf_in, const int buf_in_len,
                              const int * max_pos, const int max_num,
                              const int * min_pos, const int min_num);
    void computeBoundaryCharacteristicWaveCondition(double * buf_in, const int buf_in_len,
                              const int * max_pos, const int max_num,
                              const int * min_pos, const int min_num);
    int    *boundary_x_min;
    double *boundary_y_min;
    int    *boundary_x_max;
    double *boundary_y_max;

private:
    int     boundary_error;
    int     boundary_min_num;
    int     boundary_max_num;
    int     boundary_min_start;
    int     boundary_max_start;
};

#endif // BOUNDARYCLASS_H
