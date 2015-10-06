#ifndef WIDGET_H
#define WIDGET_H

#include <QWidget>
#include "myemd.h"
#include "plotter.h"
#include <QDebug>
//#include "spline.h"
#include <iostream>
#include "splineclass.h"
#include "emdclass.h"
#include "fileclass.h"
#include "hilbertclass.h"

#define DATA_LENGTH 200
namespace Ui {
class Widget;
}

class Widget : public QWidget
{
    Q_OBJECT

public:
    explicit Widget(QWidget *parent = 0);
    ~Widget();
//    myemd ohmy;

private slots:
    void openFile(QString filename);
    void on_startBtn_clicked();
    void on_testBtn_clicked();

    void on_forwardBtn_clicked();

    void on_backBtn_clicked();

    void on_exitBtn_clicked();

    void on_testHilbertBtn_clicked();

    void on_testForBtn_clicked();

    void on_hilbertBtn_clicked();

private:
    Ui::Widget *ui;
    Plotter* plot;
    SplineClass *mySpline;
    EmdClass *myEmd;
    FileClass *myFile;
    ExtrZeroClass *myEZ;
    BoundaryClass *myBndr;
    HilbertClass *myHilbert;
    double * x;
    double * imf;
    int imfCounter;
    int imfLength;


};

#endif // WIDGET_H
