#include "widget.h"
#include "ui_widget.h"

Widget::Widget(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::Widget)
{
    ui->setupUi(this);
    imfCounter = 0;
    imfLength = 0;
    plot = new Plotter;
    x = new double[DATA_LENGTH];
    imf = nullptr;

    ui->verticalLayout->addWidget(plot);
    ui->startBtn->setEnabled(false);
    ui->lineEdit->setText(tr("E:\\qt project\\test\\EMD_FINAL\\data1.txt"));
    connect(ui->lineEdit,SIGNAL(textChanged(QString)),this,SLOT(openFile(QString)));
    myEmd = new EmdClass;
    myFile = new FileClass;
    myEZ = new ExtrZeroClass;
    myBndr = new BoundaryClass;
    mySpline = new SplineClass;
    myHilbert = new HilbertClass(DATA_LENGTH);
//    myHilbert = new HilbertClass(DATA_LENGTH);

}

Widget::~Widget()
{
    delete [] x;
    delete [] imf;
    delete plot;
    delete ui;
}

void Widget::openFile(QString filename)
{
    qDebug() <<filename;

    int len;
    QFile filex(filename.trimmed());
    if(filex.exists())
    {
        filex.close();
        ui->startBtn->setEnabled(true);
        len = myFile->readFileToBuf(filename, x, DATA_LENGTH);
        qDebug() << "len++++++++++++++++++++++++=--" << len  << filename;
    }
    else
    {
        ui->startBtn->setEnabled(false);
    }
}

void Widget::on_startBtn_clicked()
{
//    double *re1 = new double [DATA_LENGTH];
//    double *re2 = new double [DATA_LENGTH];
    double t[DATA_LENGTH];
    for (int i = 0; i < DATA_LENGTH; i++)
    {
        t[i] = i;
//        qDebug() << x[i];
    }
/**************************************************************************************
 TEST
//    myEZ->computeExtrResult(x, DATA_LENGTH);
//    myEZ->computeZeroResult(x, DATA_LENGTH);
//    qDebug() << myEZ->getExtrTotalNum();
//    qDebug() << myEZ->getZeroNum();
//    qDebug() << myEZ->isExtrNumOK();
//    myEZ->showAllExtrData();
//    myBndr->computeBoundaryCharacteristicWaveCondition(x, DATA_LENGTH,
//                                     myEZ->extr_max_x, myEZ->getExtrMaxNum(),
//                                     myEZ->extr_min_x, myEZ->getExtrMinNum());
//    qDebug() << "Bndr_error" << myBndr->getBoundaryError();
//    qDebug() << "Bndr_maxnum" << myBndr->getBoundaryMaxNum();
//    qDebug() << "Bndr_minnum" << myBndr->getBoundaryMinNum();
//    qDebug() << "Bndr_minstart" << myBndr->getBoundaryMinStart();
//    qDebug() << "Bndr_maxstart" << myBndr->getBoundaryMaxStart();
//    myBndr->showAllBoundaryData();
//    mySpline->computeThreeSpline(myBndr->boundary_x_min, myBndr->boundary_y_min,
//                                 myBndr->getBoundaryMinNum(), myBndr->getBoundaryMinStart(),
//                                 re1, DATA_LENGTH);
//    mySpline->computeThreeSpline(myBndr->boundary_x_max, myBndr->boundary_y_max,
//                                 myBndr->getBoundaryMaxNum(), myBndr->getBoundaryMaxStart(),
//                                 re2, DATA_LENGTH);
//    qDebug() << "Spline_error" << mySpline->getSplineError();
**********************************************************************************************/
//

////    myEmd->computeMeanAndAmp(t, x, DATA_LENGTH);
////    qDebug() << myEmd->getShiftStopFlag(t, x, DATA_LENGTH);
////    myEmd->showAllEnvmoyAndAmp();
////    qDebug() << "IsStopEmd:" << myEmd->isStopEmd(x, DATA_LENGTH);

    myEmd->computeEmd(t, x, DATA_LENGTH);
    qDebug() << "层数" << myEmd->getEmdImfNum();
    qDebug() << "次数" << myEmd->getEmdIteration();
    qDebug() << "Emd error" << myEmd->getEmdError();
    imfLength = myEmd->getEmdImfNum() + 1;
//    qDebug() << "why3";
//    myEmd->showAllEmdResidueData();
    if(imf != nullptr)
    {
        delete [] imf;
        imf = nullptr;
    }
    imf = new double[imfLength * DATA_LENGTH];
//    qDebug() << "why1";
    for(int i = 0; i < imfLength * DATA_LENGTH; i++)
    {
        imf[i] = myEmd->emd_imf[i];
        if (ui->radioButton->isChecked())
        {
            qDebug() << "imf[" << i << "]=" << imf[i];
        }
    }
//    qDebug() << "why";
    plot->curve(imf, DATA_LENGTH);
//    delete [] drawH;
//    delete []re1; re1 = nullptr;
//    delete []re2; re2 = nullptr;
}


void Widget::on_testBtn_clicked()
{
//    int x[3] = {-298, 300, 896};
//    double y[3] = {12, 125.835, 250};
//    double * drawX = new double[DATA_LENGTH];
//    mySpline->computeThreeSpline(x, y, 3, 0, drawX, DATA_LENGTH);

    plot->curve(x, DATA_LENGTH);

//    delete [] drawX; drawX = nullptr;

}

void Widget::on_forwardBtn_clicked()
{
    if(imf == nullptr)
        return;
    imfCounter--;
    if(imfCounter <= 0)
    {
        imfCounter = 0;
    }
//    labNotice->setText(tr("共%1层IMF，当前第%2层").arg(imfLength).arg(imfCounter+1));
    plot->curve(imf + imfCounter * DATA_LENGTH, DATA_LENGTH);
}

void Widget::on_backBtn_clicked()
{
    if(imf == nullptr)
        return;
    imfCounter++;
    if(imfCounter >= imfLength)
    {
        imfCounter = imfLength - 1;
    }
//    labNotice->setText(tr("共%1层IMF，当前第%2层").arg(imfLength).arg(imfCounter+1));
    plot->curve(imf + imfCounter * DATA_LENGTH, DATA_LENGTH);
}

void Widget::on_exitBtn_clicked()
{
    close();
}

void Widget::on_testHilbertBtn_clicked()
{
//    double x[8] = {1,2,2,3,4,5,6,7};
//    double y[8] = {0};
    double y[DATA_LENGTH]= {0};
//    myHilbert = new HilbertClass(DATA_LENGTH);
    qDebug() << myHilbert->getHilbertLen();
//    int a = 21,b = 0;
//    if ((a & (a - 1)) != 0)
//    {
//        while (b <= 32)
//        {
//            if ((1 << b) > a)
//            {
//                break;
//            }
//            b++;
//        }
//    }
//    qDebug() << "b=" << b;
//    myHilbert->showHilbertData();
//    myHilbert->computeDIFFft(x, y, 0);
//    myHilbert->computeHilbert(y);
    myHilbert->computeInstantFrequency(x, y);
    qDebug() << "error" << myHilbert->getHilbertError();
    plot->curve(y, DATA_LENGTH);
//    for(int i = 0; i < 8; i++)
//    {
//        qDebug() << "x[" << i << "]=:" <<x[i];
//    }
    for(int i = 0; i < DATA_LENGTH; i++)
    {
        qDebug() << "y[" << i << "]=:" <<y[i];
    }
//    double a = 0.00023456789012345678901234567890;
//    qDebug() << a;
//    qDebug() << QString("%1").arg(a) << QString::number(a, 'g', 8) ;
//    delete myHilbert;
}

void Widget::on_testForBtn_clicked()
{
    myEZ->computeExtrResult(x, DATA_LENGTH);
    myEZ->computeZeroResult(x, DATA_LENGTH);
    qDebug() << myEZ->getExtrMaxNum();
    qDebug() << myEZ->getExtrMinNum();
    qDebug() << myEZ->getZeroNum();
    myEZ->showAllExtrData();
    myEZ->showAllZeroData();
}

void Widget::on_hilbertBtn_clicked()
{
    if(imf == nullptr)
        return;
    double * y = new double [DATA_LENGTH];
    for(int i = 0; i < imfLength; i++)
    {
//        qDebug() << "why";
        myHilbert->computeInstantFrequency(&imf[i * DATA_LENGTH], y);
//        qDebug() << "soga:" << i;
        for (int j = 0; j < DATA_LENGTH; j++)
        {
            imf[j + i * DATA_LENGTH] = y[j];
        }
//        qDebug() << "what";
    }
//    for (int i = 0; i < imfLength * DATA_LENGTH; i++)
//    {
//        imf[i] = y[i];
//    }
//     myHilbert->computeInstantFrequency(&imf[DATA_LENGTH], &y[DATA_LENGTH]);
    plot->curve(imf, DATA_LENGTH);
    delete [] y;
    y = nullptr;
}
