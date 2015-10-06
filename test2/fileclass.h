#ifndef FILECLASS_H
#define FILECLASS_H

#include <QString>
#include <QFile>
#include <QTextStream>
#include <QDebug>

class FileClass
{
public:
    FileClass();
    bool openFile(const QString fileName);
    void closeFile();
    bool read();
    bool save();
    int readFileToBuf(const QString fileName, double *buf_out, int buf_out_len);
    int getTotalLine();
    void showTotalLine();
    QString getStrall();
private:
    QString strall;  //��ȡ������
    int totalLine;       //����
    QFile * fp;

};

#endif // FILECLASS_H
