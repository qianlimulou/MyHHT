#include "fileclass.h"

FileClass::FileClass()
{
    this->fp = NULL;
    this->totalLine = 0;
    this->strall = "";
}

bool FileClass::openFile(const QString fileName)
{
    this->fp = new QFile(fileName.trimmed());
    if (this->fp != NULL)
    {
        return true;
    }
    else
    {
        return false;
    }
}

void FileClass::closeFile()
{
    delete this->fp;
}

bool FileClass::read() // testing passed
{
    if (this->fp->open(QIODevice::ReadOnly))
    {
        QTextStream stream(this->fp);
        strall = stream.readAll();
        this->fp->close();
    }
    else
    {
        totalLine = 0;
        return false;
    }

    int Index = 0;
    if (strall.length() == 0)
    {
        totalLine = 0;
        return false;
    }

    Index = strall.indexOf('\n');
    Index++;
    totalLine = 1;

    while(Index != 0)       //统计总行数,从1开始计数
    {
        Index = strall.indexOf('\n', Index);
        Index++;
        totalLine = totalLine + 1;
    }

    // re-style strall
    // 1. strall ends with \r\n

    if (strall.length() >= 2)
    {
        if (strall.at(strall.length() - 1) != '\r'
            || strall.at(strall.length()) != '\n')
        {
            strall[strall.length() - 1] = '\r';
            strall[strall.length()] = '\n';
        }
    }

    return true;
}

bool FileClass::save() // testing passed
{
    this->fp->resize(0);

    if (this->fp->open(QIODevice::ReadWrite))
    {
        QTextStream stream(this->fp);
        stream << strall;
        return true;
    }
    else
    {
        return false;
    }
}

int FileClass::readFileToBuf(const QString fileName, double *buf_out, int buf_out_len)
{
    openFile(fileName);
//    qDebug() << "why";
    if (this->fp->open(QIODevice::ReadOnly))
    {
        QTextStream inx(this->fp);
//        qDebug() << "what";
        for(int i = 0; i < buf_out_len; i++)
//        while(!this->fp->atEnd())

        {
            buf_out[i] = inx.readLine().toDouble();
//            if (this->fp->atEnd())
//            {
//                qDebug() << "why1";
//                this->totalLine = i;
//                break;
//            }
//            qDebug() << i << ":" <<buf_out[i];
//            this->totalLine++;
//            if (this->totalLine == buf_out_len)
//            {
//                break;
//            }
        }
        this->fp->close();
        closeFile();
        return this->totalLine;
    }
    else
    {
        totalLine = 0;
        return -1;
    }

}

void FileClass::showTotalLine()
{
    qDebug() << totalLine << endl;
}

int FileClass::getTotalLine()
{
    return this->totalLine;
}

QString FileClass::getStrall()
{
    return strall;
}
