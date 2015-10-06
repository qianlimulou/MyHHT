#include "widget.h"
#include <QApplication>
#include <QTextCodec>
int main(int argc, char *argv[])
{
    QTextCodec::setCodecForTr(QTextCodec::codecForName("gb2312"));
    QTextCodec::setCodecForCStrings(QTextCodec::codecForName("gb2312"));
    QTextCodec::setCodecForLocale(QTextCodec::codecForName("gb2312"));
    QApplication a(argc, argv);
    Widget w;
    w.show();

    return a.exec();
}
