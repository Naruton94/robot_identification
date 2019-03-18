#include "dynamicwidget.h"
#include <QApplication>


int main(int argc, char *argv[])
{

    QApplication a(argc, argv);
    DynamicWidget w;

    w.setWindowTitle("动力学参数辨识");
    w.show();

    return a.exec();
}
