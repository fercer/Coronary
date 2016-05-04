#include "coronarygui.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    coronaryGUI w;
    w.show();

    return a.exec();
}
