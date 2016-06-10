#include "coronarygui.h"
#include <QApplication>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int main(int argc, char *argv[])
{

    QApplication a(argc, argv);
    coronaryGUI w;
    w.show();

    return a.exec();
}
