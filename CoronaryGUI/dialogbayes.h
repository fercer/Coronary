#ifndef DIALOGBAYES_H
#define DIALOGBAYES_H

#include <QDialog>

namespace Ui {
    class DialogBayes;
}

class DialogBayes : public QDialog
{
        Q_OBJECT

    public:
        explicit DialogBayes(QWidget *parent = 0);
        ~DialogBayes();

    private:
        Ui::DialogBayes *ui;
};

#endif // DIALOGBAYES_H
