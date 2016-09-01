#ifndef FILTERPARS_H
#define FILTERPARS_H

#include <QDialog>

namespace Ui {
class filterpars;
}

class filterpars : public QDialog
{
    Q_OBJECT

public:
    explicit filterpars(QWidget *parent = 0);
    ~filterpars();

private slots:
    void on_chkTune_toggled(bool checked);

    void on_chkLFix_toggled(bool checked);

    void on_cmbFilter_currentIndexChanged(int index);

private:
    Ui::filterpars *ui;
};

#endif // FILTERPARS_H
