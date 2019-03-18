#ifndef DYNAMICWIDGET_H
#define DYNAMICWIDGET_H

#include <QWidget>

namespace Ui {
class DynamicWidget;
}

class DynamicWidget : public QWidget
{
    Q_OBJECT

public:
    explicit DynamicWidget(QWidget *parent = 0);
    ~DynamicWidget();
    void a_data_read();
    void alpha_data_read();
    void d_data_read();
    void result_read();
    void tau_pos_read();
    void phi_read();
    void error_read();


private slots:
    void on_pushButton_clicked();


    void on_pushButton_2_clicked();

    void on_Button_Results_clicked();


private:
    Ui::DynamicWidget *ui;
};

#endif // DYNAMICWIDGET_H
