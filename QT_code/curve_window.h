#ifndef CURVE_WINDOW_H
#define CURVE_WINDOW_H

#include <QWidget>
#include <QPainter>

namespace Ui {
class curve_window;
}

class curve_window : public QWidget
{
    Q_OBJECT

public:
    explicit curve_window(QWidget *parent = 0);
    ~curve_window();
    void tau_pos_data_get(double *tau_pos_data1, double *tau_pos_data2, double *tau_pos_data3, double *tau_pos_data4, double *tau_pos_data5,
                      double *tau_pos_data6, double *_ma1, double *_ma2,double *_ma3,double *_ma4,double *_ma5,double *_ma6);
    void tau_data_get(double *tau_data1, double *tau_data2, double *tau_data3, double *tau_data4, double *tau_data5,
                      double *tau_data6);

protected:
    void paintEvent(QPaintEvent *);
private slots:

    void on_ButtonJoint1_clicked();

    void on_ButtonJoint2_clicked();

    void on_ButtonJoint3_clicked();

    void on_ButtonJoint4_clicked();

    void on_ButtonJoint5_clicked();

    void on_ButtonJoint6_clicked();

private:
    Ui::curve_window *ui;
    QImage image;
    double *tau_pos_data, *tau_data;
    double tau_pos_data1[100];
    double tau_pos_data2[100];
    double tau_pos_data3[100];
    double tau_pos_data4[100];
    double tau_pos_data5[100];
    double tau_pos_data6[100];
    double tau_data1[100];
    double tau_data2[100];
    double tau_data3[100];
    double tau_data4[100];
    double tau_data5[100];
    double tau_data6[100];
    double _ma = 0, _ma1 = 0.0,_ma2 = 0.0,_ma3 = 0.0,_ma4 = 0.0,_ma5 = 0.0,_ma6 = 0.0;
    double _mal = 0, _mal1 = 0.0,_mal2 = 0.0,_mal3 = 0.0,_mal4 = 0.0,_mal5 = 0.0,_mal6 = 0.0;


};

#endif // CURVE_WINDOW_H
