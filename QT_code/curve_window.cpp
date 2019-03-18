#include "curve_window.h"
#include "ui_curve_window.h"
#include <QPaintEvent>
#include <QtGui>
#include <QFile>
#include <QDebug>
#include <QPen>
#define n 100

curve_window::curve_window(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::curve_window)
{
    ui->setupUi(this);


    tau_pos_data_get(tau_pos_data1,tau_pos_data2,tau_pos_data3,tau_pos_data4,tau_pos_data5,tau_pos_data6,&_ma1,&_ma2,&_ma3,&_ma4,&_ma5,&_ma6);
    tau_data_get(tau_data1,tau_data2,tau_data3,tau_data4,tau_data5,tau_data6);

    image = QImage(600, 400, QImage::Format_RGB32);
    QColor backColor = qRgb(255, 255, 255);
    image.fill(backColor);
    tau_data = tau_data1;
    tau_pos_data = tau_pos_data1;
    _ma = _ma1;



}

curve_window::~curve_window()
{
    delete ui;
}
void curve_window::tau_pos_data_get(double *tau_pos_data1, double *tau_pos_data2, double *tau_pos_data3, double *tau_pos_data4, double *tau_pos_data5,
                                double *tau_pos_data6,double *_ma1, double *_ma2,double *_ma3,double *_ma4,double *_ma5,double *_ma6)
{
    QString filename = "tau_pos.txt";
    QFile file(filename);
    bool isOk = file.open(QIODevice::ReadOnly);
    if(isOk)
    {
        int i = 0;
        QString tau_pos_array1;
        QString tau_pos_array2;
        QString tau_pos_array3;
        QString tau_pos_array4;
        QString tau_pos_array5;
        QString tau_pos_array6;
        while (file.atEnd() == false && i < n)
        {
            tau_pos_array1 = file.readLine();
            tau_pos_array2 = file.readLine();
            tau_pos_array3 = file.readLine();
            tau_pos_array4 = file.readLine();
            tau_pos_array5 = file.readLine();
            tau_pos_array6 = file.readLine();
            tau_pos_data1[i] = tau_pos_array1.toDouble();
            tau_pos_data2[i] = tau_pos_array2.toDouble();
            tau_pos_data3[i] = tau_pos_array3.toDouble();
            tau_pos_data4[i] = tau_pos_array4.toDouble();
            tau_pos_data5[i] = tau_pos_array5.toDouble();
            tau_pos_data6[i] = tau_pos_array6.toDouble();
            if(tau_pos_data1[i] > *_ma1)
            {
                 *_ma1 = tau_pos_data1[i];
            }
            if(tau_pos_data2[i] > *_ma2)
            {
                 *_ma2 = tau_pos_data2[i];
            }
            if(tau_pos_data3[i] > *_ma3)
            {
                 *_ma3 = tau_pos_data3[i];
            }
            if(tau_pos_data4[i] > *_ma4)
            {
                 *_ma4 = tau_pos_data4[i];
            }
            if(tau_pos_data5[i] > *_ma5)
            {
                 *_ma5 = tau_pos_data5[i];
            }
            if(tau_pos_data6[i] > *_ma6)
            {
                 *_ma6 = tau_pos_data6[i];
            }
            ++i;
        }
    }
    file.close();
}
void curve_window::tau_data_get(double *tau_data1, double *tau_data2, double *tau_data3, double *tau_data4, double *tau_data5,
                                double *tau_data6)
{
    QString filename = "tau_data.txt";
    QFile file(filename);
    bool isOk = file.open(QIODevice::ReadOnly);
    if(isOk)
    {
        int i = 0;
        QString tau_array1;
        QString tau_array2;
        QString tau_array3;
        QString tau_array4;
        QString tau_array5;
        QString tau_array6;
        while (file.atEnd() == false && i < n)
        {
            tau_array1 = file.readLine();
            tau_array2 = file.readLine();
            tau_array3 = file.readLine();
            tau_array4 = file.readLine();
            tau_array5 = file.readLine();
            tau_array6 = file.readLine();
            tau_data1[i] = tau_array1.toDouble();
            tau_data2[i] = tau_array2.toDouble();
            tau_data3[i] = tau_array3.toDouble();
            tau_data4[i] = tau_array4.toDouble();
            tau_data5[i] = tau_array5.toDouble();
            tau_data6[i] = tau_array6.toDouble();
            ++i;
        }
    }
    file.close();
}

void curve_window::paintEvent(QPaintEvent *){
    // QPainter painter1(this);
    // painter1.drawImage(0, 0, image);
    // QPainter painter(&image);
    // painter.begin(&image);
    QPainter painter(this);
    painter.begin(this);
    painter.drawImage(0, 0, image);

    painter.setRenderHint(QPainter::Antialiasing, true);
    int pointX = 55, pointY = 280;
    int width = 580 - pointX, height = 260;
    painter.drawRect(5,5, 600-10, 400-10);
    painter.drawLine(pointX, pointY, width+pointX, pointY);
    painter.drawLine(pointX, pointY-height, pointX, pointY);
    painter.drawLine(pointX, pointY, pointX, pointY+100);

    double kx = (double)width/(n-1);
    double ky = (double)height/_ma;
    QPen pen, pen2, penPoint, penPoint2;
    pen.setColor(Qt::blue);
    pen.setWidth(1);
    penPoint.setColor(Qt::blue);
    penPoint.setWidth(3);
    pen2.setColor(Qt::red);
    pen2.setWidth(1);
    penPoint2.setColor(Qt::red);
    penPoint2.setWidth(3);

    // draw data point
    for(int i=0;i<n-1;i++)
    {
        //由于y轴是倒着的，所以y轴坐标要pointy-a[i]*ky 其中ky为比例系数
        painter.setPen(pen);//黑色笔用于连线
        painter.drawLine(pointX+kx*i,pointY-tau_pos_data[i]*ky,pointX+kx*(i+1),pointY-tau_pos_data[i+1]*ky);
        painter.setPen(penPoint);//蓝色的笔，用于标记各个点
        painter.drawPoint(pointX+kx*i,pointY-tau_pos_data[i]*ky);
    }
    painter.drawPoint(pointX+kx*(n-1),pointY-tau_pos_data[n-1]*ky);

    for(int i=0;i<n-1;i++)
    {
        painter.setPen(pen2);
        painter.drawLine(pointX+kx*i,pointY-tau_data[i]*ky,pointX+kx*(i+1),pointY-tau_data[i+1]*ky);
        painter.setPen(penPoint2);
        painter.drawPoint(pointX+kx*i,pointY-tau_data[i]*ky);
    }
    painter.drawPoint(pointX+kx*(n-1),pointY-tau_data[n-1]*ky);



    // 绘制刻度线
    QPen penDegree;
    penDegree.setColor(Qt::black);
    penDegree.setWidth(1);
    painter.setPen(penDegree);

    //画上x轴刻度线
    for(int i=0;i<10;i++)//分成10份
    {
        //选取合适的坐标，绘制一段长度为4的直线，用于表示刻度
        painter.drawLine(pointX+(i+1)*width/10,pointY,pointX+(i+1)*width/10,pointY+4);
        painter.drawText(pointX+(i+0.65)*width/10,
                         pointY+10,QString::number((int)((i+1)*((double)n/10))));
    }
    //y轴刻度线
    double _maStep=(double)_ma/10;//y轴刻度间隔需根据最大值来表示
    for(int i=0;i<10;i++)
    {
        //主要就是确定一个位置，然后画一条短短的直线表示刻度。
        painter.drawLine(pointX,pointY-(i+1)*height/10,
                         pointX-4,pointY-(i+1)*height/10);
        painter.drawText(pointX-35,pointY-(i+0.85)*height/10,
                         QString::number((int)(_maStep*(i+1))));
    }
    for(int i=0;i<3;i++)
    {
        //主要就是确定一个位置，然后画一条短短的直线表示刻度。
        painter.drawLine(pointX,pointY+(i+1)*height/10,
                         pointX-4,pointY+(i+1)*height/10);
        painter.drawText(pointX-35,pointY+(i+1.15)*height/10,
                         QString::number((int)(-_maStep*(i+1))));
    }

    painter.end();


}



void curve_window::on_ButtonJoint1_clicked()
{
    tau_pos_data = tau_pos_data1;
    tau_data = tau_data1;
    _ma = _ma1;
    update();

}

void curve_window::on_ButtonJoint2_clicked()
{
    tau_pos_data = tau_pos_data2;
    tau_data = tau_data2;
    _ma = _ma2;
    update();
}

void curve_window::on_ButtonJoint3_clicked()
{
    tau_pos_data = tau_pos_data3;
    _ma = _ma3;
    tau_data = tau_data3;
    update();
}

void curve_window::on_ButtonJoint4_clicked()
{
    tau_pos_data = tau_pos_data4;
    tau_data = tau_data4;
    _ma = _ma4;
    update();
}

void curve_window::on_ButtonJoint5_clicked()
{
    tau_pos_data = tau_pos_data5;
    tau_data = tau_data5;
    _ma = _ma5;
    update();
}

void curve_window::on_ButtonJoint6_clicked()
{
    tau_pos_data = tau_pos_data6;
    tau_data = tau_data6;
    _ma = _ma6;
    update();
}
