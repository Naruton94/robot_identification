#include "dynamicwidget.h"
#include "ui_dynamicwidget.h"
#include <QFile>
#include <QFileDialog>
#include <QFileInfo>
#include <QDebug>
#include <QDateTime>
#include <QDataStream>
#include <QProcess>

#ifdef __cplusplus
extern "C"
{
#endif

    int algo();

#ifdef __cplusplus
}
#endif

DynamicWidget::DynamicWidget(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::DynamicWidget)
{
    ui->setupUi(this);
    ui->Button_Results->setEnabled(false);
}

DynamicWidget::~DynamicWidget()
{
    delete ui;
}

void DynamicWidget::on_pushButton_clicked()
{
    // dynamic();
    algo();
    qDebug()<< "Algorithm is ok!"<<endl;

    a_data_read();
    alpha_data_read();
    d_data_read();
    // result_read();
    // tau_pos_read();
    phi_read();
    error_read();
    ui->Button_Results->setEnabled(true);


}
void DynamicWidget::a_data_read()
{
    QString filename = "D-H data\\a_data.txt";
    QFile file(filename);
    bool isOk = file.open(QIODevice::ReadOnly);
    if(isOk)
    {
        qDebug() << "File a_data is ok";
        QString a_data1, a_data2, a_data3, a_data4, a_data5, a_data6;
        a_data1 = file.readLine();
        a_data2 = file.readLine();
        a_data3 = file.readLine();
        a_data4 = file.readLine();
        a_data5 = file.readLine();
        a_data6 = file.readLine();
        ui->a_data1->setText(a_data1);
        ui->a_data2->setText(a_data2);
        ui->a_data3->setText(a_data3);
        ui->a_data4->setText(a_data4);
        ui->a_data5->setText(a_data5);
        ui->a_data6->setText(a_data6);
    }
    file.close();
}
void DynamicWidget::alpha_data_read()
{
    QString filename = "D-H data\\alpha_data.txt";

    QFile file(filename);
    bool isOk = file.open(QIODevice::ReadOnly);
    if(isOk)
    {
        qDebug() << "File theta_dot_data is ok";
        QString alpha_data1, alpha_data2, alpha_data3, alpha_data4, alpha_data5, alpha_data6;
        alpha_data1 = file.readLine();
        alpha_data2 = file.readLine();
        alpha_data3 = file.readLine();
        alpha_data4 = file.readLine();
        alpha_data5 = file.readLine();
        alpha_data6 = file.readLine();
        ui->alpha_data1->setText(alpha_data1);
        ui->alpha_data2->setText(alpha_data2);
        ui->alpha_data3->setText(alpha_data3);
        ui->alpha_data4->setText(alpha_data4);
        ui->alpha_data5->setText(alpha_data5);
        ui->alpha_data6->setText(alpha_data6);
     }
    file.close();
}

void DynamicWidget::d_data_read()
{
    QString filename = "D-H data\\d_data.txt";

    QFile file(filename);
    bool isOk = file.open(QIODevice::ReadOnly);
    if(isOk)
    {
        qDebug() << "File d_data is ok";
        QString d_data1, d_data2, d_data3, d_data4, d_data5, d_data6;
        d_data1 = file.readLine();
        d_data2 = file.readLine();
        d_data3 = file.readLine();
        d_data4 = file.readLine();
        d_data5 = file.readLine();
        d_data6 = file.readLine();
        ui->d_data1->setText(d_data1);
        ui->d_data2->setText(d_data2);
        ui->d_data3->setText(d_data3);
        ui->d_data4->setText(d_data4);
        ui->d_data5->setText(d_data5);
        ui->d_data6->setText(d_data6);
    }
    file.close();
}
#if 0
void DynamicWidget::result_read()
{
    QString filename = "result.txt";
    QFile file(filename);

    bool isOk = file.open(QIODevice::ReadOnly);
    if(isOk)
    {
        qDebug() << "File result is ok";
        QByteArray rank, cond, rms, rel;
        while (file.atEnd() == false)
        {
            rank = file.readLine();
            cond = file.readLine();
            rms = file.readLine();
            rel = file.readLine();

            ui->rank->setText(rank);
            ui->cond->setText(cond);
            ui->rms->setText(rms);
            ui->rel->setText(rel);
        }

    }
    file.close();
}

void DynamicWidget::tau_pos_read()
{
    QString filename = "tau_pos.txt";
    QFile file(filename);
    bool isOk = file.open(QIODevice::ReadOnly);
    if(isOk)
    {
        qDebug() << "File tau_pos is ok";
        QByteArray tau_pos1, tau_pos2, tau_pos3, tau_pos4, tau_pos5, tau_pos6;
        while (file.atEnd() == false)
        {
            tau_pos1 = file.readLine();
            tau_pos2 = file.readLine();
            tau_pos3 = file.readLine();
            tau_pos4 = file.readLine();
            tau_pos5 = file.readLine();
            tau_pos6 = file.readLine();
            ui->tau_pos1->setText(tau_pos1);
            ui->tau_pos2->setText(tau_pos2);
            ui->tau_pos3->setText(tau_pos3);
            ui->tau_pos4->setText(tau_pos4);
            ui->tau_pos5->setText(tau_pos5);
            ui->tau_pos6->setText(tau_pos6);
        }

    }
    file.close();
}
#endif

void DynamicWidget::phi_read(){
    QString filename = "results/phi.txt";
    QFile file(filename);
    bool isOk = file.open(QIODevice::ReadOnly);
    if(isOk)
    {
        qDebug() << "File phi.txt is ok";
        QString phi_data = file.readAll();
        ui->phi_data->setText(phi_data);
    }
}
void DynamicWidget::error_read()
{
    QString filename = "results/error.txt";
    QFile file(filename);
    bool isOk = file.open(QIODevice::ReadOnly);
    if(isOk)
    {
        qDebug() << "File error.txt is ok";
        QString err_data1 = file.readLine();
        QString err_data2 = file.readLine();
        QString err_data3 = file.readLine();
        QString err_data4 = file.readLine();
        QString err_data5 = file.readLine();
        QString err_data6 = file.readLine();
        QString err_data7 = file.readLine();
        QString err_data8 = file.readLine();
        QString err_data9 = file.readLine();
        QString err_data10 = file.readLine();
        QString err_data11 = file.readLine();
        QString err_data12 = file.readLine();

        ui->acc1->setText(err_data1);
        ui->acc2->setText(err_data2);
        ui->acc3->setText(err_data3);
        ui->acc4->setText(err_data4);
        ui->acc5->setText(err_data5);
        ui->acc6->setText(err_data6);
        ui->acc_f1->setText(err_data7);
        ui->acc_f2->setText(err_data8);
        ui->acc_f3->setText(err_data9);
        ui->acc_f4->setText(err_data10);
        ui->acc_f5->setText(err_data11);
        ui->acc_f6->setText(err_data12);

    }
}

void DynamicWidget::on_pushButton_2_clicked()
{
    this->close();
}

void DynamicWidget::on_Button_Results_clicked()
{

    QProcess process(this);
    process.startDetached("draw_tau/draw_tau.exe");

}


