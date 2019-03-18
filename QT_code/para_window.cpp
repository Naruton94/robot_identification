#include "para_window.h"
#include "ui_para_window.h"
#include <QFile>
#include <QDebug>

para_window::para_window(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::para_window)
{
    ui->setupUi(this);
}

para_window::~para_window()
{
    delete ui;
}
void para_window:: parameters_read(QString filename)
{


    QFile file(filename);
    bool isOk = file.open(QIODevice::ReadOnly);
    if(isOk)
    {
        qDebug() << "File parameters is ok";
        QByteArray parameter1, parameter2, parameter3, parameter4, parameter5,
                parameter6, parameter7, parameter8, parameter9, parameter10, parameter11, parameter12, parameter13;
        while (file.atEnd() == false)
        {
            parameter1 = file.readLine();
            parameter2 = file.readLine();
            parameter3 = file.readLine();
            parameter4 = file.readLine();
            parameter5 = file.readLine();
            parameter6 = file.readLine();
            parameter7 = file.readLine();
            parameter8 = file.readLine();
            parameter9 = file.readLine();
            parameter10 = file.readLine();
            parameter11 = file.readLine();
            parameter12 = file.readLine();
            parameter13 = file.readLine();

            ui->parameter1->setText(parameter1);
            ui->parameter2->setText(parameter2);
            ui->parameter3->setText(parameter3);
            ui->parameter4->setText(parameter4);
            ui->parameter5->setText(parameter5);
            ui->parameter6->setText(parameter6);
            ui->parameter7->setText(parameter7);
            ui->parameter8->setText(parameter8);
            ui->parameter9->setText(parameter9);
            ui->parameter10->setText(parameter10);
            ui->parameter11->setText(parameter11);
            ui->parameter12->setText(parameter12);
            ui->parameter13->setText(parameter13);

        }

    }
    file.close();
}
void para_window::on_J1Button_clicked()
{
    QString filename = "parameters.txt";
    parameters_read(filename);
}
