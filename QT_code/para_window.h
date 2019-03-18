#ifndef PARA_WINDOW_H
#define PARA_WINDOW_H

#include <QWidget>

namespace Ui {
class para_window;
}

class para_window : public QWidget
{
    Q_OBJECT

public:
    explicit para_window(QWidget *parent = 0);
    ~para_window();
    void parameters_read(QString filename);

private slots:
    void on_J1Button_clicked();

private:
    Ui::para_window *ui;
};

#endif // PARA_WINDOW_H
