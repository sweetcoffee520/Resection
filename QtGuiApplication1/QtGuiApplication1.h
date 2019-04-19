#pragma once

#include <QtWidgets/QMainWindow>
#include "ui_QtGuiApplication1.h"

class QtGuiApplication1 : public QMainWindow
{
	Q_OBJECT

public:
	QtGuiApplication1(QWidget *parent = Q_NULLPTR);

private slots:
	void on_actionOpenFile();
	void on_calculate();
	void on_cout();

private:
	Ui::QtGuiApplication1Class ui;
};
