#include "QtGuiApplication1.h"
#include <QtWidgets/QApplication>

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);
	QtGuiApplication1 w;
	w.setWindowTitle(QStringLiteral("�󷽽���"));
	w.show();
	return a.exec();
}
