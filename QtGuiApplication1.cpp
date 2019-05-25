#include "QtGuiApplication1.h"
#include "ui_QtGuiApplication1.h"
#include <QFileDialog>
#include <QString>
#include <QTextstream>
#include <QFile>
#include <QMessagebox>
#include <Matrix.h>
#include <cmath>
#define pi 3.1415926
const int N = 6;
void Calculate_R(double fa, double w, double ka, double* a, double* b, double* c);   //旋转矩阵
void Iteration_Equation(double w, double ka, double line_num, double (*A)[N], double *l, double *a, double *b, double *c, double f, double *x, double *y, double XS0, double YS0, double ZS0, double *X, double *Y, double *Z);   //外方位元素迭代函数
int i, j, k;    //计数器
QtGuiApplication1::QtGuiApplication1(QWidget *parent)
	: QMainWindow(parent)
{
	ui.setupUi(this);
	connect(this->ui.actionOpen_File, SIGNAL(triggered()),this,SLOT(on_actionOpenFile()));
	connect(this->ui.calculate, SIGNAL(clicked()), this, SLOT(on_calculate()));
	connect(this->ui.cout, SIGNAL(clicked()), this, SLOT(on_cout()));
}

void QtGuiApplication1::on_actionOpenFile()
{
	QString filename = QFileDialog::getOpenFileName(this,QStringLiteral("打开数据文件"),"",tr("*.txt"));
	QFile file(filename);
	if (!file.open(QFile::ReadOnly | QFile::Text))
	{
		QMessageBox::warning(this, QStringLiteral("警告"), QStringLiteral("文件未打开"));
	}
	else
	{
		QTextStream in(&file);
		QString ss;
		in >> ss;
		ui.PDOC->setText(ss);  //摄影机主距
		in >> ss;
		ui.scale->setText(ss);  //摄影比例尺
		ss = in.readAll().trimmed();
		ui.textEdit->setPlainText(ss);  //坐标数据
		file.close();
	}
}
void QtGuiApplication1::on_calculate()
{

	if (ui.textEdit->document()->lineCount() <= 1)
	{
		QMessageBox::information(this, QStringLiteral("提示"), QStringLiteral("请选择文件"));
	}
	else
	{
		int line_num = ui.textEdit->document()->lineCount();   //统计坐标个数
		QString ss = ui.textEdit->toPlainText();
		QTextStream s(&ss);
		double m, f;   //摄影比例尺和摄影主距
		QString temp;
		temp = ui.scale->text();
		m = temp.toDouble();
		temp = ui.PDOC->text();
		f = temp.toDouble() / 1000;  //单位变成米
		double *x = new double[line_num], *y = new double[line_num], *X = new double[line_num], *Y = new double[line_num], *Z = new double[line_num];  //开辟空间
		//赋值
		for (i = 0;i < line_num;i++)
		{
			s >> x[i] >> y[i] >> X[i] >> Y[i] >> Z[i];
		}
		//外方位元素赋初值
		double XS0 = 0, YS0 = 0, ZS0 = 0, fa = 0, w = 0, ka = 0;
		for (i = 0;i < 4;i++)
		{
			XS0 += X[i];
			YS0 += Y[i];
			x[i] /= 1000;
			y[i] /= 1000;
		}
		XS0 /= 4;
		YS0 /= 4;
		ZS0 = m * f;
		double a[3], b[3], c[3];   //旋转矩阵
		double x0, y0, z0;
		double(*A)[N];
		A = new double[2 * line_num][N];
		double *l = new double[2 * line_num];
		Matrix V, B(2 * line_num, N), x1, ll(2 * line_num, 1);  //间接平差矩阵
		for (k = 0;k < 10;k++)
		{
			Calculate_R(fa, w, ka, a, b, c);
			Iteration_Equation(w, ka, line_num, A, l, a, b, c, f, x, y, XS0, YS0, ZS0, X, Y, Z);
			//矩阵初始化
			for (i = 0;i < 2 * line_num;i++)
			{
				for (j = 0;j < N; j++)
				{
					B.set(i, j, A[i][j]);
				}
				ll.set(i, 0, l[i]);
			}
			x1 = (B.Trans()*B).Inverse()*B.Trans()*ll;
			XS0 = XS0 + x1.get(0, 0);
			YS0 = YS0 + x1.get(1, 0);
			ZS0 = ZS0 + x1.get(2, 0);
			fa = fa + x1.get(3, 0);
			w = w + x1.get(4, 0);
			ka = ka + x1.get(5, 0);
			if (abs(x1.get(3, 0)) < 0.1 / 60 / 180 * pi|| abs(x1.get(4, 0)) < 0.1 / 60 / 180 * pi|| abs(x1.get(5, 0)) < 0.1 / 60 / 180 * pi) break;
		}
		V = B * x1 - ll;   //中误差
		ss = "XS0: " + QString::number(XS0) + '\n' + "YS0: " + QString::number(YS0) + '\n' + "ZS0: " + QString::number(ZS0) + '\n' + "fa: " + QString::number(fa) + '\n' + "w: " + QString::number(w) + '\n' + "ka: " + QString::number(ka) + '\n';
		//ss = QString::number(B.get(4,4));
		double o = sqrt((V.Trans()*V).get(0, 0) / (2 * line_num - 6));
		ui.mid_error->setText(QString::number(o));
		ui.iteration_num->setText(QString::number(k + 1));
		ui.textEdit_1->setPlainText(ss);
	}
}
//输出至文件
void QtGuiApplication1::on_cout()
{
	QString ss = ui.textEdit_1->toPlainText();
	QString filename = ui.outName->text();
	QFile File(filename);
	if (filename == "")
	{
		filename = "out.txt";
	}
	else
	{
		filename += ".txt";
	}
	QFile file(filename);
	if (!file.open(QIODevice::WriteOnly | QIODevice::Text | QIODevice::Truncate))
	{
		QMessageBox::information(this, QStringLiteral("提示"), QStringLiteral("文件未自动创建并打开"), QMessageBox::Ok);
		return;
	}
	else
	{
		QTextStream out(&file);
		out << ss;
		file.close();
	}
}

void Calculate_R(double fa, double w, double ka, double* a, double* b, double* c)
{
	a[0] = cos(fa)*cos(ka) - sin(fa)*sin(w)*sin(ka);
	a[1] = -cos(fa)*sin(ka) - sin(fa)*sin(w)*cos(ka);
	a[2] = -sin(fa)*cos(w);
	b[0] = cos(w)*sin(ka);
	b[1] = cos(w)*cos(ka);
	b[2] = -sin(w);
	c[0] = sin(fa)*cos(ka) + cos(fa)*sin(w)*sin(ka);
	c[1] = -sin(fa)*sin(ka) + cos(fa)*sin(w)*cos(ka);
	c[2] = cos(fa)*cos(w);
}

void Iteration_Equation(double w, double ka,double line_num, double (*A)[N], double *l, double *a, double *b, double *c, double f, double *x, double *y, double XS0, double YS0, double ZS0,double *X, double *Y, double *Z)
{
	double x0, y0, z0;
	for (int i = 0; i < line_num; i++)
	{
		x0 = -f * (a[0] * (X[i] - XS0) + b[0] * (Y[i] - YS0) + c[0] * (Z[i] - ZS0)) / (a[2] * (X[i] - XS0) + b[2] * (Y[i] - YS0) + c[2] * (Z[i] - ZS0));
		y0 = -f * (a[1] * (X[i] - XS0) + b[1] * (Y[i] - YS0) + c[1] * (Z[i] - ZS0)) / (a[2] * (X[i] - XS0) + b[2] * (Y[i] - YS0) + c[2] * (Z[i] - ZS0));
		z0 = a[2] * (X[i] - XS0) + b[2] * (Y[i] - YS0) + c[2] * (Z[i] - ZS0);
		A[2 * i][0] = (a[0] * f + a[2] * x[i]) / z0;
		A[2 * i][1] = (b[0] * f + b[2] * x[i]) / z0;
		A[2 * i][2] = (c[0] * f + c[2] * x[i]) / z0;
		A[2 * i][3] = y[i] * sin(w) - (x[i] * (x[i] * cos(ka) - y[i] * sin(ka)) / f + f * cos(ka))*cos(w);
		A[2 * i][4] = -f * sin(ka) - x[i] * (x[i] * sin(ka) + y[i] * cos(ka)) / f;
		A[2 * i][5] = y[i];
		A[2 * i + 1][0] = (a[1] * f + a[2] * y[i]) / z0;
		A[2 * i + 1][1] = (b[1] * f + b[2] * y[i]) / z0;
		A[2 * i + 1][2] = (c[1] * f + c[2] * y[i]) / z0;
		A[2 * i + 1][3] = -x[i] * sin(w) - (y[i] * (x[i] * cos(ka) - y[i] * sin(ka)) / f - f * sin(ka))*cos(w);
		A[2 * i + 1][4] = -f * cos(ka) - y[i] * (x[i] * sin(ka) + y[i] * cos(ka)) / f;
		A[2 * i + 1][5] = -x[i];
		l[2 * i] = x[i] - x0;
		l[2 * i + 1] = y[i] - y0;
	}
}
