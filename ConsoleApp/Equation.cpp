#include <iostream>
#include <vector>
#include <complex>
#include <iomanip>
#include <algorithm>
#include <cstdlib>
using std::cin;
using std::cout;
using std::endl;
using std::vector;
using std::complex;

vector<double>root;

auto approx = [](double a, double b) {return fabs(a - b) < DBL_EPSILON; };
auto cbrtc = [](complex<double> a) {return (a.real() < 0 && a.imag() == 0) ? -pow(-a, 1.0 / 3.0) : pow(a, 1.0 / 3.0); };

void display(complex<double> a)
{
	if (a == 0.0) {
		cout << 0;
	}
	else {
		if (a.real() != 0) {
			cout << a.real();
			if (a.imag() > 0) {
				cout << "+";
			}
		}
		if (a.imag() != 0) {
			if (a.imag() == 1) {
				cout << "i";
			}
			else if (a.imag() == -1) {
				cout << "-i";
			}
			else {
				cout << a.imag() << "i";
			}
		}
	}
	cout << endl;
}

void fun2(double a, double b, complex<double> c) {//ax^2+bx+c=0
	complex <double> t = sqrt((complex<double>)(b * b - 4 * a * c));
	complex <double> x1 = (-b + t) / (2 * a);
	complex <double> x2 = (-b - t) / (2 * a);
	display(x1);
	display(x2);
}

double fun3_subsidiary(double a, double b, double c, double d) {//返回ax^3+bx^2+cx+d=0的一个实数根
	double k = 1 / a;
	a = b * k;
	b = c * k;
	c = d * k;
	double p3 = (3 * b - a * a) / 9;
	double q2 = ((2 * a * a - 9 * b) * a / 27 + c) / 2;
	complex <double> t = sqrt((complex<double>)q2 * q2 + p3 * p3 * p3);
	complex <double> x = cbrtc(-q2 + t) + cbrtc(-q2 - t) - a / 3;//卡尔达诺公式
	return x.real();
}

void fun3(double a, double b, double c, double d) {//ax^3+bx^2+cx+d=0
	double x = (d == 0) ? 0 : fun3_subsidiary(a, b, c, d);
	display(x);
	fun2(a, a * x + b, (a * x + b) * x + c);
}

void fun4(double a, double b, double c, double d, double e) {//ax^4+bx^3+cx^2+dx+e=0	
	double y, p, q;
	if (approx(e, 0)) {
		cout << 0 << endl;
		fun3(a, b, c, d);
	}
	else {
		double k = 1 / a;
		a = b * k;
		b = c * k;
		c = d * k;
		d = e * k;
		y = fun3_subsidiary(1, -b, a * c - 4 * d, 4 * b * d - a * a * d - c * c);//费拉里法
		p = sqrt(a * a / 4 - b + y);
		if (approx(p, 0)) {
			fun2(1, a / 2, y / 2 + sqrt((complex<double>)(y * y / 4 - d)));
			fun2(1, a / 2, y / 2 - sqrt((complex<double>)(y * y / 4 - d)));
		}
		else {
			q = (a * y - 2 * c) / (4 * p * p);
			fun2(1, a / 2 - p, y / 2 - p * q);
			fun2(1, a / 2 + p, y / 2 + p * q);
		}
	}
}

void fun4_realroot(double a, double b, double c, double d, double e) {//四次方程实根
	double y, p, q, t, delta;
	if (approx(e, 0)) {
		root.push_back(0);
		root.push_back(fun3_subsidiary(a, b, c, d));
		delta = b * b - 4 * a * c - a * root[1] * (3 * a * root[1] + 2 * b);
		if (delta >= 0) {
			t = sqrt(delta);
			root.push_back(-(root[1] + (b - t) / a) / 2);
			root.push_back(-(root[1] + (b + t) / a) / 2);
		}
	}
	else {
		t = 1 / a;
		a = b * t;
		b = c * t;
		c = d * t;
		d = e * t;
		y = fun3_subsidiary(1, -b, a * c - 4 * d, 4 * b * d - a * a * d - c * c);//费拉里法
		p = sqrt(a * a / 4 - b + y);
		if (approx(p, 0)) {
			if (y * y - 4 * d >= 0) {
				delta = a * a / 16 - y / 2 + sqrt(y * y / 4 - d);
				if (delta >= 0) {
					t = sqrt(delta);
					root.push_back(-a / 4 + t);
					root.push_back(-a / 4 - t);
				}
				delta = a * a / 16 - y / 2 - sqrt(y * y / 4 - d);
				if (delta >= 0) {
					t = sqrt(delta);
					root.push_back(-a / 4 + t);
					root.push_back(-a / 4 - t);
				}
			}
		}
		else {
			q = (a * y - 2 * c) / (4 * p * p);
			delta = 4 * p * (p + 4 * q - a) + a * a - 8 * y;
			if (delta >= 0) {
				t = sqrt(delta);
				root.push_back((2 * p - a + t) / 4);
				root.push_back((2 * p - a - t) / 4);
			}
			delta = 4 * p * (p - 4 * q + a) + a * a - 8 * y;
			if (delta >= 0) {
				t = sqrt(delta);
				root.push_back((-2 * p - a + t) / 4);
				root.push_back((-2 * p - a - t) / 4);
			}
		}
	}
}

double start(double a, double b, double c, double d, double e) {//根的上界
	double q, y;
	double k[5] = { a,b,c,d,e };
	if (a > 0 && b > 0 && c > 0 && d > 0 && e > 0) {
		y = 0;
	}
	else {
		q = *std::max_element(k, k + 5);
		if (a < 0) {
			y = q;
		}
		else if (b < 0) {
			y = sqrt(q);
		}
		else if (c < 0) {
			y = cbrt(q);
		}
		else if (d < 0) {
			y = pow(q, 0.25);
		}
		else {
			y = pow(q, 1.0 / 5.0);
		}
	}
	return y;
}

double fun5_calculation(double a, double b, double c, double d, double e, double x) {//函数值
	return x * (x * (x * (x * (x + a) + b) + c) + d) + e;//x^5+ax^4+bx^3+cx^2+dx+e            
}

double fun5_derivative(double a, double b, double c, double d, double x) {//导数
	return x * (x * (x * (5 * x + 4 * a) + 3 * b) + 2 * c) + d;//5x^4+4ax^3+3bx^2+2cx+d
}

void fun5(double a, double b, double c, double d, double e) {//x^5+ax^4+bx^3+cx^2+dx+e=0
	if (approx(e, 0)) {
		cout << "x1,x2,x3,x4,x5:" << endl << 0 << endl;
		fun4(1, a, b, c, d);
	}
	else {
		unsigned int i;
		double x, t;
		bool tangents = true;//切线法可用		
		fun4_realroot(5, 4 * a, 3 * b, 2 * c, d);
		for (i = 0; i < root.size(); i++) {
			if (approx(fun5_calculation(a, b, c, d, e, root[i]), 0)) {
				tangents = false;
				x = root[i];
				break;
			}
		}
		if (tangents) {
			t = start(a, b, c, d, e);
			x = -start(-a, b, -c, d, -e);
			if (t != 0) {
				x = (x == 0) ? t : (t + x) / 2;
			}
			while (approx(fun5_derivative(a, b, c, d, x), 0)) {
				x += 0.125;
			}
			for (i = 1; i < 100000000; i++) {
				t = x;
				x -= fun5_calculation(a, b, c, d, e, x) / fun5_derivative(a, b, c, d, x);
				if (fabs(x - t) < 1e-15) {
					break;
				}
			}
			cout << "迭代次数:" << i << endl;
			t = fun5_calculation(a, b, c, d, e, x);
			if (fabs(t) > 1) {
				cout << "------误差较大------" << endl;
			}
			cout << "L-R=" << t << endl;
		}
		cout << "x1,x2,x3,x4,x5:" << endl;
		display(x);
		fun4(1, x + a, x * (x + a) + b, x * (x * (x + a) + b) + c, x * (x * (x * (x + a) + b) + c) + d);//降次
	}
}

void menu(int i) {
	vector<double>n(i + 1);
	int k;
	for (k = 0; k <= i; k++) {
		cout << (char)(k + 97) << ": ";
		cin >> n[k];
	}
	switch (i) {
	case(2):cout << "x1,x2:" << endl; fun2(n[0], n[1], n[2]); break;
	case(3):cout << "x1,x2,x3:" << endl; fun3(n[0], n[1], n[2], n[3]); break;
	case(4):cout << "x1,x2,x3,x4:" << endl; ; fun4(n[0], n[1], n[2], n[3], n[4]); break;
	case(5):
		double t = 1 / n[0];
		for (k = 1; k < 6; k++) {
			n[k] *= t;
		}
		fun5(n[1], n[2], n[3], n[4], n[5]);
		break;
	}
}

int main()
{
	int i;

	cout << std::setprecision(12);
	cout << "=============方程计算器============" << endl;
	cout << "2.ax^2+bx+c=0" << endl;
	cout << "3.ax^3+bx^2+cx+d=0" << endl;
	cout << "4.ax^4+bx^3+cx^2+dx+e=0" << endl;
	cout << "5.ax^5+bx^4+cx^3+dx^2+ex+f" << endl;
	cout << "======================By Cappuccino" << endl;
	cout << "请输入方程次数" << endl;
	cin >> i;
	cout << "===================================" << endl;

	if (i < 6 && i > 1) {
		menu(i);
		cout << "===================================" << endl;
		system("pause");
	}

	return 0;
}