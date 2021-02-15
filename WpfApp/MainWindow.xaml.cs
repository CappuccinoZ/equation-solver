using System;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Input;

namespace Solver
{
    public class Complex
    {
        public double r, i;
        public Complex(double real, double imag)
        {
            r = real;
            i = imag;
        }
        public static implicit operator Complex(double a) => new Complex(a, 0);
        public double Norm() => r * r + i * i;
        public double Abs() => Math.Sqrt(Norm());
        public double Arg() => Math.Atan2(i, r);
        public Complex Conj() => new Complex(r, -i);
        public static Complex operator +(Complex left, double right) => new Complex(left.r + right, left.i);
        public static Complex operator -(Complex left, double right) => new Complex(left.r - right, left.i);
        public static Complex operator -(Complex right) => new Complex(-right.r, -right.i);
        public static Complex operator *(Complex left, double right) => new Complex(left.r * right, left.i * right);
        public static Complex operator /(Complex left, double right) => new Complex(left.r / right, left.i / right);
        public static Complex operator +(Complex left, Complex right) => new Complex(left.r + right.r, left.i + right.i);
        public static Complex operator -(Complex left, Complex right) => new Complex(left.r - right.r, left.i - right.i);
        public static Complex operator *(Complex left, Complex right) => new Complex(left.r * right.r - left.i * right.i, left.r * right.i + left.i * right.r);
        public static Complex operator /(Complex left, Complex right) => left * right.Conj() / right.Norm();
        public Complex Log() => new Complex(Math.Log(Abs()), Arg());
        public Complex Exp() => new Complex(Math.Cos(i), Math.Sin(i)) * Math.Exp(r);
        public static Complex operator ^(Complex left, Complex right) => (left.Log() * right).Exp();
        public static Complex operator ^(Complex left, double right) => new Complex(Math.Cos(left.Arg() * right), Math.Sin(left.Arg() * right)) * Math.Pow(left.Abs(), right);
        public Complex Sqrt() => this ^ 0.5;
        public Complex Cbrt() => (i == 0 && r < 0) ? new Complex(-Math.Pow(-r, 1.0 / 3), 0) : this ^ (1.0 / 3);
        public override string ToString()
        {
            double a = Math.Round(r, 12), b = Math.Round(i, 12);
            if (a == 0 && b == 0) return "0";
            else
            {
                string str = "";
                if (a != 0)
                {
                    str += a;
                    if (b > 0) str += "+";
                }
                switch (b)
                {
                    case 0: break;
                    case -1: str += "-i"; break;
                    case 1: str += "i"; break;
                    default: str += b + "i"; break;
                }
                return str;
            }
        }
    }

    /// <summary>
    /// MainWindow.xaml 的交互逻辑
    /// </summary>
    public partial class MainWindow : Window
    {
        public bool Equal(double a, double b) => Math.Abs(a - b) <= 2.2204460492503131E-16;//判断是否约等
        public double Cbrt(double x) => (x < 0) ? -Math.Pow(-x, 1.0 / 3) : Math.Pow(x, 1.0 / 3);//立方根
        public bool Isdouble(string target) => double.TryParse(target, out _);//判断String是否可以转换成Double
        public void Show<T>(T str) => textBox7.Text += str.ToString() + "\r\n";

        public void Fun2(Complex a, Complex b, Complex c)//ax^2+bx+c=0
        {
            Complex t = (b * b - 4 * a * c).Sqrt();
            Show((-t - b) / (2 * a));
            Show((t - b) / (2 * a));
        }

        public double Fun3_one(double a, double b, double c, double d)//返回ax^3+bx^2+cx+d=0的一个实数根
        {
            double k = 1 / a;
            a = b * k;
            b = c * k;
            c = d * k;
            double e = (3 * b - a * a) / 9, f = ((2 * a * a - 9 * b) * a / 27 + c) / 2;
            Complex t = new Complex(f * f + e * e * e, 0).Sqrt();
            return ((t - f).Cbrt() + (-t - f).Cbrt()).r - a / 3;
        }

        public void Fun3(double a, double b, double c, double d)//ax^3+bx^2+cx+d=0
        {
            double x = Equal(d, 0) ? 0 : Fun3_one(a, b, c, d);
            Show(x);
            Fun2(a, a * x + b, (a * x + b) * x + c);
        }

        public void Fun4(double a, double b, double c, double d, double e)//ax^4+bx^3+cx^2+dx+e=0
        {
            if (Equal(e, 0))//x(ax^3+bx^2+cx+d)=0
            {
                Show(0);
                Fun3(a, b, c, d);
            }
            else
            {
                double k = 1 / a;
                a = b * k;
                b = c * k;
                c = d * k;
                d = e * k;
                if (Equal(a, 0) && Equal(c, 0))
                {
                    Complex t = new Complex(b * b - 4 * d, 0).Sqrt();
                    Fun2(1, 0, (b + t) / 2);
                    Fun2(1, 0, (b - t) / 2);
                }
                else
                {
                    double y = Fun3_one(1, -b, a * c - 4 * d, 4 * b * d - a * a * d - c * c);
                    double p = Math.Sqrt(a * a / 4 - b + y);
                    if (Equal(p, 0))
                    {
                        Complex t = new Complex(y * y / 4 - d, 0).Sqrt();
                        Fun2(1, a / 2, y / 2 + t);
                        Fun2(1, a / 2, y / 2 - t);
                    }
                    else
                    {
                        double q = (a * y - 2 * c) / (4 * p * p);
                        Fun2(1, a / 2 - p, y / 2 - p * q);
                        Fun2(1, a / 2 + p, y / 2 + p * q);
                    }
                }
            }
        }

        double[] root = { 0, 0, 0, 0 };
        public void Fun4_realroot(double a, double b, double c, double d, double e)//四次方程实根数量
        {
            double t, y, p, q, delta;
            t = 1 / a;
            a = b * t;
            b = c * t;
            c = d * t;
            d = e * t;
            y = Fun3_one(1, -b, a * c - 4 * d, 4 * b * d - a * a * d - c * c);
            p = Math.Sqrt(a * a / 4 - b + y);
            if (Equal(p, 0))
            {
                if (y * y - 4 * d >= 0)
                {
                    delta = a * a / 16 - y / 2 + Math.Sqrt(y * y / 4 - d);
                    if (delta >= 0)
                    {
                        t = Math.Sqrt(delta);
                        root[0] = -a / 4 + t;
                        root[1] = -a / 4 - t;
                    }
                    delta = a * a / 16 - y / 2 - Math.Sqrt(y * y / 4 - d);
                    if (delta >= 0)
                    {
                        t = Math.Sqrt(delta);
                        root[2] = -a / 4 + t;
                        root[3] = -a / 4 - t;
                    }
                }
            }
            else
            {
                q = (a * y - 2 * c) / (4 * p * p);
                delta = 4 * p * (p + 4 * q - a) + a * a - 8 * y;
                if (delta >= 0)
                {
                    t = Math.Sqrt(delta);
                    root[0] = (2 * p - a + t) / 4;
                    root[1] = (2 * p - a - t) / 4;
                }
                delta = 4 * p * (p - 4 * q + a) + a * a - 8 * y;
                if (delta >= 0)
                {
                    t = Math.Sqrt(delta);
                    root[2] = (-2 * p - a + t) / 4;
                    root[3] = (-2 * p - a - t) / 4;
                }
            }
        }

        public double Start(double a, double b, double c, double d, double e)//x^5+ax^4+bx^3+cx^2+dx+e=0根的上界
        {
            double[] list = { a, b, c, d, e };
            Array.Sort(list);
            double q = list[4], y = 0;
            if (list[0] < 0)
            {
                if (a < 0)
                {
                    y = q;
                }
                else if (b < 0)
                {
                    y = Math.Sqrt(q);
                }
                else if (c < 0)
                {
                    y = Cbrt(q);
                }
                else if (d < 0)
                {
                    y = Math.Pow(q, 0.25);
                }
                else if (e < 0)
                {
                    y = Math.Pow(q, 0.2);
                }
            }
            return y;
        }

        public double Fun5_calc(double a, double b, double c, double d, double e, double x) => x * (x * (x * (x * (x + a) + b) + c) + d) + e;//计算函数值x^5+ax^4+bx^3+cx^2+dx+e

        public double Fun5_derivative(double a, double b, double c, double d, double x) => x * (x * (x * (5 * x + 4 * a) + 3 * b) + 2 * c) + d;//计算导数5x^4+4ax^3+3bx^2+2cx+d

        public void Fun5(double a, double b, double c, double d, double e)//x^5+ax^4+bx^3+cx^2+dx+e=0
        {
            if (Equal(e, 0))
            {
                Show(0);
                Fun4(1, a, b, c, d);
            }
            else
            {
                var f = true;
                Fun4_realroot(5, 4 * a, 3 * b, 2 * c, d);
                int i = 0;
                double x = 0;
                while (i < root.Length)
                {
                    if (Equal(Fun5_calc(a, b, c, d, e, root[i]), 0))
                    {
                        f = false;
                        x = root[i];
                        break;
                    }
                    i++;
                }
                if (f)
                {
                    double t = Start(a, b, c, d, e);
                    x = -Start(-a, b, -c, d, -e);
                    if (t != 0) x = Equal(x, 0) ? t : (t + x) / 2;
                    while (Equal(Fun5_derivative(a, b, c, d, x), 0)) x += 0.125;
                    i = 0;
                    do
                    {
                        t = x;
                        x -= Fun5_calc(a, b, c, d, e, x) / Fun5_derivative(a, b, c, d, x);
                        i++;
                    } while (i < 60000000 && Math.Abs(x - t) > 1e-15);
                    Show("迭代次数:" + i);
                    t = Fun5_calc(a, b, c, d, e, x);
                    if (Math.Abs(t) > 1) Show("误差较大");
                    Show("L-R=" + t);
                }
                Show(x);
                Fun4(1, x + a, x * (x + a) + b, x * (x * (x + a) + b) + c, x * (x * (x * (x + a) + b) + c) + d);
            }
        }

        public void Judgement(double a, double b, double c, double d, double e, double f)//判断次数
        {
            if (a != 0)//五次
            {
                textBox7.Text = "x1,x2,x3,x4,x5:\r\n";
                Fun5(b / a, c / a, d / a, e / a, f / a);
            }
            else if (b != 0)//四次
            {
                textBox7.Text = "x1,x2,x3,x4:\r\n";
                Fun4(b, c, d, e, f);
            }
            else if (c != 0)//三次
            {
                textBox7.Text = "x1,x2,x3:\r\n";
                Fun3(c, d, e, f);
            }
            else if (d != 0)//二次
            {
                textBox7.Text = "x1,x2:\r\n";
                Fun2(d, e, f);
            }
            else if (e != 0)//一次
            {
                textBox7.Text = "x = " + (-f / e).ToString();
            }
            else//常值
            {
                textBox7.Text = (f == 0) ? "任意复数" : "无解";
            }
        }

        public MainWindow()
        {
            InitializeComponent();
            textBox1.Focus();
        }

        private void Button1_Click(object sender, RoutedEventArgs e)//求解
        {
            textBox7.Text = "";
            Array.Clear(root, 0, 4);
            if (Isdouble(textBox1.Text) && Isdouble(textBox2.Text) && Isdouble(textBox3.Text) && Isdouble(textBox4.Text) && Isdouble(textBox5.Text) && Isdouble(textBox6.Text))
            {
                if (comboBoxItem1.IsSelected)
                {
                    Judgement(
                        Convert.ToDouble(textBox1.Text),
                        Convert.ToDouble(textBox2.Text),
                        Convert.ToDouble(textBox3.Text),
                        Convert.ToDouble(textBox4.Text),
                        Convert.ToDouble(textBox5.Text),
                        Convert.ToDouble(textBox6.Text)
                        );
                    button2.Focus();
                }
                else if (comboBoxItem2.IsSelected)
                {
                    Judgement(
                            0,
                            Convert.ToDouble(textBox1.Text),
                            Convert.ToDouble(textBox2.Text),
                            Convert.ToDouble(textBox3.Text),
                            Convert.ToDouble(textBox4.Text),
                            Convert.ToDouble(textBox5.Text)
                            );
                }
                else if (comboBoxItem3.IsSelected)
                {
                    Judgement(
                            0,
                            0,
                            Convert.ToDouble(textBox1.Text),
                            Convert.ToDouble(textBox2.Text),
                            Convert.ToDouble(textBox3.Text),
                            Convert.ToDouble(textBox4.Text)
                            );
                }
                else
                {
                    Judgement(
                            0,
                            0,
                            0,
                            Convert.ToDouble(textBox1.Text),
                            Convert.ToDouble(textBox2.Text),
                            Convert.ToDouble(textBox3.Text)
                            );
                }
            }
            else
            {
                MessageBox.Show("请检查输入内容(Alt+F4)");
            }
        }

        private void Button2_Click(object sender, RoutedEventArgs e)//重置
        {
            Array.Clear(root, 0, 4);
            textBox1.Text = "0";
            textBox2.Text = "0";
            textBox3.Text = "0";
            textBox4.Text = "0";
            textBox5.Text = "0";
            textBox6.Text = "0";
            textBox7.Text = "";
        }

        private void FocusSet(KeyEventArgs e, TextBox t, bool c = true)
        {
            if (e.Key == Key.Enter || e.Key == Key.Down)
            {
                if (c)
                {
                    t.Focus();
                }
                else
                {
                    button1.Focus();
                }
            }
        }

        private void TextBox1_KeyUp(object sender, KeyEventArgs e) => FocusSet(e, textBox2);
        private void TextBox2_KeyUp(object sender, KeyEventArgs e) => FocusSet(e, textBox3);
        private void TextBox3_KeyUp(object sender, KeyEventArgs e) => FocusSet(e, textBox4, !comboBoxItem4.IsSelected);
        private void TextBox4_KeyUp(object sender, KeyEventArgs e) => FocusSet(e, textBox5, !comboBoxItem3.IsSelected);
        private void TextBox5_KeyUp(object sender, KeyEventArgs e) => FocusSet(e, textBox6, !comboBoxItem2.IsSelected);
        private void TextBox6_KeyUp(object sender, KeyEventArgs e) => FocusSet(e, textBox7, false);

        private void ComboBoxSet(int a)
        {
            label4.Visibility = (a == 2) ? Visibility.Hidden : Visibility.Visible;
            label5.Visibility = (a < 4) ? Visibility.Hidden : Visibility.Visible;
            label6.Visibility = (a < 5) ? Visibility.Hidden : Visibility.Visible;
            textBox4.Visibility = (a == 2) ? Visibility.Hidden : Visibility.Visible;
            textBox5.Visibility = (a < 4) ? Visibility.Hidden : Visibility.Visible;
            textBox6.Visibility = (a < 5) ? Visibility.Hidden : Visibility.Visible;
            button1.Margin = new Thickness(8, 99 + 33 * a, 0, 0);
            button2.Margin = new Thickness(8, 160 + 33 * a, 0, 0);
            textBox7.Margin = new Thickness(58, 71 + 33 * a, 0, 0);
            Height = 300 + 33 * a;
        }

        private void LabelSet(string[] s)
        {
            label7.Content = s[0];
            label1.Content = s[1];
            label2.Content = s[2];
            label3.Content = s[3];
            label4.Content = s[4];
            label5.Content = s[5];
        }

        private void ComboBoxItem1_Selected(object sender, RoutedEventArgs e)
        {
            string[] s = { "ax^5+bx^4+cx^3+dx^2+ex+f=0", "ax^5", "bx^4", "cx^3", "dx^2", "ex" };
            LabelSet(s);
            ComboBoxSet(5);
        }

        private void ComboBoxItem2_Selected(object sender, RoutedEventArgs e)
        {
            string[] s = { "ax^4+bx^3+cx^2+dx+e=0", "ax^4", "bx^3", "cx^2", "dx", "e" };
            LabelSet(s);
            ComboBoxSet(4);
        }

        private void ComboBoxItem3_Selected(object sender, RoutedEventArgs e)
        {
            string[] s = { "ax^3+bx^2+cx+d=0", "ax^3", "bx^2", "cx", "d", "" };
            LabelSet(s);
            ComboBoxSet(3);
        }

        private void ComboBoxItem4_Selected(object sender, RoutedEventArgs e)
        {
            string[] s = { "ax^2+bx+c=0", "ax^2", "bx", "c", "", "" };
            LabelSet(s);
            ComboBoxSet(2);
        }
    }
}
