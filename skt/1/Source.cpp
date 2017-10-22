#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>

#define vector_double std::vector<double>
#define M 18
#define x0 0
#define xN 1
#define u1 0
#define u2 2



vector_double progon(vector_double &a, vector_double &b, vector_double &c, vector_double &d, const int n)
{
	vector_double u(n + 1, 0);
	vector_double p(n + 1, 0), q(n + 1, 0);
	for (int i = 1; i <= n; i++)
	{
		p[i] = -b[i] / (c[i] * p[i - 1] + a[i]);
		q[i] = (d[i] - c[i] * q[i - 1]) / (c[i] * p[i - 1] + a[i]);
	}
	u[n] = (-c[n] * q[n] + d[n]) / (a[n] + c[n] * p[n]);
	for (int i = n-1; i >= 1; i--)
	{
		u[i] = p[i] * u[i + 1] + q[i];
	}
	return u;
}

int main()
{
	const int N = pow(30 - M, log(10 + M)) / 500;
	double dx = (xN - x0) / double(N);
	vector_double a(N + 1, -2), b(N + 1, 1), c(N + 1, 1), d(N + 1, 0);
	c[1] = 0;
	b[N] = 0;
	a[1] = 1;
	a[N] = 1;
	b[1] = 0;
	c[N] = 0;
	d[1] = u1;
	d[N] = u2;
	vector_double u = progon(a, b, c, d, N);
	std::ofstream out("ans.txt");
	for (int i = 1; i <= N; i++)
	{
		out << x0+(i-1)*dx << " " << u[i] << std::endl;
	}
	out.close();
	system("gnuplot ans.plt");
	return 0;
}