#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <algorithm>

using namespace std;

typedef std::vector<double> vec;

void operator+=(vec & summ, vec & summand) {
	for (size_t i = 0; i < summ.size(); i++)
		summ[i] += summand[i];
}

void operator-=(vec & summ, vec & summand) {
	for (size_t i = 0; i < summ.size(); i++)
		summ[i] -= summand[i];
}

void operator*=(vec & summ, double mult) {
	for (size_t i = 0; i < summ.size(); i++)
		summ[i] *= mult;
}

void operator<<(vec & summ, double val) {
	for (size_t i = 0; i < summ.size(); i++)
		summ[i] = val;
}

void print(vec & v) {
	for (size_t i = 0; i < v.size(); i++)
		cout << v[i] << " ";
	cout << endl;
}

int nt = 10000;
int ns = 10000;
double lt = 5, rt = 6, ls = 0, rs = M_PI, lz = 2, rz = -2;
double dt, ds, sigma, h, B;
double c_k = 0.02986315562;
double c_u = -0.035714285714;
int nt_n[3] = { 100, 1000, 10000 };
int ns_n[3] = { 100, 1000, 10000 };
double h_n[3] = { 0.005, 0.0000005, 0.000000005 };
//double h_n[3] = { 0.05, 0.05, 0.05 };
//double sigma_n[3] = { 0.05, 0.05, 0.05 };
double sigma_n[3] = { 0.005, 0.0005, 0.00005 };

ofstream res_file;

inline double sign(double val)
{
	if (val < 0)
		return -1;
	return 1;
}

inline double K(double t, double z)
{
	return  (z*z - 2 ) / (2*(t + z));
}

inline double u(double t)
{
	return (M_PI / sqrt(t*t - 4)) * pow((sqrt(t*t - 4) - t) / 2., 2);
}

inline double dK_dz(double t, double z) {
	return (z*t - 2) / pow(t + z, 2);
}

inline 

void update_mesh(int ind_s, int ind_t) {
	nt = nt_n[ind_t];
	ns = ns_n[ind_s];

	sigma = sigma_n[ind_s];
	h = h_n[ind_s];

	dt = (rt - lt) / (nt - 1.0);
	ds = (rs - ls) / (ns - 1.0);

	B = h * (h + 2 * c_k) * (rt - lt);
	B += (sigma * (h + 2 * c_k) + h * (sigma + 2 * c_u)) * (rt - lt) / (rs - ls);
	B += sigma * (sigma + 2 * c_u) / (rs - ls) / (rs - ls) * (rt - lt);
}

double calc_psi(vec & z) {
	double res = 0;
	for (size_t j = 0; j < ns - 1; j++)
		res += (z[j] * z[j] + z[j + 1] * z[j + 1]);
	res = B * (rs - ls) * (rs - ls + 1) * (1 + res * ds * 0.5);
	return res;
}

// int(Kn(t,s,z) ds) - u(t), s = a..b
double calc_inner_int(vec & z, double t)
{
	double prev = K(t, z[0]);
	double integral = 0;
	for (size_t j = 1; j < ns; j++) {
		double curr = K(t, z[j]);
		integral += prev + curr;
		prev = curr;
	}
	return integral * ds * 0.5 - u(t);
}

double calc_integral_KU2(vec & z)
{
	double prev = calc_inner_int(z, lt);
	prev *= prev;
	double res = 0;
	// т.к. уже есть prev_val
	double t = lt + dt;
	for (size_t i = 1; i < nt; i++)
	{
		double curr = calc_inner_int(z, t);
		curr *= curr;
		res += curr + prev;
		prev = curr;
		t += dt;
	}

	return res * dt * 0.5;
}

double calc_rho(vec & z)
{
	printf("Integral (K - U)^2 %.15f\n", calc_integral_KU2(z));
	return calc_integral_KU2(z) - calc_psi(z);
}

void calc_dz_dss(vec & dz_dss, vec & z)
{
	//cout << "Z IN GRAD \n";
	//print(z);
	double dss = ds * ds;
	for (int i = 1; i < dz_dss.size() - 1; i++) {
		dz_dss[i] = (z[i - 1] - 2 * z[i] + z[i + 1]) / dss;
	}
	// Производные на краях роли не играют, так как есть граничные условия
	//dz_dss[0] = 0;
	//dz_dss[dz_dss.size() - 1] = 0;
	//cout << "SUMMAND IN GRAD" << endl;
	//print(dz_dss);
}

void calc_grad(double alpha, vec & z, vec & grad)
{
	grad << 0;

	vec prev(ns);
	vec curr(ns);

	for (size_t i = 0; i < z.size(); i++)
		prev[i] = dK_dz(lt, z[i]);
	prev *= calc_inner_int(z, lt);

	double t = lt + dt;
	for (int i = 1; i < nt; i++) {
		for (size_t j = 0; j < z.size(); j++)
			curr[j] = dK_dz(t, z[j]);
		curr *= calc_inner_int(z, t);

		prev += curr;
		grad += prev;
		prev = curr;

		t += dt;
	}
	grad *= 0.5 * dt;
	//cout << "INT IN GRAD" << endl;
	//print(grad);

	vec dz_dss(ns);
	calc_dz_dss(dz_dss, z);

	dz_dss -= z;
	dz_dss *= alpha;

	grad -= dz_dss;
}

double calc_M(vec & z, double alpha)
{
	double dss = ds * ds;
	double res = calc_integral_KU2(z);

	double prev = z[0] * z[0] + (z[1] - z[0]) * (z[1] - z[0]) / dss;
	double integral = 0;
	for (int j = 1; j < ns; j++) {
		double curr = z[j] * z[j] + (z[j] - z[j - 1]) * (z[j] - z[j - 1]) / dss;
		integral += (prev + curr);
		prev = curr;
	}
	integral *= ds * 0.5 * alpha;

	return res + integral;
}

void grad_method(vec & z, double alpha)
{
	double b = 1;
	vec new_z(ns);
	vec grad(ns);
	vec temp(ns);

	double M_old = calc_M(z, alpha);
	double M_new;
		
	//printf("Starting M %.10f\n", calc_M(z, alpha));
	int iter = 0;
	do {
		calc_grad(alpha, z, grad);
		//cout << "GRAD" << endl;
		//print(grad);
		do
		{
			new_z = grad;
			new_z *= -b;
			new_z += z;
			new_z[0] = lz;
			new_z[ns - 1] = rz;
			//cout << "NEW_Z" << endl;
			//print(new_z);
			M_new = calc_M(new_z, alpha);
			//if(iter % 100 == 0)
				printf("b = %.15f  M = %.15f\n", b, M_new);
			if (M_old > M_new) 
			{
				b *= 2;
				break;
			}
			else
				b /= 2;
		} while (true);//(abs(M_old - M_new) > 0.00000001);

		z = new_z;

		//cout << "New M " << M_new << endl;

		if (abs(M_old - M_new) < 0.00000001 && iter > 100) 
		{
			M_old = M_new;
			break;
		}
		else
			M_old = M_new;

		iter++;
	} while (true);
	printf("Iter num %i\n", iter);

	//printf("Final M %.10f\n", calc_M(z, alpha));
}

void calc(string filename)
{
	vec grad(ns);
	vec z(ns);

	for (int i = 0; i < ns; i++)
		z[i] = cos(ls + ds * i);

	double alpha = 0;
	double l_alpha = 0;
	double r_alpha = 250;
	double m_alpha = -1;
	double l_rho = 1000, r_rho = 1000, m_rho;

	grad_method(z, l_alpha);
	l_rho = calc_rho(z);
	for (int i = 0; i < ns; i++)
	z[i] = cos(ls + ds * i);

	grad_method(z, r_alpha);
	r_rho = calc_rho(z);
	for (int i = 0; i < ns; i++)
	z[i] = cos(ls + ds * i);

	vec save_z(ns);

	ofstream file;
	file.open("init.txt");
	for (int i = 0; i < ns; i++)
		file << z[i] << endl;
	file.close();

	while (abs(r_rho) >= 0.00001 && abs(l_rho) >= 0.00001)
	{
		save_z = z;

		if (l_alpha != m_alpha)
		{
			grad_method(z, l_alpha);
			l_rho = calc_rho(z);
			z = save_z;
		} 

		if (r_alpha != m_alpha)
		{
			grad_method(z, r_alpha);
			r_rho = calc_rho(z);
			z = save_z;
		}

		m_alpha = (r_alpha + l_alpha) * 0.5;
		grad_method(z, m_alpha);
		m_rho = calc_rho(z);

		printf("Rho [%.10f;  %.10f;  %.10f]\n", l_rho, m_rho, r_rho);
		printf("Alpha [%.10f;  %.10f;  %.10f]\n", l_alpha, m_alpha, r_alpha);

		if (r_rho < 0)
			break;
		if (l_rho > 0)
			break;

		if (sign(m_rho) != sign(r_rho))
		{
			l_rho = m_rho;
			l_alpha = m_alpha;
		}
		else
		{
			r_rho = m_rho;
			r_alpha = m_alpha;
		}
	}

	file.open(filename);
	for (int i = 0; i < ns; i++)
		file << z[i] << endl;
	file.close();
}

double test_result(vec & res) {
	double integral = 0;
	double max_val = 0;
	for (int i = 0; i < nt; i++) {
		double val = calc_inner_int(res, i * dt);
		if (abs(val) > max_val)
			max_val = abs(val);
	}
	return max_val;
}

void main()
{
	res_file.open("res2.txt");

	update_mesh(0, 0);
	vec z(ns);
	vec z_init(ns);
	for (int i = 0; i < ns; i++)
		//z_init[i] = 0.5;
		//z_init[i] = i / (double)(ns);
		z_init[i] = 2 * cos(i * ds);
		// z_init[i] = sin(i * ds);

	z_init[0] = ls;
	z_init[ns - 1] = rs;

	z = z_init;
	// cout << calc_rho(z) << endl;

	srand(1337);

	/*grad_method(z, -1e-8);
	double rho = calc_rho(z);
	printf("Rho %lf\n", rho);
	return;*/

	grad_method(z, 0);
	cout << calc_rho(z) << endl;
	for (int i = 0; i < ns; i += 1)
		res_file << i * ds << " " << z[i] << endl;
	res_file << endl;

	res_file.open("init.txt");

	for (int i = 0; i < ns; i++) {
		res_file << i * ds << " " << 2 * cos(i * ds) << endl;
	}

	calc("res0.txt");

	return;

	double prev_rho = calc_rho(z);
	for (int i = 1; i < 50; i+=1)
	{
		// z = z_init;
		grad_method(z, i * 0.5);
		double rho = calc_rho(z);
		cout << i * 0.5 << "\t" << calc_rho(z) << endl;
		//if (abs(rho - prev_rho) < 0.000000000001)
			//break;
		prev_rho = rho;
	}

	/*double l_rho = 0;
	double r_rho = 5;
	for (int i = 0; i < 3; i++) {
		if (i % 4 == 0 && i != 0) {
			l_rho += 5;
			r_rho += 5;
			if (r_rho == 20)
				break;
		}
		grad_method(z, l_rho);
		double rho = calc_rho(z);
		cout << l_rho << "\t" << calc_rho(z) << endl;

		grad_method(z, r_rho);
		rho = calc_rho(z);
		cout << r_rho << "\t" << calc_rho(z) << endl;
	}*/

	for (int i = 0; i < ns; i += 1)
		res_file << i * ds << " " << z[i] << endl;
	res_file << endl;

	//update_mesh(100, 1);
	//calc("res100.txt");
	//update_mesh(1000, 1);
	//calc("res1000.txt");
	res_file.close();

	cout << test_result(z);
	res_file.open("init.txt");

	for (int i = 0; i < ns; i++) {
		res_file << i * ds << " " << 2*cos(i * ds) << endl;
	}
}

