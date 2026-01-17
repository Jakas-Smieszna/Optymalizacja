#include"user_funs.h"
#include "matrix.h"
#include <cmath>
#include <cstring>

//LAB0

matrix ff0T(matrix x, matrix ud1, matrix ud2)				// funkcja celu dla przypadku testowego
{
	matrix y;												// y zawiera warto�� funkcji celu
	y = pow(x(0) - ud1(0), 2) + pow(x(1) - ud1(1), 2);		// ud1 zawiera wsp�rz�dne szukanego optimum
	return y;
}

matrix ff0R(matrix x, matrix ud1, matrix ud2)				// funkcja celu dla problemu rzeczywistego
{
	matrix y;												// y zawiera warto�� funkcji celu
	matrix Y0 = matrix(2, 1),								// Y0 zawiera warunki pocz�tkowe
		MT = matrix(2, new double[2] { m2d(x), 0.5 });		// MT zawiera moment si�y dzia�aj�cy na wahad�o oraz czas dzia�ania
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, ud1, MT);	// rozwi�zujemy r�wnanie r�niczkowe
	int n = get_len(Y[0]);									// d�ugo�� rozwi�zania
	double teta_max = Y[1](0, 0);							// szukamy maksymalnego wychylenia wahad�a
	for (int i = 1; i < n; ++i)
		if (teta_max < Y[1](i, 0))
			teta_max = Y[1](i, 0);
	y = abs(teta_max - m2d(ud1));							// warto�� funkcji celu (ud1 to za�o�one maksymalne wychylenie)
	Y[0].~matrix();											// usuwamy z pami�ci rozwi�zanie RR
	Y[1].~matrix();
	return y;
}

matrix df0(double t, matrix Y, matrix ud1, matrix ud2)
{
	matrix dY(2, 1);										// definiujemy wektor pochodnych szukanych funkcji
	double m = 1, l = 0.5, b = 0.5, g = 9.81;				// definiujemy parametry modelu
	double I = m * pow(l, 2);
	dY(0) = Y(1);																// pochodna z po�o�enia to pr�dko��
	dY(1) = ((t <= ud2(1)) * ud2(0) - m * g * l * sin(Y(0)) - b * Y(1)) / I;	// pochodna z pr�dko�ci to przyspieszenie
	return dY;
}

//LAB1

matrix ff1R(matrix x, matrix ud1, matrix ud2)
{
	matrix y;

	double VA0 = 5; // m^3
	double VB0 = 1; // m^3
	double TB0 = 20; // centigrade
	// Y0 zawiera warunku poczatkowe
	matrix Y0 = matrix(3, new double[3] {VA0, VB0, TB0});
	matrix* Y = solve_ode(df1, 0, 1, 2000, Y0, x, ud2);

	int n = get_len(Y[0]);
	double teta_max = Y[1](0, 2);
	for (int i = 1; i < n; ++i)
		if (teta_max < Y[1](i, 2))
			teta_max = Y[1](i, 2);
	y = abs(teta_max - 50);
	Y[0].~matrix();
	Y[1].~matrix();
	return y;
}

// Dla DA = 50cm, temp. maks powinna byc 62 C
// LAB 1
// ud1 - DA, wielkosc otworu zbiornika A

matrix df1(double t, matrix Y, matrix ud1, matrix ud2) {
	// Wektor pochodnych funkcji po czasie; dVA/dt, dVB/dt, dTB/dt

	double DA = m2d(ud1)/(100*100);

	double PA = 2; // m^2
	double TA0 = 95; // centigrade
	double PB = 1; // m^2
	double TBin = 20; // centigrade;
	double FBin = 0.010; // m^3/s
	double DB = 0.003657; // m^2
	double a = 0.98; // wsp. lepkosci cieczy
	double b = 0.63; // wsp. zwezenia strumienia cieczy
	double g = 9.81; // m/s^2, przyspieszenie ziemskie

	double& VA = Y(0);
	double& VB = Y(1);
	double& TB = Y(2);

	auto FAout = a*b*DA*sqrt((2*g*VA)/PA);

	matrix dY(3, 1);
	dY(0) = -FAout;
	dY(1) = -a*b*DB*sqrt( (2*g*VB)/(PB) ) + FAout + FBin ;
	dY(2) = (FBin/VB)*(TBin-TB) + (FAout/VB)*(TA0 - TB);
	return dY;
}

matrix ff1T(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	y = -1.0 * std::cos(0.1 * x(0)) * std::exp(-1.0 * std::pow(0.1 * x(0) - 2.0 * M_PI, 2.0)) + 0.002 * std::pow(0.1 * x(0), 2.0);
	return y;
}

// LAB 2 (K2)

matrix ff2T(matrix x, matrix ud1, matrix ud2) {
	matrix y;
	//y = x(0)*x(0) + x(1)*x(1) - std::cos(2.5 * M_PI * x(0)) - std::cos(2.5 * M_PI * x(0)) + 2;
	y = pow(x(0), 2.0) + pow(x(1), 2.0) - cos(2.5 * M_PI * x(0)) - cos(2.5 * M_PI * x(1)) + 2;
  return y;
}
// od 0s do 100s co 0.1s

matrix df2(double t, matrix Y, matrix ud1, matrix ud2) {

	auto alphaT = Y(0);
	auto omegaT = Y(1);

	auto k1 = ud2(0);
	auto k2 = ud2(1);

	auto alphaRef = ud1(0);
	auto omegaRef = ud1(1);


	double l = 2;   // m
	double m_r = 1; // kg
	double m_c = 5; // kg

	double b = 0.25; // Nms

	const double one_third = 1.0 / 3.0;

	double I = (one_third * m_r * pow(l, 2)) +
			 	(m_c * pow(l, 2));

	double M_t = k1 * (alphaRef - alphaT) + k2 * (omegaRef - omegaT);

	matrix dY(2, 1);
	dY(0) = omegaT;
	dY(1) = (M_t - b * omegaT)/I;
	return dY;
}

matrix ff2R(matrix x, matrix ud1, matrix ud2)
{
	matrix y = {0};

	double Alfa0 = 0;
	double Omega0 = 0;

	double AlfaRef = M_PI;
	double OmegaRef = 0;
	// Y0 zawiera warunku poczatkowe
	matrix Y0 = matrix(2, new double[2] {Alfa0, Omega0});
	matrix Yref  = matrix(2, new double[2] {AlfaRef, OmegaRef});
	matrix* Y = solve_ode(df2, 0, 0.1, 100, Y0, Yref, x);

	//x = ud2;

	int n = get_len(Y[0]);
	for(int i = 0; i < n; i++) {
		y = y +
			10 * pow(Yref(0) - Y[1](i, 0), 2) +
			pow(Yref(1) - Y[1](i, 1), 2) +
			pow(
				x(0) * (Yref(0) - Y[1](i, 0)) +
				x(1) * (Yref(1) - Y[1](i, 1)), 2
			);
	}
	y = y * 0.1;
	Y[0].~matrix();
	Y[1].~matrix();
	return y;
}

//LAB3

matrix ff3T(matrix x, matrix ud1, matrix ud2)
{
	if (fabs(sqrt(pow(x(0) / M_PI, 2) + pow(x(1) / M_PI, 2))) < 1e-12)
	{
		return 1.0;
	}
	matrix y;
	y = (sin(M_PI * sqrt(pow(x(0) / M_PI, 2) + pow(x(1) / M_PI, 2))) / (M_PI * sqrt(pow(x(0) / M_PI, 2) + pow(x(1) / M_PI, 2))));
	return y;
}

matrix df3(double t, matrix Y, matrix ud1, matrix ud2) {

	// Dane
	double m = 0.6;    // kg
	double r = 0.12;   // m
	double y0 = 100.0; // m
	double C = 0.47;
	double rho = 1.2; // kg/m^3
	double S = M_PI * pow(r, 2.0);
	double g = 9.81; // m/s^2
	double omega = ud1(1);
	matrix v = matrix(2, new double[2] {Y(1), Y(3)});
	matrix D(2, new double[2] {
		0.5 * C * rho * S * v(0) * abs(v(0)),
			0.5 * C * rho * S * v(1) * abs(v(1)),
		});
	matrix FM(
		2, new double[2] {rho* v(1)* omega* S* r, rho* v(0)* omega* S* r});

	matrix dY(4, new double[4] {Y(1), (-D(0) - FM(0)) / m, Y(3),
		(-m * g - D(1) - FM(1)) / m});

	return dY;
}

matrix ff3R(matrix x, matrix ud1, matrix ud2) {
	matrix y = { 0 };

	// Y0 zawiera warunki poczatkowe
	matrix Y0 = matrix(4, new double[4] {0, x(0), 100, 0});

	matrix* Y = solve_ode(df3, 0, 0.1, 7.0, Y0, x, ud2);

	int n = get_len(Y[0]);
	int i0 = 0, i50 = 0;
	for (int i = 0; i < n; i++) {
		if (fabs(Y[1](i, 2) - 50) < fabs(Y[1](i50, 2) - 50)) {
			i50 = i;
		}
		if (fabs(Y[1](i, 2)) < fabs(Y[1](i0, 2))) {
			i0 = i;
			y = -Y[1](i0, 0);
		}
		if (fabs(x(0)) - 10 > 0) {
			y = y + ud2(0) * pow(fabs(x(0)) - 10, 2);
		}
		if (fabs(x(1)) - 10 > 0) {
			y = y + ud2(0) * pow(fabs(x(1)) - 10, 2);
		}
		if (fabs(Y[1](i50, 2) - 5) - 2 > 0) {
			y = y + ud2(0) * pow(fabs(Y[1](i50, 0) - 5) - 2, 2);
		}
	}
	//std::cout << "x(i50): " << Y[1](i50, 0) << " y :" << Y[1](i50, 2) << '\n';
	// y = y * 0.1;
	Y[0].~matrix();
	Y[1].~matrix();
	return y;
}

double g3T1(matrix x, double a) {
	return (-1 * x(0) + 1); //<= 0);
}
double g3T2(matrix x, double a) {
	return (-1 * x(1) + 1); //<= 0);
}
double g3T3(matrix x, double a) {
	return (std::sqrt(x(0) * x(0) + x(1) * x(1)) - a); //<= 0);
}

//LAB 4

matrix ff4T(matrix x, matrix ud1, matrix ud2)
{
return (1.0f/6.0f) * pow(x(0), 6) - 1.05 * pow(x(0), 4) + 2.0f * pow(x(0), 2) + pow(x(1), 2) + x(0) * x(1);
}

matrix gf4T(matrix x, matrix ud1, matrix ud2)
{
	double& x1 = x(0);
	double& x2 = x(1);
	return matrix(2, new double[2] {
		pow(x1, 5) - (21.0/5.0) * pow(x1, 3) + 4*x1 + x2,
		x1 + 2 * x2
	});
}

matrix Hf4T(matrix x, matrix ud1, matrix ud2)
{
	double& x1 = x(0);
	double& x2 = x(1);
	matrix ret = matrix(2, 2);
	ret(0, 0) = 5*pow(x1, 4) + (63.0/5.0) * pow(x1, 2) + 4;
	ret(0, 1) = 1;
	ret(1, 0) = 1;
	ret(1, 1) = 2;
	return ret;

}
#define TRAINING_DATA_AMOUNT 100
void load_data_lab4(matrix*& xD, matrix*& yD) {
    static bool initialized = false;
    static matrix dataX[TRAINING_DATA_AMOUNT];
    static matrix dataY[TRAINING_DATA_AMOUNT];

    if(initialized) goto return_data;

    { // Scope for goto
        ifstream daneX("./dane/lab4/XData.txt");
        ifstream daneY("./dane/lab4/YData.txt");

        if(!daneX.is_open() || !daneY.is_open()) {
            cerr << "Error opening data files" << endl;
            return;
        }

        // First approach: Store all X values in a 2D array first
        double xValues[3][TRAINING_DATA_AMOUNT];

        // Read all 3 lines of X data
        for(int line = 0; line < 3; line++) {
            for(int i = 0; i < TRAINING_DATA_AMOUNT; i++) {
                daneX >> xValues[line][i];
                daneX.ignore(); // Skip semicolon and space
            }
        }

        // Now create matrices
        for(int i = 0; i < TRAINING_DATA_AMOUNT; i++) {
            dataX[i] = matrix(3, 1);
            for(int row = 0; row < 3; row++) {
                dataX[i](row) = xValues[row][i];
            }
        }

        daneX.close();

        // Read Y data (single line)
        for(int i = 0; i < TRAINING_DATA_AMOUNT; i++) {
            dataY[i] = matrix(1, 1);
            double yValue;
            daneY >> yValue;
            dataY[i](0) = yValue;
            daneY.ignore(); // Skip semicolon and space
        }

        daneY.close();
    }

    initialized = true;

return_data:
    xD = dataX;
    yD = dataY;
}

double h(matrix x, matrix theta, matrix ud2) {
	return 1 / (1 + std::exp(m2d(-1 * trans(theta) * x)));
}

matrix ff4R(matrix theta, matrix ud1, matrix ud2) {
	double sum = 0.0;
	const int m = TRAINING_DATA_AMOUNT;
	static matrix *xData, *yData;
	static bool loaded;
	if(loaded) goto postLoad;
	load_data_lab4(xData, yData);
	postLoad:
	for(int i = 0; i < m; i++) {
		sum += (
			m2d(yData[i]) * std::log(h(xData[i], theta, ud2)) +
			(1 - m2d(yData[i])) * std::log(1 - h(xData[i], theta, ud2))
		);
	}
	return -sum/m;
}

matrix gf4R(matrix theta, matrix ud1, matrix ud2) {
	matrix ret = matrix(3, 1);
	const int m = TRAINING_DATA_AMOUNT;
	static matrix *xData, *yData;
	static bool loaded;
	if(loaded) goto postLoad;
	load_data_lab4(xData, yData);
	postLoad:
	for(int j = 0; j < 3; j++) {
		double deriv = 0.0;
		for(int i = 0; i < m; i++) {
			deriv += (h(theta,xData[i], ud2) - m2d(yData[i]))*xData[i](j);
		}
		ret(j) = deriv/m;
	}
	return ret;
}

double poprawne4R(matrix theta) {
    const int m = TRAINING_DATA_AMOUNT;
	matrix *xData, *yData;
	double sum = 0.0;
	load_data_lab4(xData, yData);
	for(int i = 0; i < m; i++) {
	    if(h(xData[i], theta, 0) > 0.5 && yData[i](0) == 1) sum++;
	}
	return sum/m;
}

matrix zlotf4T(matrix a, matrix d, matrix x)
{

	matrix arg = x;
	for (int i = 0; i < get_len(x); i++)
	{
		arg(i) = x(i) + d(i) * a(0);
	}
	return ff4T(arg);

}

matrix zlotf4R(matrix a, matrix d, matrix x)
{

	matrix arg = x;
	for (int i = 0; i < get_len(x); i++)
	{
		arg(i) = x(i) + d(i) * a(0);
	}
	return ff4R(arg);

}


//LAB 5

matrix ff5T1(matrix x, matrix ud1, matrix ud2)
{
	const matrix a = ud1[0];
	const matrix w = ud1[1];
	return a(0) * (
		pow(x(0) - 3.0, 2.0) +
		pow(x(1) - 3.0, 2.0)
	);
}

matrix ff5T2(matrix x, matrix ud1, matrix ud2)
{
	const matrix a = ud1[0];
	const matrix w = ud1[1];
	return (1.0 / a(0)) * (
		pow(x(0) + 3.0, 2.0) +
		pow(x(1) + 3.0, 2.0)
	);
}
matrix ff5TX(matrix x, matrix ud1, matrix ud2)
{
	const matrix a = ud1[0];
	const matrix w = ud1[1];
	return (w * ff5T1(x,ud1, ud2)) + ((1-w) * ff5T2(x, ud1, ud2));
}

matrix gg5T1(matrix a, matrix ud1, matrix ud2)
{
	const matrix d = ud2[0];
	const matrix p = ud2[1];
	return ff5T1(p + a * d, ud1, ud2);
}

matrix gg5T2(matrix a, matrix ud1, matrix ud2)
{
	const matrix d = ud2[0];
	const matrix p = ud2[1];
	return ff5T2(p + a * d, ud1, ud2);
}

matrix gg5TX(matrix a, matrix ud1, matrix ud2)
{
	const matrix d = ud2[0];
	const matrix p = ud2[1];
	return ff5TX(p + a * d, ud1, ud2);
}

// Obliczanie masy.
matrix ff5R1(matrix x, matrix ud1, matrix ud2) {
	const double l = x(0);
	const double d = x(1);
	const double Volume = M_PI * 0.25 * d * d * l; // m^3
	const double Density = 8920; // kg / m^3
	return Volume * Density;
}

// Obliczanie ugięcia
matrix ff5R2(matrix x, matrix ud1, matrix ud2) {
	const double l = x(0);
	const double d = x(1);
	const double P = 2000; // Niutonów
	const double E = 120e9; // paskali
	return (64 * P * l * l * l) / (3 * E * M_PI * d*d*d*d);
}

// Obliczanie naprężenia
matrix ff5R3(matrix x, matrix ud1, matrix ud2) {
	const double l = x(0);
	const double d = x(1);
	const double P = 2000; // Niutonów
	return (32 * P * l) / (M_PI * d * d * d);
}

matrix ff5RX(matrix x, matrix ud1, matrix ud2) {
	const double l = x(0);
	const double d = x(1);
	const matrix a = ud1[0];
	const matrix w = ud1[1];

	const double masa = m2d(ff5R1(x, ud1, ud2));
	const double ugiecie = m2d(ff5R2(x, ud1, ud2));
	const double naprezenie = m2d(ff5R3(x, ud1, ud2));

	double result = m2d(w * masa + (1-w)*ugiecie);
	double penalty = 2.137e10;
	if(l < 0.2) result += penalty * pow(0.2 - l, 2.0);
	if(l > 1) result += penalty * pow(l - 1, 2.0);
	if(d < 0.01) result += penalty * pow(0.01 - d, 2.0);
	if(d > 0.05) result += penalty * pow(d - 0.05, 2.0);

	if(ugiecie > 2.5e-3) result += penalty * pow(ugiecie - 2.5e-3, 2.0);
	if(naprezenie > 300e6) result += penalty * pow(naprezenie - 300e6, 2.0);

	return result;

}


matrix gg5RX(matrix a, matrix ud1, matrix ud2)
{
	const matrix d = ud2[0];
	const matrix p = ud2[1];
	return ff5RX(p + a * d, ud1, ud2);
}

//LAB 6

matrix ff6T(matrix x, matrix ud1, matrix ud2)
{
	return x(0) + x(1) - cos(2.5 * M_PI * x(0)) - cos(2.5 * M_PI * x(1)) + 2;
}
