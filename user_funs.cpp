#include"user_funs.h"

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

