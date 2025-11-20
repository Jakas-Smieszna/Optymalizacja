#include "opt_alg.h"
#include "exceptions.h"

solution MC(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	// Zmienne wejœciowe:
	// ff - wskaŸnik do funkcji celu
	// N - liczba zmiennych funkcji celu
	// lb, ub - dolne i górne ograniczenie
	// epslion - zak³¹dana dok³adnoœæ rozwi¹zania
	// Nmax - maksymalna liczba wywo³añ funkcji celu
	// ud1, ud2 - user data
	try
	{
		solution Xopt;
		while (true)
		{
			Xopt = rand_mat(N);									// losujemy macierz Nx1 stosuj¹c rozk³ad jednostajny na przedziale [0,1]
			for (int i = 0; i < N; ++i)
				Xopt.x(i) = (ub(i) - lb(i)) * Xopt.x(i) + lb(i);// przeskalowywujemy rozwi¹zanie do przedzia³u [lb, ub]
			Xopt.fit_fun(ff, ud1, ud2);							// obliczmy wartoœæ funkcji celu
			if (Xopt.y < epsilon)								// sprawdzmy 1. kryterium stopu
			{
				Xopt.flag = 1;									// flaga = 1 ozancza znalezienie rozwi¹zanie z zadan¹ dok³adnoœci¹
				break;
			}
			if (solution::f_calls > Nmax)						// sprawdzmy 2. kryterium stopu
			{
				Xopt.flag = 0;									// flaga = 0 ozancza przekroczenie maksymalne liczby wywo³añ funkcji celu
				break;
			}
		}
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution MC(...):\n" + ex_info);
	}
}

double* expansion(matrix(*ff)(matrix, matrix, matrix), double x0, double d, double alpha, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		double* p = new double[2] { 0, 0 };

		int i = 0;
		double x1 = x0 + d;
		double x0i, x1i, x2i;
		solution::f_calls += 2;
		if (ff(matrix(x0), ud1, ud2) == ff(matrix(x1), ud1, ud2)) {
			p[0] = x0;
			p[1] = x1;
			return p;
		}
		solution::f_calls += 2;
		if (ff(matrix(x0), ud1, ud2) < ff(matrix(x1), ud1, ud2)) {
			d = -1.0 * d;
			x1 = x0 + d;
			solution::f_calls += 2;
			if (ff(matrix(x0), ud1, ud2) <= ff(matrix(x1), ud1, ud2)) {
				p[0] = x1;
				p[1] = x0 - d;
				return p;
			}
		}
		{
			x0i = x0, x1i = x0, x2i = x1;
			do
			{
				if (solution::f_calls > Nmax) {
					throw std::runtime_error("Przekroczono maksymalna liczba poszukiwan przedzialu - niepowodzenie.");
				}
				i = i + 1;
				x0i = x1i;
				x1i = x2i;
				x2i = x0 + pow(alpha, i) * d;
				solution::f_calls += 2;
			} while (!(ff(matrix(x1i), ud1, ud2) <= ff(matrix(x2i), ud1, ud2)));
			if (d > 0) {
				p[0] = x0i;
				p[1] = x2i;
				return p;
			}
		}
		if (ud1 < x2i) p[0] = x2i;
		else p[0] = ud1(0);
		if (x0i < ud2(0)) p[1] = x0i;
		else p[1] = ud2(0);
		return p;

	}
	catch (string ex_info)
	{
		throw ("double* expansion(...):\n" + ex_info);
	}
}

double phi(int x)
{
	int first = 0, second = 1, sum = 0;

	if ((x == 1) || (x == 2)) { return 1; }
	for (int i = 2; i <= x; ++i)
	{
		sum = first + second;
		first = second;
		second = sum;
	}

	return (double)sum;
}

solution fib(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2)
{
	try
	{
		//std::cout << a << std::endl;
		int k = 1;
		while (phi(k) <= (b - a) / epsilon) { ++k; }
		//std::cout << k << " - " << phi(k) << " k-1: " << phi(k - 1) << " | " << (b - a) / epsilon << std::endl;
		solution Xopt;
		double c = b - (phi(k - 1) / phi(k)) * (b - a);
		double d = a + b - c;
		//std::cout << a << " | " << b << " | " << c << " | " << d << std::endl;
		for (int i = 0; i < k - 3; i++) {
			Xopt.x = c;
			matrix ff_from_c = Xopt.fit_fun(ff, ud1, ud2);
			Xopt.x = d;
			matrix ff_from_d = Xopt.fit_fun(ff, ud1, ud2);
			if (ff(c, ud1, ud2) < ff(d, ud1, ud2)) {
				b = d;
			}
			else {
				a = c;
			}
			//std::cout << a << " | " << b << " | " << c << " | " << d << std::endl;
			c = b - (phi(k - i - 2) / phi(k - i - 1)) * (b - a);
			d = a + b - c;
		}

		Xopt.x = c;
		Xopt.fit_fun(ff, ud1, ud2);


		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution fib(...):\n" + ex_info);
	}

}

solution lag(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;

		int i = 0;
		double l, m;
		matrix ai(2, 1, a);
		matrix bi(2, 1, b);
		matrix ci(2, 1, (a + b) / 2.0);
		matrix di(2, 1, 0.0);
		do
		{
			Xopt.f_calls = Xopt.f_calls + 3;
			l = ff(ai(0), ud1, ud2)(0) * (pow(bi(0), 2) - pow(ci(0), 2))
				+ ff(bi(0), ud1, ud2)(0) * (pow(ci(0), 2) - pow(ai(0), 2))
				+ ff(ci(0), ud1, ud2)(0) * (pow(ai(0), 2) - pow(bi(0), 2));
			Xopt.f_calls = Xopt.f_calls + 3;
			m = ff(ai(0), ud1, ud2)(0) * (bi(0) - ci(0))
				+ ff(bi(0), ud1, ud2)(0) * (ci(0) - ai(0))
				+ ff(ci(0), ud1, ud2)(0) * (ai(0) - bi(0));

			if (m <= 0) {
				throw std::invalid_argument("Niepoprawne wyniki podczas obliczen - zle dane/blad w kodzie.");
			}

			di(1) = 0.5 * l / m;
			if (ai(0) < di(1) && di(1) < ci(0)) {

				Xopt.f_calls = Xopt.f_calls + 2;
				if (ff(di(1), ud1, ud2) < ff(ci(0), ud1, ud2)) {
					ai(1) = ai(0);
					ci(1) = di(1);
					bi(1) = ci(0);
				}
				else {
					ai(1) = di(1);
					ci(1) = ci(0);
					bi(1) = bi(0);
				}

			}
			else {

				if (ci(0) < di(1) && di(1) < bi(0)) {

					Xopt.f_calls = Xopt.f_calls + 2;
					if (ff(di(1), ud1, ud2) < ff(ci(0), ud1, ud2)) {
						ai(1) = ci(0);
						ci(1) = di(1);
						bi(1) = bi(0);
					}
					else {
						ai(1) = ai(1);
						ci(1) = ci(0);
						bi(1) = di(0);
					}

				}
				else {
					throw std::invalid_argument("Niepoprawne wyniki podczas obliczen - zle dane/blad w kodzie.");
				}

			}
			//std::cout << ai(0) << " | " << bi(0) << std::endl;

			i = i + 1;
			ai(0) = ai(1);
			bi(0) = bi(1);
			ci(0) = ci(1);
			di(0) = di(1);

			if (Xopt.f_calls > Nmax) {
				Xopt.flag = 0;
				throw std::runtime_error("Przekroczono limit wywolan funkcji.");
			}



		} while (!(bi(0) - ai(0) < epsilon || abs(di(1) - di(0)) < gamma));
		Xopt.x(0) = di(1);
		Xopt.fit_fun(ff, ud1, ud2);
		return Xopt;

	}
	catch (string ex_info)
	{
		throw ("solution lag(...):\n" + ex_info);
	}
}

solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2)
{

	solution Xopt;

	try
	{

		Xopt.x = x0;
		if (alpha - 1.0 > -TOL || alpha < TOL) {
			Xopt.flag = 0;
			throw BadArguments("Wspolczynnik ekspansi poza przedzialem (0, 1).");
		}

		do {
			solution XB;
			XB.x = Xopt.x;
			XB.fit_fun(ff, ud1, ud2);
			Xopt = HJ_trial(ff, XB, s, ud1, ud2);
			Xopt.fit_fun(ff, ud1, ud2);
			if (Xopt.y - XB.y < -TOL){

				do {

					solution XB_ = XB;
					XB = Xopt;
					for (int i = 0; i < get_len(Xopt.x); i++) {
						Xopt.x(i) = 2 * XB.x(i) - XB_.x(i);
					}
					Xopt = HJ_trial(ff, Xopt, s, ud1, ud2);
					XB.fit_fun(ff, ud1, ud2);
					Xopt.fit_fun(ff, ud1, ud2);
					if (solution::f_calls > Nmax){
						Xopt.flag = 0;
						throw ToManyCalls("Etap roboczy.\nLiczba wywolan = " + to_string(solution::f_calls) + "\nLimit wywyolan = " + to_string(Nmax));
					}

				} while (!(Xopt.y - XB.y > -TOL));
				Xopt = XB;

			}
			else {
				s = alpha * s;
			}

			if (solution::f_calls > Nmax) {
				Xopt.flag = 0;
				throw ToManyCalls("Koniec iteracji algorytmu.\nLiczba wywolan = " + to_string(solution::f_calls) + "\nLimit wywyolan = " + to_string(Nmax));
			}

		} while (!(s < epsilon));
		Xopt.flag = 1;

	}
	catch (const std::exception& ex)
	{
		std::cerr << "PRZECHWYCONO WYJATEK - Hook-Jeeves: " << ex.what() << "\n";
	}
	catch (string ex_info)
	{
		throw ("WYJATEK - Hook-Jeeves:\n" + ex_info);
	}

	Xopt.fit_fun(ff, ud1, ud2);
	return Xopt;

}

solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2)
{
	try
	{

		int XD = get_len(XB.x);//Wymiar przestrzeni
		matrix* e = new matrix[XD];
		for (int i = 0; i < XD; i++) {
			e[i] = matrix(XD, 1, 0.0);
			e[i](i) = 1.0;
		}
		for (int j = 0; j < XD; j++) {

			XB.fit_fun(ff, ud1, ud2);
			solution::f_calls += 1;
			if (XB.x(j) + s - ud1(j) > -TOL 
				&& XB.x(j) + s - ud2(j) < TOL 
				&& (*ff)(XB.x + e[j] * s, ud1, ud2) - XB.y < -TOL) {
				
				XB.x = XB.x + s * e[j];

			}
			else {

				solution::f_calls += 1;
				if (XB.x(j) - s - ud1(j) > -TOL
					&& XB.x(j) - s - ud2(j) < TOL
					&& (*ff)(XB.x - e[j] * s, ud1, ud2) - XB.y < -TOL) {

					XB.x = XB.x - s * e[j];

				}

			}

		}
		delete[] e;
		
	}
	catch (const std::exception& ex)
	{
		std::cerr << "PRZECHWYCONO WYJATEK - Hook-Jeeves - etap probny: " << ex.what() << "\n";
	}
	catch (string ex_info)
	{
		throw ("WYJATEK - Hook-Jeeves - etap probny:\n" + ex_info);
	}

	return XB;

}

solution Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{

	solution Xopt;

	try
	{

		if (alpha - 1.0 < TOL) {
			Xopt.flag = 0;
			throw BadArguments("Wspolczynnik ekspansi mniejszy lub rowny 1.");
		}
		if (beta - 1.0 > -TOL || beta < TOL) {
			Xopt.flag = 0;
			throw BadArguments("Wspolczynnik korekcji poza przedzialem (0, 1).");
		}

		int i = 0;
		int n = get_len(x0);
		matrix d(n, n, 0.0);
		for (int i = 0; i < n; i++) {
			d(i,i) = exp(i + 1);
		}
		matrix lambda(n, 1, 0.0);
		matrix p(n, 1, 0.0);
		matrix xb = x0;
		matrix s = s0;

		do {

			for (int j = 0; j < n; j++)
			{
				solution::f_calls += 2;
				if (ff(xb + s(j) * d(j), ud1, ud2)(0) < ff(xb, ud1, ud2)(0))
				{
					xb = xb + (s(j) * d(j));
					lambda(j) = lambda(j) + s(j);
					s(j) = alpha * s(j);
				}
				else
				{
					s(j) = -beta * s(j);
					p(j) = p(j) + 1;
				}
			}
			i++;
			Xopt.x = xb;
			//NOWA WERSJA
			bool warunek = 1;
			for (int j = 0; warunek && j < n; j++) {
				if (abs(lambda(j)) < TOL) warunek = 0;
				if (abs(p(j)) < TOL) warunek = 0;
			}
			if (warunek) {

				//LICZENIE dj(i+1):
				{

					matrix Q(n, n, 0.0);
					for (int j = 0; j < n; j++) {
						for (int i = 0; i < n; i++) {
							if (j >= i) {
								Q(j, i) = lambda(j);
							}
							else break;
						}
					}
					matrix D = ident_mat(n);
					Q = D * Q;

					matrix v = Q[0];
					d[0] = v / norm(v);

					for (int j = 1; j < n; j++) {

						matrix pom = Q[j];
						for (int k = 0; k < j; k++) {

							pom = pom - ((trans(Q) * d[k])() * d[k]);

						}
						v = pom;
						d[j] = v / norm(v);

					}
					
					lambda = matrix(n, 1, 0.0);
					p = matrix(n, 1, 0.0);
					s = s0;

				}
				//KONIEC LICZENIA dj(i+1)

			}
			if (solution::f_calls > Nmax)
			{
				Xopt.flag = 0;
				throw ToManyCalls("Koniec iteracji.\nLiczba wywolan = " + to_string(solution::f_calls) + "\nLimit wywyolan = " + to_string(Nmax));
			}

		} while (norm(s) > epsilon);

		Xopt.flag = 1;

	}
	catch (const std::exception& ex)
	{
		std::cerr << "PRZECHWYCONO WYJATEK - Rosen: " << ex.what() << "\n";
	}
	catch (string ex_info)
	{
		throw ("WYJATEK - ROSEN:\n" + ex_info);
	}

	Xopt.fit_fun(ff, ud1, ud2);
	return Xopt;

}

solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try {
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution pen(...):\n" + ex_info);
	}
}

solution sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution sym_NM(...):\n" + ex_info);
	}
}

solution SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution SD(...):\n" + ex_info);
	}
}

solution CG(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution CG(...):\n" + ex_info);
	}
}

solution Newton(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix),
	matrix(*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Newton(...):\n" + ex_info);
	}
}

solution golden(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution golden(...):\n" + ex_info);
	}
}

solution Powell(matrix(*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Powell(...):\n" + ex_info);
	}
}

solution EA(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, int mi, int lambda, matrix sigma0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution EA(...):\n" + ex_info);
	}
}
