/*********************************************
Kod stanowi uzupe�nienie materia��w do �wicze�
w ramach przedmiotu metody optymalizacji.
Kod udost�pniony na licencji CC BY-SA 3.0
Autor: dr in�. �ukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia G�rniczo-Hutnicza
Data ostatniej modyfikacji: 30.09.2025
*********************************************/


//Pliki z wynikami typu csv znajdują się w folderze out, przykładowo dla debugowania x64:
//Optymalizacja\out\build\x64-debug\nazwa_pliku.csv
//Wykonanie zadania: MF, JG, MG, AG


#include "matrix.h"
#include "user_funs.h"
#include <cstdlib>
#include <ctime>
#define _USE_MATH_DEFINES
#include"opt_alg.h"
#include <cmath>


void lab0();
void lab1();
void lab2();
void lab3();
void lab4();
void lab5();
void lab6();

int main()
{
	try
	{
		lab5();
	}
	catch (string EX_INFO)
	{
		cerr << "ERROR:\n";
		cerr << EX_INFO << endl << endl;
	}
	return 0;
}

void lab0()
{
	//Funkcja testowa
	double epsilon = 1e-2;									// dok�adno��
	int Nmax = 10000;										// maksymalna liczba wywo�a� funkcji celu
	matrix lb(2, 1, -5), ub(2, 1, 5),						// dolne oraz g�rne ograniczenie
		a(2, 1);											// dok�adne rozwi�zanie optymalne
	solution opt;											// rozwi�zanie optymalne znalezione przez algorytm
	a(0) = -1;
	a(1) = 2;
	opt = MC(ff0T, 2, lb, ub, epsilon, Nmax, a);			// wywo�anie procedury optymalizacji
	cout << opt << endl << endl;							// wypisanie wyniku
	solution::clear_calls();								// wyzerowanie licznik�w

	//Wahadlo
	Nmax = 1000;											// dok�adno��
	epsilon = 1e-2;											// maksymalna liczba wywo�a� funkcji celu
	lb = 0, ub = 5;											// dolne oraz g�rne ograniczenie
	double teta_opt = 1;									// maksymalne wychylenie wahad�a
	opt = MC(ff0R, 1, lb, ub, epsilon, Nmax, teta_opt);		// wywo�anie procedury optymalizacji
	cout << opt << endl << endl;							// wypisanie wyniku
	solution::clear_calls();								// wyzerowanie licznik�w

	//Zapis symulacji do pliku csv
	matrix Y0 = matrix(2, 1),								// Y0 zawiera warunki pocz�tkowe
		MT = matrix(2, new double[2] { m2d(opt.x), 0.5 });	// MT zawiera moment si�y dzia�aj�cy na wahad�o oraz czas dzia�ania
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, NAN, MT);	// rozwi�zujemy r�wnanie r�niczkowe
	ofstream Sout("symulacja_lab0.csv");					// definiujemy strumie� do pliku .csv
	Sout << hcat(Y[0], Y[1]);								// zapisyjemy wyniki w pliku
	Sout.close();											// zamykamy strumie�
	Y[0].~matrix();											// usuwamy z pami�ci rozwi�zanie RR
	Y[1].~matrix();
}

void lab1()
{
	srand(time(NULL));
	//Funkcja testowa
	double alfa = 1.5;										//wspolczynnik ekspansji
	double krok_d = 0.002;									//krok/odleglosc do ekspansji
	double gamma = 1e-2;									//kolejna dokladnosc
	double epsilon = 1e-2;									// dok�adno��
	int Nmax = 10000;										// maksymalna liczba wywo�a� funkcji celu
	matrix lb(1, 1, -100), ub(1, 1, 100),						// dolne oraz g�rne ograniczenie
		ps(1, 1, rand() % 201 - 100);						// punkt startowy
	solution opt;											// rozwi�zanie optymalne znalezione przez algorytm

	//-----FUNKCJA TESTOWA-----------------------------------------------------------

	char kont = '1';
	fstream Sout;
	Sout.open("testy_lab1.csv", std::ios::out);
	while (kont == '1') {
		for (int i = 0; i < 1; i++) {							//JG:mozna wybrac liczbe powtorzen

			ps(0) = rand() % 201 - 100;
			double* obszar = expansion(*ff1T, ps(0), krok_d, alfa, Nmax, lb, ub);

			cout << "\nKrok d = " << krok_d << "\tWspolczynnik ekspansji alfa = " << alfa << ".\n";
			cout << "Punkt startowy = " << ps(0) << "\tUzyskany przedzial = [" << obszar[0] << ", " << obszar[1] << "].\n";			//JG:pozwala przywrocic liczbe po ekspansji, gdy zostanie wycyszczona.
			cout << "EKSPANSJA: fcalls = " << solution::f_calls << ".\n\n";
			if (Sout.good() == true) Sout << ps(0) << "\t" << obszar[0] << "\t" << obszar[1] << "\t" << solution::f_calls << "\t";
			solution::clear_calls();

			cout << "FIBBONACI:\n";
			opt = fib(ff1T, obszar[0], obszar[1], epsilon, lb, ub);								// wywołanie procedury optymalizacji
			cout << opt << endl << endl;														// wypisanie wyniku
			if (Sout.good() == true) Sout << opt.x(0) << "\t" << opt.y(0) << "\t" << solution::f_calls << "\tlokalne\t";
			solution::clear_calls();

			cout << "LAGRANGE:\n";
			opt = lag(ff1T, obszar[0], obszar[1], epsilon, gamma, Nmax, lb, ub);				// wywołanie procedury optymalizacji
			delete[] obszar;
			cout << opt << endl << endl;														// wypisanie wyniku
			if (Sout.good() == true) Sout << opt.x(0) << "\t" << opt.y(0) << "\t" << solution::f_calls << "\tlokalne\n";
			solution::clear_calls();

		}

		cin >> kont;

	}

	if (Sout.good() == true) Sout << "\nBEZ EKSPANSJI:\n\n";
	kont = '1';
	while (kont == '1') {

		cout << "\nBEZ EKSPANSJI:\n";
		cout << "FIBBONACI:\n";
		opt = fib(ff1T, -100.0, 100.0, epsilon, lb, ub);											// wywołanie procedury optymalizacji
		cout << opt << endl << endl;																// wypisanie wyniku
		if (Sout.good() == true) Sout << opt.x(0) << "\t" << opt.y(0) << "\t" << solution::f_calls << "\tlokalne\t";
		solution::clear_calls();

		cout << "LAGRANGE:\n";
		opt = lag(ff1T, -100.0, 100.0, epsilon, gamma, Nmax, lb, ub);								// wywołanie procedury optymalizacji
		cout << opt << endl << endl;																// wypisanie wyniku
		if (Sout.good() == true) Sout << opt.x(0) << "\t" << opt.y(0) << "\n" << solution::f_calls << "\tlokalne\n";
		solution::clear_calls();

		cin >> kont;

	}
	Sout.close();

	//-----PROBLEM RZECZYWISTY-------------------------------------------------------

	kont = '1';
	Sout.open("problem_przebieg_lab1.csv", std::ios::out);

	// NOWY NAGŁÓWEK z wymaganymi kolumnami
	if (Sout.good() == true) {
		Sout << "Metoda\tDA*\ty*\tLiczba_wywolan_funkcji_celu\n";
	}

	auto zapiszSymulacje = [](matrix* Y, const string& nazwa_pliku) {
		ofstream plik(nazwa_pliku);
		if (plik.good()) {
			plik << "Czas[s]\tVA[m3]\tVB[m3]\tTB[C]\n";
			int n = get_len(Y[0]);
			for (int i = 0; i < n; i++) {
				plik << Y[0](i) << "\t"
					<< Y[1](i, 0) << "\t"
					<< Y[1](i, 1) << "\t"
					<< Y[1](i, 2) << "\n";
			}
			cout << "Zapisano symulacje: " << nazwa_pliku << " (" << n << " punktow czasowych)\n";
		}
		plik.close();
		Y[0].~matrix();
		Y[1].~matrix();
		};

	while (kont == '1') {
		lb = matrix(1, 1, 1);
		ub = matrix(1, 1, 100);
		ps = matrix(1, 1, rand() % 101);

		for (int i = 0; i < 1; i++) {
			ps(0) = rand() % 101;
			double* obszar = expansion(*ff1R, ps(0), krok_d, alfa, Nmax, lb, ub);

			cout << "\nKrok d = " << krok_d << "\tWspolczynnik ekspansji alfa = " << alfa << ".\n";
			cout << "Punkt startowy = " << ps(0) << "\tUzyskany przedzial = [" << obszar[0] << ", " << obszar[1] << "].\n";
			int fcalls_ekspansja = solution::f_calls;
			cout << "EKSPANSJA: fcalls = " << solution::f_calls << ".\n\n";
			solution::clear_calls();

			cout << "LAGRANGE:\n";
			solution opt_lag = lag(ff1R, obszar[0], obszar[1], epsilon, gamma, Nmax, lb, ub);
			cout << opt_lag << endl << endl;

			int fcalls_lagrange = solution::f_calls;

			if (Sout.good() == true) {
				Sout << "Lagrange\t" << opt_lag.x(0) << "\t" << opt_lag.y(0) << "\t" << fcalls_lagrange << "\n";
			}
			solution::clear_calls();

			cout << "FIBBONACI:\n";
			solution opt_fib = fib(ff1R, obszar[0], obszar[1], epsilon, lb, ub);
			cout << opt_fib << endl << endl;

			int fcalls_fibonacci = solution::f_calls;

			if (Sout.good() == true) {
				Sout << "Fibonacci\t" << opt_fib.x(0) << "\t" << opt_fib.y(0) << "\t" << fcalls_fibonacci << "\n";
			}
			solution::clear_calls();

			cout << "\n=== ZAPIS PEŁNYCH SYMULACJI ===\n";

			// Warunki początkowe
			matrix Y0 = matrix(3, new double[3] {5, 1, 20});

			// SYMULACJA DLA LAGRANGE'A
			cout << "Symulacja dla Lagrange'a (DA = " << opt_lag.x(0) << ")...\n";
			matrix DA_lag = matrix(1, 1, opt_lag.x(0));
			matrix* Y_lag = solve_ode(df1, 0, 1, 2000, Y0, DA_lag, NAN);
			zapiszSymulacje(Y_lag, "symulacja_Lagrange_DA_" + to_string((int)opt_lag.x(0)) + ".csv");

			// SYMULACJA DLA FIBONACCI'EGO
			cout << "Symulacja dla Fibonacci'ego (DA = " << opt_fib.x(0) << ")...\n";
			matrix DA_fib = matrix(1, 1, opt_fib.x(0));
			matrix* Y_fib = solve_ode(df1, 0, 1, 2000, Y0, DA_fib, NAN);
			zapiszSymulacje(Y_fib, "symulacja_Fibonacci_DA_" + to_string((int)opt_fib.x(0)) + ".csv");



			delete[] obszar;
		}

		cin >> kont;
	}

	Sout.close();


	//Zapis symulacji do pliku csv
	//matrix Y0 = matrix(2, 1),								// Y0 zawiera warunki pocz�tkowe
	//	MT = matrix(2, new double[2] { m2d(opt.x), 0.5 });	// MT zawiera moment si�y dzia�aj�cy na wahad�o oraz czas dzia�ania
	//matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, NAN, MT);	// rozwi�zujemy r�wnanie r�niczkowe
	//ofstream Sout("testy_lab1.csv");						// definiujemy strumie� do pliku .csv
	//Sout << hcat(Y[0], Y[1]);								// zapisyjemy wyniki w pliku
	//Sout.close();											// zamykamy strumie�
	//Y[0].~matrix();											// usuwamy z pami�ci rozwi�zanie RR
	//Y[1].~matrix();

}

void lab2()
{

	srand(time(NULL));
	//Funkcja testowa
	double alfa = 0.0;										// wspolczynnik ekspansji (0.0 do 1.0 dla HJ i > 1 dla R - zmiana w petli, nie tutaj)
	double krok_s = 1.5;									// krok
	double beta = 0.5;										// wspolczynnik kontrakcji
	double epsilon = 1e-8;									// dokladnosc
	int Nmax = 10000;										// maksymalna liczba wywolan funkcji celu
	matrix lb(2, 1, -1.0), ub(2, 1, 1.0),					// dolne oraz g�rne ograniczenie
		ps(2, 1, double(rand() % 20001 - 10000)/10000.0);	// punkt startowy
	solution opt;											// rozwiazanie optymalne znalezione przez algorytm
	solution::clear_calls();

	//-----FUNKCJA TESTOWA-----------------------------------------------------------

  char kont = '1';
  fstream Sout;
  Sout.open("testy_lab2.csv", std::ios::out);
  while (kont == '1') {
	  for (int i = 0; i < 1; i++) {							//JG:mozna wybrac liczbe powtorzen

  		ps(0) = double(rand() % 20001 - 10000) / 10000.0;
  		ps(1) = double(rand() % 20001 - 10000) / 10000.0;
  		cout << "Punkt startowy = [" << ps(0) << ", " << ps(1) << "].\n";
  		cout << "Krok startowy = " << krok_s << ".\n\n";

	  	alfa = 0.5;

	  	cout << "HOOK-JEEVES:\n";
		if (Sout.good() == true) Sout << ps(0) << "\t" << ps(1) << '\t';
	  	opt = HJ(ff2T, ps, krok_s, alfa, epsilon, Nmax, lb, ub);							// wywołanie procedury optymalizacji
	  	cout << opt << endl << endl;

	  	if (Sout.good() == true) {
			Sout << opt.x(0) << "\t" << opt.x(1) << '\t' << opt.y(0)
			<< "\t" << solution::f_calls;
			if(fabs(opt.x(0)) < 10*epsilon && fabs(opt.x(1)) < 10*epsilon) {
			Sout << "\tTAK\t";
			} else Sout << "\tNIE\t";
		}

	  	solution::clear_calls();

  		alfa = 1.5;

        cout << "ROSENBROCK:\n";
	    opt = Rosen(ff2T, ps, matrix(2, 1, krok_s), alfa, beta, epsilon, Nmax, lb, ub);	// wywołanie procedury optymalizacji
	    cout << opt << endl << endl;
	  	if (Sout.good() == true) {
			Sout << opt.x(0) << "\t" << opt.x(1) << '\t' << opt.y(0)
			<< "\t" << solution::f_calls;
			if(fabs(opt.x(0)) < 10*epsilon && fabs(opt.x(1)) < 10*epsilon) {
			Sout << "\tTAK\n";
			} else Sout << "\tNIE\n";
		}
		solution::clear_calls();

  	}

  	cin >> kont;

  }
  Sout.close();

  //-----PROBLEM RZECZYWISTY-------------------------------------------------------

	matrix Yref(2, new double[2] {
	    M_PI, // Alpha ref
		0 // Omega ref
	}); // ud1
	matrix k(2, new double[2] {
    	5, // k1
    	5 // k2
	});
	cout << ff2R(k, Yref) << endl << endl;
	lb = matrix(2, 1, 0.0);
	ub = matrix(2, 1, 20.0),
	ps = matrix(2, 1, double(rand() % 200001)/10000.0);
	krok_s = 0.01;
	kont = '1';
    Sout.open("real_lab2.csv", std::ios::out);
    while (kont == '1') {
        for (int i = 0; i < 1; i++) {							//JG:mozna wybrac liczbe powtorzen
      		ps(0) = double(rand() % 200001)/10000.0;
      		ps(1) = double(rand() % 200001)/10000.0;
      		cout << "Punkt startowy = [" << ps(0) << ", " << ps(1) << "].\n";
      		cout << "Krok startowy = " << krok_s << ".\n\n";

           	alfa = 0.5;

           	cout << "HOOK-JEEVES:\n";
           	if (Sout.good() == true) Sout << ps(0) << "\t" << ps(1) << '\t';
           	opt = HJ(ff2R, ps, krok_s, alfa, epsilon, Nmax, Yref, ub);							// wywołanie procedury optymalizacji
           	cout << opt << endl << endl;

           	if (Sout.good() == true) {
          		Sout << opt.x(0) << "\t" << opt.x(1) << '\t' << opt.y(0)
          		<< "\t" << solution::f_calls << '\t';
           	}

       	solution::clear_calls();

          		alfa = 1.5;

                cout << "ROSENBROCK:\n";
            opt = Rosen(ff2R, ps, matrix(2, 1, krok_s), alfa, beta, epsilon, Nmax, Yref, ub);	// wywołanie procedury optymalizacji
            cout << opt << endl << endl;
       	if (Sout.good() == true) {
      		Sout << opt.x(0) << "\t" << opt.x(1) << '\t' << opt.y(0)
      		<< "\t" << solution::f_calls << '\n';
       	}
       	solution::clear_calls();

   	}

   	cin >> kont;

    }
  Sout.close();

 	auto zapiszSymulacje = [](matrix* Y, const string& nazwa_pliku) {
		ofstream plik(nazwa_pliku);
		if (plik.good()) {
			plik << "Czas[s]\tAlfa[rad]\tomega[rad/s]\n";
			int n = get_len(Y[0]);
			for (int i = 0; i < n; i++) {
				plik << Y[0](i) << "\t"
					<< Y[1](i, 0) << "\t"
					<< Y[1](i, 1) << "\n";
			}
			cout << "Zapisano symulacje: " << nazwa_pliku << " (" << n << " punktow czasowych)\n";
		}
		plik.close();
		Y[0].~matrix();
		Y[1].~matrix();
		};

		solution::clear_calls();

		cout << "\n=== ZAPIS PEŁNYCH SYMULACJI ===\n";

		// Warunki początkowe
		matrix k_HJ = matrix(2, new double[2] {2.8664, 10.3436});
		matrix k_Rosen = matrix(2, new double[2] {3.14159, 11.3324});

		double Alfa0 = 0;
		double Omega0 = 0;
		// Y0 zawiera warunku poczatkowe
		matrix Y0 = matrix(2, new double[2] {Alfa0, Omega0});
		// matrix Yref = matrix(2, new double[2] {AlfaRef, OmegaRef});
		// Yref - wyżej w pliku

		// SYMULACJA DLA HJ'A
		cout << "Symulacja dla HJ'a...\n";
		matrix* Y_HJ = solve_ode(df2, 0, 0.1, 100, Y0, Yref, k_HJ);
		zapiszSymulacje(Y_HJ, "symulacja_HJ.csv");

		// SYMULACJA DLA ROSENA
		cout << "Symulacja dla Rosenbrocka...\n";
		matrix* Y_Rosen = solve_ode(df2, 0, 0.1, 100, Y0, Yref, k_HJ);
		zapiszSymulacje(Y_Rosen, "symulacja_Rosen.csv");

	//Zapis symulacji do pliku csv
	//matrix Y0 = matrix(2, 1),								// Y0 zawiera warunki pocz�tkowe
	//	MT = matrix(2, new double[2] { m2d(opt.x), 0.5 });	// MT zawiera moment si�y dzia�aj�cy na wahad�o oraz czas dzia�ania
	//matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, NAN, MT);	// rozwi�zujemy r�wnanie r�niczkowe
	//ofstream Sout("testy_lab1.csv");						// definiujemy strumie� do pliku .csv
	//Sout << hcat(Y[0], Y[1]);								// zapisyjemy wyniki w pliku
	//Sout.close();											// zamykamy strumie�
	//Y[0].~matrix();											// usuwamy z pami�ci rozwi�zanie RR
	//Y[1].~matrix();

}

void lab3()
{

	auto zapiszSymulacje = [](matrix* Y, const string& nazwa_pliku) {
		ofstream plik(nazwa_pliku);
		if (plik.good()) {
			plik << "Czas[s]\tx[m]\tv_x[m/s]\ty[m]\tv_y[m/s]\n";
			int n = get_len(Y[0]);
			for (int i = 0; i < n; i++) {
				plik << Y[0](i) << "\t" << Y[1](i, 0) << "\t" << Y[1](i, 2) << "\n";
			}
			cout << "Zapisano symulacje: " << nazwa_pliku << " (" << n
				<< " punktow czasowych)\n";
		}
		plik.close();
		Y[0].~matrix();
		Y[1].~matrix();
		};

	srand(time(NULL));
	double s = 0.5;
	double alpha = 1, beta = 0.5, gamma = 1.5, delta = 0.5;
	double epsilon = 1e-5;
	solution::clear_calls();
	// std::cout << penSol << norm(penSol.x) << " g3:" << g3T3(x0, S_data(0)) << "
	// | " << g3T2(penSol.x,0) << " | " << g3T1(penSol.x,0);
	char kont = '1';
	fstream Sout;
	Sout.open("testy_lab3.csv", std::ios::out);
	matrix const_data = matrix(5, new double[5] {alpha, beta, gamma, delta, s});
	matrix ps(2, new double[2] {1.2, 1.2});

	// matrix x0(2, new double[2] {2, 5});
	// for (int i = 0; i < 1; i++) {
	// 	ps(0) = double(rand() % 20001 - 10000) / 1000.0;
	// 	ps(1) = double(rand() % 20001 - 10000) / 1000.0;
	// 	std::cout << "p0: " << ps << '\n';

	// 	solution solvedReal =
	// 		pen(ff3R, ps, 1, 1.5, epsilon, 1e5, const_data, const_data);
	// 	std::cout << solvedReal;
	// 	std::cout << ff3R(ps, ps, 0);
	// 	matrix Y0 = matrix(4, new double[4] {0, solvedReal.x(0), 100, 0});
	// 	matrix ud2(0);
	// 	matrix* Y = solve_ode(df3, 0, 0.01, 7.0, Y0, solvedReal.x, ud2);
	// 	zapiszSymulacje(Y, "pls_work_solved.csv");
	// 	int n = get_len(Y[0]);
	// 	int i0 = 0, i50 = 0;
	// 	for (int i = 0; i < n; i++) {
	// 		if (fabs(Y[1](i, 2) - 50) < fabs(Y[1](i50, 2) - 50)) {
	// 			i50 = i;
	// 		}
	// 		if (fabs(Y[1](i, 2)) < fabs(Y[1](i0, 2))) {
	// 			i0 = i;
	// 		}
	// 	}
	// 	std::cout << "x(i50): " << Y[1](i50, 0) << " y :" << Y[1](i50, 2) << '\n';
	// }

	// while (kont != '1') {
	// 	for (auto a : { 4.0, 4.4934, 5.0 }) {
	// 		matrix S_data = matrix(2, new double[2] {a, 0});
	// 		for (int i = 0; i < 100; i++) {
	// 			ps(0) = double(rand() % 10001 + 13000) / 10000.0;
	// 			ps(1) = double(rand() % 10001 + 13000) / 10000.0;

	// 			S_data(1) = 0; // ustawianie na zewn. testową kare
	// 			std::cout << "x0: \n" << ps << '\n';
	// 			std::cout << "\nFUNKCJA TESTOWA METODA SYN_NM Z KARĄ - TESTOWA ZEWN.\n";
	// 			if (Sout.good() == true)
	// 				Sout << ps(0) << "\t" << ps(1) << '\t';
	// 			solution penSol =
	// 				pen(ff3T_zewn, ps, 1, 1.5, epsilon, 1000000, S_data, const_data);
	// 			std::cout << penSol << norm(penSol.x) << '\n';
	// 			if (Sout.good() == true) {
	// 				Sout << penSol.x(0) << "\t" << penSol.x(1) << '\t' << norm(penSol.x)
	// 					<< '\t' << penSol.y(0) << '\t' << solution::f_calls << '\t';
	// 			}
	// 			solution::clear_calls();
	// 			std::cout << "\nFUNKCJA TESTOWA METODA SYN_NM Z KARĄ - TESTOWA WEWN.\n";
	// 			S_data(1) = 1; // ustawianie na wewn. testową kare
	// 			penSol =
	// 				pen(ff3T_wewn, ps, 1000, 0.1, epsilon, 1000000, S_data, const_data);
	// 			std::cout << penSol << norm(penSol.x) << '\n';
	// 			if (Sout.good() == true) {
	// 				Sout << penSol.x(0) << "\t" << penSol.x(1) << '\t' << norm(penSol.x)
	// 					<< '\t' << penSol.y(0) << '\t' << solution::f_calls << '\n';
	// 			}
	// 		}
	// 	}
	// 	Sout.close();

	// 	cin >> kont;
	// }

}

void lab4()
{
	//std::cout << ff4R(matrix(3, new double[3]{0.1, 0.1, 0.1})) << '\n';
	//std::cout << gf4R(matrix(3, new double[3]{0.1, 0.1, 0.1})) << '\n';
	//-----FUNKCJA TESTOWA-----------------------------------------------------------
	srand(time(NULL));
	char kont = '1';
	fstream Sout;
	matrix ps(2, 1);
	Sout.open("testy_lab4.csv", std::ios::out);
	solution opt;
	double epsilon = 1e-4;
	int limit = 1e6;
	double h0 = 0.05; // start step

	while (kont == '1') {
		for (int i = 0; i < 1; i++) {							//JG:mozna wybrac liczbe powtorzen

			ps = matrix(2, new double[2] {-0.544, -1.7704});

			//opt = SD(ff4T, gf4T, zlotf4T, ps, h0, epsilon, limit, 0, 0);
			//CG(ff4T, gf4T, ps, h0, epsilon, limit, 0, 0);
			//Newton(ff4T, gf4T, Hf4T, ps, h0, epsilon, limit, 0, 0);
			goto test;
			return;
		test:
			while (kont == '1') {
				for (int i = 0; i < 100; i++) {							//JG:mozna wybrac liczbe powtorzen

					ps(0) = double(rand() % 10001 - 20000) / 10000.0;
					ps(1) = double(rand() % 40001 - 20000) / 10000.0;
					cout << "Punkt startowy = [" << ps(0) << ", " << ps(1) << "].\n";
					cout << "Krok startowy = " << h0 << ".\n\n";

					cout << "SD:\n";
					if (Sout.good() == true) Sout << ps(0) << "\t" << ps(1) << '\t';
					try {
						opt = SD(ff4T, gf4T, zlotf4T, ps, h0, epsilon, limit, 0, 0);
					}
					catch (...) {
						Sout << "nan\tnan\tnan\t" << solution::f_calls << '\t' << solution::g_calls;
						Sout << "\tNIE\t"; solution::clear_calls(); goto cg;
					}
					cout << opt << endl << endl;

					if (Sout.good() == true) {
						Sout << opt.x(0) << "\t" << opt.x(1) << '\t' << opt.y(0)
							<< "\t" << solution::f_calls << "\t" << solution::g_calls;
						if (fabs(opt.y(0)) < 10 * epsilon) {
							Sout << "\tTAK\t";
						}
						else Sout << "\tNIE\t";
					}
					solution::clear_calls();
				cg:
					cout << "CG:\n";
					opt = CG(ff4T, gf4T, zlotf4T, ps, h0, epsilon, limit, 0, 0);
					cout << opt << endl << endl;

					if (Sout.good() == true) {
						Sout << opt.x(0) << "\t" << opt.x(1) << '\t' << opt.y(0)
							<< "\t" << solution::f_calls << "\t" << solution::g_calls;
						if (fabs(opt.y(0)) < 10 * epsilon) {
							Sout << "\tTAK\t";
						}
						else Sout << "\tNIE\t";
					}
					else throw "aaaa";
					solution::clear_calls();

					cout << "Newton:\n";
					opt = Newton(ff4T, gf4T, Hf4T, zlotf4T, ps, h0, epsilon, limit, 0, 0);
					cout << opt << endl << endl;

					if (Sout.good() == true) {
						Sout << opt.x(0) << "\t" << opt.x(1) << '\t' << opt.y(0)
							<< "\t" << solution::f_calls << "\t" << solution::g_calls << "\t" << solution::H_calls;
						if (fabs(opt.y(0)) < 10 * epsilon) {
							Sout << "\tTAK\n";
						}
						else Sout << "\tNIE\n";
					}
					else throw "aaaa";
					solution::clear_calls();

					std::cout << "loop: " << i << '\n';

				}

				std::cout << "koniec petli\n";
				cin >> kont;
			}
			return;
		real:
			Sout.close();
			Sout.open("real_lab4.csv", std::ios::out);
			while (kont == '1') {
				limit = 100000; // bo to powolne jest idk why
				for (auto h : { 0.01, 0.001, 0.0001 }) {
					solution::clear_calls();
					ps = matrix(3, new double[3] {
						double(rand() % 100001) / 10000.0,
							double(rand() % 100001) / 10000.0,
							double(rand() % 100001) / 10000.0
						});
					cout << "Punkt startowy = [" << ps(0) << ", " << ps(1) << ", " << ps(2) << "].\n";
					cout << "Krok startowy = " << h << ".\n\n";
					try {
						opt = SD(ff4R, gf4R, zlotf4R, ps, h, epsilon, limit, 0, 0);
						std::cout << opt;
						Sout << opt.x(0) << '\t'
							<< opt.x(1) << '\t'
							<< opt.x(2) << '\t'
							<< opt.y(0) << '\t'
							<< poprawne4R(opt.x) << '\t'
							<< solution::g_calls << '\n';
					}
					catch (...) {
						std::cout << "Blad dla kroku: " << h << '\n';
					}
				}
				std::cin >> kont;
			}
			Sout.close();
		}
	}
}

void lab5()
{


	srand(time(NULL));
	//Funkcja testowa
	double alfa = 0.0;										// wspolczynnik ekspansji (0.0 do 1.0 dla HJ i > 1 dla R - zmiana w petli, nie tutaj)
	double krok_s = 1.5;									// krok
	double beta = 0.5;										// wspolczynnik kontrakcji
	double epsilon = 1e-8;									// dokladnosc
	int Nmax = 100000;										// maksymalna liczba wywolan funkcji celu
	matrix lb(2, 1, -1.0), ub(2, 1, 1.0),					// dolne oraz g�rne ograniczenie
		ps(2, 1, double(rand() % 20001 - 10000) / 10000.0);	// punkt startowy
	solution opt;											// rozwiazanie optymalne znalezione przez algorytm
	solution::clear_calls();

	//-----FUNKCJA TESTOWA-----------------------------------------------------------

	char kont = '1';
	fstream Sout;
	Sout.open("testy_lab5.csv", std::ios::out);
	goto real5; // POMIJA TESTOWE; IDZIE DO RZECZYWISTEJ

	while (kont == '1') {
		for(double w = 0.00; w <= 1.01; w += 0.01) {
			ps(0) = double(rand() % 200001 - 100000) / 10000.0;
			ps(1) = double(rand() % 200001 - 100000) / 10000.0;
			matrix a = (matrix)(1.0);
			matrix wm = (matrix)(w);
			cout << "Punkt startowy = [" << ps(0) << ", " << ps(1) << "].\n";
			if (Sout.good() == true) {
				Sout << ps(0) << "\t" << ps(1) << '\t';
			}
			for(double a : {1.0, 10.0, 100.0}) {
				cout << "Dla współczynnika a = " << a << ":\n";
				matrix ud1(1, 2, 0.0);
				ud1(0) = a; ud1(0,1) = w;
				opt = Powell(gg5TX, ps, epsilon, Nmax, ud1, 0.0);
				opt.fit_fun(ff5TX, ud1, 0.0);
				cout << opt << endl << endl;
				exit(0);
				if (Sout.good() == true) {
					Sout << opt.x(0) << '\t' << opt.x(1) << '\t'
					<< ff5T1(opt.x, ud1, 0) << '\t' << ff5T2(opt.x, ud1, 0) << '\t'
					<< solution::f_calls << '\t';
				}
				solution::clear_calls();
			}
			if(Sout.good()) Sout << '\n';
		}
		cin >> kont;
	}
real5:
	Sout.close();
	Sout.open("rzeczywista_lab5.csv", std::ios::out);
	while (kont == '1') {
		for(double w = 0.00; w <= 1.01; w += 0.01) {
			ps(0) = double(rand() % 8000 + 2000) / 10000.0; // l
			ps(1) = double(rand() % 400 + 100) / 10000.0; // d
			cout << "Punkt startowy = [" << ps(0) << ", " << ps(1) << "].\n";
			if (Sout.good() == true) {
				Sout << 1000 * ps(0) << "\t" << 1000 * ps(1) << '\t';
			}
			for(double a : {1.0}) { // Nie chciało mi sie usuwać for-a.
				matrix ud1(1, 2, 0.0);
				ud1(0) = a; ud1(0,1) = w;
				opt = Powell(gg5RX, ps, epsilon, Nmax, ud1, 0.0);
				opt.fit_fun(ff5RX, ud1, 0.0);
				cout << opt << endl << endl;
				if (Sout.good() == true) {
					Sout << 1000 * opt.x(0) << '\t' << 1000 * opt.x(1) << '\t'
					<< ff5R1(opt.x, ud1, 0) << '\t' // masa
					<< 1000 * ff5R2(opt.x, ud1, 0) << '\t' // ugięcie
					<< 1e-6 * ff5R3(opt.x, ud1, 0) << '\t' // naprężenie
					<< solution::f_calls << '\t';
				}
				solution::clear_calls();
			}
			if(Sout.good()) Sout << '\n';
		}
		cin >> kont;
	}
	Sout.close();


}

void lab6()
{

}
