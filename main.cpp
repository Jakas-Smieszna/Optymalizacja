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

#include"opt_alg.h"

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
		lab2();
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
	double krok_d = 1.0;									//krok/odleglosc do ekspansji
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
	double alfa = 0.1;										//wspolczynnik ekspansji (0.0 do 1.0)
	double krok_s = 1.0;									//krok
	double beta = 1e-6;										//kolejna dokladnosc
	double epsilon = 1e-2;									// dokladnosc
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

			cout << "HOOK-JEEVES:\n";
			opt = HJ(ff2T, ps, krok_s, alfa, epsilon, Nmax, lb, ub);							// wywołanie procedury optymalizacji
			cout << opt << endl << endl;														// wypisanie wyniku
			//if (Sout.good() == true) Sout << opt.x(0) << "\t" << opt.y(0) << "\t" << solution::f_calls << "\tlokalne\t";
			solution::clear_calls();

			cout << "ROSENBROCK:\n";
			opt = Rosen(ff2T, ps, matrix(2, 1, krok_s), alfa, beta, epsilon, Nmax, lb, ub);						// wywołanie procedury optymalizacji
			cout << opt << endl << endl;														// wypisanie wyniku
			//if (Sout.good() == true) Sout << opt.x(0) << "\t" << opt.y(0) << "\t" << solution::f_calls << "\tlokalne\n";
			solution::clear_calls();

		}

		cin >> kont;

	}
	Sout.close();

	//-----PROBLEM RZECZYWISTY-------------------------------------------------------



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

}

void lab4()
{

}

void lab5()
{

}

void lab6()
{

}
