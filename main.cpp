/*********************************************
Kod stanowi uzupe³nienie materia³ów do æwiczeñ
w ramach przedmiotu metody optymalizacji.
Kod udostêpniony na licencji CC BY-SA 3.0
Autor: dr in¿. £ukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia Górniczo-Hutnicza
Data ostatniej modyfikacji: 30.09.2025
*********************************************/

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
		lab1();
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
	double epsilon = 1e-2;									// dok³adnoœæ
	int Nmax = 10000;										// maksymalna liczba wywo³añ funkcji celu
	matrix lb(2, 1, -5), ub(2, 1, 5),						// dolne oraz górne ograniczenie
		a(2, 1);											// dok³adne rozwi¹zanie optymalne
	solution opt;											// rozwi¹zanie optymalne znalezione przez algorytm
	a(0) = -1;
	a(1) = 2;
	opt = MC(ff0T, 2, lb, ub, epsilon, Nmax, a);			// wywo³anie procedury optymalizacji
	cout << opt << endl << endl;							// wypisanie wyniku
	solution::clear_calls();								// wyzerowanie liczników

	//Wahadlo
	Nmax = 1000;											// dok³adnoœæ
	epsilon = 1e-2;											// maksymalna liczba wywo³añ funkcji celu
	lb = 0, ub = 5;											// dolne oraz górne ograniczenie
	double teta_opt = 1;									// maksymalne wychylenie wahad³a
	opt = MC(ff0R, 1, lb, ub, epsilon, Nmax, teta_opt);		// wywo³anie procedury optymalizacji
	cout << opt << endl << endl;							// wypisanie wyniku
	solution::clear_calls();								// wyzerowanie liczników

	//Zapis symulacji do pliku csv
	matrix Y0 = matrix(2, 1),								// Y0 zawiera warunki pocz¹tkowe
		MT = matrix(2, new double[2] { m2d(opt.x), 0.5 });	// MT zawiera moment si³y dzia³aj¹cy na wahad³o oraz czas dzia³ania
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, NAN, MT);	// rozwi¹zujemy równanie ró¿niczkowe
	ofstream Sout("symulacja_lab0.csv");					// definiujemy strumieñ do pliku .csv
	Sout << hcat(Y[0], Y[1]);								// zapisyjemy wyniki w pliku
	Sout.close();											// zamykamy strumieñ
	Y[0].~matrix();											// usuwamy z pamiêci rozwi¹zanie RR
	Y[1].~matrix();
}

void lab1()
{
	srand(time(NULL));
	//Funkcja testowa
	double alfa = 1.5;										//wspolczynnik ekspansji
	double krok_d = 1.0;									//krok/odleglosc do ekspansji
	double gamma = 1e-2;									//kolejna dokladnosc
	double epsilon = 1e-2;									// dok³adnoœæ
	int Nmax = 10000;										// maksymalna liczba wywo³añ funkcji celu
	matrix lb(1, 1, -100), ub(1, 1, 100),					// dolne oraz górne ograniczenie
		ps(1, 1, rand()%201 - 100);							// punkt startowy
	solution opt;											// rozwi¹zanie optymalne znalezione przez algorytm
	
	for (int i = 0; i < 5; i++) {							//JG:mozna wybrac liczbe powtorzen
		
		ps(0) = rand() % 201 - 100;
		double* obszar = expansion(*ff1T, ps(0), krok_d, alfa, Nmax, lb, ub);
		
		cout << "\nKrok d = " << krok_d << "\tWspolczynnik ekspansji alfa = " << alfa << ".\n";
		cout << "Punkt startowy = " << ps(0) << "\tUzyskany przedzial = [" << obszar[0] << ", " << obszar[1] << "].\n";
		int pamiec_fcalls = solution::f_calls;				//JG:pozwala przywrocic liczbe po ekspansji, gdy zostanie wycyszczona.
		cout << "EKSPANSJA: fcalls = " << solution::f_calls << ".\n\n";

		cout << "LAGRANGE:\n";
		opt = lag(ff1T, obszar[0], obszar[1], epsilon, gamma, Nmax, lb, ub);				// wywo³anie procedury optymalizacji
		cout << opt << endl << endl;							// wypisanie wyniku
		solution::clear_calls();
		solution::f_calls = pamiec_fcalls;

		cout << "FIBBONACI:\n";
		opt = fib(ff1T, obszar[0], obszar[1], epsilon, lb, ub);								// wywo³anie procedury optymalizacji
		delete[] obszar;
		cout << opt << endl << endl;							// wypisanie wyniku
		solution::clear_calls();
	}

	//Zapis symulacji do pliku csv
	//matrix Y0 = matrix(2, 1),								// Y0 zawiera warunki pocz¹tkowe
	//	MT = matrix(2, new double[2] { m2d(opt.x), 0.5 });	// MT zawiera moment si³y dzia³aj¹cy na wahad³o oraz czas dzia³ania
	//matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, NAN, MT);	// rozwi¹zujemy równanie ró¿niczkowe
	//ofstream Sout("testy_lab1.csv");						// definiujemy strumieñ do pliku .csv
	//Sout << hcat(Y[0], Y[1]);								// zapisyjemy wyniki w pliku
	//Sout.close();											// zamykamy strumieñ
	//Y[0].~matrix();											// usuwamy z pamiêci rozwi¹zanie RR
	//Y[1].~matrix();

}

void lab2()
{

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
