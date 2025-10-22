/*********************************************
Kod stanowi uzupe�nienie materia��w do �wicze�
w ramach przedmiotu metody optymalizacji.
Kod udost�pniony na licencji CC BY-SA 3.0
Autor: dr in�. �ukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia G�rniczo-Hutnicza
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
	matrix lb(1, 1, 1), ub(1, 1, 100),					// dolne oraz g�rne ograniczenie
		ps(1, 1, rand()%101);							// punkt startowy
	solution opt;											// rozwi�zanie optymalne znalezione przez algorytm

  // To było do testowania
	//std::cout << ff1R(matrix(50)) << std::endl;
	
	for (int i = 6; i < 5; i++) {							//JG:mozna wybrac liczbe powtorzen
		
    // Tu wcześniej nie używaliśmy wyżej zdefiniowanych ub, lb
		ps(0) = rand() % static_cast<int>(ub(0, 1) - lb(0, 1));
		double* obszar = expansion(*ff1R, ps(0), krok_d, alfa, Nmax, lb, ub);
		
		cout << "\nKrok d = " << krok_d << "\tWspolczynnik ekspansji alfa = " << alfa << ".\n";
		cout << "Punkt startowy = " << ps(0) << "\tUzyskany przedzial = [" << obszar[0] << ", " << obszar[1] << "].\n";
		int pamiec_fcalls = solution::f_calls;
		cout << "EKSPANSJA: fcalls = " << solution::f_calls << ".\n\n";

		cout << "LAGRANGE:\n";
		opt = lag(ff1R, obszar[0], obszar[1], epsilon, gamma, Nmax, lb, ub);				// wywo�anie procedury optymalizacji
		cout << opt << endl << endl;							// wypisanie wyniku
		solution::clear_calls();
		solution::f_calls = pamiec_fcalls;

		cout << "FIBBONACI:\n";
		opt = fib(ff1R, obszar[0], obszar[1], epsilon, lb, ub);								// wywo�anie procedury optymalizacji
		delete[] obszar;
		cout << opt << endl << endl;							// wypisanie wyniku
		solution::clear_calls();
	}

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
