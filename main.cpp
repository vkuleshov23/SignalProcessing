#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <fstream>
#include <math.h>

using namespace std;

double T = 2.0;
double f1 = 15.50;
// double f1 = 1.50;
double f2 = -2.2;
double pi = 3.14159265359;

double eps = 1e-6;
bool isZero(double num){
	return abs(num) < eps;
}
double u(double t, double f){
	return sin(2*pi*f*t);
}
vector<pair<double, double>> sampling(double t1, double t2, double dt, double (*func)(double, double), double f){
	vector<pair<double, double>> vals((t2-t1)/dt);
	for (int i = 0; i < vals.size(); i++){
		vals[i].first = t1 + dt*i;
		vals[i].second = func(vals[i].first, f);
	}
	return vals;
}
vector<pair<double, double>> normalizing(double t1, double t2, double dt, double (*func)(double, double), double f, double norm){
	vector<pair<double, double>> vals((t2-t1)/dt);
	for (int i = 0; i < vals.size(); i++){
		vals[i].first = t1 + dt*i;
		vals[i].second = func(vals[i].first, f) / norm;
	}
	return vals;
}
double norm(vector<pair<double, double>> vals, double period, double dt){
	double norm = 0;
	for (pair<double, double> val : vals){
		norm += val.second * val.second;
	}
	// cout << "\n\tnorm: " << norm << '\n';
	norm = (norm / period) * dt ;
	return sqrt(norm);
}
bool isOrtNorm(double t1, double t2, double dt, double (*func)(double, double), double fr1, double fr2){
	
	vector<pair<double, double>> a = sampling(t1, t2, dt, func, fr1);
	vector<pair<double, double>> b = sampling(t1, t2, dt, func, fr2);

	double aNORM = norm(a, t2-t1, 0.001);
	cout << "\n\tOLD 1 Norm: " << aNORM;
	
	double bNORM = norm(b, t2-t1, 0.001);
	cout << "\n\tOLD 2 Norm: " << bNORM;

	cout << "\n\t\t Normalizing...";
	a = normalizing(t1, t2, dt, func, fr1, aNORM);
	b = normalizing(t1, t2, dt, func, fr2, bNORM);

	aNORM = norm(a, t2-t1, 0.001);
	cout << "\n\t\t1 Norm: " << aNORM;
	
	bNORM =  norm(b, t2-t1, 0.001);
	cout << "\n\t\t2 Norm: " << bNORM;

	double ort = 0;
	for(int i =0; i < a.size(); i++ ){
		ort += a[i].second * b[i].second;
	}
	cout << "\n\tort: " << ort << "\n";
	return (isZero(abs(aNORM)-1) && isZero(abs(bNORM)-1) && isZero(ort));
}

vector<pair<double, double>> kotelnikovRow(vector<pair<double, double>> func, double dt1, double dt2){
	vector<pair<double, double>> vals(func.size()*dt1/dt2);
	for(int i = 0; i < vals.size(); i++){
		double sum = 0;
		double t = i * dt2;
		for(int j = 0; j < func.size(); j++){
			double k = pi * (t/dt1 - j);
			sum += (k == 0 ? func[j].second : func[j].second * (sin(k)/k));
		}
		vals[i].first = t;
		vals[i].second = sum;
	}
	return vals;
}

int main(){
	ofstream fout;
	cout << "input frequency & T(f1 f2 T): ";
	cin >> f1 >> f2 >> T;
	cout << "f1: " << f1 << " f2: " << f2 << " T:"<< T << "\n";
	fout << fixed << setprecision(22);
	// cout << fixed << setprecision(22);
	cout << "first function counting...\n";
	vector<pair<double, double>> a = sampling(-T/2, T/2, 0.001, u, f1);
	fout.open("u1.txt");
	for(pair<double, double> val : a){
		fout << val.first << " " << val.second << '\n';
	} fout.close();
	cout << "end first.\n\n";
	system("./plot.sh u1.txt u-f1 2> /tmp/plotlog.log");

	cout << "second function counting...\n";
	vector<pair<double, double>> b = sampling(-T/2, T/2, 0.001, u, f2);
	fout.open("u2.txt");
	for(pair<double, double> val : b){
		fout << val.first << " " << val.second << '\n';
	} fout.close();
	cout << "end second.\n\n";
	system("./plot.sh u2.txt u-f2 2> /tmp/plotlog.log");

	cout << "approximate scalar product: ";
	double appScalProd = 0;
	for (int i = 0; i < a.size(); i++) {
		appScalProd += a[i].second * b[i].second;
	}
	cout << appScalProd << "\n\n";

	cout << "the norm of the 1 function: ";
	double aNorm = norm(a, T, 0.001);
	cout << aNorm << "\n\n";

	cout << "the norm of the 2 function: ";
	double bNorm = norm(b, T, 0.001);
	cout << bNorm << "\n\n";

	cout << "are these functions orthogonal?: ";
	
	cout << (isZero(appScalProd) ? "Yes" : "No") << "; appScalProd = " << appScalProd <<"\n\n";
	// cout 	<< "are the functions an orthonormal basis?: (-T/2, T/2, 0.001, u, 0.0065, 0.0065)" 
	// 		<< (isOrtNorm(-T/2, T/2, 0.001, u, 0.0065, 0.0065) ? "Yes" : "No") << "\n\n";
	cout 	<< "are the functions an orthonormal basis?: (-T/2, T/2, 0.001, u/norm, " << f1 << ", "<< f2 <<")" 
			<< (isOrtNorm(-T/2, T/2, 0.001, u, f1, f2) ? "Yes" : "No") << "\n\n";

	cout 	<< "are the functions an orthonormal basis?: (-T/2, T/2, 0.001, u/norm, " << f1*2 << ", " << f2*2 <<")" 
			<< (isOrtNorm(-T/2, T/2, 0.001, u, f1*2, f2*2) ? "Yes" : "No") << "\n\n";

	cout 	<< "are the functions an orthonormal basis?: (-T/4, T/4, 0.001, u/norm, " << f1<< ", " << f2 <<")" 
			<< (isOrtNorm(-T/4, T/4, 0.001, u, f1, f2) ? "Yes" : "No") << "\n\n";

	cout 	<< "are the functions an orthonormal basis?: (T, T, 0.001, u/norm, " << f1 << ", "<< f2 <<")" 
			<< (isOrtNorm(-T, T, 0.001, u, f1, f2) ? "Yes" : "No") << "\n\n";

	
	cout << "\n------------------------------\nKotelnikov Row\n------------------------------\n";
	cout << "input F: ";
	cin >> f1;
	cin.get();
	vector<double> freqs = {1.5, 1.75, 2, 3, 1000};

	for(int i = 0; i < freqs.size(); i++){
		vector<pair<double, double>> src = sampling(0, 10/f1, 1/(freqs[i]*f1), u, f1);
		fout.open("u1.txt");
		for(pair<double, double> val : src){
			fout << val.first << " " << val.second << "\n";
		} fout.close();
		system("./plot.sh u1.txt src 2> /tmp/plotlog.log");
		vector<pair<double, double>> dst = kotelnikovRow(src, 1/(freqs[i]*f1), 1/(3000*f1));

		fout.open("u2.txt");
		for(pair<double, double> val : dst){
			fout << val.first << " " << val.second << "\n";
		} fout.close();
		system("./plot.sh u2.txt dst 2> /tmp/plotlog.log");
		cin.get();
	}
}