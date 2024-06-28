#include <iostream>
#include <cmath>
#include "Random64.h"

using namespace std;

// Constantes Globales
const double Gamma = 16;
const double Deltat = 0.001;
const double alpha = 1 - exp(-Gamma * Deltat);
const double alphaAux = alpha * (2 - alpha);
const int N = 1;

// Constantes para fuerzas
const double Cd = 0.47; // Coeficiente de arrastre (esfera)
const double rho = 1.225; // Densidad del aire (kg/m^3)
const double S = 1.0; // Constante de Magnus


class Cuerpo {
private:
  double x, y, Vx, Vy, Fx, Fy, omega; 
  double m, R;
public:
  void Inicie(double x0, double y0, double Vx0, double Vy0, double omega0, double m0, double R0);
  void CalculeFuerza(bool Mag);
  void Arranque(double dt);
  void Muevase(double dt, double kT, Crandom &ran64);
  void Dibujese(void);
  double Getx(void) { return x; } 
  double GetVx(void) { return Vx; }
  double GetVy(void) { return Vy; }
};

// Implementación de las funciones
void Cuerpo::Inicie(double x0, double y0, double Vx0, double Vy0, double omega0, double m0, double R0) {
  x = x0; y = y0; Vx = Vx0; Vy = Vy0; omega = omega0; m = m0; R = R0;
}

void Cuerpo::CalculeFuerza(bool Mag) {
  double v = sqrt(Vx * Vx + Vy * Vy);
  double area = M_PI * R * R;

  // Fuerza de arrastre (drag)
  double Fdragx = -0.5 * Cd * rho * area * v * Vx;
  double Fdragy = -0.5 * Cd * rho * area * v * Vy;

  // Fuerza de Magnus
  double Fmagnusx = S * omega * Vy;
  double Fmagnusy = -S * omega * Vx;

  // Sumando las fuerzas
  Fx = Fdragx; //+ Fmagnusx;
  Fy = Fdragy; // + Fmagnuy;
  if (Mag){
     Fx = Fx+Fmagnusx;
     Fy = Fy+Fmagnusy;
  }
}

void Cuerpo::Arranque(double dt) {
  // Algoritmo LeapFrog
  Vx -= Fx * (0.5 * dt / m);
  Vy -= Fy * (0.5 * dt / m);
}

void Cuerpo::Muevase(double dt, double kT, Crandom &ran64) {
  // Algoritmo Browniano de Leap Frog (van Gunsteren)
  double Vxprime = Vx + Fx * dt;
  double Vyprime = Vy + Fy * dt;
  double DeltaVx = -alpha * Vxprime + sqrt(alphaAux * kT / m) * ran64.gauss(0, 1);
  double DeltaVy = -alpha * Vyprime + sqrt(alphaAux * kT / m) * ran64.gauss(0, 1);
  x += (Vxprime + DeltaVx / 2) * dt;
  y += (Vyprime + DeltaVy / 2) * dt;
  Vx = Vxprime + DeltaVx;
  Vy = Vyprime + DeltaVy;
}

void Cuerpo::Dibujese(void) {
  cout << " , " << x << "+" << R << "*cos(t)," << y << "+" << R << "*sin(t)";
}

// Funciones Globales
void InicieAnimacion(const std::string& str) {
  cout << "set terminal gif animate" << endl;
  cout << "set output '"<< str <<"'" << endl;
  cout << "unset key" << endl;
  cout << "set xrange[-5:5]" << endl;
  cout << "set yrange[-5:5]" << endl;
  cout << "set size ratio -1" << endl;
  cout << "set parametric" << endl;
  cout << "set trange [0:7]" << endl;
  cout << "set isosamples 12" << endl;
}

void InicieCuadro(void) {
  cout << "plot 0,0 ";
}

void TermineCuadro(void) {
  cout << endl;
}

double Sigma2(Cuerpo *Polen) {
  double suma = 0, xprom = 0, sigma2 = 0;
  // Calculo xprom
  for (int i = 0; i < N; i++)
    suma += Polen[i].Getx();
  xprom = suma / N;
  // Calculo sigma2
  suma = 0;
  for (int i = 0; i < N; i++)
    suma += pow(Polen[i].Getx() - xprom, 2.0);
  //sigma2 = suma / (N - 1);
  sigma2=0;
  return sigma2;
}

bool stringToBool(const std::string& str) {
    return (str == "true" || str == "1");
}

int main(int argc, char *argv[]) {
  std::string input = argv[1];
  std::string str = argv[2];
  bool Mag = stringToBool(input);
  
  Cuerpo Polen[N];
  Crandom ran64(1);
  double x0 = 0, y0 = 0, Vx0 = 0, Vy0 = 0, omega0 = 100, R0 = 1, m0 = 1;
  double kT = 4.0;
  double t, ttotal = 20;
  int Ncuadros = 1000;
  double tdibujo, tcuadro = ttotal / Ncuadros;

  InicieAnimacion(str);

  // Inicialización de cuerpos
  for (int i = 0; i < N; i++) Polen[i].Inicie(x0, y0, Vx0, Vy0, omega0, m0, R0);
  for (t = tdibujo = 0; t < ttotal; t += Deltat, tdibujo += Deltat) {
    for (int i = 0; i < N; i++) Polen[i].CalculeFuerza(Mag);
    for (int i = 0; i < N; i++) Polen[i].Arranque(Deltat);
    if (tdibujo > tcuadro) {
      InicieCuadro();
      for (int i = 0; i < N; i++) Polen[i].Dibujese();
      TermineCuadro();
      tdibujo = 0;
    }
    for (int i = 0; i < N; i++) Polen[i].CalculeFuerza(Mag);
    for (int i = 0; i < N; i++) Polen[i].Muevase(Deltat, kT, ran64);
  }

  return 0;
}

