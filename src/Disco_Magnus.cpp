#include <iostream>
#include <cmath>
#include "vector.h"
#include <random>
using namespace std;

//* ------------------------------ Constantes globales ------------------------------
// Constantes Globales
const double Gamma = 16.0;
const double dt = 1e-3;
const double alpha = 1 - exp(-Gamma * dt);
const int N = 1;

// Constantes para fuerzas
const double Cd = 0.47; // Coeficiente de arrastre (esfera)
const double rho = 1.225; // Densidad del aire (kg/m^3)
const double S = 1.0; // Constante de Magnus

//* ------------------------------ Declaraion de clases, metodos de clases y funciones globales ------------------------------
// Declaracion de clases

class Cuerpo;

// Implementacion de las clases

class Cuerpo
{
private:
    vector3D r, v, F;
    double m, R, omega;

public:
    void Inicie(double x0, double y0, double vx0, double vy0, double omega0, double m0, double R0);
    void CalculateFuerzas(void);
    void Muevase(double dt, double KT, mt19937 &gen, normal_distribution<> &dis);
    void Arranque(double dt);
    void Dibujese(void);
    double get_x() { return r.x(); }
    double get_y() { return r.y(); }
    double get_z() { return r.z(); }
    double get_vx() { return v.x(); }
    double get_vy() { return v.y(); }
    double get_vz() { return v.z(); }
    double get_R() { return R; }
};

// Declaracion de funciones globales

void InicieAnimacion(const std::string& str);
void InicieCuadro(void);
void TermineCuadro(void);
bool stringToBool(const std::string& str);
void Imprimase(double t, Cuerpo *Polen, bool gnuplot);

//* ------------------------------ main ------------------------------

int main(int argc, char *argv[]) 
{
    std::string input = argv[1];
    bool gnu = stringToBool(input);

    double r0 = 0.0;
    double KT = 4.0;
    double v0 = r0;
    double m0 = 1.0;
    double R0 = 1.0;
    double omega0 = 100.0;

    double Ncuadros = 5000; //TODO: Revisar impresion, a 5000 se etsan generadndo menos datos (4167 en realidad)
    double dt = 1e-3, T = 25, tcuadro = T / Ncuadros;

    random_device rd;
    mt19937 rng(rd()); // Generador con semilla aleatoria
    // mt19937 rng(1); // Generador con semilla fija
    normal_distribution<> dis(0.0, 1.0);

    Cuerpo Polen[N];

    for (int i = 0; i < N; i++) Polen[i].Inicie(r0, 0.0, v0, 0.0, 0.0, m0, R0);

    for (int i = 0; i < N; i++) Polen[i].CalculateFuerzas();

    for (int i = 0; i < N; i++) Polen[i].Arranque(dt);

    for (double t = 0.0, tdibujo = 0.0; t < T + dt; t += dt, tdibujo += dt)
    {
        if (tdibujo > tcuadro || t == 0)
        {
            Imprimase(t, Polen, gnu);
            tdibujo = 0;
        }

        for (int i = 0; i < N; i++) Polen[i].CalculateFuerzas();
        
        for (int i = 0; i < N; i++) Polen[i].Muevase(dt, KT, rng, dis);
        
    }
    return 0;
}

//* ------------------------------ Implementacion de metodos de clases y funciones globales ------------------------------
// implementacion de los metodos de clases

void Cuerpo::Inicie(double x0, double y0, double vx0, double vy0, double omega0, double m0, double R0)
{
    r.load(x0, y0, 0.0);
    v.load(vx0, vy0, 0.0);
    omega = omega0;
    m = m0;
    R = R0;
}

void Cuerpo::CalculateFuerzas(void)
{
    double A = M_PI * R * R;
    vector3D w;
    w.load(0, 0, omega);

    // borrar fuerza 
    F.load(0, 0, 0);

    // Agrergar fuerzas

    // Fuerza de arrastre
    // esta fuerza esta incluida implicitamente en la ecuacion de movimiento

    // Fuerza de Magnus
    vector3D Fm = -0.5 * Cd * rho * A * R * (w ^ v);

    // Sumando las fuerzas
    F += Fm;

}

void Cuerpo::Muevase(double dt, double KT, mt19937 &gen, normal_distribution<> &dis)
{
    // algoritmo de leapfrog
    vector3D vp = v + F * dt;
    vector3D epsilon;
    epsilon.load(dis(gen), dis(gen), 0);
    vector3D deltav = -alpha * vp + sqrt(alpha * (2 - alpha) * (KT / m)) * epsilon;
    r += (vp + deltav * 0.5) * dt;
    v = vp + deltav;
}

void Cuerpo::Arranque(double dt)
{
    v -= F * dt / (2 * m);
}

void Cuerpo::Dibujese(void) 
{
    cout << " , " << r.x() << "+" << R << "*cos(t)," << r.y() << "+" << R << "*sin(t)";
}

// implementacion de funciones globales

void InicieAnimacion(const std::string& str) 
{
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

void InicieCuadro(void) 
{
    cout << "plot 0,0 ";
}

void TermineCuadro(void) 
{
    cout << endl;
}

bool stringToBool(const std::string& str) 
{
    return (str == "true" || str == "1");
}

void Imprimase(double t, Cuerpo *Polen, bool gnuplot) 
{
    ////cout << gnuplot << endl;
    if (gnuplot) 
    {
        if (t == 0) InicieAnimacion("mag.gif"); //! Cambiar nombre del gif cuando se vaya a ejecutar
        InicieCuadro();
        for (int i = 0; i < N; i++) Polen[i].Dibujese();
        TermineCuadro();
    }
    else
    {
        for (int i = 0; i < N; i++) 
        {
            cout << t << " " << i << " " << Polen[i].get_x() << " " << Polen[i].get_y() << " " << Polen[i].get_vx() << " " << Polen[i].get_vy() << endl;
        }
    }
}
