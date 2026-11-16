#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

// Precalcula las probabilidades de giro p_j para las 10 clases (índices 0 a 9)
std::array<double, 10> calcular_probabilidades_glauber(double J, double H, double kBT) {
    std::array<double, 10> p{};
    const std::array<int, 5> vecinos_arriba = {4, 3, 2, 1, 0};

    for (int idx = 0; idx < 5; ++idx) {
        int u = vecinos_arriba[idx];
        int suma_vecinos = 2 * u - 4;

        // Clases 0..4: espín UP (+1) -> al girar pasa a DOWN (-1)
        double dE_up = 2.0 * (J * suma_vecinos + H);
        p[idx] = 1.0 / (1.0 + std::exp(dE_up / kBT));

        // Clases 5..9: espín DOWN (-1) -> al girar pasa a UP (+1)
        double dE_down = -2.0 * (J * suma_vecinos + H);
        p[idx + 5] = 1.0 / (1.0 + std::exp(dE_down / kBT));
    }
    return p;
}

// Precalcula los 4 vecinos de cada sitio i en una red LxL periódica (arreglo plano N*4)
std::vector<int> construir_vecinos(int L) {
    int N = L * L;
    std::vector<int> vecinos(N * 4);
    for (int r = 0; r < L; ++r) {
        for (int c = 0; c < L; ++c) {
            int i = r * L + c;
            vecinos[i * 4 + 0] = ((r - 1 + L) % L) * L + c; // Arriba
            vecinos[i * 4 + 1] = ((r + 1) % L) * L + c;     // Abajo
            vecinos[i * 4 + 2] = r * L + ((c - 1 + L) % L); // Izquierda
            vecinos[i * 4 + 3] = r * L + ((c + 1) % L);     // Derecha
        }
    }
    return vecinos;
}

// Ejecuta una simulación hasta M = 0 y devuelve el tiempo de escape en MCSS
double simular_escape_nfold(int L,
                            const std::array<double, 10>& p,
                            const std::vector<int>& vecinos,
                            std::mt19937_64& rng) {
    const int N = L * L;

    // Inicialización: todos los N espines UP en la clase 0 (clase 1 del texto)
    std::array<int, 10> nncls{};
    nncls[0] = N;

    std::vector<int> loc(10 * N, 0);
    std::iota(loc.begin(), loc.begin() + N, 0);

    std::vector<int> looc(N, 0);
    std::vector<int> look(N);
    std::iota(look.begin(), look.end(), 0);

    int M = N;
    int64_t t = 0; // Tiempo acumulado en intentos de Monte Carlo (mcs)
    std::uniform_real_distribution<double> dist_uniforme(0.0, 1.0);

    auto mover_spin = [&](int spin_i, int nueva_clase) {
        int clase_vieja = looc[spin_i];
        int pos_vieja = look[spin_i];

        int ultima_pos = nncls[clase_vieja] - 1;
        int ultimo_spin = loc[clase_vieja * N + ultima_pos];
        loc[clase_vieja * N + pos_vieja] = ultimo_spin;
        look[ultimo_spin] = pos_vieja;
        nncls[clase_vieja]--;

        int nueva_pos = nncls[nueva_clase];
        loc[nueva_clase * N + nueva_pos] = spin_i;
        looc[spin_i] = nueva_clase;
        look[spin_i] = nueva_pos;
        nncls[nueva_clase]++;
    };

    std::array<double, 10> Q{};
    while (M > 0) {
        double r_tilde = dist_uniforme(rng);
        double r_bar = 1.0 - dist_uniforme(rng); // Evita ln(0)

        // Probabilidades acumuladas Q_i (Ec. 42)
        double acumulado = 0.0;
        for (int j = 0; j < 10; ++j) {
            acumulado += (static_cast<double>(nncls[j]) * p[j]) / N;
            Q[j] = acumulado;
        }
        double Q10 = Q[9];

        // Tiempo de escape m de la configuración actual (Ec. 44)
        int64_t m = static_cast<int64_t>(std::floor(std::log(r_bar) / std::log1p(-Q10))) + 1;

        // Seleccionar clase k tal que Q_{k-1} <= r_tilde * Q10 < Q_k (Ec. 45)
        double objetivo = r_tilde * Q10;
        int k = 0;
        for (int j = 0; j < 10; ++j) {
            if (nncls[j] > 0) {
                k = j;
                if (objetivo < Q[j]) break;
            }
        }

        // Elegir al azar uno de los n_k espines de la clase k
        std::uniform_int_distribution<int> dist_spin(0, nncls[k] - 1);
        int spin_elegido = loc[k * N + dist_spin(rng)];

        // Actualizar clases del espín elegido y de sus 4 vecinos
        if (k < 5) {
            mover_spin(spin_elegido, k + 5);
            for (int v = 0; v < 4; ++v) {
                int nb = vecinos[spin_elegido * 4 + v];
                mover_spin(nb, looc[nb] + 1);
            }
            M -= 2;
        } else {
            mover_spin(spin_elegido, k - 5);
            for (int v = 0; v < 4; ++v) {
                int nb = vecinos[spin_elegido * 4 + v];
                mover_spin(nb, looc[nb] - 1);
            }
            M += 2;
        }

        t += m;
    }

    // Tiempo físico en Monte Carlo Steps per Spin (MCSS = t / N)
    return static_cast<double>(t) / N;
}

int main() {
    int L = 32;        // Red de 32x32 (N = 1024 espines)
    double J = 1.0;    // Constante de acoplamiento
    double H = -3.0;   // Campo magnético externo (|H| > 2J)
    double kBT = 0.5;  // Baja temperatura
    int K = 1000;      // Número de escapes para promediar

    // Precalcular tablas estáticas una sola vez fuera del bucle de K escapes
    const auto p = calcular_probabilidades_glauber(J, H, kBT);
    const auto vecinos = construir_vecinos(L);

    std::mt19937_64 rng(42);
    double suma_tiempos = 0.0;

    for (int i = 0; i < K; ++i) {
        suma_tiempos += simular_escape_nfold(L, p, vecinos, rng);
    }

    std::cout << "Tiempo medio de escape (K=" << K << "): "
              << (suma_tiempos / K) << " MCSS\n";

    return 0;
}