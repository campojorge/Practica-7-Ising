
from pathlib import Path
import re
import numpy as np
import matplotlib.pyplot as plt


# ============================================================
# CONFIGURACIÓN
# ============================================================

BASE = Path(__file__).resolve().parent

CARPETA_DATOS = BASE / "../data_results"
CARPETA_PLOTS = BASE / "../plot_results"

N = 64
J = 1.0
H = -0.3

# None: espacios o tabulaciones. Usa "," si separas con comas.
DELIMITADOR = None

# Pon 1 si hay una cabecera que no comienza por #.
SALTAR_FILAS = 0

NUM_BINS = 60


# ============================================================
# DETECTAR TEMPERATURAS DISPONIBLES
# ============================================================

def buscar_temperaturas():

    patron = (
        f"N={N}_J={J:.6f}_H={H:.6f}_kBT=*.csv"
    )

    temperaturas = []

    for archivo in CARPETA_DATOS.glob(patron):

        # Excluimos los archivos N-fold Way
        if archivo.name.startswith("N-FW_"):
            continue

        coincidencia = re.search(
            r"kBT=([-+]?\d+(?:\.\d+)?)\.csv$",
            archivo.name
        )

        if coincidencia:
            temperaturas.append(float(coincidencia.group(1)))

    return sorted(set(temperaturas))


# ============================================================
# LECTURA DE DATOS
# ============================================================

def leer_tiempos(archivo):

    datos = np.loadtxt(
        archivo,
        comments="#",
        delimiter=DELIMITADOR,
        skiprows=SALTAR_FILAS,
        ndmin=2,
    )

    if datos.size == 0 or datos.shape[1] != 4:
        raise ValueError(
            f"{archivo}: se esperan datos con cuatro columnas."
        )

    # Cuarta columna: tiempo en MCSS.
    tiempos = datos[:, 3]

    if not np.all(np.isfinite(tiempos)) or np.any(tiempos < 0):
        raise ValueError(
            f"{archivo}: hay tiempos MCSS no válidos."
        )

    return tiempos


# ============================================================
# GRÁFICA PARA UNA TEMPERATURA
# ============================================================

def plot_distribucion(kBT):

    NOMBRE = (
        f"N={N}_J={J:.6f}_H={H:.6f}_kBT={kBT:.6f}"
    )

    ARCHIVO_GLAUBER = CARPETA_DATOS / f"{NOMBRE}.csv"
    ARCHIVO_NFOLD = CARPETA_DATOS / f"N-FW_{NOMBRE}.csv"

    SALIDA = CARPETA_PLOTS / f"distribucion_MCSS_{NOMBRE}.png"

    if not ARCHIVO_GLAUBER.exists():
        print(f"No existe: {ARCHIVO_GLAUBER}")
        return

    if not ARCHIVO_NFOLD.exists():
        print(f"No existe: {ARCHIVO_NFOLD}")
        return

    tiempos_glauber = leer_tiempos(ARCHIVO_GLAUBER)
    tiempos_nfold = leer_tiempos(ARCHIVO_NFOLD)

    # Mismos intervalos para ambos histogramas.
    bins = np.histogram_bin_edges(
        np.concatenate([tiempos_glauber, tiempos_nfold]),
        bins=NUM_BINS,
    )

    fig, ejes = plt.subplots(
        1, 2,
        figsize=(12, 5),
        sharex=True,
        sharey=True,
        constrained_layout=True,
    )

    conjuntos = [
        ("Glauber", tiempos_glauber, "tab:blue"),
        ("N-fold way", tiempos_nfold, "tab:orange"),
    ]

    print(f"\n========== kBT = {kBT:.6f} ==========")

    for ax, (nombre, tiempos, color) in zip(ejes, conjuntos):

        media = np.mean(tiempos)

        # Varianza muestral.
        varianza = (
            np.var(tiempos, ddof=1)
            if len(tiempos) > 1
            else None
        )

        ax.hist(
            tiempos,
            bins=bins,
            density=True,
            color=color,
            edgecolor="white",
            linewidth=0.7,
            alpha=0.85,
        )

        # Recuadro con media y varianza.
        texto = f"Media = {media:.2f} MCSS\n"

        if varianza is not None:
            texto += f"Varianza = {varianza:.2f} MCSS²"
        else:
            texto += "Varianza: no disponible (n < 2)"

        ax.text(
            0.97, 0.95,
            texto,
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=10,
            bbox=dict(
                boxstyle="round",
                facecolor="white",
                edgecolor="gray",
                alpha=0.9,
            ),
        )

        ax.set_title(
            f"{nombre} — {len(tiempos)} simulaciones",
            fontsize=14,
            fontweight="bold",
        )

        ax.set_xlabel(
            r"Tiempo hasta alcanzar $M=0$, $\tau$ (MCSS)",
            fontsize=12,
        )

        ax.tick_params(axis="both", labelsize=11)
        ax.grid(axis="y", alpha=0.25)
        ax.set_axisbelow(True)
        ax.set_ylim(bottom=0)

        if varianza is not None:

            error_media = np.sqrt(varianza / len(tiempos))

            print(
                f"{nombre}: media = {media:.6g} "
                f"± {error_media:.3g} MCSS (error estándar)"
            )

            print(
                f"Varianza muestral = {varianza:.6g} MCSS²\n"
            )

        else:
            print(
                f"{nombre}: solo hay una simulación; "
                "no permite estimar una distribución."
            )

    ejes[0].set_ylabel(
        r"Densidad de probabilidad (MCSS$^{-1}$)",
        fontsize=12,
    )

    fig.suptitle(
        rf"Distribución de tiempos de llegada a $M=0$"
        rf" — $k_BT={kBT:g}$",
        fontsize=17,
        fontweight="bold",
    )

    fig.savefig(SALIDA, dpi=300)
    plt.close(fig)

    print(f"Gráfica guardada en: {SALIDA.resolve()}")


# ============================================================
# MAIN
# ============================================================

def main():

    CARPETA_PLOTS.mkdir(parents=True, exist_ok=True)

    temperaturas = buscar_temperaturas()

    if not temperaturas:
        print("No se encontraron archivos para los parámetros indicados.")
        return

    print("Temperaturas encontradas:", temperaturas)

    for kBT in temperaturas:
        plot_distribucion(kBT)

    print("\nTodas las gráficas han sido generadas.")


if __name__ == "__main__":
    main()
