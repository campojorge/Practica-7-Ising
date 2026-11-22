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

# None para espacios/tabulaciones; "," para comas.
DELIMITADOR = None

# Usa 1 si la cabecera no empieza por #.
SALTAR_FILAS = 0

MOSTRAR_GRAFICAS = False


# ============================================================
# PARÁMETROS DEL NOMBRE DEL FICHERO
# ============================================================

def parametros(nombre):

    numero = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"

    pares = re.findall(
        rf"(?:^|_)(N|J|H|kBT)=({numero})(?=_|$)",
        nombre
    )

    p = {
        clave: float(valor)
        for clave, valor in pares
    }

    if not all(clave in p for clave in ("N", "J", "H", "kBT")):
        raise ValueError(
            "El nombre debe incluir N, J, H y kBT"
        )

    if p["N"] < 2 or not p["N"].is_integer():
        raise ValueError("N debe ser entero >= 2")

    return p


# ============================================================
# LECTURA, RECORTE Y PROMEDIO
# ============================================================

def leer_correlacion(archivo):

    datos = np.loadtxt(
        archivo,
        comments="#",
        delimiter=DELIMITADOR,
        skiprows=SALTAR_FILAS,
        ndmin=2
    )

    if datos.size == 0 or datos.shape[1] != 3:
        raise ValueError(
            "Se esperan 3 columnas: simulacion, t_MCSS, Q"
        )

    if not np.isfinite(datos).all():
        raise ValueError("Hay valores NaN o infinitos")

    if np.any(datos[:, 1] < 0):
        raise ValueError("Hay tiempos negativos")

    if np.any(np.abs(datos[:, 2]) > 1 + 1e-8):
        raise ValueError(
            "La correlación debe estar entre -1 y 1"
        )

    simulaciones = []

    for sim in np.unique(datos[:, 0]):

        filas = datos[datos[:, 0] == sim]
        filas = filas[np.argsort(filas[:, 1])]

        tiempos = np.round(filas[:, 1], decimals=9)

        if np.any(np.diff(tiempos) <= 0):
            raise ValueError(
                f"Tiempos repetidos o demasiado próximos "
                f"en la simulación {sim:g}"
            )

        simulaciones.append((tiempos, filas[:, 2]))

    # Menor tiempo final de las simulaciones del fichero.
    t_final = min(t[-1] for t, _ in simulaciones)

    # Conservar únicamente tiempos presentes en todas.
    comunes = simulaciones[0][0]
    comunes = comunes[comunes <= t_final]

    for t, _ in simulaciones[1:]:
        comunes = np.intersect1d(
            comunes,
            t[t <= t_final]
        )

    if comunes.size == 0:
        raise ValueError(
            "No hay tiempos comunes a todas las simulaciones"
        )

    curvas = np.asarray([
        q[np.searchsorted(t, comunes)]
        for t, q in simulaciones
    ])

    n = len(curvas)
    media = curvas.mean(axis=0)

    if n > 1:
        error = curvas.std(axis=0, ddof=1) / np.sqrt(n)
    else:
        error = None

    print(
        f"  {archivo.name}: "
        f"{n} simulaciones, {len(comunes)} tiempos comunes, "
        f"hasta {comunes[-1]:g} MCSS"
    )

    return comunes, media, error, n


# ============================================================
# DIBUJAR UN MÉTODO EN SU PANEL
# ============================================================

def dibujar(ax, archivo, metodo, p):

    ax.set_title(metodo)
    ax.set_xlabel("Tiempo (MCSS)")
    ax.grid(alpha=0.25)

    if not archivo.is_file():

        ax.text(
            0.5, 0.5,
            "Fichero no encontrado",
            ha="center",
            transform=ax.transAxes
        )

        print(f"Falta: {archivo.name}")
        return False

    try:
        t, media, error, n = leer_correlacion(archivo)

    except (ValueError, OSError) as exc:

        ax.text(
            0.5, 0.5,
            "No se pudieron leer los datos\nConsulta la terminal",
            ha="center",
            transform=ax.transAxes
        )

        print(f"ERROR en {archivo.name}: {exc}")
        return False

    color = (
        "tab:blue"
        if metodo == "Glauber"
        else "tab:orange"
    )

    ax.plot(
        t,
        media,
        "o-",
        markersize=3,
        color=color,
        label=f"Simulación: {n} realizaciones"
    )

    if error is not None:

        ax.fill_between(
            t,
            media - error,
            media + error,
            color=color,
            alpha=0.2,
            label="Media ± error estándar"
        )

    # Teoría solo para J=H=0 y temperatura positiva.
    if p["J"] == 0 and p["H"] == 0 and p["kBT"] > 0:

        if metodo == "n-fold Way":

            # Teoría de tiempo continuo.
            tt = np.linspace(0, t.max(), 600)

            ax.plot(
                tt,
                np.exp(-tt),
                "k--",
                label=r"Teórica: $e^{-t}$"
            )

        else:

            # Teoría discreta de Glauber.
            pasos = p["N"] * t

            pasos = np.where(
                np.isclose(
                    pasos,
                    np.rint(pasos),
                    rtol=0,
                    atol=1e-8
                ),
                np.rint(pasos),
                np.floor(pasos)
            )

            teorica = np.exp(
                pasos * np.log1p(-1 / p["N"])
            )

            ax.plot(
                t,
                teorica,
                "k--",
                label=r"Teórica: $(1-1/N)^{\lfloor Nt\rfloor}$"
            )

    if len(t) == 1:

        ax.text(
            0.5, 0.5,
            "Solo hay un tiempo común",
            ha="center",
            transform=ax.transAxes
        )

    ax.legend(fontsize=9)

    return True


# ============================================================
# EMPAREJAR FICHEROS Y GUARDAR LOS PNG
# ============================================================

def main():

    grupos = {}

    # Nombres esperados:
    # correlacion_N=..._J=..._H=..._kBT=....csv
    # correlacion_N-FW_N=..._J=..._H=..._kBT=....csv

    for metodo, prefijo in (
        ("Glauber", "correlacion_"),
        ("n-fold Way", "correlacion_N-FW_")
    ):

        archivos = sorted(
            CARPETA_DATOS.glob(prefijo + "N=*.csv")
        )

        for archivo in archivos:

            nombre = archivo.stem[len(prefijo):]

            grupos.setdefault(nombre, {})[metodo] = archivo

    if not grupos:
        raise FileNotFoundError(
            "No hay ficheros de correlación en "
            f"{CARPETA_DATOS.resolve()}"
        )

    CARPETA_PLOTS.mkdir(parents=True, exist_ok=True)

    for nombre, archivos in sorted(grupos.items()):

        try:
            p = parametros(nombre)

        except ValueError as exc:
            print(f"ERROR en {nombre}: {exc}")
            continue

        fig, axes = plt.subplots(
            1, 2,
            figsize=(13, 5),
            sharex=True,
            sharey=True,
            constrained_layout=True
        )

        validos = []

        for ax, metodo, prefijo in zip(
            axes,
            ("Glauber", "n-fold Way"),
            ("correlacion_", "correlacion_N-FW_")
        ):

            archivo = archivos.get(
                metodo,
                CARPETA_DATOS / f"{prefijo}{nombre}.csv"
            )

            validos.append(
                dibujar(ax, archivo, metodo, p)
            )

        axes[0].set_ylabel(
            r"$C(t)=\langle N^{-1}\sum_i"
            r"\sigma_i(0)\sigma_i(t)\rangle$"
        )

        titulo = ", ".join(
            f"{clave}={p[clave]:g}"
            for clave in ("N", "J", "H", "kBT")
        )

        fig.suptitle(
            f"Correlación de espines\n{titulo}"
        )

        if any(validos):

            salida = (
                CARPETA_PLOTS
                / f"comparacion_correlacion_{nombre}.png"
            )

            fig.savefig(salida, dpi=300)
            print(f"Guardado: {salida}")

            if MOSTRAR_GRAFICAS:
                plt.show()

        plt.close(fig)


if __name__ == "__main__":
    main()