"""Test della coda della binomiale negativa di gene_level_drivers_classes.py contro mpmath (60 cifre).
Uso: python test_nb_tail.py (gene_level_drivers_classes.py nella stessa cartella, ambiente con mpmath)"""
import ast, pathlib, numpy as np, mpmath as mp
from scipy import stats, special
src = pathlib.Path(__file__).with_name("gene_level_drivers_classes.py").read_text()
node = next(n for n in ast.parse(src).body if isinstance(n, ast.FunctionDef) and n.name == "nb_upper_score")
ns = {"np": np, "stats": stats, "special": special}; exec(compile(ast.Module([node], []), "f", "exec"), ns); f = ns["nb_upper_score"]
mp.mp.dps = 60
def ref(O, E, th):
    # somma DIRETTA della coda (k da O a O+3000): niente sottrazione 1 - CDF, che perde precisione per code minuscole
    if th == np.inf:
        tail = mp.fsum(mp.e ** (-mp.mpf(E)) * mp.mpf(E) ** k / mp.factorial(k) for k in range(O, O + 3000))
    else:
        th, E = mp.mpf(th), mp.mpf(E); p = th / (th + E)
        tail = mp.fsum(mp.gamma(k + th) / (mp.factorial(k) * mp.gamma(th)) * p ** th * (1 - p) ** k for k in range(O, O + 3000))
    return float(-mp.log10(tail))
cases = [(369, 5.4, 2.0), (150, 10.0, 5.0), (2, 0.05, 1.0), (10, 3.0, np.inf), (1, 0.5, 3.0), (60, 2.0, np.inf), (0, 2.0, 2.0)]
for O, E, th in cases:
    got = f(np.array([O], float), np.array([E]), th)[0]; exp = 0.0 if O == 0 else ref(O, E, th)
    print(f"O={O:4d} E={E:5.2f} theta={th}: script {got:10.4f}  mpmath {exp:10.4f}")
    assert abs(got - exp) < 1e-6 * max(1, exp), (O, E, th)
print("TEST SUPERATO")