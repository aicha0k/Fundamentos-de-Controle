import control as ctrl
import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

def criar_sistema(G_num, G_den, feedback=1):
    G = ctrl.TransferFunction(G_num, G_den)
    T = ctrl.feedback(G, feedback)
    return G,T

def resposta_ao_degrau(T, tempo_final=10):
    t, y = ctrl.step_response(T, np.linspace(0, tempo_final, 1000))
    return t, y

def resposta_ao_degrau_funcao(T):
    s,t = sp.symbols('s t', real=True, positive=True)
    # num, den = T.num[0][0], T.den[0][0]
    # G_s = sp.Poly(num, s) / sp.Poly(den, s).as_expr()
    # Ys = G_s / s
    # yt = sp.inverse_laplace_transform(Ys, s, t)
    # return sp.simplify(yt)
    num = sp.Poly(T.num[0][0], s).as_expr()
    den = sp.Poly(T.den[0][0], s).as_expr()

    G_s = num / den
    Ys = G_s / s
    yt = sp.inverse_laplace_transform(Ys, s, t)
    yt_simplified = sp.simplify(sp.expand(yt))

    latex_expr = sp.latex(yt_simplified)
    return yt_simplified, latex_expr

def desempenho_sistema(T):
    info = ctrl.step_info(T)
    Mp = info['Overshoot']
    ts = info['SettlingTime']
    return Mp, ts, info

def sistema_P43(K):
    num = [5*K]
    den = [1, 15, K]
    T = ctrl.TransferFunction(num, den)
    return T

def solucao_43():
    K_values = [10, 200, 500]
    cores = ['b', 'g', 'r']

    results = []
    plt.figure()
    for K, cor in zip(K_values, cores):
        T = sistema_P43(K)
        y, y = resposta_ao_degrau(T, tempo_final=5)
        Mp, ts, info = desempenho_sistema(T)
        ess = erro_estacionario(T, tipo="degrau")
        plt.plot(y, y, cor, label=f'K={K}')
        results.append({
            "K": K,
            "Mp": Mp,
            "ts": ts,
            "erro_estacionario": erro_estacionario(T, tipo="degrau")
        })

    plt.title('Resposta ao Degrau para Diferentes Valores de K')
    plt.xlabel('Tempo (s)')
    plt.ylabel('Saída y(t)')
    plt.grid(True)
    plt.legend()
    plt.show()
    print("Resultados para diferentes valores de K:")
    for res in results:
        print("--------")
        print("K =", res["K"])
        print("  Overshoot (Mp):", res["Mp"])
        print("  Tempo de acomodação (ts):", res["ts"])
        print("  Erro em regime estacionário:", res["erro_estacionario"])
        print("--------")

def erro_estacionario(G, tipo="degrau"):
    s = ctrl.TransferFunction.s
    Kp = ctrl.dcgain(G)
    if tipo == "degrau":
        ess = 1 / (1 + Kp)
    elif tipo == "rampa":
        Kv = ctrl.dcgain(G * s)
        ess = 1 / Kv if Kv != 0 else float('inf')
    else:
        ess = np.nan
    return ess

def simular_resposta(G_num, G_den, titulo= "resposta ao degrau"):
    G = ctrl.TransferFunction(G_num, G_den)
    t, y = ctrl.step_response(G, T=np.linspace(0, 5, 1000))
    ss_value = y[-1] # valor em regime estacionário aproximado
    ess = 1 - ss_value  #erro estacionário para degrau unitário
    plt.figure()
    plt.plot(t, y, label='Resposta ao Degrau')
    plt.title(titulo)
    plt.axhline(1, color='r', linestyle='--', label='Entrada degrau (1)')
    plt.axhline(ss_value, color='g', linestyle='--',
                label=f'Valor em regime estacionário y(∞)={ss_value:.2f}')
    plt.xlabel('Tempo (s)')
    plt.ylabel('Saída y(t)')
    plt.grid(True)
    plt.legend()
    plt.show()
    print(f"Erro em regime estacionário: {ess}")
    print("valor final da saída y(∞): ", ss_value)

if __name__ == "__main__":
    # 4.1
    # G_num = [12]
    # G_den = [1, 2, 10]
    # 4.2
    G_num = [4]
    G_den = [1, 2, 20]
    G, T = criar_sistema(G_num, G_den, feedback=1)
    # 4.3
    print("Solução do exercício 4.3:")
    solucao_43()
    # resposta ao degrau
    t, y = resposta_ao_degrau(T, tempo_final=10)
    Mp, ts, info = desempenho_sistema(T)
    yt_func, yt_latex = resposta_ao_degrau_funcao(T)

    print(f"Overshoot (máxima ultrapassagem): {Mp}%")
    print(f"Tempo de acomodação: {ts}s")

    erro = erro_estacionario(G, tipo="degrau")
    print(f"Erro em regime estacionario para entrada degrau: {erro}")
    print("resposta ao degrau: y(t) =", yt_func)
    print("Expressão LaTeX da resposta ao degrau:", yt_latex)
    # simular_resposta(G_num, G_den, titulo="Resposta ao Degrau do Sistema em Malha Fechada")

