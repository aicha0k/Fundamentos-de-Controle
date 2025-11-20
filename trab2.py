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
        t, y = resposta_ao_degrau(T, tempo_final=5)
        Mp, ts, info = desempenho_sistema(T)
        ess = erro_estacionario(T, tipo="degrau")
        plt.plot(t, y, cor, label=f'K={K}')
        results.append({
            "K": K,
            "Mp": Mp,
            "ts": ts,
            "erro_estacionario": ess
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
        # Para sistema em malha fechada T(s): ess = 1 - T(0) = 1 - dcgain(T)
        # ess = 1 - Kp
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

def PC45(k):
    C = ctrl.TransferFunction([10], [1, 0])
    P = ctrl.TransferFunction([1], [1, k])
    G = P * C # malha aberta
    T = ctrl.feedback(G, 1) # malha fechada
    return T, G

#     ks = np.linspace(k_min, k_max, passos)
#     candidatos = []
#     for k in ks:
#         T, G = PC45(k)
#         Mp, ts, info = desempenho_sistema(T)
#         # Verificar se Mp está no intervalo desejado (1% < Mp < 10%)
#         # Mp pode ser None se não houver overshoot
#         if Mp is not None and not np.isnan(Mp) and 1 < Mp < 10:
#             candidatos.append({
#                 "k": k,
#                 "Mp": Mp,
#                 "ts": ts,
#             })
    
#     if not candidatos:
#         print("Nenhum valor de k encontrado com 1% < Mp < 10%")
#         return
    
#     # Selecionar o primeiro candidato (ou pode escolher outro critério, como menor ts)
#     k_selecionado = candidatos[0]["k"]
#     print(f"Valor de k selecionado: {k_selecionado:.4f}")
#     print(f"Mp = {candidatos[0]['Mp']:.2f}%")
#     print(f"ts = {candidatos[0]['ts']:.4f}s")
    
#     # Calcular T(s) para o k selecionado
#     T, G = PC45(k_selecionado)
    
#     # Mostrar a função de transferência T(s) = Y(s)/R(s)
#     print("\nFunção de transferência em malha fechada T(s) = Y(s)/R(s):")
#     print(f"Numerador: {T.num[0][0]}")
#     print(f"Denominador: {T.den[0][0]}")
#     print(f"T(s) = {T}")
    
#     # Gerar resposta ao degrau
#     t, y = resposta_ao_degrau(T, tempo_final=10)
    
#     # Calcular valor em regime estacionário
#     ss_value = y[-1]  # valor final
#     ess = 1 - ss_value  # erro estacionário (deve ser próximo de zero)
    
#     # Plotar resposta ao degrau
#     plt.figure(figsize=(10, 6))
#     plt.plot(t, y, 'b-', linewidth=2, label=f'Resposta ao degrau (k={k_selecionado:.4f})')
#     plt.axhline(1, color='r', linestyle='--', linewidth=2, label='Entrada degrau unitário (1)')
#     plt.axhline(ss_value, color='g', linestyle='--', linewidth=2, 
#                 label=f'Valor em regime estacionário y(∞)={ss_value:.6f}')
#     plt.title(f'Resposta ao Degrau Unitário - T(s) com k={k_selecionado:.4f}')
#     plt.xlabel('Tempo (s)')
#     plt.ylabel('Saída y(t)')
#     plt.grid(True, alpha=0.3)
#     plt.legend()
#     plt.xlim([0, max(t)])
#     plt.ylim([0, max(1.2, max(y) * 1.1)])
#     plt.show()
    
#     # Verificar erro estacionário
#     print(f"\nVerificação do erro em regime estacionário:")
#     print(f"Valor final da saída y(∞) = {ss_value:.10f}")
#     print(f"Erro estacionário ess = 1 - y(∞) = {ess:.10f}")
#     if abs(ess) < 1e-6:
#         print("✓ Erro estacionário é praticamente nulo (ess ≈ 0)")
#     else:
#         print(f"⚠ Erro estacionário não é exatamente nulo, mas é muito pequeno: {ess:.10f}")
    
#     return T, k_selecionado, candidatos
def encontrar_k_P45(k_min=0.1, k_max=100, passos=500):
    ks = np.linspace(k_min, k_max, passos)
    candidatos = []
    for k in ks:
        T, G = PC45(k)
        Mp, ts, info = desempenho_sistema(T)
        if Mp is not None and not np.isnan(Mp) and 1 < Mp < 10:
            candidatos.append({
                "k": k,
                "Mp": Mp,
                "ts": ts,
            })
    return candidatos

def plotar_P45(k):
    T, G = PC45(k)
    t, y = resposta_ao_degrau(T, tempo_final=10)
    plt.figure()
    plt.plot(t, y, label=f'k={k}')
    plt.axhline(1, color='r', linestyle='--', label='Entrada degrau (1)')
    plt.title(f'Resposta ao Degrau para k={k}')
    plt.xlabel('Tempo (s)')
    plt.ylabel('Saída y(t)')
    plt.grid(True)
    plt.legend()
    plt.show()
    
    Mp, ts, info = desempenho_sistema(T)
    print("Informações do sistema:")
    print("k = ", k)
    print(f"Overshoot (máxima ultrapassagem): {Mp}%")
    print(f"Tempo de acomodação: {ts}s")
    ess = erro_estacionario(T, tipo="degrau")
    print(f"Erro em regime estacionario para entrada degrau: {ess}")

def PC48(K):
    # retorna o sistema em malha fechada com controlador proporcional
    # P (s) = 10 / (s + 10)
    Gc = ctrl.TransferFunction([K], [1]) 
    P = ctrl.TransferFunction([10], [1, 10])
    G = Gc * P
    T = ctrl.feedback(G, 1)
    return T, G

def sistema_malha_fechada_PC48(K0=2, K1=20):
    Gc = ctrl.TransferFunction([K0,K1], [1, 0]) # (K0s + K1)/s
    P = ctrl.TransferFunction([10], [1, 10])
    G = Gc * P
    T = ctrl.feedback(G, 1)
    return T, G


def executar_PC48():
    print("\nA - Controlador proporcional k = 2")
    Tprop, Gprop = PC48(2)
    t1, y1 = resposta_ao_degrau(Tprop, tempo_final=10)
    Mp1, ts1, info1 = desempenho_sistema(Tprop)
    ess1 = erro_estacionario(Gprop, tipo="degrau")
    
    print(f"Mp = {Mp1:.3f}%")
    print(f"ts = {ts1:.3f}s")
    print(f"ess = {ess1:.3f}")
    plt.figure()
    plt.plot(t1, y1, label='Controlador Proporcional k = 2')
    plt.title('Resposta ao Degrau para k = 2')
    plt.xlabel('Tempo (s)')
    plt.ylabel('Saída y(t)')
    plt.grid(True)
    plt.legend()
    plt.show()
    # B 
    print("\nB - Controlador proporcional PI K0 = 2, K1 = 20")
    Tpi, Gpi = sistema_malha_fechada_PC48(2, 20)
    t2, y2 = resposta_ao_degrau(Tpi, tempo_final=10)
    Mp2, ts2, info2 = desempenho_sistema(Tpi)
    ess2 = erro_estacionario(Gpi, tipo="degrau")
    
    print(f"Mp = {Mp2:.3f}%")
    print(f"ts = {ts2:.3f}s")
    print(f"ess = {ess2:.3f}")
    plt.figure()
    plt.plot(t2, y2, label='Controlador Proporcional PI K0 = 2, K1 = 20')
    plt.title('Resposta ao Degrau para K0 = 2, K1 = 20')
    plt.xlabel('Tempo (s)')
    plt.ylabel('Saída y(t)')
    plt.grid(True)
    
    # C - comparacao entre A e B
    plt.figure()
    plt.plot(t1, y1, label='Proporcional k = 2')
    plt.plot(t2, y2, label='Proporcional PI K0 = 2, K1 = 20')
    plt.axhline(1, color='r', linestyle='--')
    plt.title("Comparação entre Proporcional e PI")
    plt.grid(True)
    plt.legend()
    plt.show()
    print("resumo:")
    print(f"ess do proporcional: {ess1:.2f} -> erro nao nulo (sistema tipo 0)")
    print(f"ess do PI: {ess2:.2f} -> erro nulo (sistema tipo 1)")
    print("conclusao:\n O controlador proporcional PI reduz, mas nao elimina o erro em regime estacionario.")
    print("O controlador PI adiciona um polo na origem, o que aumenta o tipo do sistema de 0 para 1 e consequente o erro em regime estacionario é nulo.")
    
def sistema_PC411(K):
    C = ctrl.TransferFunction([K], [1]) # controlador
    P = ctrl.TransferFunction([20], [1, 4.5, 64]) # processo
    H = ctrl.TransferFunction([1], [1, 1]) # sensor
    G = C * P # malha aberta
    T = ctrl.feedback(G, H) # malha fechada com sensor
    return T, G, H
    
def executar_PC411(K_values=[10,12,15], tempo_final=10):
    plt.figure()
    results = []
    for K in K_values:
        T, G, H = sistema_PC411(K)
        t, y = resposta_ao_degrau(T, tempo_final=tempo_final)
        
        plt.plot(t, y, label=f'K={K}')
        Mp, ts, info = desempenho_sistema(T)
        ess = erro_estacionario(G, tipo="degrau")
        results.append({
            "K": K,
            "Mp": Mp,
            "ts": ts,
            "ess": ess
        })
    plt.title('Resposta ao Degrau para Diferentes Valores de K')
    plt.xlabel('tempo (s)')
    plt.ylabel('saída y(t)')
    plt.grid(True)
    plt.legend()
    plt.show()
    print("resumo de resultados PC411:")
    for res in results:
        print(f"K={res['K']}")
        print(f"Mp={res['Mp']:.3f}%")
        print(f"Tempo de acomodação: {res['ts']:.3f}s")
        print(f"Erro em regime estacionario: {res['ess']:.3f}")
        print("-----------------")
        
    return results
   
def criar_L():
    s = ctrl.TransferFunction.s
    L = (s + 10) / (s**2 * (s + 15))
    return L

def resposta_rampa(T, tempo_final=50, n_pontos=2000):
    t = np.linspace(0, tempo_final, n_pontos)
    r = t # rampa unitaria r(t) = t
    t_out, y = ctrl.forced_response(T, T=t, U=r)
    e = r - y
    ess_apox = e[-1]
    return t_out, r, y, e, ess_apox

def executar_rampa_L():
    L = criar_L()
    T = ctrl.feedback(L, 1) # malha fechada, realimentacao unitaria
    t,r, y, e, ess_apox = resposta_rampa(T)
    
    plt.figure()
    plt.plot(t, y, label='Resposta à Rampa')
    plt.axhline(1, color='r', linestyle='--', label='Entrada rampa (t)')
    plt.axhline(ess_apox, color='g', linestyle='--', label=f'Erro em regime estacionario: {ess_apox:.3f}')
    plt.title('Resposta à Rampa do Sistema em Malha Fechada')
    plt.xlabel('Tempo (s)')
    plt.ylabel('Saída y(t)')
    plt.grid(True)
    plt.legend()
    
    # plot erro
    plt.figure()
    plt.plot(t, e, label='Erro')
    # plt.axhline(0, color='b', linestyle='--', label='Erro nulo')
    plt.title('Erro em Regime Estacionário PC5.2')
    plt.xlabel('Tempo (s)')
    plt.ylabel('Erro e(t)')
    plt.grid(True)
    plt.legend()
    plt.show()
    print(f"Erro em regime estacionario: {ess_apox:.6f}")
    
def PC54_sistema():
    s = ctrl.TransferFunction.s
    Gc = 21/s
    Gp = 1/(s+2)
    L = Gc * Gp
    T = ctrl.feedback(L, 1)
    return T
def calcular_MUP(T):
    info = ctrl.step_info(T)
    Mp = info['Overshoot']
    return Mp

def plotar_resposta_PC54():
    T = PC54_sistema()
    t, y = ctrl.step_response(T, np.linspace(0, 10, 2000))

    plt.figure()
    plt.plot(t, y, label="Resposta ao degrau")
    plt.axhline(1, color='r', linestyle='--', label="Valor final = 1")
    plt.title("Resposta ao Degrau do PC5.4")
    plt.xlabel("Tempo (s)")
    plt.ylabel("Saída y(t)")
    plt.grid(True)
    plt.legend()
    plt.show()

    # cálculo aproximado visual
    y_max = max(y)
    MUP_aprox = (y_max - 1) * 100
    print(f"y_max = {y_max:.3f}")
    print(f"MUP (estimada pelo gráfico) = {MUP_aprox:.2f}%")

def PC58_sistema():
    s = ctrl.TransferFunction.s
    Gc = (0.1*s + 5) / s
    Gp = 100*(s+1)/ (s**2 + 2*s + 100)
    
    L = Gc * Gp
    T = ctrl.feedback(L, 1)
    return T, L

def parametros_PC58(T):
    den = T.den[0][0]
    polos = np.roots(den)
    polos_ordenados = np.sort(polos)
    p1, p2 = polos_ordenados[:2] # polos dominantes
    print(f"Polos dominantes: {p1:.3f}, {p2:.3f}")
    # parâmetros equivalentes
    wn = np.sqrt(p1.real**2 + p1.imag**2)
    zeta = -p1.real / wn
    print(f"Frequência natural: {wn:.3f}")
    print(f"Coeficiente de amortecimento: {zeta:.3f}")
    return wn, zeta

def calcular_PC58(zeta, wn):
    Mp=np.exp(-np.pi*zeta/np.sqrt(1-zeta**2))*100
    Tp = np.pi / (wn * np.sqrt(1-zeta**2))
    Ts = 4 / (zeta * wn)
    return Mp, Tp, Ts

def simular_PC58():
    T,L = PC58_sistema()
    t,y = ctrl.step_response(T, np.linspace(0, 5, 2000))
    info = ctrl.step_info(T)
    
    Mp = info['Overshoot']
    Tp = info['PeakTime']
    Ts = info['SettlingTime']
    
    wn, zeta = parametros_PC58(T)
    Mp_analitico, Tp_analitico, Ts_analitico = calcular_PC58(zeta, wn)
    
    
    plt.figure()
    plt.plot(t, y, label="Resposta ao degrau")
    plt.axhline(1, linestyle='--', color='r', label="Valor final = 1")
    plt.title("Resposta ao Degrau do PC5.8")
    plt.xlabel("Tempo (s)")
    plt.ylabel("Saída y(t)")
    plt.grid(True)
    plt.legend()
    plt.show()
    print("Resultados do PC5.8:")
    print(f"Mp = {Mp:.3f}%")
    print(f"Tp = {Tp:.3f}s")
    print(f"Ts = {Ts:.3f}s")
    print("--------------------------------")
    print(f"Mp analitico = {Mp_analitico:.3f}%")
    print(f"Tp analitico = {Tp_analitico:.3f}s")
    print(f"Ts analitico = {Ts_analitico:.3f}s")
    
def sistema_PC59():
    G = ctrl.TransferFunction([10], [1, 10])
    H = ctrl.TransferFunction([0.5], [10,0.5])
    T = ctrl.feedback(G, H)
    return G, H, T

def analisar_PC59():
    G, H, T = sistema_PC59()
    t, y = ctrl.step_response(T, np.linspace(0, 30, 2000))
    info = ctrl.step_info(T)
    Mp = info['Overshoot']
    Tp = info['PeakTime']
    Ts = info['SettlingTime']
    print("Resultados do PC5.9:")
    print(f"Máxima ultrapassagem: Mp = {Mp:.3f}%")
    print(f"Tp = {Tp:.3f}s")
    print(f"Tempo de acomodação: Ts = {Ts:.3f}s")
    plt.figure()
    plt.plot(t, y, label="Resposta ao degrau")
    plt.axhline(1, linestyle='--', color='r', label="Valor final = 1")
    plt.title("Resposta ao Degrau do PC5.9")
    plt.xlabel("Tempo (s)")
    plt.ylabel("Saída y(t)")
    plt.grid(True)
    plt.legend()
    plt.show()
    
def sistema_PC510():
    s = ctrl.TransferFunction.s
    G = 10 / (s * (s + 15)* (s + 5))
    T = ctrl.feedback(G, 1)
    return G, T

def analisar_PC510():
    G, T = sistema_PC510()
    t = np.linspace(0, 20, 2000)
    r = t
    t_out, y= ctrl.forced_response(T, T=t, U=r)
    e = r - y
    ess_apox = e[-1]
    
    t_out, r, y, e, ess_apox = resposta_rampa(T, tempo_final=40)
    # ganho de velocidade Kv, pra calcular o erro em regime estacionario
    s = ctrl.TransferFunction.s
    # Kv = ctrl.dcgain(G*s)
    # Kv = ctrl.dcgain((G.minreal()) * s)
    Kv = 10 / (15* 5)
    ess_teorico = 1/Kv
    print("Resultados do PC5.10:")
    print(f"Erro em regime estacionario teorico: {ess_teorico:.3f}")
    print(f"Erro em regime estacionario aproximado: {ess_apox:.3f}")
    print("--------------------------------")
    plt.figure()
    plt.plot(t_out, y, label="Resposta à rampa")
    plt.axhline(1, linestyle='--', color='r', label="Entrada rampa (t)")
    # plt.axhline(ess_apox, linestyle='--', color='g', label=f"Erro em regime estacionario: {ess_apox:.3f}")
    plt.title("Resposta à Rampa do PC5.10")
    plt.xlabel("Tempo (s)")
    plt.ylabel("Saída y(t)")
    plt.grid(True)
    plt.legend()
    plt.show()
    

    

if __name__ == "__main__":
    # 4.1
    # G_num = [12]
    # G_den = [1, 2, 10]
    # 4.2
    G_num = [4]
    G_den = [1, 2, 20]
    G, T = criar_sistema(G_num, G_den, feedback=1)
    # 4.3
    # print("Solução do exercício 4.3:")
    # solucao_43()
    # resposta ao degrau
    t, y = resposta_ao_degrau(T, tempo_final=10)
    Mp, ts, info = desempenho_sistema(T)
    yt_func, yt_latex = resposta_ao_degrau_funcao(T)

    print(f"Overshoot (máxima ultrapassagem): {Mp}%")
    print(f"Tempo de acomodação: {ts}s")

    erro = erro_estacionario(G, tipo="degrau")
    print(f"Erro em regime estacionario para entrada degrau: {erro}")
    #print("resposta ao degrau: y(t) =", yt_func)
    #print("Expressão LaTeX da resposta ao degrau:", yt_latex)
    # simular_resposta(G_num, G_den, titulo="Resposta ao Degrau do Sistema em Malha Fechada")

    # 4.5
    #print("Solução do exercício 4.5:")
    # print("Encontrando o valor de k para 1% < Mp < 10%:")
    #candidatos = encontrar_k_P45()
    # if len(candidatos) == 0:
    #     print("Nenhum valor de k encontrado com 1% < Mp < 10%")
    # else:
    #     print("Candidatos encontrados:")
    #     for candidato in candidatos:
    #         print(f"k = {candidato['k']}, Mp = {candidato['Mp']}%")
    #     print("Plotando a resposta ao degrau para o primeiro candidato:")
    #     plotar_P45(candidatos[0]['k'])
    #print("Solução do exercício 4.8:")
    #executar_PC48()
    # print("Solução do exercício 4.11:")
    #executar_PC411()
    
    # Resposta à rampa
    # print("\n" + "="*50)
    # executar_resposta_rampa()
    
    #executar_rampa_L()
    # print("Solução do exercício 5.4:")
    # T = PC54_sistema()
    # Mp = calcular_MUP(T)
    # print(f"Mp calculado = {Mp:.2f}%")
    # plotar_resposta_PC54()
    # print("Solução do exercício 5.8:")
    # simular_PC58()
    # print("Solução do exercício 5.9:")
    # analisar_PC59()
    print("Solução do exercício 5.10:")
    analisar_PC510()