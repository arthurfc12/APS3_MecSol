from math import *
from cmath import *

from numpy import linalg
import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg
from tqdm import tqdm
"""
A funcao 'plota' produz um gráfico da estrutura definida pela matriz de nos N 
e pela incidencia Inc.

Sugestao de uso:

from funcoesTermosol import plota
plota(N,Inc)
-------------------------------------------------------------------------------
A funcao 'importa' retorna o numero de nos [nn], a matriz dos nos [N], o numero
de membros [nm], a matriz de incidencia [Inc], o numero de cargas [nc], o vetor
carregamento [F], o numero de restricoes [nr] e o vetor de restricoes [R] 
contidos no arquivo de entrada.

Sugestao de uso:
    
from funcoesTermosol import importa
[nn,N,nm,Inc,nc,F,nr,R] = importa('entrada.xlsx')
-------------------------------------------------------------------------------
A funcao 'geraSaida' cria um arquivo nome.txt contendo as reacoes de apoio Ft, 
deslocamentos Ut, forcas Fi e tensoes Ti internas. As entradas devem ser 
vetores coluna.

Sugestao de uso:
    
from funcoesTermosol import geraSaida
geraSaida(nome,Ft,Ut,Epsi,Fi,Ti)
-------------------------------------------------------------------------------

"""
def plota(N,Inc):
    # Numero de membros
    nm = len(Inc[:,0])
    
    import matplotlib as mpl
    import matplotlib.pyplot as plt

    #plt.show()
    fig = plt.figure()
    # Passa por todos os membros
    for i in range(nm):
        
        # encontra no inicial [n1] e final [n2] 
        n1 = int(Inc[i,0])
        n2 = int(Inc[i,1])        

        plt.plot([N[0,n1-1],N[0,n2-1]],[N[1,n1-1],N[1,n2-1]],color='r',linewidth=3)


    plt.xlabel('x [m]')
    plt.ylabel('y [m]')
    plt.grid(True)
    plt.axis('equal')
    plt.show()
    
def importa(entradaNome):
    
    import numpy as np
    import xlrd
    
    arquivo = xlrd.open_workbook(entradaNome)
    
    ################################################## Ler os nos
    nos = arquivo.sheet_by_name('Nos')
    
    # Numero de nos
    nn = int(nos.cell(1,3).value)
                 
    # Matriz dos nós
    N = np.zeros((2,nn))
    
    for c in range(nn):
        N[0,c] = nos.cell(c+1,0).value
        N[1,c] = nos.cell(c+1,1).value
    
    ################################################## Ler a incidencia
    incid = arquivo.sheet_by_name('Incidencia')
    
    # Numero de membros
    nm = int(incid.cell(1,5).value)
                 
    # Matriz de incidencia
    Inc = np.zeros((nm,4))
    
    for c in range(nm):
        Inc[c,0] = int(incid.cell(c+1,0).value)
        Inc[c,1] = int(incid.cell(c+1,1).value)
        Inc[c,2] = incid.cell(c+1,2).value
        Inc[c,3] = incid.cell(c+1,3).value
    
    ################################################## Ler as cargas
    carg = arquivo.sheet_by_name('Carregamento')
    
    # Numero de cargas
    nc = int(carg.cell(1,4).value)
                 
    # Vetor carregamento
    F = np.zeros((nn*2,1))
    
    for c in range(nc):
        no = carg.cell(c+1,0).value
        xouy = carg.cell(c+1,1).value
        GDL = int(no*2-(2-xouy)) 
        F[GDL-1,0] = carg.cell(c+1,2).value
         
    ################################################## Ler restricoes
    restr = arquivo.sheet_by_name('Restricao')
    
    # Numero de restricoes
    nr = int(restr.cell(1,3).value)
                 
    # Vetor com os graus de liberdade restritos
    R = np.zeros((nr,1))
    
    for c in range(nr):
        no = restr.cell(c+1,0).value
        xouy = restr.cell(c+1,1).value
        GDL = no*2-(2-xouy) 
        R[c,0] = GDL-1


    return nn,N,nm,Inc,nc,F,nr,R

def geraSaida(nome,Ft,Ut,Epsi,Fi,Ti):
    nome = nome + '.txt'
    f = open(nome,"w+")
    f.write('Reacoes de apoio [N]\n')
    f.write(str(Ft))
    f.write('\n\nDeslocamentos [m]\n')
    f.write(str(Ut))
    f.write('\n\nDeformacoes []\n')
    f.write(str(Epsi))
    f.write('\n\nForcas internas [N]\n')
    f.write(str(Fi))
    f.write('\n\nTensoes internas [Pa]\n')
    f.write(str(Ti))
    f.close()
    
#################################################################################################
def faz_conectividade(m_incidencia, n_membros, n_nos):
  
    C = []
    for i in range(n_membros):
        C.append(conectividade(i+1, m_incidencia, n_nos))
    return np.array(C).T

def conectividade(n_elemento, matriz, n_nos):
    
    conectividade = n_nos*[0]
    
    no_1 = int(matriz[n_elemento-1, 0])
    no_2 = int(matriz[n_elemento-1, 1])
    
    conectividade[no_1-1] = -1
    conectividade[no_2-1] = 1

    return conectividade

def soma_mKs(n_nos, n_membros, m_incidencia, m_nos, m_membros):
   
    get_shape = calculate_K(1, m_incidencia, m_nos, m_membros, n_nos) 
    
    x = get_shape.shape[0] 
    y = get_shape.shape[1] 
    
    kg = np.zeros((x, y)) 
    for i in range(n_membros):
        kg = kg + calculate_K(i, m_incidencia, m_nos, m_membros, n_nos)        
    return kg

def vetor_forcas(v_restricoes, v_carregamento):
    
    matriz_resp = v_carregamento.copy() 
    
    v_restricoes = v_restricoes[:,0].tolist()
    v_restricoes_int = [int(item) for item in v_restricoes]
    
    matriz_resp = np.delete(matriz_resp, v_restricoes_int, axis=0) 
        
    return matriz_resp

def MS(matriz, v_rest):

    v_rest = v_rest[:,0].tolist() 
    v_rest_int = [int(item) for item in v_rest] 
    
    matriz_resp = matriz.copy() 
    matriz_resp = np.delete(matriz_resp, v_rest_int, axis=0) 
    matriz_resp = np.delete(matriz_resp, v_rest_int, axis=1) 
    
    return matriz_resp

def calcula_deslocamentos(matriz_rigidez, matriz_força):
    L,U = scipy.linalg.lu(matriz_rigidez, permute_l=True)
    y = scipy.linalg.solve(L, matriz_força)
    x = scipy.linalg.solve_triangular(U, y)
    return x

def forca(matriz_k, u, linha_number):
    
    linha = matriz_k[linha_number,:]
    return np.dot(linha, u) 

def forca_total(somaKs, u2, vet_restricoes):
    lista_forcas = []
    for i in vet_restricoes[:,0]:
        lista_forcas.append(calculate_force(somaKs, u2, int(i)))
        
    return lista_forcas

def tensao_forca_deformacao(area, n_elemento, n_membros, matriz_U, matriz_incidencia, matriz_nos):
   
    n1 = int(matriz_incidencia[n_elemento-1, 0])
    n2 = int(matriz_incidencia[n_elemento-1, 1])     
    
    matriz_aux = np.array((
            matriz_U[(n1-1)*2], 
            matriz_U[(n1-1)*2 +1], 
            matriz_U[(n2-1)*2], 
            matriz_U[(n2-1)*2 +1]))
    
    E =  matriz_incidencia[n_elemento-1, 2]  
    
    n1 = int(matriz_incidencia[n_elemento-1, 0])
    n2 = int(matriz_incidencia[n_elemento-1, 1])
    
    X_n1 = matriz_nos[0, no_1-1]
    Y_n1 = matriz_nos[1, no_1-1]
    x_n2 = matriz_nos[0, no_2-1]
    y_n2 = matriz_nos[1, no_2-1]
    
    l = sqrt((x_no2-x_no1)**2+(y_no2-y_no1)**2)
    sen = (y_no2-y_no1)/l
    cos = (x_no2-x_no1)/l
    
    c = np.array(([-cos, -sen, cos, sen]))
    
    deformacao = (1/l) * np.dot(c, matriz_aux)
    tensao = (E/l) * np.dot(c, matriz_aux)
    forca = (E/l) * np.dot(c, matriz_aux)*area
    
    return deformacao[0], tensao[0], forca[0]


def jacobi(ite,tol,K,F):
    if type(K) != list:
        U = np.zeros((K.shape[0],1))
        n = K.shape[0]
    else:
        U = np.zeros((len(K),1))
        n = len(K)
    U_ = U.copy()
    count = 0
    while count < ite:
        for i in range(n):
            U_[i][0] = F[i]
            for j in range(n):
                if j != i:
                    U_[i][0] -= K[i][j] * U[j][0]
            U_[i][0] /= K[i][i]

        if 0 not in U:     # Cálculo do erro
            ERROR_list = []
            for i,j in zip(U_,U):
                ERROR_list.append(abs((i-j)/j))
            ERRO = max(ERROR_list)
            if ERRO < tol:
                print("Iterações: {0}".format(count))
                return U
        U = U_.copy() 
        count += 1
    print("Iterações: {0}".format(count))
    return U

def gauss_seidel(ite,tol,K,F):
    if type(K) != list:
        U = np.zeros((K.shape[0],1))
        n = K.shape[0]
    else:
        U = np.zeros((len(K),1))
        n = len(K)
    U_ = U.copy()
    count = 0
    while count < ite:
        for i in range(n):
            U_[i][0] = F[i]
            for j in range(n):
                if j != i:
                    U_[i][0] -= K[i][j] * U_[j][0]
            U_[i][0] /= K[i][i]

        if 0 not in U:     # Cálculo do erro
            ERROR_list = []
            for i,j in zip(U_,U):
                ERROR_list.append(abs((i-j)/j))
            ERRO = max(ERROR_list)
            if ERRO < tol:
                print("Iterações: {0}".format(count))
                return U
        U = U_.copy() 
        count += 1
    print("Iterações: {0}".format(count))
    return U












def dist(x1, y1, x2, y2):
    
    dist = sqrt((x2 - x1)**2 + (y2 - y1)**2)
    return dist

def compara_solucoes(a1, a2):
    return max( abs( (s2 - s1)/s2) for s1,s2 in zip(a1, a2))


def  calculate_K(num_membro, matriz_incidencia, matriz_nos, matriz_membros, num_nos):
    """
    função responsável por calcular a matriz K para um elemento
    recebe: número do membro [inteiro], matriz incidência, matriz de nós e matriz de membros
    retorna: matriz K para o dado elemento
    """
    
    MC = np.array([matriz_conectividade(num_membro, matriz_incidencia, num_nos)]).T #vetor (9,1)
    MC_T = MC.T #vetor (1,9)
    Se = calculate_Se(num_membro, matriz_incidencia, matriz_nos, matriz_membros)#(2,2)
    dot = np.dot(MC, MC_T) 

    return np.kron(dot, Se)

def calculate_Se(num_membro, matriz_incidencia, matriz_nos, matriz_membros):
    """ 
    Função responsável por calcular o valor de Se para um dado elemento
    recebe: número do membro [inteiro], matriz de incidência [matriz], matriz dos nós [matriz], matriz dos membros [matriz]
    retorna: valor de Se para o elemento escolhido [matriz/inteiro]
    """
    
    E = matriz_incidencia[num_membro-1, 2] #elemento linear elástico
    A = matriz_incidencia[num_membro-1, 3] #área da senção transversal
    
    no_1 = int(matriz_incidencia[num_membro-1, 0])
    no_2 = int(matriz_incidencia[num_membro-1, 1])
    
    # Pegar a coordenada do no_1
    x_no1 = matriz_nos[0, no_1-1]
    y_no1 = matriz_nos[1, no_1-1]
    
    # Pegar a coordenada do no_2
    x_no2 = matriz_nos[0, no_2-1]
    y_no2 = matriz_nos[1, no_2-1]
    
    l = distancia_entre_pontos(x_no2, x_no1, y_no2, y_no1)
    
    cordeenadas_membro = matriz_membros[:, num_membro-1] #pega as coordenadas do membro
    
    coornadas_matriz = np.array([cordeenadas_membro]) #transforma em matriz
    coornadas_matriz_T = coornadas_matriz.T #calcula a transposta
    
    me = sum(i**2 for i in cordeenadas_membro)
    
    segunda_parte = (np.dot(coornadas_matriz_T, coornadas_matriz))/me
    
    return ((E*A)/l)*segunda_parte


