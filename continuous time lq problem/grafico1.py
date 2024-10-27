import numpy as np
import pandas as pd
from scipy.linalg import expm
import seaborn as sns
import matplotlib.pyplot as plt

sns.set_theme(style="darkgrid", font="Calibri")
plt.figure(figsize=(8.49, 5.18))

dados = pd.read_excel(r"C:\Users\DESKTOP\OneDrive\Documentos\Python\lq continuo.xlsx") # local onde a base foi salva
dados = dados.set_index(dados.index)

def create_and_save_plot(data, column_name):
    sns.set_theme(style="darkgrid", font="Calibri")
    plt.figure(figsize=(8.49, 5.18))
    sns.lineplot(x=data.index/100, y=column_name, data=data, color="#002060", linewidth=1.7)
    plt.ylabel(f"Valor do componente")
    plt.xlabel('Tempo (t)')
    plt.savefig(f"C:\\Users\\DESKTOP\\OneDrive\\Documentos\\Python\\{column_name}cont.png")

# Exemplo de uso:
create_and_save_plot(dados, "serie_1")
create_and_save_plot(dados, "serie_2")
create_and_save_plot(dados, "serie_3")

# sns.scatterplot(x=dados.index, y="serie_3", data=dados, color="#002060", linewidth=1.7)
# plt.ylabel("Valor do Componente 3")
# plt.savefig(r"C:\Users\DESKTOP\OneDrive\Documentos\Python\componente3.png") # local onde vai ser salva a imagem do gráfico.