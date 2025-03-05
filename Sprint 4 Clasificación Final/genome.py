import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import MinMaxScaler
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.ndimage import gaussian_filter1d
import csv
import argparse


def process_genomic_data(csv_file, visualizar):
    """
    Carga un archivo CSV, filtra los datos, elimina valores atípicos usando MAD en 'Depth' y genera visualizaciones.

    Parámetros:
        csv_file (str): Ruta del archivo CSV.

    Retorna:
        pd.DataFrame: DataFrame procesado.
    """
    # Configurar estilo de gráficos
    plt.style.use('ggplot')
    sns.set_style("whitegrid")

    # Cargar datos
    with open(csv_file, "r", newline="") as f:
      sample = f.read(1024)  # Read a small portion of the file
      dialect = csv.Sniffer().sniff(sample)

    df = pd.read_csv(csv_file, delimiter=dialect.delimiter)
    #print(df)
    antes = len(df)
    # Filtrar por longitud mínima
    df = df[df['Length'] >= 1000]

    # Eliminar valores atípicos en 'Depth' usando el método MAD
    median_depth = df['Depth'].median()
    MAD = np.median(np.abs(df['Depth'] - median_depth))
    threshold_mad = 4 * MAD
    lower_bound_mad = median_depth - threshold_mad
    upper_bound_mad = median_depth + threshold_mad
    df = df[(df['Depth'] >= lower_bound_mad) & (df['Depth'] <= upper_bound_mad)]
    despues = len(df)
    # Graficar distribuciones
    if visualizar:
      fig, axes = plt.subplots(2, 2, figsize=(12, 10))

      sns.histplot(df["Length"], bins=20, kde=True, color="seagreen", ax=axes[0, 0])
      axes[0, 0].set_title("Distribución de la Longitud de los Scaffolds")

      sns.histplot(df["GC"], bins=20, kde=True, ax=axes[0, 1])
      axes[0, 1].set_title("Distribución del contenido GC")

      sns.histplot(df["Depth"], bins=20, kde=True, color="orange", ax=axes[1, 0])
      axes[1, 0].set_title("Distribución de la Profundidad de los Scaffolds")

      sns.histplot(df["AltAllels"], bins=10, kde=True, color="cyan", ax=axes[1, 1])
      axes[1, 1].set_title("Distribución de los Alelos Alternativos")

      plt.tight_layout()
      plt.savefig('distribuciones.png', dpi=300, bbox_inches='tight')
      plt.close()  

      plt.figure(figsize=(8, 6))  
      #plt.show()

      # Ajustar el número de grupos para evitar desajustes con la paleta
      num_groups = min(len(df), 5)  # Máximo 6 grupos, pero limitado al tamaño de df
      df = df.sort_values(by='Length', ascending=False)
      df['Scaffold group'] = np.repeat(range(1, num_groups + 1), np.ceil(len(df) / num_groups))[:len(df)]
      custom_palette = sns.color_palette("viridis", num_groups)
      '''
      # Graficar pairplot
      g = sns.pairplot(
          df,
          diag_kind="kde",
          vars=['Length', 'GC', 'Depth', 'AltAllels'],
          hue='Scaffold group',
          palette=custom_palette
      )

      g._legend.remove()
      '''
      #plt.show()

      # Graficar Length vs Depth
      sns.scatterplot(x='Depth', y='Length', data=df, hue='Scaffold group', palette=custom_palette)
      plt.xlabel('Profundidad')
      plt.ylabel('Longitud')
      plt.title('Profundidad vs Longitud')
      plt.savefig('Longitud_Profundidad.png', dpi=300, bbox_inches='tight')

      # Calcular el total de AltAllels
      #print(f"Total de Alelos Alternativos: {total_altallels}")
    total_altallels = df["AltAllels"].sum()
    return df, antes-despues, total_altallels


def calcular_n50(df):
    df_sorted = df.sort_values(by='Length', ascending=False)
    total_length = df_sorted['Length'].sum()
    half_length = total_length / 2

    cumulative_length = 0
    for length in df_sorted['Length']:
        cumulative_length += length
        if cumulative_length >= half_length:
            return length


def mezcla_de_gaussianas(df, n, visualizar):
    """
    Genera visualizaciones de distribuciones gaussianas basadas solo en Depth.
    Normaliza los datos antes del ajuste del modelo y los devuelve en su escala original al final.

    Parámetros:
        df (pandas.DataFrame): DataFrame del experimento
        n (int): Número de gaussianas a ajustar en el modelo.

    Retorna:
        pandas.DataFrame: DataFrame con la clasificación de gaussianas y datos en su escala original.
    """

    if not isinstance(n, int) or n < 1:
        raise ValueError("El número de gaussianas 'n' debe ser un entero positivo.")

    df_clean = df.dropna(subset=['Depth']).copy()

    scaler = StandardScaler()
    df_clean[['Depth']] = scaler.fit_transform(df_clean[['Depth']])

    X = df_clean[['Depth']].values

    gmm = GaussianMixture(n_components=n, covariance_type='tied', random_state=13)
    gmm.fit(X)

    df_clean['gaussian'] = gmm.predict(X)
    if visualizar:
      plt.figure(figsize=(8, 6))
      sns.scatterplot(x=df_clean["Depth"], y=df_clean["Length"], hue=df_clean["gaussian"], palette="viridis", alpha=0.9)
      plt.xlabel("Depth")
      plt.ylabel("Length")
      plt.title("Fitting Gaussian Mixture Model to Depth vs Length")

      x_min, x_max = df_clean["Depth"].min(), df_clean["Depth"].max()
      x = np.linspace(x_min, x_max, 100).reshape(-1, 1)

      z = -gmm.score_samples(x).reshape(-1)

      plt.show()

    df_clean[['Depth']] = scaler.inverse_transform(df_clean[['Depth']])

    return df_clean

def max_Depth(df):
  return df["Length"].max()

def segunda_maxima(df):
  return df.nlargest(2, 'Length')['Length'].iloc[1]

def tercera_maxima(df):
  return df.nlargest(3, 'Length')['Length'].iloc[2]

def depth_difference(df, length1, length2):
    depth = 'Depth'
    if len(df.columns) < 4:
        depth = 'Depth_bin'

    depth1 = df.loc[df['Length'] == length1, depth].values[0]
    depth2 = df.loc[df['Length'] == length2, depth].values[0]

    return abs(depth1 - depth2)

def peaks(df, n):
    return df.nlargest(n, 'Length')['Length'].iloc[n - 1]

def two_bacterias(df):
    P = 0
    L_1 = peaks(df, 1)
    for i in range(2, 8):
      L_j = peaks(df, i)
      P += np.sqrt(L_1 * L_j) * abs(depth_difference(df, L_1, L_j))

    return P*100

def acumular_longitudes(df, bins, visualizar):

  num_bins = bins
  depth_bins = np.linspace(df["Depth"].min(), df["Depth"].max(), num_bins)
  bin_centers = (depth_bins[:-1] + depth_bins[1:]) / 2
  df['Depth_bin'] = pd.cut(df['Depth'], bins=depth_bins, labels=bin_centers)

  # Calcular el total de Length en cada bin de Depth
  binned_data = df.groupby(['Depth_bin', 'gaussian'], observed=False)['Length'].sum().reset_index()
  binned_data['Depth_bin'] = binned_data['Depth_bin'].astype(float)
  if visualizar:
    # Graficar los datos binned como un bar plot con colores según la clase
    plt.figure(figsize=(8,6))
    classes = binned_data['gaussian'].unique()
    colors = sns.color_palette("viridis", len(classes))
    for g, color in zip(classes, colors):
      subset = binned_data[binned_data['gaussian'] == g]
      plt.bar(subset['Depth_bin'], subset['Length'], width=np.diff(depth_bins).mean(), color=color, alpha=0.6, label=f"Class {g}")

    # Etiquetas y título
    plt.xlabel("Profundidad")
    plt.ylabel("Suma de Longitudes")
    plt.title("Longitudes acumuladas")
    plt.legend()
    plt.grid(True)
    plt.savefig('longitud_acumulada.png', dpi=300, bbox_inches='tight')
  binned_data = binned_data.groupby("Depth_bin", as_index=False)["Length"].sum()
  scaler = MinMaxScaler()
  df_scaled = pd.DataFrame(scaler.fit_transform(binned_data), columns=binned_data.columns)
  return df_scaled


def GC_peaks(df, visualizar):

  # Extract GC content column (Replace 'gc_content_column' with the actual column name)
  gc_content = df["GC"].dropna()  # Remove NaN values if any
  min = gc_content.min()
  max = gc_content.max()

  bins = np.linspace(min, max, 40)  # Adjust bin range
  hist, bin_edges = np.histogram(gc_content, bins=bins)
  smoothed_hist = gaussian_filter1d(hist, sigma=2)
  # Find local maxima
  peaks, properties = find_peaks(smoothed_hist)

  # Find global maximum
  global_max_index = np.argmax(smoothed_hist)

  bins = []
  bins.append(global_max_index)
  for p in bin_edges[peaks]:
    if smoothed_hist[np.where(bin_edges == p)[0][0]] != smoothed_hist[global_max_index]:
      bins.append(np.where(bin_edges == p)[0][0])

  bins = np.array(bins)
  if visualizar:
    # Plot histogram
    plt.figure(figsize=(8, 5))
    plt.plot(bin_edges[:-1], smoothed_hist, label="Histograma", color="blue")
    plt.plot(bin_edges[:-1][bins], smoothed_hist[bins], "ro", label="Maximo Local")
    plt.plot(bin_edges[:-1][global_max_index], smoothed_hist[global_max_index], "go", markersize=10, label="Global Maximum")
    plt.xlabel("Contenido GC")
    plt.ylabel("Frecuencia")
    plt.legend()
    plt.savefig('picos_de_GC.png', dpi=300, bbox_inches='tight')
  return len(bins)

def clasificacion(final):
  if final > 1000:
    return 'La muestra está contaminada'
  elif final > 323 and final <= 1000:
    return '2 bacterias detectadas'
  else:
    return '1 bacteria detectada'

def s1(x):
    a = np.log(999) / (1000 - 323)
    return 1 / (1 + np.exp(a * (x - 323)))

def s2(x):
    a = np.log(999) / (1000 - 323)
    return 1 / (1 + np.exp(-a * (x - 323)))

def s3(x):
    a = np.log(999) / 1000
    return 1 / (1 + np.exp(-a * (x - 1000)))

def main():
    '''
    Usamos un parser de argumentos para recibir los argumentos de entrada.
    '''
    parser = argparse.ArgumentParser(description="genome.")
    parser.add_argument("-input", "--archivo_entrada", type=str, required=True, help="Archivo de entrada csv con los resultados del experimento")
    parser.add_argument("-output", "--gráficos", type=bool, required=False, help="Generar gráficos de los datos", default=False)
    args  = parser.parse_args() # initialize argument parser and store in args
    
    # Abrir el archivo de entrada y leer las líneas
    try:
        # Procesar los datos del experimento
        if args.gráficos:
            visualizar = True
        else:
           visualizar = False
        processing = process_genomic_data(args.archivo_entrada,visualizar)
    except FileNotFoundError:
        print(f"Error: Archivo no encontrado - {args.archivo_entrada}")
        return
    except Exception as e:
        print(f"Error inesperado al abrir el archivo de entrada: {e}")
        return
    
    # Crear un diccionario con los datos de ejemplo
    data = {
        "número de datos": 0,
        "número de datos eliminados": 0,
        "alelos alternativos": 0,
        "N50": 0,
        "longitud máxima": 0,
        "segunda longitud máxima": 0,
        "tercera longitud máxima": 0,
        "diferencia de profundidades longitudes 1-2": 0,
        "diferencia de profundidades 1-3": 0,
        "suma de longitudes": 0,
        "contaminación": 0,
        "2 bacterias" : 0,
        "GC peaks": 0,
    }

    # dataframe preprocesado
    df = processing[0]

    # Análisis básico
    avg_depth = df["Depth"].mean()
    avg_gc = df["GC"].mean()
    num_contigs = df.shape[0]

    # Aplicar la mezcla de gaussianas y acumular longitudes
    dfg = mezcla_de_gaussianas(df,2, False)
    dfa = acumular_longitudes(dfg,50, False)

    # Calcular medidas para las métricas
    data["número de datos"] = len(df)
    data["número de datos eliminados"] = processing[1]
    data["alelos alternativos"] = processing[2]
    data["N50"] = int(calcular_n50(dfg))
    data["longitud máxima"] = int(max_Depth(dfg))
    data["segunda longitud máxima"] = int(segunda_maxima(dfg))
    data["tercera longitud máxima"] = int(tercera_maxima(dfg))
    data["diferencia de profundidades longitudes 1-2"] = depth_difference(dfg, max_Depth(dfg), segunda_maxima(dfg))
    data["diferencia de profundidades 1-3"] = depth_difference(dfg, max_Depth(dfg), tercera_maxima(dfg))
    data["suma de longitudes"] = int(df['Length'].sum())
    data["contaminación"] = (int(df['Length'].sum())*((processing[2]/len(df))+0.1)/calcular_n50(df))
    data["2 bacterias"] = two_bacterias(dfa)
    data["GC peaks"] = GC_peaks(df, visualizar)

    # Métrica de dos bacterias enhanced
    metrica1 = data["2 bacterias"] * data["GC peaks"]

    # Métrica de contaminación  
    metrica2 = data["contaminación"] 

    # Métrica final
    final = metrica1 * 3 + metrica2 * 2

    if clasificacion(final) == '2 bacterias detectadas':
        dfg = mezcla_de_gaussianas(df,2, False)
        dfa = acumular_longitudes(dfg,50, visualizar)
        certeza = s2(final)
    elif clasificacion(final) == '1 bacteria detectada':
        dfg = mezcla_de_gaussianas(df,1, False)
        dfa = acumular_longitudes(dfg,50, visualizar)
        certeza = s1(final)
    else:
        dfg = mezcla_de_gaussianas(df,3, False)
        dfa = acumular_longitudes(dfg,50, visualizar)
        certeza = s3(final)

    print(f"📂 Procesando archivo: {args.archivo_entrada}")
    print(f"- Contigs analizados: {num_contigs}")
    print(f"- Profundidad media de cobertura: {avg_depth:.2f}")
    print(f"- Porcentaje medio de GC: {avg_gc:.4f}")
    print(f"🔬 Resultado: {clasificacion(final)}")
    print(f"🦠 Métricas:")
    print(f"   - Alelos alternativos: {data['alelos alternativos']}")
    print(f"   - N50: {data['N50']}")
    print(f"   - Longitud máxima: {data['longitud máxima']}")
    print(f"   - Longitud total: {data['suma de longitudes']}")
    print(f"   - Contaminación: {metrica2}")
    print(f"   - Picos de GC: {data['GC peaks']}")
    print(f"   - Certeza: {certeza:.2f}")


if __name__ == "__main__":
    main()
