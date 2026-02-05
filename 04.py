# Varredura de Biossegurança Profissional (Batch Processing)
import subprocess
import pandas as pd
from Bio import SeqIO
import os
import math

# --- 1. CONFIGURAÇÕES GERAIS ---
# Caminho do arquivo com os genes do fungo
ARQUIVO_QUERY = "c_abscissum_data/ncbi_dataset/data/GCA_023376855.1/cds_from_genomic.fna"

# Nome do arquivo onde o relatório será salvo
ARQUIVO_SAIDA_TXT = "relatorio_biosseguranca.txt"
ARQUIVO_SAIDA_CSV = "dados_biosseguranca.csv"

# Número de processadores a utilizar (ajuste conforme seu PC)
NUM_THREADS = 32 

# Tamanho do lote de genes para processar por vez 
# (Quanto menor, mais atualizações de progresso, mas ligeiramente mais lento)
TAMANHO_LOTE = 500 

# Configuração dos Bancos de Dados (Nome Exibição : Nome Arquivo BLAST DB)
BANCOS_ALVO = {
    "Citrus": "db_citrus",
    "Human":  "db_human",
    "Bee":    "db_bee"
}

# Parâmetros de Risco (RNAi)
MIN_LENGTH = 21   # Tamanho mínimo do match (pb)
MIN_IDENT = 100   # Identidade necessária (%)

# --- 2. FUNÇÕES AUXILIARES ---

def rodar_blast_chunk(query_fasta_str, db_name, threads):
    """
    Roda o BLAST para um pedaço (chunk) de sequências.
    Retorna um DataFrame com os resultados.
    """
    # Cria arquivo temporário para o lote atual
    temp_query = "temp_batch_query.fasta"
    with open(temp_query, "w") as f:
        f.write(query_fasta_str)
        
    cmd = [
        "blastn",
        "-task", "blastn-short",
        "-query", temp_query,
        "-db", db_name,
        "-outfmt", "6 qseqid length pident", # Formato tabular simplificado
        "-word_size", "7",
        "-evalue", "1000",
        "-perc_identity", str(MIN_IDENT),
        "-num_threads", str(threads)
    ]
    
    # Executa o BLAST
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    # Limpa arquivo temporário
    if os.path.exists(temp_query):
        os.remove(temp_query)
        
    if result.returncode != 0:
        return None # Erro na execução
        
    # Se não houve output (nenhum match), retorna vazio
    if not result.stdout:
        return pd.DataFrame(columns=["qseqid", "length", "pident"])
        
    # Processa os dados
    from io import StringIO
    df = pd.read_csv(StringIO(result.stdout), sep="\t", names=["qseqid", "length", "pident"])
    return df

# --- 3. EXECUÇÃO PRINCIPAL ---

print(f"--- INICIANDO PIPELINE DE BIOSSEGURANCA ---")
print(f"Threads: {NUM_THREADS}")
print(f"Lendo arquivo de entrada: {ARQUIVO_QUERY}...")

# 3.1 Carregar todos os genes na memória
todos_registros = list(SeqIO.parse(ARQUIVO_QUERY, "fasta"))
total_genes = len(todos_registros)
print(f"Total de genes carregados: {total_genes}")

# Inicializa dicionário de resultados
# Estrutura: { "ID_GENE": {"Citrus": False, "Human": False, "Bee": False} }
matriz_risco = {rec.id: {sp: False for sp in BANCOS_ALVO} for rec in todos_registros}

# 3.2 Loop por Espécie (Banco de Dados)
for nome_especie, db_arquivo in BANCOS_ALVO.items():
    print(f"\nAnalise iniciada para: {nome_especie.upper()}")
    
    genes_processados = 0
    total_lotes = math.ceil(total_genes / TAMANHO_LOTE)
    
    # Loop por Lotes (Chunks)
    for i in range(0, total_genes, TAMANHO_LOTE):
        # Pega fatia da lista
        lote = todos_registros[i : i + TAMANHO_LOTE]
        
        # Converte para string FASTA
        fasta_str = ""
        for rec in lote:
            fasta_str += f">{rec.id}\n{str(rec.seq)}\n"
            
        # Roda BLAST
        df_resultado = rodar_blast_chunk(fasta_str, db_arquivo, NUM_THREADS)
        
        if df_resultado is not None:
            # Filtra hits perigosos
            df_perigoso = df_resultado[df_resultado['length'] >= MIN_LENGTH]
            genes_com_match = df_perigoso['qseqid'].unique()
            
            # Atualiza a matriz de risco
            for gene_id in genes_com_match:
                if gene_id in matriz_risco:
                    matriz_risco[gene_id][nome_especie] = True
        else:
            # Marca como erro se o BLAST falhou
            for rec in lote:
                matriz_risco[rec.id][nome_especie] = "ERRO"

        # Atualiza progresso
        genes_processados += len(lote)
        progresso = (genes_processados / total_genes) * 100
        print(f"\rProcessando {nome_especie}: {progresso:.1f}% concluido ({genes_processados}/{total_genes})", end="")

    print("") # Nova linha após terminar a espécie

# --- 4. GERAÇÃO DO RELATÓRIO ---
print(f"\n--- GERANDO ARQUIVOS DE SAIDA ---")

aprovados = 0
reprovados = 0

# Prepara listas para o Pandas (CSV)
dados_csv = []

with open(ARQUIVO_SAIDA_TXT, "w") as f:
    # Cabeçalho
    header = f"{'GENE ID':<40} | {'CITRUS':<10} | {'HUMAN':<10} | {'BEE':<10}\n"
    divider = "="*80 + "\n"
    
    f.write(divider)
    f.write(f"RELATORIO DE BIOSSEGURANCA - RNAi\n")
    f.write(divider)
    f.write(header)
    f.write(divider)
    
    for rec in todos_registros:
        gene_id = rec.id
        status_line = f"{gene_id[:37]:<40} | "
        is_safe_global = True
        
        row_csv = {"GeneID": gene_id}
        
        for sp in BANCOS_ALVO:
            risco = matriz_risco[gene_id][sp]
            
            if risco == True:
                tag = "RISK"
                is_safe_global = False
            elif risco == "ERRO":
                tag = "ERROR"
                is_safe_global = False
            else:
                tag = "SAFE"
            
            row_csv[sp] = tag
            status_line += f"{'['+tag+']':<10} | "
        
        f.write(status_line + "\n")
        
        # Define status final para o CSV
        row_csv["Global_Status"] = "APPROVED" if is_safe_global else "REJECTED"
        dados_csv.append(row_csv)
        
        if is_safe_global:
            aprovados += 1
        else:
            reprovados += 1
            
    # Resumo Final no arquivo
    f.write(divider)
    f.write(f"TOTAL ANALISADO: {total_genes}\n")
    f.write(f"APROVADOS (100% SEGURO): {aprovados}\n")
    f.write(f"REPROVADOS (COM RISCO):  {reprovados}\n")

# Salva CSV
df_final = pd.DataFrame(dados_csv)
df_final.to_csv(ARQUIVO_SAIDA_CSV, index=False, sep=";")

print(f"Sucesso!")
print(f"Relatorio de texto salvo em: {ARQUIVO_SAIDA_TXT}")
print(f"Tabela de dados salva em: {ARQUIVO_SAIDA_CSV}")
print(f"Total de Genes Aprovados: {aprovados}")