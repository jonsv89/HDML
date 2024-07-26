import distanceclosure as dc
import pandas as pd
import networkx as nx
import numpy as np

import os
import multiprocessing as mp
import time

def ensurePathExists(filepath):
    """Given a file with path, ensures a path exists by creating a directory if needed. """
    if not os.path.exists(os.path.dirname(filepath)):
        os.makedirs(os.path.dirname(filepath))


def _data_list():
    
    # Layers that have FS/IFS and weighted by cosine and jaccards
    layer_type = ['drug', 'gene', 'metabolome', 'methylome1e3', 'methylome1e7', 'microbiome', 'mirna', 'symptom']
    weight_type = ['fs', 'ifs']
    normalization_type = ['cosinesimilarities', 'jaccardindex', 'weightedjaccardindex']
    
    triples = []
    for layer in layer_type:
        for weight in weight_type:
            for normalization in normalization_type:
                triples.append((layer, weight, normalization))
    
    weight = 'tfidf'
    for layer in ['mirna', 'symptom']:
        for normalization in normalization_type:
            triples.append((layer, weight, normalization))
    
    normalization = 'cosinesimilarities'
    weight = 'tvalue'
    for layer in ['methylome1e3', 'methylome1e7']:
        triples.append((layer, weight, normalization))
    
    layer = 'microbiome'
    normalization = 'cosinesimilarities'
    for weight in ['fsd','ifsd']:
        triples.append((layer, weight, normalization))
    
    layer = 'transcriptome'
    weight = 'fc'
    normalization = 'cosinesimilarities'
    triples.append((layer, weight, normalization))
    
    return triples


def closure_evolution(layer, weight, normalization):
    
    rGfile = f'Networks/Intralayer/ICD10_{layer}_{weight}_{normalization}.txt'
    if layer == 'symptom':
        rGfile += '.zip'
            
    df_data = pd.read_csv(rGfile, sep='\t')
    
    df_graph = df_data[df_data[normalization]>0]
    df_graph.loc[:, 'distance'] = 1./df_graph[normalization] - 1.
    
    G = nx.from_pandas_edgelist(df_graph, source='disease1', target='disease2', edge_attr='distance')
    
    T = 5
    
    D = np.empty(T, dtype=object)
    D[0] = G.copy()
    
    for t in range(1, T):
        #print(t)
        
        D[t] = D[t-1].copy()
                
        for u in G.nodes():
            d = dict()
            
            for k in D[t-1].neighbors(u):
                if k == u:
                    d[k] = 0
                else:
                    d[k] = D[t-1][u][k]['distance']
                
                for v in G.neighbors(k):
                    if v in d.keys():
                        d[v] = min(d[v], d[k] + G[k][v]['distance'])
                    else:
                        d[v] = d[k] + G[k][v]['distance']
            
            for v, dval in d.items():
                if (v == u) and D[t-1].has_edge(u, v):
                    dval = min(dval, D[t-1][u][v]['distance'])
                
                if dval > 0.0:
                    if D[t].has_edge(u, v):
                        D[t][u][v]['distance'] = dval
                    else:
                        D[t].add_edge(u, v, distance=dval)
                
        dc_df = nx.to_pandas_edgelist(D[t])
        
        wGfile = f'Closure/Step{t+1}/' + rGfile
        ensurePathExists(wGfile)
        dc_df.to_csv(wGfile)


def closure_evolution_epidemiology():
    
    rGfile = 'Networks/Multilayer/epidemiology.txt'
    df_data = pd.read_csv(rGfile, sep='\t')
    
    maxRR = 1.01+df_data['RR'].max()
    df_data['weight'] = df_data['RR']/maxRR 
    df_data.loc[:, 'distance'] = 1./df_data['weight'] - 1.
    
    G = nx.from_pandas_edgelist(df_data, source='disease1', target='disease2', edge_attr='distance')
    
    T = 5
    
    D = np.empty(T, dtype=object)
    D[0] = G.copy()
    
    for t in range(1, T):
        #print(t)
        
        D[t] = D[t-1].copy()
                
        for u in G.nodes():
            d = dict()
            
            for k in D[t-1].neighbors(u):
                if k == u:
                    d[k] = 0
                else:
                    d[k] = D[t-1][u][k]['distance']
                
                for v in G.neighbors(k):
                    if v in d.keys():
                        d[v] = min(d[v], d[k] + G[k][v]['distance'])
                    else:
                        d[v] = d[k] + G[k][v]['distance']
            
            for v, dval in d.items():
                if (v == u) and D[t-1].has_edge(u, v):
                    dval = min(dval, D[t-1][u][v]['distance'])
                
                if dval > 0.0:
                    if D[t].has_edge(u, v):
                        D[t][u][v]['distance'] = dval
                    else:
                        D[t].add_edge(u, v, distance=dval)
                
        dc_df = nx.to_pandas_edgelist(D[t])
        
        wGfile = f'Closure/Step{t+1}/' + rGfile
        ensurePathExists(wGfile)
        dc_df.to_csv(wGfile)



if __name__ == "__main__":
    
    triples = _data_list()
    
    pool = mp.Pool(processes=10)
    t0 = time.time()
    mp_out = pool.starmap(closure_evolution, triples)
    t1 = time.time()
    pool.close()
    print(f'Execution time {(t1 - t0)//60}m{round((t1 - t0)%60)}s')
    print(mp_out) 
    
    closure_evolution_epidemiology()