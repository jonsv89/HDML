import distanceclosure as dc
import pandas as pd
import networkx as nx
import numpy as np

import os


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


def unidentifiable_edge_graph(G: nx.Graph, H: nx.Graph, weight:str = 'weight'):
    '''
    Compute missing edges in the overlap of both graphs G and H which are closer via connections outside of the overlap.
    
    weight is the edge distance attribute
    
    '''
    
    #graph = sorted([G, H], key=len) # Sort graphs by number of nodes
    graph = [G, H]
    overlap_nodes = [u for u in graph[0].nodes if graph[1].has_node(u)]
    
    new_edges = [[], []]
    for idx in range(2):
        
        subgraph = graph[idx].subgraph(overlap_nodes).copy()
        
        graph[idx] = dc.metric_backbone(graph[idx], weight=weight)
        subgraph = dc.metric_backbone(subgraph, weight=weight)
        
        for i in range(len(overlap_nodes)):
            for j in range(i, len(overlap_nodes)):
                
                u = overlap_nodes[i]
                v = overlap_nodes[j]
                
                if not graph[idx].has_edge(u, v):
                    if subgraph.has_edge(u, v):
                        new_edges[idx].append((u, v))
                    else:
                        try:
                            duv_orig = nx.shortest_path_length(graph[idx], source=u, target=v, weight=weight)
                        except:
                            duv_orig = np.inf
                        
                        try:
                            duv_over = nx.shortest_path_length(subgraph, source=u, target=v, weight=weight)
                        except:
                            duv_over = np.inf
                        
                        if duv_orig < duv_over:
                            new_edges[idx].append((u, v))
             
                    
    #print(new_edges)    
    return new_edges
    


def overlap_with_epidemiology(G, layer, weight, normalization):
    
    rGfile = f'Networks/Intralayer/ICD10_{layer}_{weight}_{normalization}.txt'
    if layer == 'symptom':
        rGfile += '.zip'
    
    df_data = pd.read_csv(rGfile, sep='\t')
        
    df_data = pd.read_csv(rGfile, sep='\t')
    df_data = df_data[df_data[normalization]>0]
    df_data.loc[:, 'distance'] = 1./df_data[normalization] - 1.
    
    H = nx.from_pandas_edgelist(df_data, source='disease1', target='disease2', edge_attr='distance')
    
    new_edges = unidentifiable_edge_graph(G, H, weight='distance')
    
    print(new_edges)

if __name__ == "__main__":
    
    '''
    # Test
    G = nx.Graph()
    G.add_weighted_edges_from([(1, 2, 2), (1, 3, 2), (2, 4, 1), (3, 4, 1)])
    
    H = nx.Graph()
    H.add_weighted_edges_from([(0, 1, 1), (1, 2, 1), (1, 3, 1), (2, 3, 1)])
    
    new_edges = unidentifiable_edge_graph(G, H) 
    '''
    

    rGfile = f'Networks/Multilayer/epidemiology.txt'
    df_data = pd.read_csv(rGfile, sep='\t')
    maxRR = 1.01+df_data['RR'].max()
    df_data['weight'] = df_data['RR']/maxRR 
    df_data.loc[:, 'distance'] = 1./df_data['weight'] - 1.
    
    G = nx.from_pandas_edgelist(df_data, source='disease1', target='disease2', edge_attr='distance')
    
    triples = _data_list()
    for layer, weight, normalization in triples:
        overlap_with_epidemiology(G, layer, weight, normalization)
        
        break