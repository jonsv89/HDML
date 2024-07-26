import distanceclosure as dc
import pandas as pd
import networkx as nx

import os


def ensurePathExists(filepath):
    """Given a file with path, ensures a path exists by creating a directory if needed. """
    if not os.path.exists(os.path.dirname(filepath)):
        os.makedirs(os.path.dirname(filepath))


def _order_dataframe(df):
    
    for idx, row in df.iterrows():
        d1, d2, _ = row
        if d2 < d1:   
            df.loc[idx, 'disease1'] = d2
            df.loc[idx, 'disease2'] = d1
    
    return df


def compute_backbone(layer, weight, normalization):
    
    rGfile = f'Networks/Intralayer/ICD10_{layer}_{weight}_{normalization}.txt'
    if layer == 'symptom':
        rGfile += '.zip'
    
    df_data = pd.read_csv(rGfile, sep='\t')
    
    df_graph = df_data[df_data[normalization] > 0]    
    df_graph.loc[:, 'distance'] = 1./df_graph[normalization] - 1.
    
    # No Filter
    G = nx.from_pandas_edgelist(df_graph, source='disease1', target='disease2', edge_attr='distance')
    B = dc.metric_backbone(G, weight='distance')                

    df_backbone = nx.to_pandas_edgelist(B, source='disease1', target='disease2')
    df_backbone = _order_dataframe(df_backbone)
    df_backbone['all_backbone'] = 1
    df_data = pd.merge(df_data, df_backbone, how='outer', on=['disease1', 'disease2'])
    
    # FDR Filter
    df_graph = df_graph[df_graph['fdr_up'] <= 0.05]
    
    G = nx.from_pandas_edgelist(df_graph, source='disease1', target='disease2', edge_attr='distance')
    B = dc.metric_backbone(G, weight='distance')                

    df_backbone = nx.to_pandas_edgelist(B, source='disease1', target='disease2')
    df_backbone = _order_dataframe(df_backbone)
    df_backbone['fdr_backbone'] = 1
    df_data = pd.merge(df_data, df_backbone, how='outer', on=['disease1', 'disease2'])
    
    df_data.fillna(0, inplace=True)
    df_data.drop(columns=['distance_x', 'distance_y'], inplace=True)                
    
    #print(df_data.sort_values('all_backbone'))
    
    wGfile = f'Backbone/Networks/Intralayer/ICD10_{layer}_{weight}_{normalization}.txt'
    ensurePathExists(wGfile)
    df_data.to_csv(wGfile, sep='\t')


def compute_backbone_transcriptome():
    
    layer = 'transcriptome'
    weight = 'fc'
    normalization = 'cosinesimilarities'
    rGfile = f'Networks/Intralayer/ICD10_{layer}_{weight}_{normalization}.txt'
    if layer == 'symptom':
        rGfile += '.zip'
    
    df_data = pd.read_csv(rGfile, sep='\t')
    
    df_graph = df_data[df_data[normalization] > 0]    
    df_graph.loc[:, 'distance'] = 1./df_graph[normalization] - 1.
    
    # No Filter
    G = nx.from_pandas_edgelist(df_graph, source='disease1', target='disease2', edge_attr='distance')
    B = dc.metric_backbone(G, weight='distance')                

    df_backbone = nx.to_pandas_edgelist(B, source='disease1', target='disease2')
    df_backbone = _order_dataframe(df_backbone)
    df_backbone['all_backbone'] = 1
    df_data = pd.merge(df_data, df_backbone, how='outer', on=['disease1', 'disease2'])
    
    # FDR Filter
    df_graph = df_graph[df_graph['fdr'] <= 0.05]
    
    G = nx.from_pandas_edgelist(df_graph, source='disease1', target='disease2', edge_attr='distance')
    B = dc.metric_backbone(G, weight='distance')                

    df_backbone = nx.to_pandas_edgelist(B, source='disease1', target='disease2')
    df_backbone = _order_dataframe(df_backbone)
    df_backbone['fdr_backbone'] = 1
    df_data = pd.merge(df_data, df_backbone, how='outer', on=['disease1', 'disease2'])
    
    df_data.fillna(0, inplace=True)
    df_data.drop(columns=['distance_x', 'distance_y'], inplace=True)                
    
    #print(df_data.sort_values('all_backbone'))
    
    wGfile = f'Backbone/Networks/Intralayer/ICD10_{layer}_{weight}_{normalization}.txt'
    ensurePathExists(wGfile)
    df_data.to_csv(wGfile, sep='\t')    
 
 
def compute_backbone_epidemiology():
    
    layer = 'epidemiology'
    normalization = 'relativerisk'
    rGfile = f'Networks/Multilayer/{layer}.txt'
    
    df_data = pd.read_csv(rGfile, sep='\t')
    
    maxRR = 1.01+df_data['RR'].max()
    df_data['weight'] = df_data['RR']/maxRR 
    df_data.loc[:, 'distance'] = 1./df_data['weight'] - 1.
    
    # No Filter
    G = nx.from_pandas_edgelist(df_data, source='disease1', target='disease2', edge_attr='distance')
    B = dc.metric_backbone(G, weight='distance')                

    df_backbone = nx.to_pandas_edgelist(B, source='disease1', target='disease2')
    df_backbone = _order_dataframe(df_backbone)
    df_backbone['all_backbone'] = 1
    df_data = pd.merge(df_data, df_backbone, how='outer', on=['disease1', 'disease2'])
    
    df_data.fillna(0, inplace=True)    
    df_data.drop(columns=['distance_x', 'distance_y'], inplace=True)                
    
    #print(df_data.sort_values('all_backbone'))
    
    wGfile = f'Backbone/Networks/Multilayer/{layer}.txt'
    ensurePathExists(wGfile)
    df_data.to_csv(wGfile, sep='\t')    


if __name__ == "__main__":
    
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
    
    #compute_backbone_transcriptome()
    #compute_backbone_epidemiology()
     
    data = {'layer': [], 'weight': [], 'normalization': [], 'all_edges': [], 'all_backbone': [], 'fdr_edges': [], 'fdr_backbone': []}
    for layer, weight, normalization in triples:
        rGfile = f'Backbone/Networks/Intralayer/ICD10_{layer}_{weight}_{normalization}.txt'
        df_data = pd.read_csv(rGfile, sep='\t')
        
        #edges = len(df_data[df_data[normalization]>0])
        #backbone_edges = len(df_data[df_data['all_backbone']>0])
        
        data['layer'].append(layer)
        data['weight'].append(weight)
        data['normalization'].append(normalization)
        
        data['all_edges'].append(len(df_data[df_data[normalization]>0]))
        data['all_backbone'].append(len(df_data[df_data['all_backbone']>0]))
        
        if layer == 'transcriptome':
            data['fdr_edges'].append(len(df_data[df_data['fdr']<=0.05]))
        else:
            data['fdr_edges'].append(len(df_data[df_data['fdr_up']<=0.05]))
        data['fdr_backbone'].append(len(df_data[df_data['fdr_backbone']>0]))
    
    df = pd.DataFrame(data)
    df['all_backbone_size'] = df['all_backbone']/df['all_edges']
    df['fdr_backbone_size'] = df['fdr_backbone']/df['fdr_edges']
    
    df.to_csv('Backbone_Summary.csv')