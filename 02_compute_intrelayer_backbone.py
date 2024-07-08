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

if __name__ == "__main__":
    
    # Layers that have FS/IFS and weighted by cosine and jaccards
    layer_type = ['drug', 'gene', 'metabolome', 'methylome1e3', 'methylome1e7', 'microbiome', 'mirna', 'symptom']
    weight_type = ['fs', 'ifs']
    normalization_type = ['cosinesimilarities', 'jaccardindex', 'weightedjaccardindex']
    
    for layer in layer_type:
        for weight in weight_type:
            for normalization in normalization_type:
                
                print(layer, weight, normalization)
                
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
                
                #break
            #break
        #break