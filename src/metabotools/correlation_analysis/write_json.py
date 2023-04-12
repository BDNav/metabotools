from IPython import embed
import networkx as nx



def convert_graph_for_cytoscapejs(G):
    from networkx.readwrite import json_graph
    data = json_graph.node_link_data(G)
    data['edges'] = data['links']
    new_nodes = []
    for n in data['nodes']:
        n['label'] = n['id']
        # embed()
        new_nodes.append({"data":n})
    new_edges = []
    for e in data['edges']:
        source = e['source']
        target = e['target']
        source = G.nodes()[source]
        target = G.nodes()[target]
        e['source'] = source
        e['target'] = target
        new_edges.append({"data":e})
    data['nodes'] = new_nodes
    data['edges'] = new_edges
    return data
    
def convert_graph_for_d3(G):
    from networkx.readwrite import json_graph
    data = json_graph.node_link_data(G)
    return data
    

    

G=nx.Graph()

G.add_node('1')
G.add_node('b')
G.add_node('c')
G.add_node('d')
G.add_node('e')
G.add_node('f')
G.add_node('n')
G.add_edge('1','b')
G.add_edge('1','1')
G.add_edge('e','f')
G.add_edge('b','c')
G.add_edge('d','e')

    
    
data = convert_graph_for_cytoscapejs(G)

import json

s = json.dumps(data, sort_keys=True, indent=4, separators=(',', ': '))
with open('networkx_test.json', 'w') as outfile:
    outfile.write(s)