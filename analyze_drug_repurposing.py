import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from sklearn.cluster import KMeans
from sklearn.manifold import TSNE
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.decomposition import PCA
from wordcloud import WordCloud
import networkx as nx
from community import community_louvain
import spacy
from transformers import pipeline
import dash
from dash import dcc, html
from dash.dependencies import Input, Output
import dash_bootstrap_components as dbc

# Load spaCy model for NLP
nlp = spacy.load('en_core_web_sm')
sentiment_analyzer = pipeline('sentiment-analysis')

# Read and preprocess data
df = pd.read_csv('drug_repurposing_data.csv')

# Advanced Statistics
print("\nAdvanced Statistics:")
print(f"Total number of drug-disease pairs: {len(df)}")
print(f"Number of unique diseases: {df['disease_name'].nunique()}")
print(f"Number of unique drugs: {df['drug_name'].nunique()}")

# Calculate disease-drug network metrics
G = nx.Graph()
for _, row in df.iterrows():
    G.add_edge(row['disease_name'], row['drug_name'], weight=row['confidence'])

print("\nNetwork Analysis:")
print(f"Number of nodes: {G.number_of_nodes()}")
print(f"Number of edges: {G.number_of_edges()}")
print(f"Average clustering coefficient: {nx.average_clustering(G):.3f}")
print(f"Network density: {nx.density(G):.3f}")

# Community detection
partition = community_louvain.best_partition(G)
df['community'] = [partition[node] for node in df['disease_name']]

# Sentiment Analysis of Reasons
df['sentiment'] = [sentiment_analyzer(reason)[0]['label'] for reason in df['reason']]
df['sentiment_score'] = [sentiment_analyzer(reason)[0]['score'] for reason in df['reason']]

# Create visualizations
plt.style.use('seaborn')

# 1. Advanced Distribution Analysis
fig = make_subplots(rows=2, cols=2, 
                   subplot_titles=('Confidence Distribution', 'Sentiment Distribution',
                                 'Disease-Drug Network', 'Community Structure'))

# Confidence distribution with KDE
fig.add_trace(go.Histogram(x=df['confidence'], nbinsx=20, name='Confidence'), row=1, col=1)
fig.add_trace(go.Scatter(x=np.linspace(0, 1, 100), 
                        y=np.exp(-0.5*((np.linspace(0, 1, 100)-df['confidence'].mean())/df['confidence'].std())**2),
                        name='Normal Distribution'), row=1, col=1)

# Sentiment distribution
fig.add_trace(go.Pie(labels=df['sentiment'].value_counts().index,
                    values=df['sentiment'].value_counts().values), row=1, col=2)

# Network visualization
pos = nx.spring_layout(G)
edge_trace = go.Scatter(x=[], y=[], line=dict(width=0.5, color='#888'), hoverinfo='none', mode='lines')
node_trace = go.Scatter(x=[], y=[], text=[], mode='markers', hoverinfo='text',
                       marker=dict(showscale=True, colorscale='YlGnBu', size=10))

for edge in G.edges():
    x0, y0 = pos[edge[0]]
    x1, y1 = pos[edge[1]]
    edge_trace['x'] += (x0, x1, None)
    edge_trace['y'] += (y0, y1, None)

for node in G.nodes():
    x, y = pos[node]
    node_trace['x'] += (x,)
    node_trace['y'] += (y,)
    node_trace['text'] += (node,)

fig.add_trace(edge_trace, row=2, col=1)
fig.add_trace(node_trace, row=2, col=1)

# Community structure
community_sizes = pd.Series(partition).value_counts()
fig.add_trace(go.Bar(x=community_sizes.index, y=community_sizes.values), row=2, col=2)

fig.update_layout(height=800, showlegend=False)
fig.write_html('advanced_analysis.html')

# 2. Interactive Dashboard
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

app.layout = dbc.Container([
    dbc.Row([
        dbc.Col(html.H1("Drug Repurposing Analysis Dashboard"), width=12)
    ]),
    dbc.Row([
        dbc.Col([
            html.H3("Disease Selection"),
            dcc.Dropdown(
                id='disease-dropdown',
                options=[{'label': disease, 'value': disease} for disease in df['disease_name'].unique()],
                value=df['disease_name'].unique()[0]
            )
        ], width=4),
        dbc.Col([
            html.H3("Confidence Threshold"),
            dcc.Slider(
                id='confidence-slider',
                min=0,
                max=1,
                step=0.1,
                value=0.5,
                marks={i/10: str(i/10) for i in range(11)}
            )
        ], width=4)
    ]),
    dbc.Row([
        dbc.Col(dcc.Graph(id='disease-drug-network'), width=6),
        dbc.Col(dcc.Graph(id='confidence-distribution'), width=6)
    ]),
    dbc.Row([
        dbc.Col(dcc.Graph(id='sentiment-analysis'), width=6),
        dbc.Col(dcc.Graph(id='community-structure'), width=6)
    ])
])

@app.callback(
    [Output('disease-drug-network', 'figure'),
     Output('confidence-distribution', 'figure'),
     Output('sentiment-analysis', 'figure'),
     Output('community-structure', 'figure')],
    [Input('disease-dropdown', 'value'),
     Input('confidence-slider', 'value')]
)
def update_graphs(selected_disease, confidence_threshold):
    filtered_df = df[df['confidence'] >= confidence_threshold]
    
    # Disease-drug network
    network_fig = px.scatter(filtered_df, x='disease_name', y='drug_name',
                           color='confidence', size='confidence',
                           hover_data=['reason', 'pmid'])
    
    # Confidence distribution
    dist_fig = px.histogram(filtered_df, x='confidence', nbins=20,
                          title='Confidence Score Distribution')
    
    # Sentiment analysis
    sentiment_fig = px.pie(filtered_df, names='sentiment',
                         title='Sentiment Analysis of Drug Reasons')
    
    # Community structure
    community_fig = px.bar(pd.Series(partition).value_counts(),
                          title='Community Structure')
    
    return network_fig, dist_fig, sentiment_fig, community_fig

if __name__ == '__main__':
    app.run_server(debug=True)

# 3. Machine Learning Analysis
# TF-IDF Vectorization of reasons
vectorizer = TfidfVectorizer(max_features=100)
X = vectorizer.fit_transform(df['reason'])

# Dimensionality reduction
pca = PCA(n_components=2)
X_pca = pca.fit_transform(X.toarray())

# Clustering
kmeans = KMeans(n_clusters=5)
clusters = kmeans.fit_predict(X_pca)

# Create cluster visualization
cluster_fig = px.scatter(x=X_pca[:, 0], y=X_pca[:, 1], color=clusters,
                        hover_name=df['disease_name'],
                        title='Drug-Disease Clusters')
cluster_fig.write_html('drug_clusters.html')

# 4. Word Cloud of Drug Reasons
wordcloud = WordCloud(width=800, height=400, background_color='white').generate(' '.join(df['reason']))
plt.figure(figsize=(10, 5))
plt.imshow(wordcloud, interpolation='bilinear')
plt.axis('off')
plt.savefig('reason_wordcloud.png')
plt.close()

print("\nAdvanced visualizations have been created:")
print("- advanced_analysis.html")
print("- drug_clusters.html")
print("- reason_wordcloud.png")
print("\nInteractive dashboard is running at http://localhost:8050") 