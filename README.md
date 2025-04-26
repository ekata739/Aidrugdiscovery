# Repurposing Hope: AI for Rare Disease Treatment Discovery

An advanced AI-powered platform for analyzing and visualizing drug repurposing data for rare diseases, featuring machine learning, natural language processing, and interactive dashboards.

## Project Structure

- `drug_repurposing_data.csv`: Contains the dataset with disease-drug pairs, confidence scores, and supporting information
- `analyze_drug_repurposing.py`: Advanced Python script for data analysis, machine learning, and visualization
- `requirements.txt`: List of required Python packages

## Advanced Features

### 1. Interactive Dashboard
- Real-time filtering of diseases and confidence scores
- Dynamic network visualization of disease-drug relationships
- Sentiment analysis of drug repurposing reasons
- Community detection in the disease-drug network

### 2. Machine Learning Analysis
- TF-IDF vectorization of drug repurposing reasons
- Dimensionality reduction using PCA
- K-means clustering of drug-disease pairs
- Network analysis and community detection

### 3. Natural Language Processing
- Sentiment analysis of drug repurposing reasons
- Word cloud visualization of common terms
- Text mining for pattern discovery

### 4. Advanced Visualizations
- Interactive network graphs
- Community structure analysis
- Confidence score distributions
- Drug-disease clusters
- Sentiment analysis pie charts

## Data Format

The dataset contains the following columns:
- `disease_name`: Name of the rare disease
- `drug_name`: Name of the repurposed drug
- `reason`: Explanation of the drug's potential use
- `confidence`: Confidence score (0-1) for the drug-disease pair
- `pmid`: PubMed ID of the supporting research

## Setup

1. Install the required packages:
```bash
pip install -r requirements.txt
python -m spacy download en_core_web_sm
```

2. Run the analysis script:
```bash
python analyze_drug_repurposing.py
```

3. Access the interactive dashboard at http://localhost:8050

## Output

The script generates several advanced visualizations:
- `advanced_analysis.html`: Comprehensive analysis dashboard
- `drug_clusters.html`: Interactive cluster visualization
- `reason_wordcloud.png`: Word cloud of drug repurposing reasons
- Interactive dashboard with real-time filtering

## Technical Features

- Network analysis using NetworkX
- Community detection using Louvain algorithm
- Sentiment analysis using Transformers
- Interactive visualizations using Plotly and Dash
- Machine learning clustering using scikit-learn
- Natural language processing using spaCy

## Analysis Features

The script provides:
- Basic statistics about the dataset
- Distribution analysis of confidence scores
- Disease-specific drug repurposing patterns
- Interactive visualizations for exploring relationships
- Network analysis of disease-drug connections

## Contributing

Feel free to contribute to this project by:
1. Adding new drug-disease pairs
2. Improving the machine learning models
3. Enhancing the visualizations
4. Adding new features
5. Optimizing the performance

## License

This project is open source and available under the MIT License. 

## Running the API

To run the API, use the following command:
```bash
python api.py 
```

## Testing the API

You can test the API using the following command:
```bash
curl -X POST "http://localhost:8000/search_diseases" \
     -H "Content-Type: application/json" \
     -d '{"query": "Duchenne"}'
```

To get all diseases:
```bash
curl "http://localhost:8000/all_diseases"