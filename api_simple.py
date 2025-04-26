from flask import Flask, jsonify, request
import pandas as pd
import json
import os
from flask_cors import CORS

app = Flask(__name__)
CORS(app)  # Enable CORS for all routes

# Load the data
try:
    print("Attempting to load data...")
    if not os.path.exists('drug_repurposing_data.csv'):
        print("Error: drug_repurposing_data.csv not found in current directory")
        print("Current directory:", os.getcwd())
        raise FileNotFoundError("drug_repurposing_data.csv not found")
    
    df = pd.read_csv('drug_repurposing_data.csv')
    print("Data loaded successfully")
    print("Available diseases:", df['disease_name'].unique().tolist())
except Exception as e:
    print(f"Error loading data: {str(e)}")
    raise

@app.route('/search_diseases', methods=['POST'])
def search_diseases():
    try:
        print("Received search request")
        data = request.get_json()
        if not data or 'query' not in data:
            return jsonify({"error": "Missing query parameter"}), 400
        
        query = data['query'].lower()
        print(f"Searching for: {query}")
        
        # Search in disease names
        disease_matches = df[df['disease_name'].str.lower().str.contains(query)]
        print(f"Found {len(disease_matches)} matches")
        
        results = []
        for _, row in disease_matches.iterrows():
            results.append({
                'name': row['disease_name'],
                'source': 'Custom',
                'drugs': [{
                    'drug_name': row['drug_name'],
                    'confidence': float(row['confidence']),
                    'reason': row['reason'],
                    'pmid': row['pmid']
                }]
            })
        
        return jsonify(results)
    except Exception as e:
        print(f"Error in search_diseases: {str(e)}")
        return jsonify({"error": str(e)}), 500

@app.route('/disease/<disease_name>', methods=['GET'])
def get_disease_details(disease_name):
    try:
        print(f"Getting details for: {disease_name}")
        disease_data = df[df['disease_name'] == disease_name]
        if disease_data.empty:
            return jsonify({'error': 'Disease not found'}), 404
        
        drugs = []
        for _, row in disease_data.iterrows():
            drugs.append({
                'drug_name': row['drug_name'],
                'confidence': float(row['confidence']),
                'reason': row['reason'],
                'pmid': row['pmid']
            })
        
        return jsonify({
            'name': disease_name,
            'source': 'Custom',
            'drugs': drugs
        })
    except Exception as e:
        print(f"Error in get_disease_details: {str(e)}")
        return jsonify({"error": str(e)}), 500

@app.route('/all_diseases', methods=['GET'])
def get_all_diseases():
    try:
        print("Getting all diseases")
        return jsonify(df['disease_name'].unique().tolist())
    except Exception as e:
        print(f"Error in get_all_diseases: {str(e)}")
        return jsonify({"error": str(e)}), 500

if __name__ == '__main__':
    print("Starting API server...")
    app.run(host='0.0.0.0', port=8000, debug=True) 