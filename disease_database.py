import pandas as pd
from typing import Dict, List
import requests
from mygene import MyGeneInfo
import json

class DiseaseDatabase:
    def __init__(self):
        self.mg = MyGeneInfo()
        self.disease_data = {}
        
    def load_orphanet_data(self):
        """Load disease data from Orphanet"""
        try:
            # Orphanet API endpoint for rare diseases
            response = requests.get('https://www.orphadata.com/api/data/diseases')
            orphanet_data = response.json()
            
            for disease in orphanet_data:
                self.disease_data[disease['name']] = {
                    'source': 'Orphanet',
                    'orpha_code': disease['orpha_code'],
                    'genes': self._get_disease_genes(disease['name']),
                    'prevalence': disease.get('prevalence', 'Unknown'),
                    'inheritance': disease.get('inheritance', 'Unknown')
                }
        except Exception as e:
            print(f"Error loading Orphanet data: {e}")
            
    def load_nord_data(self):
        """Load disease data from NORD (National Organization for Rare Disorders)"""
        try:
            # NORD API endpoint
            response = requests.get('https://rarediseases.org/api/diseases')
            nord_data = response.json()
            
            for disease in nord_data:
                if disease['name'] not in self.disease_data:
                    self.disease_data[disease['name']] = {
                        'source': 'NORD',
                        'nord_id': disease['id'],
                        'genes': self._get_disease_genes(disease['name']),
                        'description': disease.get('description', ''),
                        'symptoms': disease.get('symptoms', [])
                    }
        except Exception as e:
            print(f"Error loading NORD data: {e}")
            
    def load_custom_data(self, csv_path: str):
        """Load custom disease data from CSV"""
        try:
            df = pd.read_csv(csv_path)
            for _, row in df.iterrows():
                if row['disease_name'] not in self.disease_data:
                    self.disease_data[row['disease_name']] = {
                        'source': 'Custom',
                        'genes': self._get_disease_genes(row['disease_name']),
                        'drugs': self._get_disease_drugs(row['disease_name'], df),
                        'confidence': row.get('confidence', 0.0),
                        'pmid': row.get('pmid', '')
                    }
        except Exception as e:
            print(f"Error loading custom data: {e}")
            
    def _get_disease_genes(self, disease_name: str) -> List[str]:
        """Get genes associated with a disease using MyGene"""
        try:
            results = self.mg.query(disease_name, 
                                  fields='entrezgene,symbol',
                                  species='human')
            return [gene['symbol'] for gene in results if 'symbol' in gene]
        except Exception as e:
            print(f"Error getting genes for {disease_name}: {e}")
            return []
            
    def _get_disease_drugs(self, disease_name: str, df: pd.DataFrame) -> List[Dict]:
        """Get drugs associated with a disease from the custom dataset"""
        drugs = df[df['disease_name'] == disease_name]
        return [{
            'drug_name': row['drug_name'],
            'confidence': row['confidence'],
            'reason': row['reason'],
            'pmid': row['pmid']
        } for _, row in drugs.iterrows()]
        
    def search_diseases(self, query: str) -> List[Dict]:
        """Search diseases by name or gene"""
        results = []
        query = query.lower()
        
        for name, data in self.disease_data.items():
            if query in name.lower():
                results.append({
                    'name': name,
                    'source': data['source'],
                    'genes': data['genes'],
                    'drugs': data.get('drugs', []),
                    'confidence': data.get('confidence', 0.0)
                })
            elif 'genes' in data and any(query in gene.lower() for gene in data['genes']):
                results.append({
                    'name': name,
                    'source': data['source'],
                    'genes': data['genes'],
                    'drugs': data.get('drugs', []),
                    'confidence': data.get('confidence', 0.0)
                })
                
        return sorted(results, key=lambda x: x['confidence'], reverse=True)
        
    def get_disease_details(self, disease_name: str) -> Dict:
        """Get detailed information about a specific disease"""
        return self.disease_data.get(disease_name, {})
        
    def get_all_diseases(self) -> List[str]:
        """Get list of all diseases in the database"""
        return list(self.disease_data.keys()) 