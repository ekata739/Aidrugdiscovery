import pandas as pd
import torch
from knowledge_graph import KnowledgeGraphBuilder, HeteroGNN
from model_trainer import DrugRepurposingTrainer, optimize_hyperparameters
from api import app
import uvicorn
from typing import Dict, List
import json

def load_data() -> Dict[str, pd.DataFrame]:
    """Load and preprocess data from various sources"""
    # Load drug repurposing data
    drug_data = pd.read_csv('drug_repurposing_data.csv')
    
    # Add placeholder data for demonstration
    drug_data['smiles'] = ['CC(=O)OC1=CC=CC=C1C(=O)O'] * len(drug_data)  # Aspirin SMILES
    
    # Create placeholder data for other sources
    disease_data = pd.DataFrame({
        'disease_name': drug_data['disease_name'].unique(),
        'gene_list': ['GENE1,GENE2'] * len(drug_data['disease_name'].unique())
    })
    
    protein_data = pd.DataFrame({
        'uniprot_id': ['P12345'] * 100,
        'sequence': ['MALWMRLLPLL'] * 100
    })
    
    return {
        'drug_data': drug_data,
        'disease_data': disease_data,
        'protein_data': protein_data
    }

def build_knowledge_graph(data: Dict[str, pd.DataFrame]) -> HeteroData:
    """Build the knowledge graph from multiple data sources"""
    kg_builder = KnowledgeGraphBuilder()
    
    # Add nodes
    kg_builder.add_drug_nodes(data['drug_data'])
    kg_builder.add_disease_nodes(data['disease_data'])
    kg_builder.add_protein_nodes(data['protein_data'])
    
    # Add edges (placeholder for demonstration)
    edge_data = {
        ('drug', 'interacts_with', 'protein'): torch.tensor([[0, 1], [1, 2]]),
        ('disease', 'associated_with', 'protein'): torch.tensor([[0, 1], [1, 2]])
    }
    kg_builder.add_edges(edge_data)
    
    return kg_builder.build_graph()

def train_model(data: HeteroData) -> HeteroGNN:
    """Train the model with hyperparameter optimization"""
    # Optimize hyperparameters
    best_params = optimize_hyperparameters(data)
    
    # Create model with best parameters
    model = HeteroGNN(
        hidden_channels=best_params['best_params']['hidden_channels'],
        out_channels=64,
        num_layers=best_params['best_params']['num_layers'],
        metadata=data.metadata()
    )
    
    # Train model
    trainer = DrugRepurposingTrainer(model, lr=best_params['best_params']['lr'])
    
    # Save trained model
    torch.save(model, 'trained_model.pt')
    
    return model

def main():
    # Load data
    print("Loading data...")
    data = load_data()
    
    # Build knowledge graph
    print("Building knowledge graph...")
    kg_data = build_knowledge_graph(data)
    
    # Train model
    print("Training model...")
    model = train_model(kg_data)
    
    # Start API server
    print("Starting API server...")
    uvicorn.run(app, host="0.0.0.0", port=8000)

if __name__ == "__main__":
    main() 