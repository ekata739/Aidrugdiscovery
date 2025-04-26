import torch
from torch_geometric.data import HeteroData
from torch_geometric.nn import HeteroConv, GCNConv
import networkx as nx
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
from typing import Dict, List, Tuple
import pandas as pd
from mygene import MyGeneInfo
import gseapy as gp

class KnowledgeGraphBuilder:
    def __init__(self):
        self.mg = MyGeneInfo()
        self.hetero_data = HeteroData()
        
    def add_drug_nodes(self, drug_data: pd.DataFrame):
        """Add drug nodes with molecular graph features"""
        drug_features = []
        for _, row in drug_data.iterrows():
            mol = Chem.MolFromSmiles(row['smiles'])
            if mol:
                # Generate Morgan fingerprints
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
                drug_features.append(np.array(fp))
            else:
                drug_features.append(np.zeros(1024))
        
        self.hetero_data['drug'].x = torch.tensor(drug_features, dtype=torch.float)
        
    def add_disease_nodes(self, disease_data: pd.DataFrame):
        """Add disease nodes with gene expression features"""
        disease_features = []
        for _, row in disease_data.iterrows():
            # Get associated genes from MyGene
            genes = self.mg.query(row['disease_name'], 
                                fields='entrezgene,symbol',
                                species='human')
            # Convert to gene expression features
            features = self._get_gene_expression_features(genes)
            disease_features.append(features)
            
        self.hetero_data['disease'].x = torch.tensor(disease_features, dtype=torch.float)
        
    def add_protein_nodes(self, protein_data: pd.DataFrame):
        """Add protein nodes with sequence features"""
        protein_features = []
        for _, row in protein_data.iterrows():
            # Get protein sequence features
            features = self._get_protein_features(row['uniprot_id'])
            protein_features.append(features)
            
        self.hetero_data['protein'].x = torch.tensor(protein_features, dtype=torch.float)
        
    def add_edges(self, edge_data: Dict[str, torch.Tensor]):
        """Add various types of edges to the graph"""
        for edge_type, edge_index in edge_data.items():
            self.hetero_data[edge_type].edge_index = edge_index
            
    def _get_gene_expression_features(self, genes: List[Dict]) -> np.ndarray:
        """Convert gene list to expression features"""
        # Implement gene expression feature extraction
        return np.random.rand(1000)  # Placeholder
        
    def _get_protein_features(self, uniprot_id: str) -> np.ndarray:
        """Get protein sequence features"""
        # Implement protein feature extraction
        return np.random.rand(1000)  # Placeholder
        
    def build_graph(self) -> HeteroData:
        """Build the final heterogeneous graph"""
        return self.hetero_data

class HeteroGNN(torch.nn.Module):
    def __init__(self, hidden_channels: int, out_channels: int, 
                 num_layers: int, metadata: Tuple[List[str], List[Tuple[str, str, str]]]):
        super().__init__()
        
        self.convs = torch.nn.ModuleList()
        for _ in range(num_layers):
            conv = HeteroConv({
                edge_type: GCNConv(-1, hidden_channels)
                for edge_type in metadata[1]
            })
            self.convs.append(conv)
            
        self.lin = torch.nn.Linear(hidden_channels, out_channels)
        
    def forward(self, x_dict, edge_index_dict):
        for conv in self.convs:
            x_dict = conv(x_dict, edge_index_dict)
            x_dict = {key: torch.relu(x) for key, x in x_dict.items()}
            
        return self.lin(x_dict['drug']), self.lin(x_dict['disease'])
        
    def predict_drug_disease(self, drug_features: torch.Tensor, 
                           disease_features: torch.Tensor) -> torch.Tensor:
        """Predict drug-disease relationships"""
        drug_emb = self.lin(drug_features)
        disease_emb = self.lin(disease_features)
        return torch.sigmoid(torch.sum(drug_emb * disease_emb, dim=1)) 