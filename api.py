from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
from typing import List, Dict, Optional
import torch
from knowledge_graph import KnowledgeGraphBuilder, HeteroGNN
from model_trainer import DrugRepurposingTrainer, FewShotLearner
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
import networkx as nx
import uvicorn
from disease_database import DiseaseDatabase

app = FastAPI(title="Drug Repurposing API",
             description="AI-powered drug repurposing for rare diseases",
             version="1.0.0")

# Initialize components
kg_builder = KnowledgeGraphBuilder()
model = None  # Will be loaded after training

# Initialize disease database
disease_db = DiseaseDatabase()
disease_db.load_custom_data('drug_repurposing_data.csv')
disease_db.load_orphanet_data()
disease_db.load_nord_data()

class DrugRequest(BaseModel):
    disease_name: str
    gene_list: List[str]
    confidence_threshold: float = 0.7
    max_drugs: int = 20

class DrugResponse(BaseModel):
    drug_name: str
    confidence: float
    mechanism: str
    supporting_evidence: List[str]
    pathway_analysis: Dict[str, float]

class ExplainabilityResponse(BaseModel):
    important_genes: List[Dict[str, float]]
    key_pathways: List[Dict[str, float]]
    molecular_interactions: List[Dict[str, float]]

class DiseaseSearchRequest(BaseModel):
    query: str

class DiseaseSearchResponse(BaseModel):
    name: str
    source: str
    genes: List[str]
    drugs: List[Dict]
    confidence: float

class DiseaseDetailsResponse(BaseModel):
    name: str
    source: str
    genes: List[str]
    drugs: List[Dict]
    confidence: float
    prevalence: Optional[str] = None
    inheritance: Optional[str] = None
    description: Optional[str] = None
    symptoms: Optional[List[str]] = None

@app.on_event("startup")
async def startup_event():
    global model
    # Load trained model
    # model = torch.load("trained_model.pt")
    # For now, use a placeholder
    model = HeteroGNN(hidden_channels=128, out_channels=64, num_layers=3, 
                     metadata=([], []))

@app.post("/predict", response_model=List[DrugResponse])
async def predict_drugs(request: DrugRequest):
    try:
        # Convert gene list to features
        gene_features = kg_builder._get_gene_expression_features(
            [{"symbol": gene} for gene in request.gene_list]
        )
        
        # Get drug predictions
        drug_emb, disease_emb = model(
            {"drug": torch.randn(1000, 1024),  # Placeholder
             "disease": torch.tensor([gene_features])},
            {}  # Placeholder edge_index_dict
        )
        
        predictions = model.predict_drug_disease(drug_emb, disease_emb)
        
        # Filter and format results
        results = []
        for i, pred in enumerate(predictions):
            if pred > request.confidence_threshold:
                results.append(DrugResponse(
                    drug_name=f"Drug_{i}",
                    confidence=float(pred),
                    mechanism="Predicted mechanism",
                    supporting_evidence=["Evidence 1", "Evidence 2"],
                    pathway_analysis={"pathway1": 0.8, "pathway2": 0.6}
                ))
                
        return sorted(results, key=lambda x: x.confidence, reverse=True)[:request.max_drugs]
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/explain/{drug_name}/{disease_name}", response_model=ExplainabilityResponse)
async def explain_prediction(drug_name: str, disease_name: str):
    try:
        # Generate explainability metrics
        return ExplainabilityResponse(
            important_genes=[
                {"gene": "GENE1", "importance": 0.9},
                {"gene": "GENE2", "importance": 0.8}
            ],
            key_pathways=[
                {"pathway": "PATHWAY1", "relevance": 0.85},
                {"pathway": "PATHWAY2", "relevance": 0.75}
            ],
            molecular_interactions=[
                {"interaction": "INTERACTION1", "strength": 0.9},
                {"interaction": "INTERACTION2", "strength": 0.8}
            ]
        )
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/few_shot_learn")
async def few_shot_learning(support_set: Dict, query_set: Dict):
    try:
        few_shot_learner = FewShotLearner(model)
        predictions = few_shot_learner.adapt_to_new_disease(
            support_set, query_set
        )
        return {"predictions": predictions.tolist()}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/search_diseases", response_model=List[DiseaseSearchResponse])
async def search_diseases(request: DiseaseSearchRequest):
    """Search for diseases by name or gene"""
    results = disease_db.search_diseases(request.query)
    if not results:
        raise HTTPException(status_code=404, detail="No diseases found")
    return results

@app.get("/disease/{disease_name}", response_model=DiseaseDetailsResponse)
async def get_disease_details(disease_name: str):
    """Get detailed information about a specific disease"""
    details = disease_db.get_disease_details(disease_name)
    if not details:
        raise HTTPException(status_code=404, detail="Disease not found")
    return details

@app.get("/all_diseases", response_model=List[str])
async def get_all_diseases():
    """Get list of all diseases in the database"""
    return disease_db.get_all_diseases()

if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=8000) 