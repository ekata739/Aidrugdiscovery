import torch
import torch.nn.functional as F
from torch_geometric.data import HeteroData
import optuna
from typing import Dict, Tuple
import numpy as np
from knowledge_graph import HeteroGNN
from torch_geometric.loader import NeighborLoader
import pytorch_lightning as pl
from pytorch_lightning.callbacks import EarlyStopping, ModelCheckpoint

class DrugRepurposingTrainer(pl.LightningModule):
    def __init__(self, model: HeteroGNN, lr: float = 0.001):
        super().__init__()
        self.model = model
        self.lr = lr
        
    def forward(self, x_dict, edge_index_dict):
        return self.model(x_dict, edge_index_dict)
        
    def training_step(self, batch, batch_idx):
        drug_emb, disease_emb = self(batch.x_dict, batch.edge_index_dict)
        pred = torch.sigmoid(torch.sum(drug_emb * disease_emb, dim=1))
        loss = F.binary_cross_entropy(pred, batch.y)
        self.log('train_loss', loss)
        return loss
        
    def validation_step(self, batch, batch_idx):
        drug_emb, disease_emb = self(batch.x_dict, batch.edge_index_dict)
        pred = torch.sigmoid(torch.sum(drug_emb * disease_emb, dim=1))
        loss = F.binary_cross_entropy(pred, batch.y)
        self.log('val_loss', loss)
        return loss
        
    def configure_optimizers(self):
        return torch.optim.Adam(self.parameters(), lr=self.lr)

def objective(trial: optuna.Trial, data: HeteroData) -> float:
    # Hyperparameter optimization
    hidden_channels = trial.suggest_int('hidden_channels', 32, 256)
    num_layers = trial.suggest_int('num_layers', 2, 6)
    lr = trial.suggest_float('lr', 1e-4, 1e-2, log=True)
    
    # Create model
    model = HeteroGNN(
        hidden_channels=hidden_channels,
        out_channels=64,
        num_layers=num_layers,
        metadata=data.metadata()
    )
    
    # Create trainer
    trainer = DrugRepurposingTrainer(model, lr=lr)
    
    # Create data loaders
    train_loader = NeighborLoader(
        data,
        num_neighbors=[10, 10],
        batch_size=32,
        input_nodes=('drug', data['drug'].train_mask)
    )
    
    val_loader = NeighborLoader(
        data,
        num_neighbors=[10, 10],
        batch_size=32,
        input_nodes=('drug', data['drug'].val_mask)
    )
    
    # Train model
    pl_trainer = pl.Trainer(
        max_epochs=100,
        callbacks=[
            EarlyStopping(monitor='val_loss', patience=10),
            ModelCheckpoint(monitor='val_loss')
        ]
    )
    
    pl_trainer.fit(trainer, train_loader, val_loader)
    
    return pl_trainer.callback_metrics['val_loss'].item()

def optimize_hyperparameters(data: HeteroData, n_trials: int = 50) -> Dict:
    study = optuna.create_study(direction='minimize')
    study.optimize(lambda trial: objective(trial, data), n_trials=n_trials)
    
    return {
        'best_params': study.best_params,
        'best_value': study.best_value
    }

class FewShotLearner:
    def __init__(self, base_model: HeteroGNN):
        self.base_model = base_model
        
    def adapt_to_new_disease(self, support_set: HeteroData, 
                           query_set: HeteroData, n_epochs: int = 10):
        """Adapt the model to a new disease with few examples"""
        # Freeze base model parameters
        for param in self.base_model.parameters():
            param.requires_grad = False
            
        # Create adaptation layer
        self.adaptation_layer = torch.nn.Linear(64, 64)
        
        optimizer = torch.optim.Adam(self.adaptation_layer.parameters(), lr=0.001)
        
        for epoch in range(n_epochs):
            # Adaptation phase
            drug_emb, disease_emb = self.base_model(
                support_set.x_dict, 
                support_set.edge_index_dict
            )
            
            adapted_emb = self.adaptation_layer(disease_emb)
            loss = F.mse_loss(adapted_emb, disease_emb)
            
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            
            # Evaluation phase
            with torch.no_grad():
                drug_emb, disease_emb = self.base_model(
                    query_set.x_dict,
                    query_set.edge_index_dict
                )
                adapted_emb = self.adaptation_layer(disease_emb)
                pred = torch.sigmoid(torch.sum(drug_emb * adapted_emb, dim=1))
                
        return pred 