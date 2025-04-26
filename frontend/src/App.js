import React, { useState, useEffect } from 'react';
import {
  Container,
  Grid,
  Paper,
  Typography,
  TextField,
  Button,
  Card,
  CardContent,
  Chip,
  CircularProgress,
  Box,
  Divider,
  List,
  ListItem,
  ListItemText,
  Avatar,
  Snackbar,
  Alert,
} from '@mui/material';
import { Search as SearchIcon, Science as ScienceIcon, LocalHospital as HospitalIcon } from '@mui/icons-material';
import { LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, Legend, ResponsiveContainer } from 'recharts';

const App = () => {
  const [searchQuery, setSearchQuery] = useState('');
  const [results, setResults] = useState([]);
  const [loading, setLoading] = useState(false);
  const [selectedDisease, setSelectedDisease] = useState(null);
  const [allDiseases, setAllDiseases] = useState([]);
  const [error, setError] = useState(null);

  useEffect(() => {
    fetchAllDiseases();
  }, []);

  const fetchAllDiseases = async () => {
    try {
      console.log('Fetching all diseases...');
      const response = await fetch('http://localhost:8000/all_diseases');
      const data = await response.json();
      
      if (!response.ok) {
        throw new Error(data.error || `HTTP error! status: ${response.status}`);
      }
      
      console.log('All diseases:', data);
      setAllDiseases(data);
    } catch (error) {
      console.error('Error fetching diseases:', error);
      setError(`Failed to fetch diseases: ${error.message}. Please check if the backend is running at http://localhost:8000`);
    }
  };

  const handleSearch = async () => {
    if (!searchQuery) return;
    setLoading(true);
    setError(null);
    try {
      console.log('Searching for:', searchQuery);
      const response = await fetch('http://localhost:8000/search_diseases', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({ query: searchQuery }),
      });
      
      const data = await response.json();
      
      if (!response.ok) {
        throw new Error(data.error || `HTTP error! status: ${response.status}`);
      }
      
      console.log('Search results:', data);
      setResults(data);
      if (data.length > 0) {
        setSelectedDisease(data[0]);
      } else {
        setError('No results found. Try a different search term.');
      }
    } catch (error) {
      console.error('Error searching:', error);
      setError(`Search failed: ${error.message}. Please check if the backend is running at http://localhost:8000`);
    }
    setLoading(false);
  };

  const DiseaseCard = ({ disease }) => (
    <Card 
      sx={{ 
        mb: 2, 
        cursor: 'pointer',
        transition: 'transform 0.2s',
        '&:hover': {
          transform: 'scale(1.02)',
          boxShadow: 3,
        }
      }}
      onClick={() => setSelectedDisease(disease)}
    >
      <CardContent>
        <Typography variant="h6" gutterBottom>
          {disease.name}
        </Typography>
        <Box sx={{ display: 'flex', gap: 1, flexWrap: 'wrap' }}>
          {disease.drugs.map((drug, index) => (
            <Chip
              key={index}
              label={drug.drug_name}
              color="primary"
              variant="outlined"
              size="small"
            />
          ))}
        </Box>
      </CardContent>
    </Card>
  );

  const DrugDetails = ({ drug }) => (
    <Paper sx={{ p: 2, mb: 2 }}>
      <Typography variant="h6" gutterBottom>
        {drug.drug_name}
      </Typography>
      <Typography variant="body2" color="text.secondary" gutterBottom>
        Confidence: {Math.round(drug.confidence * 100)}%
      </Typography>
      <Typography variant="body1" paragraph>
        {drug.reason}
      </Typography>
      <Button
        variant="outlined"
        size="small"
        href={`https://pubmed.ncbi.nlm.nih.gov/${drug.pmid}`}
        target="_blank"
      >
        View Research
      </Button>
    </Paper>
  );

  return (
    <Container maxWidth="lg" sx={{ py: 4 }}>
      <Paper sx={{ p: 4, mb: 4, background: 'linear-gradient(45deg, #2196F3 30%, #21CBF3 90%)' }}>
        <Typography variant="h3" component="h1" gutterBottom color="white">
          AI Drug Repurposing Platform
        </Typography>
        <Typography variant="h6" color="white" paragraph>
          Discover potential treatments for rare diseases using advanced AI
        </Typography>
      </Paper>

      <Grid container spacing={3}>
        <Grid item xs={12} md={4}>
          <Paper sx={{ p: 2 }}>
            <Box sx={{ display: 'flex', mb: 2 }}>
              <TextField
                fullWidth
                label="Search Diseases"
                variant="outlined"
                value={searchQuery}
                onChange={(e) => setSearchQuery(e.target.value)}
                onKeyPress={(e) => e.key === 'Enter' && handleSearch()}
              />
              <Button
                variant="contained"
                onClick={handleSearch}
                sx={{ ml: 1 }}
                disabled={loading}
              >
                <SearchIcon />
              </Button>
            </Box>

            {loading ? (
              <Box sx={{ display: 'flex', justifyContent: 'center', p: 3 }}>
                <CircularProgress />
              </Box>
            ) : (
              <List>
                {results.map((disease, index) => (
                  <ListItem key={index}>
                    <DiseaseCard disease={disease} />
                  </ListItem>
                ))}
              </List>
            )}
          </Paper>
        </Grid>

        <Grid item xs={12} md={8}>
          {selectedDisease ? (
            <Paper sx={{ p: 3 }}>
              <Typography variant="h4" gutterBottom>
                {selectedDisease.name}
              </Typography>
              <Divider sx={{ mb: 3 }} />

              <Grid container spacing={2}>
                <Grid item xs={12}>
                  <Typography variant="h6" gutterBottom>
                    Potential Treatments
                  </Typography>
                  {selectedDisease.drugs.map((drug, index) => (
                    <DrugDetails key={index} drug={drug} />
                  ))}
                </Grid>

                <Grid item xs={12}>
                  <Typography variant="h6" gutterBottom>
                    Confidence Distribution
                  </Typography>
                  <Box sx={{ height: 300 }}>
                    <ResponsiveContainer width="100%" height="100%">
                      <LineChart
                        data={selectedDisease.drugs.map(drug => ({
                          name: drug.drug_name,
                          confidence: drug.confidence * 100
                        }))}
                      >
                        <CartesianGrid strokeDasharray="3 3" />
                        <XAxis dataKey="name" />
                        <YAxis />
                        <Tooltip />
                        <Legend />
                        <Line
                          type="monotone"
                          dataKey="confidence"
                          stroke="#8884d8"
                          activeDot={{ r: 8 }}
                        />
                      </LineChart>
                    </ResponsiveContainer>
                  </Box>
                </Grid>
              </Grid>
            </Paper>
          ) : (
            <Paper sx={{ p: 3, textAlign: 'center' }}>
              <ScienceIcon sx={{ fontSize: 60, color: 'text.secondary', mb: 2 }} />
              <Typography variant="h6" color="text.secondary">
                Search for a disease to see potential treatments
              </Typography>
            </Paper>
          )}
        </Grid>
      </Grid>

      <Snackbar 
        open={!!error} 
        autoHideDuration={6000} 
        onClose={() => setError(null)}
        anchorOrigin={{ vertical: 'bottom', horizontal: 'center' }}
      >
        <Alert onClose={() => setError(null)} severity="error" sx={{ width: '100%' }}>
          {error}
        </Alert>
      </Snackbar>
    </Container>
  );
};

export default App; 