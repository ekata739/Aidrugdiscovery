import React, { useState } from 'react';
import {
  Box,
  Container,
  TextField,
  Button,
  Typography,
  Card,
  CardContent,
  CircularProgress,
  Alert,
  Chip,
  Divider,
  Grid,
} from '@mui/material';
import SearchIcon from '@mui/icons-material/Search';
import ScienceIcon from '@mui/icons-material/Science';

function Search() {
  const [searchTerm, setSearchTerm] = useState('');
  const [results, setResults] = useState(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);

  const handleSearch = async (e) => {
    e.preventDefault();
    if (!searchTerm.trim()) return;

    setLoading(true);
    setError(null);
    
    try {
      const response = await fetch('http://localhost:5001/search_diseases', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify({ query: searchTerm }),
      });

      if (!response.ok) {
        throw new Error('Search failed. Please try again.');
      }

      const data = await response.json();
      // Handle the results array structure
      setResults(data.results[0]);
    } catch (err) {
      setError(err.message);
      setResults(null);
    } finally {
      setLoading(false);
    }
  };

  const renderConfidenceChip = (confidence) => {
    let color = 'default';
    if (confidence >= 80) color = 'success';
    else if (confidence >= 60) color = 'primary';
    else if (confidence >= 40) color = 'warning';
    else color = 'error';

    return (
      <Chip
        label={`${confidence}% Confidence`}
        color={color}
        size="small"
        sx={{ ml: 1 }}
      />
    );
  };

  return (
    <Container maxWidth="lg" sx={{ py: 4 }}>
      <Box component="form" onSubmit={handleSearch} sx={{ mb: 4 }}>
        <Typography variant="h2" gutterBottom align="center">
          Search Rare Diseases
        </Typography>
        <Typography variant="body1" gutterBottom align="center" color="text.secondary" sx={{ mb: 4 }}>
          Enter a disease name to find potential treatments and research insights
        </Typography>
        <Box sx={{ display: 'flex', gap: 2, maxWidth: 600, mx: 'auto' }}>
          <TextField
            fullWidth
            variant="outlined"
            placeholder="Enter disease name..."
            value={searchTerm}
            onChange={(e) => setSearchTerm(e.target.value)}
          />
          <Button
            type="submit"
            variant="contained"
            disabled={loading}
            startIcon={loading ? <CircularProgress size={20} /> : <SearchIcon />}
          >
            Search
          </Button>
        </Box>
      </Box>

      {error && (
        <Alert severity="error" sx={{ mb: 4 }}>
          {error}
        </Alert>
      )}

      {results && (
        <Box>
          <Card sx={{ mb: 4 }}>
            <CardContent>
              <Typography variant="h4" gutterBottom>
                {results.disease_name}
              </Typography>
              <Typography variant="body1" paragraph>
                {results.description}
              </Typography>
            </CardContent>
          </Card>

          {results.potential_treatments?.candidates && results.potential_treatments.candidates.length > 0 ? (
            <Box>
              <Typography variant="h5" gutterBottom sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
                <ScienceIcon color="primary" />
                Potential Treatments
              </Typography>
              <Grid container spacing={3}>
                {results.potential_treatments.candidates.map((treatment, index) => (
                  <Grid item xs={12} key={index}>
                    <Card>
                      <CardContent>
                        <Box sx={{ display: 'flex', alignItems: 'center', mb: 1 }}>
                          <Typography variant="h6">
                            {treatment.name}
                          </Typography>
                          {renderConfidenceChip(treatment.confidence)}
                        </Box>
                        <Typography variant="body2" color="text.secondary" paragraph>
                          {treatment.reason}
                        </Typography>
                        <Divider sx={{ my: 1 }} />
                        <Typography variant="body2">
                          <strong>Mechanism:</strong> {treatment.mechanism}
                        </Typography>
                        {treatment.current_uses && treatment.current_uses.length > 0 && (
                          <Typography variant="body2" sx={{ mt: 1 }}>
                            <strong>Current Uses:</strong> {treatment.current_uses.join(', ')}
                          </Typography>
                        )}
                      </CardContent>
                    </Card>
                  </Grid>
                ))}
              </Grid>
            </Box>
          ) : (
            <Alert severity="info">
              No potential treatments found for this disease.
            </Alert>
          )}
        </Box>
      )}
    </Container>
  );
}

export default Search; 