import React from 'react';
import { useNavigate } from 'react-router-dom';
import {
  Box,
  Button,
  Container,
  Grid,
  Typography,
  Card,
  CardContent,
  CardMedia,
} from '@mui/material';
import SearchIcon from '@mui/icons-material/Search';
import BiotechIcon from '@mui/icons-material/Biotech';
import PrecisionManufacturingIcon from '@mui/icons-material/PrecisionManufacturing';
import VerifiedIcon from '@mui/icons-material/Verified';

const features = [
  {
    title: 'Advanced Disease Search',
    description: 'Search through our comprehensive database of rare diseases with detailed information and descriptions.',
    icon: <SearchIcon sx={{ fontSize: 40 }} />,
  },
  {
    title: 'Drug Discovery AI',
    description: 'Leverage AI-powered algorithms to identify potential treatments and drug candidates.',
    icon: <PrecisionManufacturingIcon sx={{ fontSize: 40 }} />,
  },
  {
    title: 'Scientific Research',
    description: 'Access peer-reviewed research and clinical studies supporting treatment recommendations.',
    icon: <BiotechIcon sx={{ fontSize: 40 }} />,
  },
  {
    title: 'Validated Results',
    description: 'Get confidence scores and validation metrics for all suggested treatments and matches.',
    icon: <VerifiedIcon sx={{ fontSize: 40 }} />,
  },
];

function Home() {
  const navigate = useNavigate();

  return (
    <Box component="main">
      {/* Hero Section */}
      <Box
        sx={{
          background: 'linear-gradient(45deg, #2196f3 30%, #21CBF3 90%)',
          color: 'white',
          py: 8,
          position: 'relative',
          overflow: 'hidden',
        }}
      >
        <Container maxWidth="lg">
          <Grid container spacing={4} alignItems="center">
            <Grid item xs={12} md={6}>
              <Typography variant="h1" gutterBottom>
                AI-Powered Drug Discovery
              </Typography>
              <Typography variant="h5" paragraph sx={{ mb: 4 }}>
                Accelerating rare disease treatment discovery through artificial intelligence and machine learning.
              </Typography>
              <Button
                variant="contained"
                size="large"
                onClick={() => navigate('/search')}
                sx={{
                  bgcolor: 'white',
                  color: 'primary.main',
                  '&:hover': {
                    bgcolor: 'grey.100',
                  },
                }}
              >
                Start Searching
              </Button>
            </Grid>
            <Grid item xs={12} md={6}>
              <Box
                component="img"
                src="/hero-image.svg"
                alt="AI Drug Discovery Illustration"
                sx={{
                  width: '100%',
                  maxWidth: 600,
                  height: 'auto',
                }}
              />
            </Grid>
          </Grid>
        </Container>
      </Box>

      {/* Features Section */}
      <Container maxWidth="lg" sx={{ py: 8 }}>
        <Typography variant="h2" align="center" gutterBottom>
          Features
        </Typography>
        <Grid container spacing={4} sx={{ mt: 2 }}>
          {features.map((feature, index) => (
            <Grid item xs={12} sm={6} md={3} key={index}>
              <Card
                sx={{
                  height: '100%',
                  display: 'flex',
                  flexDirection: 'column',
                  alignItems: 'center',
                  textAlign: 'center',
                  p: 2,
                }}
              >
                <Box
                  sx={{
                    display: 'flex',
                    alignItems: 'center',
                    justifyContent: 'center',
                    width: 80,
                    height: 80,
                    borderRadius: '50%',
                    bgcolor: 'primary.light',
                    color: 'primary.main',
                    mb: 2,
                  }}
                >
                  {feature.icon}
                </Box>
                <Typography variant="h6" gutterBottom>
                  {feature.title}
                </Typography>
                <Typography variant="body2" color="text.secondary">
                  {feature.description}
                </Typography>
              </Card>
            </Grid>
          ))}
        </Grid>
      </Container>
    </Box>
  );
}

export default Home; 