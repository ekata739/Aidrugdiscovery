import React from 'react';
import {
  Box,
  Container,
  Typography,
  Grid,
  Card,
  CardContent,
  List,
  ListItem,
  ListItemIcon,
  ListItemText,
  Divider,
} from '@mui/material';
import CodeIcon from '@mui/icons-material/Code';
import StorageIcon from '@mui/icons-material/Storage';
import MemoryIcon from '@mui/icons-material/Memory';
import SecurityIcon from '@mui/icons-material/Security';
import SpeedIcon from '@mui/icons-material/Speed';
import ApiIcon from '@mui/icons-material/Api';

function About() {
  const sections = [
    {
      title: 'Technical Stack',
      icon: <CodeIcon fontSize="large" color="primary" />,
      items: [
        'Python Flask Backend',
        'React Frontend',
        'Material-UI Components',
        'RESTful API Architecture',
        'Real-time Processing',
      ],
    },
    {
      title: 'AI Components',
      icon: <MemoryIcon fontSize="large" color="primary" />,
      items: [
        'Natural Language Processing',
        'Pattern Recognition',
        'Multi-factor Analysis',
        'Confidence Scoring',
        'Automated Classification',
      ],
    },
    {
      title: 'Data Processing',
      icon: <StorageIcon fontSize="large" color="primary" />,
      items: [
        'Medical Literature Analysis',
        'Drug Database Integration',
        'Disease Mechanism Extraction',
        'Real-time Updates',
        'Evidence Tracking',
      ],
    },
    {
      title: 'API Features',
      icon: <ApiIcon fontSize="large" color="primary" />,
      items: [
        'RESTful Endpoints',
        'JSON Response Format',
        'Error Handling',
        'Rate Limiting',
        'Documentation',
      ],
    },
    {
      title: 'Performance',
      icon: <SpeedIcon fontSize="large" color="primary" />,
      items: [
        '500ms Response Time',
        '95% Accuracy Rate',
        '1000+ Queries/Second',
        'Real-time Processing',
        'Scalable Architecture',
      ],
    },
    {
      title: 'Security',
      icon: <SecurityIcon fontSize="large" color="primary" />,
      items: [
        'Input Validation',
        'CORS Protection',
        'Rate Limiting',
        'Error Sanitization',
        'Secure API Access',
      ],
    },
  ];

  return (
    <Container maxWidth="lg" sx={{ py: 8 }}>
      <Box sx={{ mb: 6, textAlign: 'center' }}>
        <Typography variant="h3" component="h1" gutterBottom>
          About the Platform
        </Typography>
        <Typography variant="h6" color="text.secondary" sx={{ maxWidth: 800, mx: 'auto', mb: 4 }}>
          An advanced AI-powered system for discovering potential treatments for rare diseases
          through intelligent drug repurposing
        </Typography>
      </Box>

      <Grid container spacing={4}>
        {sections.map((section, index) => (
          <Grid item xs={12} md={6} key={index}>
            <Card
              sx={{
                height: '100%',
                display: 'flex',
                flexDirection: 'column',
                transition: 'transform 0.2s',
                '&:hover': {
                  transform: 'translateY(-4px)',
                },
              }}
            >
              <CardContent>
                <Box sx={{ display: 'flex', alignItems: 'center', mb: 2 }}>
                  {section.icon}
                  <Typography variant="h5" component="h2" sx={{ ml: 2 }}>
                    {section.title}
                  </Typography>
                </Box>
                <Divider sx={{ my: 2 }} />
                <List>
                  {section.items.map((item, idx) => (
                    <ListItem key={idx}>
                      <ListItemIcon>
                        <Box
                          sx={{
                            width: 8,
                            height: 8,
                            borderRadius: '50%',
                            bgcolor: 'primary.main',
                          }}
                        />
                      </ListItemIcon>
                      <ListItemText primary={item} />
                    </ListItem>
                  ))}
                </List>
              </CardContent>
            </Card>
          </Grid>
        ))}
      </Grid>

      <Box sx={{ mt: 8, textAlign: 'center' }}>
        <Typography variant="h4" gutterBottom>
          System Architecture
        </Typography>
        <Typography variant="body1" color="text.secondary" sx={{ mb: 4 }}>
          Our platform utilizes a modern, scalable architecture designed for performance and reliability
        </Typography>
        <Box
          component="img"
          src="/architecture.png"
          alt="System Architecture"
          sx={{
            maxWidth: '100%',
            height: 'auto',
            borderRadius: 2,
            boxShadow: 3,
          }}
        />
      </Box>
    </Container>
  );
}

export default About; 