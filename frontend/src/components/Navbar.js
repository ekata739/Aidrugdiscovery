import React from 'react';
import { Link as RouterLink } from 'react-router-dom';
import {
  AppBar,
  Toolbar,
  Typography,
  Button,
  Box,
  Container,
} from '@mui/material';
import ScienceIcon from '@mui/icons-material/Science';

function Navbar() {
  return (
    <AppBar position="static" elevation={0} sx={{ backgroundColor: 'white' }}>
      <Container maxWidth="lg">
        <Toolbar disableGutters>
          <Box sx={{ display: 'flex', alignItems: 'center', flexGrow: 1 }}>
            <ScienceIcon sx={{ color: 'primary.main', fontSize: 32, mr: 1 }} />
            <Typography
              variant="h6"
              component={RouterLink}
              to="/"
              sx={{
                color: 'primary.main',
                textDecoration: 'none',
                fontWeight: 700,
                letterSpacing: 0.5,
              }}
            >
              AI Drug Discovery
            </Typography>
          </Box>
          
          <Box sx={{ display: 'flex', gap: 2 }}>
            <Button
              component={RouterLink}
              to="/"
              color="primary"
              sx={{ fontWeight: 500 }}
            >
              Home
            </Button>
            <Button
              component={RouterLink}
              to="/search"
              color="primary"
              sx={{ fontWeight: 500 }}
            >
              Search
            </Button>
            <Button
              component={RouterLink}
              to="/about"
              color="primary"
              sx={{ fontWeight: 500 }}
            >
              About
            </Button>
          </Box>
        </Toolbar>
      </Container>
    </AppBar>
  );
}

export default Navbar; 