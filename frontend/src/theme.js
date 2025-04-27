import { extendTheme } from '@chakra-ui/react';

const theme = extendTheme({
  colors: {
    brand: {
      50: '#e5f0ff',
      100: '#b8d4ff',
      200: '#8ab8ff',
      300: '#5c9cff',
      400: '#2e80ff',
      500: '#0064ff',
      600: '#004ecc',
      700: '#003999',
      800: '#002466',
      900: '#000f33',
    },
  },
  fonts: {
    heading: 'Inter, sans-serif',
    body: 'Inter, sans-serif',
  },
  styles: {
    global: {
      body: {
        bg: 'gray.50',
      },
    },
  },
  components: {
    Button: {
      defaultProps: {
        colorScheme: 'brand',
      },
    },
  },
});

export default theme; 