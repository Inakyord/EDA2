# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <omp.h>

int main ( );
void ccopy ( int n, double x[], double y[] );
void cfft2 ( int n, double x[], double y[], double w[], double sgn );
void cffti ( int n, double w[] );
double ggl ( double *ds );
void step ( int n, int mj, double a[], double b[], double c[], double d[],
double w[], double sgn );
void timestamp ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Propósito:

    MAIN es el programa principal para FFT_OPENMP.

  Discusion:

    El vector "complejo" A es almacenado como un vector double B.

    Las entradas del vector "complejo" A[I] son almacenadas como:
        B[I*2+0], la parte real,
        B[I*2+1], la parte imaginaria.

  Modificado: (Original)

    20 Marzo 2009

  Autor:

    Version C original por Wesley Petersen.
    Ultima versión por John Burkardt.

  Referencia:

    Wesley Petersen, Peter Arbenz,
    Introduction to Parallel Computing - A practical guide with examples in C,
    Oxford University Press,
    ISBN: 0-19-851576-6,
    LC: QA76.58.P47.
*/
{
  double error;
  int first;
  double flops;
  double fnm1;
  int i;
  int icase;
  int it;
  int ln2;
  int ln2_max = 20;
  double mflops;
  int n;
  int nits = 1000;
  static double seed;
  double sgn;
  double *w;
  double wtime;
  double *x;
  double *y;
  double *z;
  double z0;
  double z1;

  timestamp ( );
  printf ( "\n" );
  printf ( "FFT_OPENMP\n" );
  printf ( "  C/OpenMP version\n" );
  printf ( "\n" );
  printf ( "  Muestra la implementacion de la Transformada Rapida de Fourier\n");
  printf ( "  de un vector complejo, utilizando OpenMP para su ejecucion en paralelo.\n");

  printf ( "\n" );
  printf ( "  Numero de procesadores disponibles = %d\n", omp_get_num_procs ( ) );
  printf ( "  Numero de hilos =              %d\n", omp_get_max_threads ( ) );
/*
  Preparando las pruebas.
*/
  printf ( "\n" );
  printf ( "  Revision de precision:\n" );
  printf ( "\n" );
  printf ( "    FFT ( FFT ( X(1:N) ) ) == N * X(1:N)\n" );
  printf ( "\n" );
  printf ( "             N      NITS    Error         Tiempo        Tiempo/Llamada     MFLOPS\n" );
  printf ( "\n" );

  seed  = 331.0;
  n = 15;
/*
  LN2 es el logaritmo base 2 de N. Cada incremento de LN2 duplica N.
*/
  for ( ln2 = 1; ln2 <= ln2_max; ln2++ )
  {
    n = 2 * n;
/*
  Asigna almacenamiento para los arreglos complejos W, X, Y, Z.

  Manejamos la aritmética compleja,
  y se almacena un numero complejo como una pareja de doubles, un vector complejo como un arreglo
  de doble dimensión del cual su segunda dimension es 2.

*/
    w = ( double * ) malloc (     n * sizeof ( double ) );
    x = ( double * ) malloc ( 2 * n * sizeof ( double ) );
    y = ( double * ) malloc ( 2 * n * sizeof ( double ) );
    z = ( double * ) malloc ( 2 * n * sizeof ( double ) );

    first = 1;

    for ( icase = 0; icase < 2; icase++ )
    {
      if ( first )
      {
        for ( i = 0; i < 2 * n; i = i + 2 )
        {
          z0 = ggl ( &seed );
          z1 = ggl ( &seed );
          x[i] = z0;
          z[i] = z0;
          x[i+1] = z1;
          z[i+1] = z1;
        }
      }
      else
      {
# pragma omp parallel \
    shared ( n, x, z ) \
    private ( i, z0, z1 )

# pragma omp for nowait

        for ( i = 0; i < 2 * n; i = i + 2 )
        {
          z0 = 0.0;              /* parte real del arreglo */
          z1 = 0.0;              /* parte imaginaria del arreglo */
          x[i] = z0;
          z[i] = z0;           /* copia del dato real inicial */
          x[i+1] = z1;
          z[i+1] = z1;         /* copia del dato imaginario inicial */
        }
      }
/*
  Inicializa las tablas de seno y coseno.
*/
      cffti ( n, w );
/*
  Transforma hacia adelante, atras
*/
      if ( first )
      {
        sgn = + 1.0;
        cfft2 ( n, x, y, w, sgn );
        sgn = - 1.0;
        cfft2 ( n, y, x, w, sgn );
/*
  Los resultados deberían ser los mismos que los datos iniciales multiplicados por N.
*/
        fnm1 = 1.0 / ( double ) n;
        error = 0.0;
        for ( i = 0; i < 2 * n; i = i + 2 )
        {
          error = error
          + pow ( z[i]   - fnm1 * x[i], 2 )
          + pow ( z[i+1] - fnm1 * x[i+1], 2 );
        }
        error = sqrt ( fnm1 * error );
        printf ( "  %12d  %8d  %12e", n, nits, error );
        first = 0;
      }
      else
      {
        wtime = omp_get_wtime ( );
        for ( it = 0; it < nits; it++ )
        {
          sgn = + 1.0;
          cfft2 ( n, x, y, w, sgn );
          sgn = - 1.0;
          cfft2 ( n, y, x, w, sgn );
        }
        wtime = omp_get_wtime ( ) - wtime;

        flops = 2.0 * ( double ) nits
          * ( 5.0 * ( double ) n * ( double ) ln2 );

        mflops = flops / 1.0E+06 / wtime;

        printf ( "  %12e  %12e  %12f\n", wtime, wtime / ( double ) ( 2 * nits ), mflops );
      }
    }
    if ( ( ln2 % 4 ) == 0 )
    {
      nits = nits / 10;
    }
    if ( nits < 1 )
    {
      nits = 1;
    }
    free ( w );
    free ( x );
    free ( y );
    free ( z );
  }
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "FFT_OPENMP:\n" );
  printf ( "  Fin normal de ejecucion.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void ccopy ( int n, double x[], double y[] )

/******************************************************************************/
/*
  Propósito:

    CCOPY copia un vector complejo.

  Discusion:

    El vector "complejo" A[N] es almacenado como un vector double B[2*N].

    La entrada al vector "complejo" A[I] es almacenado como:

      B[I*2+0], la parte real,
      B[I*2+1], la parte imaginaria.

  Modificado: (Original)

    20 Marzo 2009

  Autor:

    Versión original en C por Wesley Petersen.
    Última versión en C por John Burkardt.

  Referencia:

    Wesley Petersen, Peter Arbenz,
    Introduction to Parallel Computing - A practical guide with examples in C,
    Oxford University Press,
    ISBN: 0-19-851576-6,
    LC: QA76.58.P47.

  Parametros:

    Entrada, int N, la longitud del vector.

    Entrada, double X[2*N], el vector a ser copiado.

    Salida, double Y[2*N], una copia de X.
*/
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    y[i*2+0] = x[i*2+0];
    y[i*2+1] = x[i*2+1];
   }
  return;
}
/******************************************************************************/

void cfft2 ( int n, double x[], double y[], double w[], double sgn )

/******************************************************************************/
/*
  Proposito:

    CFFT2 realiza la Transformada Rápida de Fourier compleja.

  Modificado: (Original)

    20 Marzo 2009

  Autor:

    Version C original por Wesley Petersen.
    Esta versión C por John Burkardt.

  Referencia:

    Wesley Petersen, Peter Arbenz,
    Introduction to Parallel Computing - A practical guide with examples in C,
    Oxford University Press,
    ISBN: 0-19-851576-6,
    LC: QA76.58.P47.

  Parametros:

    Entrada, int N, el tamaño del arreglo a ser transformado.

    Entrada/Salida, double X[2*N], los datos a ser transformados.
    En salida, los contenidos de X han sido sobre-escritos por información de trabajo.

    Salida, double Y[2*N], la FFT de X hacia adelante o hacia atrás.

    Entrada, double W[N], una tabla de senos y cosenos.

    Entrada, double SGN, es +1 para un FFT hacia adelante y -1 para FFT hacia atrás.
*/
{
  int j;
  int m;
  int mj;
  int tgle;

   m = ( int ) ( log ( ( double ) n ) / log ( 1.99 ) );
   mj   = 1;
/*
  Habilitando switch para el arreglo de trabajo.
  Toggling switch for work array.
*/
  tgle = 1;
  step ( n, mj, &x[0*2+0], &x[(n/2)*2+0], &y[0*2+0], &y[mj*2+0], w, sgn );

  if ( n == 2 )
  {
    return;
  }

  for ( j = 0; j < m - 2; j++ )
  {
    mj = mj * 2;
    if ( tgle )
    {
      step ( n, mj, &y[0*2+0], &y[(n/2)*2+0], &x[0*2+0], &x[mj*2+0], w, sgn );
      tgle = 0;
    }
    else
    {
      step ( n, mj, &x[0*2+0], &x[(n/2)*2+0], &y[0*2+0], &y[mj*2+0], w, sgn );
      tgle = 1;
    }
  }
/*
  Ultima pasada por los datos: mueve Y a X si es necesario.
*/
  if ( tgle )
  {
    ccopy ( n, y, x );
  }

  mj = n / 2;
  step ( n, mj, &x[0*2+0], &x[(n/2)*2+0], &y[0*2+0], &y[mj*2+0], w, sgn );

  return;
}
/******************************************************************************/

void cffti ( int n, double w[] )

/******************************************************************************/
/*
  Propósito:

    CFFTI prepara las tablas de senos y cosenos necesarias para el cálculo de la FFT.

  Modificado: (Original)

    20 Marzo 2009

  Autor:

    Version C original por Wesley Petersen.
    Esta versión C por John Burkardt.

  Referencia:

    Wesley Petersen, Peter Arbenz,
    Introduction to Parallel Computing - A practical guide with examples in C,
    Oxford University Press,
    ISBN: 0-19-851576-6,
    LC: QA76.58.P47.

  Parametros:

    Entrada, int N, el tamaño del arreglo a ser transformado.

    Salida, double W[N], una tabla de senos y cosenos.
*/
{
  double arg;
  double aw;
  int i;
  int n2;
  const double pi = 3.141592653589793;

  n2 = n / 2;
  aw = 2.0 * pi / ( ( double ) n );

# pragma omp parallel \
    shared ( aw, n, w ) \
    private ( arg, i )

# pragma omp for nowait

  for ( i = 0; i < n2; i++ )
  {
    arg = aw * ( ( double ) i );
    w[i*2+0] = cos ( arg );
    w[i*2+1] = sin ( arg );
  }
  return;
}
/******************************************************************************/

double ggl ( double *seed )

/******************************************************************************/
/*
  Propósito:

    GGL genera números reales pseudo-aleatorios uniformemente distribuidos entre [0,1].

  Modificado: (Original)

    20 Marzo 2009

  Autor:

    Version C original por Wesley Petersen, M Troyer, I Vattulainen.
    Esta versión C por John Burkardt.

  Referencia:

    Wesley Petersen, Peter Arbenz,
    Introduction to Parallel Computing - A practical guide with examples in C,
    Oxford University Press,
    ISBN: 0-19-851576-6,
    LC: QA76.58.P47.

  Parámetros:

    Entrada/Salida, double *SEED, utilizada como una semilla para la secuencia.

    Salida, double GGL, el siguiente valor pseudo-aleatorio.
*/
{
  double d2 = 0.2147483647e10;
  double t;
  double value;

  t = ( double ) *seed;
  t = fmod ( 16807.0 * t, d2 );
  *seed = ( double ) t;
  value = ( double ) ( ( t - 1.0 ) / ( d2 - 1.0 ) );

  return value;
}
/******************************************************************************/

void step ( int n, int mj, double a[], double b[], double c[],
  double d[], double w[], double sgn )

/******************************************************************************/
/*
  Propósito:

    STEP lleva a cabo un paso de la versión de trabajo de CFFT2.

  Modificado: (Original)

    20 Marzo 2009

  Autor:

    Version C original por Wesley Petersen.
    Esta versión C por John Burkardt.

  Referencia:

    Wesley Petersen, Peter Arbenz,
    Introduction to Parallel Computing - A practical guide with examples in C,
    Oxford University Press,
    ISBN: 0-19-851576-6,
    LC: QA76.58.P47.

  Parámetros:

*/
{
  double ambr;
  double ambu;
  int j;
  int ja;
  int jb;
  int jc;
  int jd;
  int jw;
  int k;
  int lj;
  int mj2;
  double wjw[2];

  mj2 = 2 * mj;
  lj  = n / mj2;

# pragma omp parallel \
    shared ( a, b, c, d, lj, mj, mj2, sgn, w ) \
    private ( ambr, ambu, j, ja, jb, jc, jd, jw, k, wjw )

# pragma omp for nowait

  for ( j = 0; j < lj; j++ )
  {
    jw = j * mj;
    ja  = jw;
    jb  = ja;
    jc  = j * mj2;
    jd  = jc;

    wjw[0] = w[jw*2+0];
    wjw[1] = w[jw*2+1];

    if ( sgn < 0.0 )
    {
      wjw[1] = - wjw[1];
    }

    for ( k = 0; k < mj; k++ )
    {
      c[(jc+k)*2+0] = a[(ja+k)*2+0] + b[(jb+k)*2+0];
      c[(jc+k)*2+1] = a[(ja+k)*2+1] + b[(jb+k)*2+1];

      ambr = a[(ja+k)*2+0] - b[(jb+k)*2+0];
      ambu = a[(ja+k)*2+1] - b[(jb+k)*2+1];

      d[(jd+k)*2+0] = wjw[0] * ambr - wjw[1] * ambu;
      d[(jd+k)*2+1] = wjw[1] * ambr + wjw[0] * ambu;
    }
  }
  return;
}
/******************************************************************************/

void timestamp ( void )

/******************************************************************************/
/*
  Propósito:

    TIMESTAMP imprime la fecha YMDHMS actual como un sello de tiempo.

  Ejemplo:

    31 Mayo 2001 09:45:54 AM

  Licensia:

    Este código es distribuid bajo la licensia GNU LGPL.

  Modificado:

    24 Septiembre 2003

  Autor:

    John Burkardt

  Parámetros:

    Ninguno
*/
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  printf ( "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
