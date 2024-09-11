#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int main()
{
  int i, j;
  double Re;
  int m = 100,n = 100;
  printf("\nEnter the value of Reynolds Number(Re):");
  scanf("%lf", &Re);
  double dx = 1.0/(m-1);
  double dy = 1.0/(n-1);

  double B = dx / dy;
  double N=(m-2)*(n-2);//number of interior points

  double u[m][n], v[m][n];
  double psi[m][n], psi_old[m][n];
  double W[m][n], W_old[m][n];
  double psi_Err, W_Err;
  int iteration = 0;

  // Boundary conditions
  for (j = 0; j < n; j++)
  {
    for (i = 0; i < m; i++)
    {
      if (i == 0) // Left BC
      {
        u[i][j] = 0.0, v[i][j] = 0.0, psi[i][j] = 0.0;
      }
      else if (i == (m - 1)) // Right BC
      {
        u[i][j] = 0.0, v[i][j] = 0.0, psi[i][j] = 0.0;
      }
      else if (j == 0) // Bottom BC
      {
        u[i][j] = 0.0, v[i][j] = 0.0, psi[i][j] = 0.0;
      }
      else if (j == (n - 1)) // Top BC
      {
        u[i][j] = 1.0, v[i][j] = 0.0, psi[i][j] = 0.0;
      }
      else // Interior pts
      {
        u[i][j] = 0.0, v[i][j] = 0.0, psi[i][j] = 0.0;
      }
    }
  }
  // W intialisation
   for(j=0;j<n;j++)
    {
        for(i=0;i<m;i++)
        {
            if(j==0)
            {
                W[i][j]=(2.0/pow(dy,2))*(psi[i][j]-psi[i][j+1]);

            }
            else if(j==(n-1))
            {
                W[i][j]=((2.0/pow(dy,2))*(psi[i][j]-psi[i][j-1]))-((2.0/dy)*u[i][j]);
            }
            else if(i==0)
            {
                W[i][j]=(2.0/pow(dx,2))*(psi[i][j]-psi[i+1][j]);
            }
            else if(i==(m-1))
            {
                W[i][j]=(2.0/pow(dx,2))*(psi[i][j]-psi[i-1][j]);
            }
             else 
            {
                W[i][j]=0.0;
            }
        }
    }
  do
  {
    for (j = 0; j < n; j++)
    {
      for (i = 0; i < m; i++) 
      {
        psi_old[i][j] = psi[i][j];
        W_old[i][j] = W[i][j];
      }
    }
    // stream function solution 
    for (j = 1; j < (n - 1); j++)
    {
      for (i = 1; i < (m - 1); i++)
      {
        psi[i][j] = ((1.0 / (2.0 * (1 + (B*B)))) * (psi[i + 1][j] + psi[i - 1][j] + (pow(B, 2) * (psi[i][j + 1] + psi[i][j - 1])) + ((dx*dx) * W[i][j])));
      }
    }
    // vorticity solution 
    for (j = 1; j < (n - 1); j++)
    {
      for (i = 1; i < (m - 1); i++)
      {
        W[i][j] = ((1.0 / (2.0 * (1.0 + (B*B)))) * (((1.0 - ((psi[i][j + 1] - psi[i][j - 1]) * ((B * Re) / 4.0))) * W[i + 1][j])
        + ((1.0 + ((psi[i][j + 1] - psi[i][j - 1]) * ((B * Re) / 4.0))) * W[i - 1][j])
        + ((1.0 + ((psi[i + 1][j] - psi[i - 1][j]) * (Re / (4.0 * B)))) * ((B*B) * W[i][j + 1]))
        + ((1.0 - ((psi[i + 1][j] - psi[i - 1][j]) * (Re / (4.0 * B)))) * ((B*B) * W[i][j - 1]))));
      }
    }
    
    // to update vorticity at boundaries
    for(j=0;j<n;j++)
       {
        for(i=0;i<m;i++)
        {
            if(j==0)
            {
                W[i][j]=(2.0/(dy*dy))*(psi[i][j]-psi[i][j+1]);
            }
            else if(j==(n-1))
            {
                W[i][j]=((2.0/(dy*dy))*(psi[i][j]-psi[i][j-1]))-((2.0/dy)*u[i][j]);
            }
            else if (i==0)
            {
                W[i][j]=(2.0/(dx*dx))*(psi[i][j]-psi[i+1][j]);
            }
            else if (i==(m-1))
            {
                W[i][j]=(2.0/(dx*dx))*(psi[i][j]-psi[i-1][j]);
            }
        }
       }
    // Error calculation
    psi_Err = 0.0, W_Err = 0.0;

    for (j = 1; j < (n - 1); j++)
    {
      for (i = 1; i < (m - 1); i++)
      {
        psi_Err += (pow((psi[i][j] - psi_old[i][j]), 2.0));
        W_Err += (pow((W[i][j] - W_old[i][j]), 2.0));
      }
    }
    psi_Err = sqrt(psi_Err / N); 
    W_Err = sqrt(W_Err /N); 

    printf("iteration=%d\t", iteration);
    printf("error=%.20lf\n", psi_Err);
    printf("error=%.20lf\n", W_Err);
    printf("\n\n");

    iteration++;
  } 
while (psi_Err > 1.0e-6 || W_Err > 1.0e-6);

  // To Update velocity u and v
  for (j = 1; j < (n - 1); j++)
  {
    for (i = 1; i < (m - 1); i++)
    {
      
      u[i][j] = (0.5/dy)*(psi[i][j + 1] - psi[i][j - 1]); 
      v[i][j] = (-0.5/dx)*(psi[i + 1][j] - psi[i - 1][j]); 
    }
  }
  FILE *file1 ;
  double x,y;
  file1 = fopen("All DATA.dat", "w");

  fprintf(file1, "Grid_Size I=%d, J=%d\n",m,n);
  for (j = 0; j < n; j++)
  {
    y = j*dy;
    for (i = 0; i < m; i++)
    {
      x= i*dx;
      fprintf(file1, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", x,y,u[i][j],v[i][j], psi[i][j],W[i][j]);
    }
  }
 FILE *file2;
	file2 = fopen("u_along_vert_centerline.dat", "w");
	for (j = 0; j < n; j++)
		fprintf(file2, "%lf \t%lf \n", u[50][j], j * dy);

	FILE *file3;
	file3 = fopen("v_along_hzt_centerline.dat", "w");
	for (i = 0; i < m; i++)
		fprintf(file3, "%lf \t%lf \n", v[i][50], i * dx);

  FILE *file4;
file4=fopen("Stream.plt","w");
fprintf(file4,"VARIABLES = \"X\", \"Y\", \"PHI\"\n");
fprintf(file4,"ZONE T = \"BLOCK1\", I = 100, J = 100, F = POINT\n\n");
for(int i=1;i<=m;i++){
	for(int j=1;j<=n;j++){
		fprintf(file4,"%1f \t %1f \t %1f \n",i*dx,j*dy,psi[i][j]);}
		}

FILE *file5;
file5=fopen("Vorticity.plt","w");
fprintf(file5,"VARIABLES = \"X\", \"Y\", \"W\"\n");
fprintf(file5,"ZONE T = \"BLOCK1\", I = 100, J = 100, F = POINT\n\n");
for(int i=1;i<=m;i++){
	for(int j=1;j<=n;j++){
		fprintf(file5,"%1f \t %1f \t %1f \n",i*dx,j*dy,W[i][j]);}
		}
    
  fclose(file1);
  fclose(file2);
  fclose(file3);
  return 0;
}

